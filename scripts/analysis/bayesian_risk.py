"""Bayesian Risk Assessment Module (BRAM).

Implements a Beta-Binomial conjugate Bayesian model to estimate
the posterior probability that pipeline performance has degraded
below a clinically meaningful threshold. Augments deterministic
acceptance criteria with calibrated uncertainty estimates.

Uses exact conjugate updates (no normal approximation). The Beta
prior is subjective and the baseline is treated as fixed (its
uncertainty is ignored because baseline validations typically use
large N, making their variance negligible relative to verification
runs). The posterior computation itself is exact given the model
assumptions.
"""

from __future__ import annotations

import logging
import math

import numpy as np
from scipy.stats import beta as beta_dist

from analysis.verdicts import (  # noqa: F401
    BRAM_FLAG,
    BRAM_NOT_RUN,
    BRAM_PASS,
    FAIL,
    REVIEW_REQUIRED,
)

log = logging.getLogger(__name__)

# Default configuration values
DEFAULT_ALERT_THRESHOLD = 0.80
DEFAULT_PRIOR_CONCENTRATION = 20
DEFAULT_AGGREGATION_METHOD = "max"
DEFAULT_MIN_HISTORICAL_RUNS = 5
DEFAULT_DEGRADATION_THRESHOLDS = {
    "tier_1": 0.005,
    "tier_2": 0.01,
    "tier_3": 0.02,
}

# Minimum concentration to prevent degenerate priors
MIN_CONCENTRATION = 2.0


def beta_std(alpha: float, beta_param: float) -> float:
    """Compute the standard deviation of a Beta(alpha, beta) distribution.

    Args:
        alpha: Alpha parameter (> 0)
        beta_param: Beta parameter (> 0)

    Returns:
        Standard deviation of the Beta distribution
    """
    ab = alpha + beta_param
    return math.sqrt(alpha * beta_param / (ab * ab * (ab + 1)))


def compute_beta_prior(
    baseline_value: float,
    concentration: float,
) -> tuple[float, float]:
    """Compute Beta prior parameters centered at a baseline value.

    The prior encodes "we expect the true proportion to be near the
    baseline, with confidence proportional to concentration."

    Args:
        baseline_value: Baseline metric value (e.g., 0.999 for 99.9% sensitivity)
        concentration: Effective prior sample size (kappa). Higher values
            produce tighter priors.

    Returns:
        (alpha_prior, beta_prior)
    """
    concentration = max(concentration, MIN_CONCENTRATION)
    p = np.clip(baseline_value, 1e-6, 1.0 - 1e-6)
    alpha = p * concentration
    beta_param = (1.0 - p) * concentration
    return float(alpha), float(beta_param)


def compute_beta_posterior(
    alpha_prior: float,
    beta_prior: float,
    successes: int,
    trials: int,
) -> tuple[float, float]:
    """Compute Beta-Binomial conjugate posterior.

    posterior = Beta(alpha_prior + successes, beta_prior + trials - successes)

    Args:
        alpha_prior: Prior alpha parameter
        beta_prior: Prior beta parameter
        successes: Number of successes (e.g., TP count)
        trials: Total trials (e.g., TP + FN for sensitivity)

    Returns:
        (alpha_post, beta_post)
    """
    successes = max(0, successes)
    trials = max(0, trials)
    failures = trials - successes
    return alpha_prior + successes, beta_prior + failures


def compute_beta_tail(
    alpha_post: float,
    beta_post: float,
    baseline_value: float,
    threshold: float,
) -> float:
    """Compute P(p_true < baseline - threshold | data).

    This is the posterior probability that the true metric value has
    degraded by more than the threshold from the baseline.

    Args:
        alpha_post: Posterior alpha parameter
        beta_post: Posterior beta parameter
        baseline_value: Baseline metric value
        threshold: Clinically meaningful degradation magnitude

    Returns:
        Tail probability (0 to 1)
    """
    cutoff = baseline_value - threshold
    if cutoff <= 0.0:
        return 0.0
    if cutoff >= 1.0:
        return 1.0
    return float(beta_dist.cdf(cutoff, alpha_post, beta_post))


def aggregate_risks(
    weighted_risks: list[float],
    method: str = "max",
) -> float:
    """Aggregate per-stratum risk scores.

    Args:
        weighted_risks: List of per-stratum tail probabilities
        method: Aggregation method ("max", "weighted_mean", "product")

    Returns:
        Aggregate risk score
    """
    if not weighted_risks:
        return 0.0

    arr = np.array(weighted_risks, dtype=np.float64)

    if method == "max":
        return float(arr.max())
    elif method == "weighted_mean":
        return float(arr.mean())
    elif method == "product":
        # Joint probability of at least one failure: 1 - product(1 - risk_i)
        # Uses numpy.prod for numerical stability over iterative multiplication.
        return float(1.0 - np.prod(1.0 - arr))
    else:
        log.warning("Unknown aggregation method '%s', falling back to max", method)
        return float(arr.max())


class BRAMEngine:
    """Bayesian Risk Assessment Module engine.

    Computes per-stratum posterior risk estimates using a Beta-Binomial
    conjugate model and aggregates them for decision-making.

    Args:
        config: BRAM configuration dict (from bram_config.yaml)
        acceptance_engine: Optional AcceptanceEngine for tier assignment
    """

    def __init__(self, config: dict | None = None, acceptance_engine=None):
        self.config = config or {}
        self.acceptance_engine = acceptance_engine

        self.alert_threshold = self.config.get(
            "alert_threshold", DEFAULT_ALERT_THRESHOLD
        )
        self.aggregation_method = self.config.get(
            "aggregation_method", DEFAULT_AGGREGATION_METHOD
        )

        default_prior = self.config.get("default_prior", {})
        self.default_concentration = default_prior.get(
            "concentration", DEFAULT_PRIOR_CONCENTRATION
        )

        self.degradation_thresholds = self.config.get(
            "degradation_thresholds", DEFAULT_DEGRADATION_THRESHOLDS
        )

        self.use_historical = self.config.get("use_historical_priors", True)
        self.min_historical_runs = self.config.get(
            "min_historical_runs", DEFAULT_MIN_HISTORICAL_RUNS
        )
        self.per_stratum_overrides = self.config.get("per_stratum_overrides", {})

    def get_prior(
        self,
        stratum: str,
        baseline_value: float,
        historical_deltas: list[float] | None = None,
    ) -> tuple[float, float]:
        """Get Beta prior parameters for a stratum.

        Priority: per-stratum override > empirical > default.
        All priors are centered at the baseline value (i.e., "expect no change").

        Args:
            stratum: Stratum name
            baseline_value: Baseline metric value (prior center)
            historical_deltas: Historical delta values (converted to
                proportions internally as baseline + delta)

        Returns:
            (alpha_prior, beta_prior)
        """
        # Per-stratum concentration override
        if stratum in self.per_stratum_overrides:
            override = self.per_stratum_overrides[stratum]
            concentration = override.get("prior_concentration", self.default_concentration)
            return compute_beta_prior(baseline_value, concentration)

        # Empirical prior from historical data
        if (
            self.use_historical
            and historical_deltas
            and len(historical_deltas) >= self.min_historical_runs
        ):
            # Convert deltas to proportions
            proportions = np.array(
                [baseline_value + d for d in historical_deltas], dtype=np.float64
            )
            proportions = np.clip(proportions, 1e-6, 1.0 - 1e-6)

            mu_emp = float(proportions.mean())
            # Use sample variance (ddof=1, Bessel's correction) consistent
            # with qc_trending.py and spc.py.
            var_emp = float(proportions.var(ddof=1))

            if var_emp > 0:
                # Method of moments: kappa = mu(1-mu)/var - 1
                kappa = mu_emp * (1.0 - mu_emp) / var_emp - 1.0
                kappa = max(kappa, MIN_CONCENTRATION)
                alpha = mu_emp * kappa
                beta_param = (1.0 - mu_emp) * kappa
                return alpha, beta_param

        # Default: centered at baseline with default concentration
        return compute_beta_prior(baseline_value, self.default_concentration)

    def get_degradation_threshold(self, stratum: str, tier: str) -> float:
        """Get degradation threshold for a stratum and tier.

        Per-stratum override takes precedence over per-tier defaults.

        Args:
            stratum: Stratum name
            tier: Acceptance tier ("tier_1", "tier_2", "tier_3")

        Returns:
            Degradation threshold value
        """
        if stratum in self.per_stratum_overrides:
            override = self.per_stratum_overrides[stratum]
            if "degradation_threshold" in override:
                return override["degradation_threshold"]
        return self.degradation_thresholds.get(
            tier, self.degradation_thresholds.get("tier_2", 0.01)
        )

    def assess_stratum(
        self,
        dimension: str,
        stratum: str,
        metric_name: str,
        baseline_value: float,
        verification_value: float,
        sample_size: int,
        tier: str = "tier_3",
        historical_deltas: list[float] | None = None,
        successes: int | None = None,
        trials: int | None = None,
    ) -> dict:
        """Assess risk for a single metric in a single stratum.

        Uses a Beta-Binomial conjugate model: the verification proportion
        is modeled as Binomial(n, p) with a Beta prior centered at the
        baseline value.

        Args:
            dimension: Stratification dimension (e.g. "variant_class")
            stratum: Stratum value (e.g. "SNP")
            metric_name: Metric name ("sensitivity" or "precision")
            baseline_value: Baseline metric value
            verification_value: Verification metric value
            sample_size: Sample size (TP+FN for sensitivity, TP+FP for precision)
            tier: Acceptance tier
            historical_deltas: Historical deltas for empirical prior
            successes: Raw success count (e.g. TP). When provided with trials,
                used directly instead of rounding from proportion.
            trials: Raw trial count (e.g. TP+FN). When provided with successes,
                used directly instead of rounding from proportion.

        Returns:
            Dict with full posterior computation details
        """
        delta_observed = verification_value - baseline_value

        # Use raw counts when available; fall back to rounding from proportion
        if successes is not None and trials is not None:
            trials = max(0, trials)
            successes = max(0, min(trials, successes))
        else:
            trials = max(0, sample_size)
            successes = max(0, min(trials, round(verification_value * trials)))

        # Get Beta prior centered at baseline
        alpha_prior, beta_prior = self.get_prior(
            stratum, baseline_value, historical_deltas
        )

        # Conjugate update
        alpha_post, beta_post = compute_beta_posterior(
            alpha_prior, beta_prior, successes, trials
        )

        # Tier-aware degradation threshold
        threshold = self.get_degradation_threshold(stratum, tier)

        # Tail probability: P(p_true < baseline - threshold | data)
        tail_prob = compute_beta_tail(alpha_post, beta_post, baseline_value, threshold)

        flagged = tail_prob > self.alert_threshold

        # Compute posterior summary statistics for storage/reporting
        post_mean = alpha_post / (alpha_post + beta_post)
        post_std = beta_std(alpha_post, beta_post)
        prior_mean = alpha_prior / (alpha_prior + beta_prior)
        prior_std = beta_std(alpha_prior, beta_prior)

        return {
            "dimension": dimension,
            "stratum": stratum,
            "metric_name": metric_name,
            "delta_observed": delta_observed,
            "prior_mu": prior_mean,
            "prior_sigma": prior_std,
            "sigma_obs": None,
            "posterior_mu": post_mean,
            "posterior_sigma": post_std,
            "tail_probability": tail_prob,
            "risk_weight": 1.0,
            "weighted_risk": tail_prob,
            "flagged": flagged,
            "tier": tier,
            "degradation_threshold": threshold,
        }

    def assess_all(
        self,
        strata_data: list[dict],
        historical_data: dict[str, list[float]] | None = None,
    ) -> dict:
        """Assess risk across all strata.

        Args:
            strata_data: List of dicts, each with:
                - dimension, stratum, metric_name
                - baseline_value, verification_value
                - sample_size, tier
                - Optional: successes, trials (raw counts; avoids rounding)
            historical_data: Optional dict mapping stratum name to historical deltas

        Returns:
            Dict with:
            - verdict: BRAM_PASS or BRAM_FLAG
            - aggregate_risk_score: Global aggregate risk
            - mean_risk_score: Mean risk across strata
            - flagged_stratum_count: Number of flagged strata
            - per_stratum: List of per-stratum results
            - alert_threshold, degradation_threshold, aggregation_method
        """
        historical_data = historical_data or {}
        per_stratum_results = []

        for sd in strata_data:
            result = self.assess_stratum(
                dimension=sd["dimension"],
                stratum=sd["stratum"],
                metric_name=sd["metric_name"],
                baseline_value=sd["baseline_value"],
                verification_value=sd["verification_value"],
                sample_size=sd["sample_size"],
                tier=sd.get("tier", "tier_3"),
                historical_deltas=historical_data.get(sd["stratum"]),
                successes=sd.get("successes"),
                trials=sd.get("trials"),
            )
            per_stratum_results.append(result)

        tail_probs = [r["tail_probability"] for r in per_stratum_results]
        aggregate = aggregate_risks(tail_probs, self.aggregation_method)
        mean_risk = float(np.mean(tail_probs)) if tail_probs else 0.0
        flagged_count = sum(1 for r in per_stratum_results if r["flagged"])

        verdict = BRAM_FLAG if flagged_count > 0 else BRAM_PASS

        log.info(
            "BRAM assessment: %s (aggregate=%.3f, flagged=%d/%d)",
            verdict,
            aggregate,
            flagged_count,
            len(per_stratum_results),
        )

        return {
            "verdict": verdict,
            "aggregate_risk_score": aggregate,
            "mean_risk_score": mean_risk,
            "flagged_stratum_count": flagged_count,
            "per_stratum": per_stratum_results,
            "alert_threshold": self.alert_threshold,
            "degradation_threshold": self.degradation_thresholds.get("tier_2", 0.01),
            "aggregation_method": self.aggregation_method,
        }


def combine_verdicts(deterministic_verdict: str, bram_verdict: str) -> str:
    """Combine deterministic and BRAM verdicts.

    Two-stage decision table:
    | Deterministic   | BRAM       | Final              |
    |-----------------|------------|--------------------|
    | PASS            | BRAM_PASS  | PASS               |
    | PASS            | BRAM_FLAG  | REVIEW_REQUIRED    |
    | CONDITIONAL_PASS| BRAM_PASS  | CONDITIONAL_PASS   |
    | CONDITIONAL_PASS| BRAM_FLAG  | REVIEW_REQUIRED    |
    | REVIEW_REQUIRED | any        | REVIEW_REQUIRED    |
    | FAIL            | any        | FAIL               |

    Args:
        deterministic_verdict: Verdict from AcceptanceEngine
        bram_verdict: Verdict from BRAM (BRAM_PASS, BRAM_FLAG, BRAM_NOT_RUN)

    Returns:
        Final combined verdict
    """
    # BRAM not run — deterministic stands
    if bram_verdict == BRAM_NOT_RUN:
        return deterministic_verdict

    # Deterministic FAIL always stands
    if deterministic_verdict == FAIL:
        return FAIL

    # REVIEW_REQUIRED always stands
    if deterministic_verdict == REVIEW_REQUIRED:
        return REVIEW_REQUIRED

    # BRAM can escalate PASS or CONDITIONAL_PASS
    if bram_verdict == BRAM_FLAG:
        return REVIEW_REQUIRED

    # BRAM_PASS — deterministic stands
    return deterministic_verdict
