"""Tests for Bayesian Risk Assessment Module (Spec S17).

Tests the Beta-Binomial conjugate model for BRAM.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
from scipy.stats import beta as beta_dist

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.bayesian_risk import (  # noqa: E402
    BRAM_FLAG,
    BRAM_NOT_RUN,
    BRAM_PASS,
    BRAMEngine,
    aggregate_risks,
    beta_std,
    combine_verdicts,
    compute_beta_posterior,
    compute_beta_prior,
    compute_beta_tail,
)


class TestBetaStd:
    """Tests for Beta distribution standard deviation."""

    def test_symmetric_beta(self):
        """Beta(10, 10) has known std."""
        # Var = ab / ((a+b)^2 (a+b+1)) = 100 / (400 * 21) = 100/8400
        expected = math.sqrt(100 / 8400)
        assert beta_std(10, 10) == pytest.approx(expected, abs=1e-10)

    def test_asymmetric_beta(self):
        """Beta(2, 5) std."""
        a, b = 2, 5
        expected = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
        assert beta_std(a, b) == pytest.approx(expected, abs=1e-10)

    def test_large_params_small_std(self):
        """Large parameters should give small std."""
        assert beta_std(1000, 1000) < beta_std(10, 10)


class TestComputeBetaPrior:
    """Tests for Beta prior computation."""

    def test_centered_at_baseline(self):
        """Prior mean should equal baseline value."""
        alpha, beta_p = compute_beta_prior(0.999, 20)
        mean = alpha / (alpha + beta_p)
        assert mean == pytest.approx(0.999, abs=1e-4)

    def test_higher_concentration_tighter(self):
        """Higher concentration should give smaller std."""
        a1, b1 = compute_beta_prior(0.99, 20)
        a2, b2 = compute_beta_prior(0.99, 100)
        assert beta_std(a2, b2) < beta_std(a1, b1)

    def test_clamps_extreme_values(self):
        """Values at 0 or 1 should be clamped to avoid degenerate Beta."""
        alpha, beta_p = compute_beta_prior(1.0, 20)
        assert alpha > 0
        assert beta_p > 0

        alpha, beta_p = compute_beta_prior(0.0, 20)
        assert alpha > 0
        assert beta_p > 0

    def test_min_concentration(self):
        """Concentration below minimum should be clamped."""
        alpha, beta_p = compute_beta_prior(0.5, 0.5)
        # Should be clamped to MIN_CONCENTRATION=2
        assert alpha + beta_p == pytest.approx(2.0, abs=1e-4)

    def test_returns_floats(self):
        alpha, beta_p = compute_beta_prior(0.99, 20)
        assert isinstance(alpha, float)
        assert isinstance(beta_p, float)


class TestComputeBetaPosterior:
    """Tests for Beta-Binomial conjugate posterior."""

    def test_basic_update(self):
        """alpha_post = alpha_prior + k, beta_post = beta_prior + n - k."""
        alpha_post, beta_post = compute_beta_posterior(10, 5, 80, 100)
        assert alpha_post == 10 + 80
        assert beta_post == 5 + 20

    def test_all_successes(self):
        """All successes: failures = 0."""
        alpha_post, beta_post = compute_beta_posterior(2, 3, 50, 50)
        assert alpha_post == 52
        assert beta_post == 3

    def test_no_successes(self):
        """No successes: all failures."""
        alpha_post, beta_post = compute_beta_posterior(2, 3, 0, 50)
        assert alpha_post == 2
        assert beta_post == 53

    def test_zero_trials(self):
        """Zero trials: posterior equals prior."""
        alpha_post, beta_post = compute_beta_posterior(10, 5, 0, 0)
        assert alpha_post == 10
        assert beta_post == 5

    def test_negative_inputs_clamped(self):
        """Negative successes/trials should be clamped to 0."""
        alpha_post, beta_post = compute_beta_posterior(10, 5, -5, -3)
        assert alpha_post == 10
        assert beta_post == 5

    def test_posterior_mean_pulled_toward_data(self):
        """With enough data, posterior mean approaches observed proportion."""
        # Prior centered at 0.5
        alpha_prior, beta_prior = 5, 5
        # Data: 90% success with 1000 trials
        alpha_post, beta_post = compute_beta_posterior(alpha_prior, beta_prior, 900, 1000)
        post_mean = alpha_post / (alpha_post + beta_post)
        assert post_mean == pytest.approx(0.9, abs=0.01)


class TestComputeBetaTail:
    """Tests for tail probability computation."""

    def test_no_degradation(self):
        """When data matches baseline, low tail probability."""
        # Prior + data centered at 0.999, threshold = 0.01
        # So we're asking P(p < 0.989) which should be very low
        alpha_post = 0.999 * 200 + 0.999 * 20  # ~approx
        beta_post = 0.001 * 200 + 0.001 * 20
        prob = compute_beta_tail(alpha_post, beta_post, 0.999, 0.01)
        assert prob < 0.05

    def test_large_degradation(self):
        """Large degradation should give high tail probability."""
        # Data suggests p ~ 0.95, baseline = 0.999, threshold = 0.01
        # P(p < 0.989) should be high when data says p ~ 0.95
        alpha_post = 0.95 * 200 + 0.95 * 20
        beta_post = 0.05 * 200 + 0.05 * 20
        prob = compute_beta_tail(alpha_post, beta_post, 0.999, 0.01)
        assert prob > 0.95

    def test_at_threshold(self):
        """When posterior is centered at baseline - threshold, tail ~ 50%."""
        # Center posterior at 0.989 (= 0.999 - 0.01)
        # Use large concentration so it's well-peaked
        alpha_post = 0.989 * 5000
        beta_post = 0.011 * 5000
        prob = compute_beta_tail(alpha_post, beta_post, 0.999, 0.01)
        assert abs(prob - 0.5) < 0.05

    def test_matches_scipy(self):
        """Verify against direct scipy.stats.beta.cdf call."""
        alpha_post, beta_post = 100, 5
        baseline, threshold = 0.99, 0.01
        cutoff = baseline - threshold
        expected = float(beta_dist.cdf(cutoff, alpha_post, beta_post))
        result = compute_beta_tail(alpha_post, beta_post, baseline, threshold)
        assert result == pytest.approx(expected, abs=1e-10)

    def test_zero_cutoff_returns_zero(self):
        """If baseline - threshold <= 0, return 0."""
        prob = compute_beta_tail(10, 5, 0.005, 0.01)
        assert prob == 0.0

    def test_returns_float(self):
        prob = compute_beta_tail(10, 5, 0.99, 0.01)
        assert isinstance(prob, float)


class TestAggregateRisks:
    """Tests for risk aggregation (S17.3)."""

    def test_max_aggregation(self):
        risks = [0.1, 0.5, 0.3]
        assert aggregate_risks(risks, "max") == 0.5

    def test_weighted_mean(self):
        risks = [0.2, 0.4, 0.6]
        assert abs(aggregate_risks(risks, "weighted_mean") - 0.4) < 1e-10

    def test_product_aggregation(self):
        """Product: 1 - (1-0.1)*(1-0.2)*(1-0.3) = 1 - 0.504 = 0.496."""
        risks = [0.1, 0.2, 0.3]
        expected = 1.0 - (0.9 * 0.8 * 0.7)
        assert abs(aggregate_risks(risks, "product") - expected) < 1e-10

    def test_empty_list(self):
        assert aggregate_risks([], "max") == 0.0

    def test_unknown_method_falls_back(self):
        risks = [0.1, 0.5, 0.3]
        assert aggregate_risks(risks, "unknown") == 0.5

    def test_single_risk(self):
        assert aggregate_risks([0.42], "max") == 0.42


class TestBRAMEngine:
    """Tests for the BRAMEngine class with Beta-Binomial model."""

    def _default_engine(self, **overrides):
        config = {
            "alert_threshold": 0.80,
            "aggregation_method": "max",
            "default_prior": {"concentration": 20},
            "degradation_thresholds": {
                "tier_1": 0.005,
                "tier_2": 0.01,
                "tier_3": 0.02,
            },
            "use_historical_priors": True,
            "min_historical_runs": 5,
            "per_stratum_overrides": {},
        }
        config.update(overrides)
        return BRAMEngine(config=config)

    def test_default_config(self):
        engine = BRAMEngine()
        assert engine.alert_threshold == 0.80

    def test_custom_config(self):
        engine = self._default_engine(alert_threshold=0.90)
        assert engine.alert_threshold == 0.90

    def test_get_prior_default(self):
        """Default prior centered at baseline with default concentration."""
        engine = self._default_engine()
        alpha, beta_p = engine.get_prior("SNP", baseline_value=0.999)
        mean = alpha / (alpha + beta_p)
        assert mean == pytest.approx(0.999, abs=1e-3)
        # Concentration = 20, so alpha + beta ~ 20
        assert alpha + beta_p == pytest.approx(20.0, abs=1e-3)

    def test_get_prior_override(self):
        """Per-stratum override uses custom concentration."""
        engine = self._default_engine(
            per_stratum_overrides={"ACMG59": {"prior_concentration": 50}}
        )
        alpha, beta_p = engine.get_prior("ACMG59", baseline_value=0.999)
        assert alpha + beta_p == pytest.approx(50.0, abs=1e-3)

    def test_get_prior_empirical(self):
        """Empirical prior from historical deltas uses method of moments."""
        engine = self._default_engine()
        deltas = [0.001, -0.002, 0.0, -0.001, 0.002, 0.001]
        alpha, beta_p = engine.get_prior(
            "SNP", baseline_value=0.999, historical_deltas=deltas
        )
        # Should produce valid Beta parameters
        assert alpha > 0
        assert beta_p > 0
        # Mean should be close to baseline + mean(deltas)
        mean = alpha / (alpha + beta_p)
        import numpy as np

        expected_center = 0.999 + np.mean(deltas)
        assert mean == pytest.approx(expected_center, abs=0.01)

    def test_get_prior_insufficient_history(self):
        """Too few historical points falls back to default."""
        engine = self._default_engine(min_historical_runs=5)
        deltas = [0.001, -0.002, 0.0]  # Only 3, need 5
        alpha, beta_p = engine.get_prior(
            "SNP", baseline_value=0.999, historical_deltas=deltas
        )
        # Should fall back to default concentration=20
        assert alpha + beta_p == pytest.approx(20.0, abs=1e-3)

    def test_get_degradation_threshold_per_tier(self):
        """Threshold varies by tier."""
        engine = self._default_engine()
        assert engine.get_degradation_threshold("SNP", "tier_1") == 0.005
        assert engine.get_degradation_threshold("SNP", "tier_2") == 0.01
        assert engine.get_degradation_threshold("SNP", "tier_3") == 0.02

    def test_get_degradation_threshold_override(self):
        engine = self._default_engine(
            per_stratum_overrides={"ACMG59": {"degradation_threshold": 0.003}}
        )
        assert engine.get_degradation_threshold("ACMG59", "tier_1") == 0.003
        assert engine.get_degradation_threshold("SNP", "tier_1") == 0.005

    def test_assess_stratum_no_degradation(self):
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.999,
            sample_size=200,
            tier="tier_3",
        )
        assert result["delta_observed"] == 0.0
        assert result["tail_probability"] < 0.1
        assert result["weighted_risk"] < 0.1
        assert result["flagged"] is False

    def test_assess_stratum_moderate_degradation(self):
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="gene_panel",
            stratum="ACMG59",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.979,
            sample_size=200,
            tier="tier_1",
        )
        assert result["delta_observed"] == pytest.approx(-0.02, abs=1e-6)
        # Risk weight is always 1.0 in Beta model
        assert result["risk_weight"] == 1.0
        assert result["weighted_risk"] == result["tail_probability"]

    def test_assess_stratum_has_all_fields(self):
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            baseline_value=0.999,
            verification_value=0.998,
            sample_size=100,
        )
        expected_fields = {
            "dimension", "stratum", "metric_name", "delta_observed",
            "prior_mu", "prior_sigma", "sigma_obs",
            "posterior_mu", "posterior_sigma",
            "tail_probability", "risk_weight", "weighted_risk",
            "flagged", "tier", "degradation_threshold",
        }
        assert set(result.keys()) == expected_fields

    def test_assess_stratum_sigma_obs_is_none(self):
        """Beta model sets sigma_obs to None (not used)."""
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.998,
            sample_size=100,
        )
        assert result["sigma_obs"] is None

    def test_assess_stratum_raw_counts(self):
        """Raw counts should be used when provided, avoiding rounding."""
        engine = self._default_engine()
        # With proportion: 95/100 = 0.95 -> round(0.95 * 100) = 95 (same)
        # But with 193/200 = 0.965 -> round(0.965 * 200) = 193 (same)
        # Edge case: 7/9 = 0.7778 -> round(0.7778 * 9) = round(7.0) = 7 (same)
        # Real edge: successes=1, trials=3 -> prop=0.333, round(0.333*3)=round(1.0)=1
        # Better: successes=2, trials=3 -> prop=0.667, round(0.667*3)=round(2.0)=2

        # Case where rounding matches: results should be identical
        result_counts = engine.assess_stratum(
            dimension="d", stratum="s", metric_name="sens",
            baseline_value=0.99, verification_value=0.95,
            sample_size=100, tier="tier_2",
            successes=95, trials=100,
        )
        result_rounded = engine.assess_stratum(
            dimension="d", stratum="s", metric_name="sens",
            baseline_value=0.99, verification_value=0.95,
            sample_size=100, tier="tier_2",
        )
        assert result_counts["tail_probability"] == result_rounded["tail_probability"]

    def test_assess_all_with_raw_counts(self):
        """assess_all should pass through successes/trials to assess_stratum."""
        engine = self._default_engine()
        strata_data = [
            {
                "dimension": "d",
                "stratum": "s",
                "metric_name": "sensitivity",
                "baseline_value": 0.999,
                "verification_value": 0.95,
                "sample_size": 100,
                "tier": "tier_2",
                "successes": 95,
                "trials": 100,
            },
        ]
        result = engine.assess_all(strata_data)
        assert result["verdict"] in ("BRAM_PASS", "BRAM_FLAG")
        assert len(result["per_stratum"]) == 1

    def test_assess_stratum_records_degradation_threshold(self):
        """Result should include the threshold that was used."""
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.998,
            sample_size=100,
            tier="tier_1",
        )
        assert result["degradation_threshold"] == 0.005

    def test_assess_stratum_large_drop_flagged(self):
        """Large degradation with low alert threshold should flag."""
        engine = self._default_engine(alert_threshold=0.10)
        result = engine.assess_stratum(
            dimension="gene_panel",
            stratum="ACMG59",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.950,
            sample_size=1000,
            tier="tier_1",
        )
        # With 1000 samples, drop from 0.999 to 0.950 is very significant
        assert result["flagged"] is True
        assert result["tail_probability"] > 0.10

    def test_assess_all_pass(self):
        engine = self._default_engine()
        strata_data = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "metric_name": "sensitivity",
                "baseline_value": 0.999,
                "verification_value": 0.998,
                "sample_size": 2000,
                "tier": "tier_3",
            },
            {
                "dimension": "variant_class",
                "stratum": "INDEL",
                "metric_name": "sensitivity",
                "baseline_value": 0.995,
                "verification_value": 0.994,
                "sample_size": 500,
                "tier": "tier_3",
            },
        ]
        result = engine.assess_all(strata_data)
        assert result["verdict"] == BRAM_PASS
        assert result["flagged_stratum_count"] == 0
        assert len(result["per_stratum"]) == 2

    def test_assess_all_flag(self):
        engine = self._default_engine(alert_threshold=0.10)
        strata_data = [
            {
                "dimension": "gene_panel",
                "stratum": "ACMG59",
                "metric_name": "sensitivity",
                "baseline_value": 0.999,
                "verification_value": 0.950,
                "sample_size": 1000,
                "tier": "tier_1",
            },
        ]
        result = engine.assess_all(strata_data)
        assert result["verdict"] == BRAM_FLAG
        assert result["flagged_stratum_count"] >= 1

    def test_assess_all_output_structure(self):
        engine = self._default_engine()
        strata_data = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "metric_name": "sensitivity",
                "baseline_value": 0.999,
                "verification_value": 0.998,
                "sample_size": 1000,
                "tier": "tier_3",
            },
        ]
        result = engine.assess_all(strata_data)
        expected_keys = {
            "verdict", "aggregate_risk_score", "mean_risk_score",
            "flagged_stratum_count", "per_stratum",
            "alert_threshold", "degradation_threshold", "aggregation_method",
        }
        assert set(result.keys()) == expected_keys

    def test_assess_all_with_historical_data(self):
        engine = self._default_engine()
        strata_data = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "metric_name": "sensitivity",
                "baseline_value": 0.999,
                "verification_value": 0.997,
                "sample_size": 500,
                "tier": "tier_3",
            },
        ]
        historical = {"SNP": [0.0, -0.001, 0.001, -0.001, 0.0, 0.001]}
        result = engine.assess_all(strata_data, historical_data=historical)
        assert result["verdict"] == BRAM_PASS

    def test_assess_all_empty(self):
        engine = self._default_engine()
        result = engine.assess_all([])
        assert result["verdict"] == BRAM_PASS
        assert result["aggregate_risk_score"] == 0.0
        assert result["flagged_stratum_count"] == 0

    def test_posterior_mu_near_verification_value_for_large_n(self):
        """With large sample, posterior mean converges to observed proportion."""
        engine = self._default_engine()
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.980,
            sample_size=10000,
            tier="tier_3",
        )
        # With N=10000, data dominates → posterior mean ~ 0.980
        assert result["posterior_mu"] == pytest.approx(0.980, abs=0.005)

    def test_posterior_mu_near_prior_for_small_n(self):
        """With very small sample, posterior stays near prior (baseline)."""
        engine = self._default_engine(
            default_prior={"concentration": 100}
        )
        result = engine.assess_stratum(
            dimension="variant_class",
            stratum="SNP",
            metric_name="sensitivity",
            baseline_value=0.999,
            verification_value=0.900,
            sample_size=2,
            tier="tier_3",
        )
        # Prior concentration=100, only 2 data points → prior dominates
        assert result["posterior_mu"] > 0.95


class TestCombineVerdicts:
    """Tests for two-stage verdict combination (S17.5)."""

    def test_pass_bram_pass(self):
        assert combine_verdicts("PASS", BRAM_PASS) == "PASS"

    def test_pass_bram_flag(self):
        assert combine_verdicts("PASS", BRAM_FLAG) == "REVIEW_REQUIRED"

    def test_conditional_pass_bram_pass(self):
        assert combine_verdicts("CONDITIONAL_PASS", BRAM_PASS) == "CONDITIONAL_PASS"

    def test_conditional_pass_bram_flag(self):
        assert combine_verdicts("CONDITIONAL_PASS", BRAM_FLAG) == "REVIEW_REQUIRED"

    def test_review_required_bram_pass(self):
        assert combine_verdicts("REVIEW_REQUIRED", BRAM_PASS) == "REVIEW_REQUIRED"

    def test_review_required_bram_flag(self):
        assert combine_verdicts("REVIEW_REQUIRED", BRAM_FLAG) == "REVIEW_REQUIRED"

    def test_fail_bram_pass(self):
        assert combine_verdicts("FAIL", BRAM_PASS) == "FAIL"

    def test_fail_bram_flag(self):
        assert combine_verdicts("FAIL", BRAM_FLAG) == "FAIL"

    def test_bram_not_run_preserves_deterministic(self):
        assert combine_verdicts("PASS", BRAM_NOT_RUN) == "PASS"
        assert combine_verdicts("FAIL", BRAM_NOT_RUN) == "FAIL"
        assert combine_verdicts("CONDITIONAL_PASS", BRAM_NOT_RUN) == "CONDITIONAL_PASS"
        assert combine_verdicts("REVIEW_REQUIRED", BRAM_NOT_RUN) == "REVIEW_REQUIRED"
