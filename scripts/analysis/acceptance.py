"""Acceptance criteria engine (Spec S15).

Evaluates stratified metrics against mode-specific, risk-based
acceptance thresholds. Produces per-stratum and overall verdicts.
"""

from __future__ import annotations

import logging

from analysis.verdicts import (  # noqa: F401
    CONDITIONAL_PASS,
    EXCLUDED,
    FAIL,
    PASS,
    REVIEW_REQUIRED,
    TIER_1,
    TIER_2,
    TIER_3,
    WARN,
)

log = logging.getLogger(__name__)

# Default thresholds
DEFAULT_MIN_SENSITIVITY = 0.99
DEFAULT_MIN_PRECISION = 0.999
DEFAULT_LOW_COUNT_THRESHOLD = 20
DEFAULT_LOW_COUNT_POLICY = "warn"


class AcceptanceEngine:
    """Evaluates metrics against acceptance criteria.

    Args:
        config: Acceptance criteria configuration dict
    """

    def __init__(self, config: dict | None = None):
        self.config = config or {}
        self.global_thresholds = self.config.get("global", {})
        self.per_stratum = self.config.get("per_stratum", {})
        self.per_panel = self.config.get("per_panel", {})
        self.tiers_config = self.config.get("tiers", {})
        self.low_count_threshold = self.config.get(
            "low_count_threshold", DEFAULT_LOW_COUNT_THRESHOLD
        )
        self.low_count_policy = self.config.get(
            "low_count_policy", DEFAULT_LOW_COUNT_POLICY
        )

    def assign_tier(self, dimension: str, stratum: str) -> str:
        """Assign a risk tier to a stratum (S15.2).

        Args:
            dimension: Stratification dimension (e.g. "variant_class")
            stratum: Stratum value (e.g. "SNP", "HIGH", "ACMG59")

        Returns:
            Tier string: "tier_1", "tier_2", or "tier_3"
        """
        tier_1 = self.tiers_config.get("tier_1", {})
        tier_2 = self.tiers_config.get("tier_2", {})

        # Check tier 1: regulated panels and high impact
        if dimension == "gene_panel" and stratum in tier_1.get("panels", []):
            return TIER_1
        if dimension == "functional_impact" and stratum in tier_1.get("impacts", []):
            return TIER_1

        # Check tier 2: standard coding variants
        if dimension == "functional_impact" and stratum in tier_2.get("impacts", []):
            return TIER_2

        # Everything else is tier 3
        return TIER_3

    def get_thresholds(self, dimension: str, stratum: str) -> dict:
        """Get acceptance thresholds for a specific stratum.

        Checks per-panel, per-stratum, then falls back to global.

        Args:
            dimension: Stratification dimension
            stratum: Stratum value

        Returns:
            Dict with min_sensitivity and min_precision
        """
        # Per-panel overrides
        if dimension == "gene_panel" and stratum in self.per_panel:
            thresholds = self.per_panel[stratum].copy()
        # Per-stratum overrides
        elif stratum in self.per_stratum:
            thresholds = self.per_stratum[stratum].copy()
        else:
            thresholds = {}

        # Fill in defaults from global
        thresholds.setdefault(
            "min_sensitivity",
            self.global_thresholds.get("min_sensitivity", DEFAULT_MIN_SENSITIVITY),
        )
        thresholds.setdefault(
            "min_precision",
            self.global_thresholds.get("min_precision", DEFAULT_MIN_PRECISION),
        )
        return thresholds

    def evaluate_stratum(self, metric: dict) -> dict:
        """Evaluate a single stratum against its acceptance criteria.

        Args:
            metric: Dict with dimension, stratum, precision, recall,
                    total_variants, low_confidence

        Returns:
            Dict with dimension, stratum, tier, result (PASS/FAIL/WARN/EXCLUDED),
            details
        """
        dimension = metric.get("dimension", "")
        stratum = metric.get("stratum", "")
        total = metric.get("total_variants", 0)
        low_conf = metric.get("low_confidence", False)

        tier = self.assign_tier(dimension, stratum)
        thresholds = self.get_thresholds(dimension, stratum)

        result_entry = {
            "dimension": dimension,
            "stratum": stratum,
            "tier": tier,
            "thresholds": thresholds,
            "total_variants": total,
            "violations": [],
        }

        # Handle low count
        if low_conf or total < self.low_count_threshold:
            if self.low_count_policy == "exclude":
                result_entry["result"] = EXCLUDED
                return result_entry
            elif self.low_count_policy == "fail":
                result_entry["result"] = FAIL
                result_entry["violations"].append("low_variant_count")
                return result_entry
            # "warn" - continue evaluation but flag

        # Check sensitivity (recall)
        recall = metric.get("recall")
        min_sens = thresholds["min_sensitivity"]
        if recall is not None and recall < min_sens:
            result_entry["violations"].append(
                f"sensitivity {recall:.4f} < {min_sens:.4f}"
            )

        # Check precision
        precision = metric.get("precision")
        min_prec = thresholds["min_precision"]
        if precision is not None and precision < min_prec:
            result_entry["violations"].append(
                f"precision {precision:.4f} < {min_prec:.4f}"
            )

        if result_entry["violations"]:
            if low_conf and self.low_count_policy == "warn":
                result_entry["result"] = WARN
            else:
                result_entry["result"] = FAIL
        elif low_conf and self.low_count_policy == "warn":
            result_entry["result"] = WARN
            result_entry["violations"].append("low_variant_count")
        else:
            result_entry["result"] = PASS

        return result_entry

    def evaluate_all(self, metrics: list[dict]) -> dict:
        """Evaluate all strata and determine overall verdict (S15.3).

        Args:
            metrics: List of metric dicts from StratificationEngine.compute_metrics()

        Returns:
            Dict with:
            - verdict: PASS, CONDITIONAL_PASS, FAIL, or REVIEW_REQUIRED
            - strata_results: list of per-stratum results
            - summary: human-readable summary
        """
        strata_results = [self.evaluate_stratum(m) for m in metrics]
        verdict = self._determine_overall_verdict(strata_results)

        fail_count = sum(1 for r in strata_results if r["result"] == FAIL)
        warn_count = sum(1 for r in strata_results if r["result"] == WARN)
        pass_count = sum(1 for r in strata_results if r["result"] == PASS)

        summary = (
            f"{verdict}: {pass_count} passed, {fail_count} failed, "
            f"{warn_count} warnings out of {len(strata_results)} strata"
        )

        log.info("Acceptance evaluation: %s", summary)

        return {
            "verdict": verdict,
            "strata_results": strata_results,
            "summary": summary,
        }

    def _determine_overall_verdict(self, strata_results: list[dict]) -> str:
        """Determine overall verdict from per-stratum results (S15.3).

        - PASS: all Tier 1 and Tier 2 strata pass
        - CONDITIONAL_PASS: all Tier 1 pass; Tier 2 has warnings only
        - FAIL: any Tier 1 stratum fails
        - REVIEW_REQUIRED: Tier 2 failures or unresolved warnings
        """
        tier1_results = [r for r in strata_results if r["tier"] == TIER_1]
        tier2_results = [r for r in strata_results if r["tier"] == TIER_2]

        tier1_failed = any(r["result"] == FAIL for r in tier1_results)
        tier2_failed = any(r["result"] == FAIL for r in tier2_results)
        tier2_warned = any(r["result"] == WARN for r in tier2_results)

        if tier1_failed:
            return FAIL
        if tier2_failed:
            return REVIEW_REQUIRED
        if tier2_warned:
            return CONDITIONAL_PASS
        return PASS

    def evaluate_with_bram(
        self,
        metrics: list[dict],
        bram_result: dict | None = None,
    ) -> dict:
        """Evaluate acceptance with optional BRAM augmentation (S17.5).

        Runs deterministic evaluation first, then combines with BRAM
        verdict using the two-stage decision table.

        Args:
            metrics: List of metric dicts from StratificationEngine
            bram_result: Optional BRAM assessment result dict

        Returns:
            Dict with combined verdict, deterministic and BRAM sub-results
        """
        from analysis.bayesian_risk import BRAM_NOT_RUN, combine_verdicts

        deterministic = self.evaluate_all(metrics)
        det_verdict = deterministic["verdict"]

        bram_verdict = BRAM_NOT_RUN
        if bram_result and "verdict" in bram_result:
            bram_verdict = bram_result["verdict"]

        final_verdict = combine_verdicts(det_verdict, bram_verdict)

        result = {
            "verdict": final_verdict,
            "deterministic_verdict": det_verdict,
            "bram_verdict": bram_verdict,
            "strata_results": deterministic["strata_results"],
            "summary": deterministic["summary"],
        }

        if bram_result:
            result["bram_result"] = bram_result

        if final_verdict != det_verdict:
            log.info(
                "BRAM escalated verdict: %s -> %s", det_verdict, final_verdict
            )

        return result
