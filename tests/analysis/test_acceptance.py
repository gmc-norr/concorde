"""Tests for acceptance criteria engine (Spec S15)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.acceptance import AcceptanceEngine  # noqa: E402


def _default_config():
    return {
        "global": {"min_sensitivity": 0.99, "min_precision": 0.999},
        "per_stratum": {
            "SNP": {"min_sensitivity": 0.995, "min_precision": 0.999},
            "INDEL": {"min_sensitivity": 0.98, "min_precision": 0.995},
        },
        "per_panel": {
            "ACMG59": {"min_sensitivity": 0.999, "min_precision": 0.9999},
        },
        "tiers": {
            "tier_1": {"panels": ["ACMG59"], "impacts": ["HIGH"]},
            "tier_2": {"impacts": ["MODERATE", "LOW"]},
        },
        "low_count_threshold": 20,
        "low_count_policy": "warn",
    }


def _make_metric(dimension="variant_class", stratum="SNP",
                 precision=0.9995, recall=0.998,
                 total_variants=100, low_confidence=False):
    return {
        "dimension": dimension,
        "stratum": stratum,
        "precision": precision,
        "recall": recall,
        "total_variants": total_variants,
        "low_confidence": low_confidence,
    }


class TestAssignTier:
    def test_tier1_acmg_panel(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("gene_panel", "ACMG59") == "tier_1"

    def test_tier1_high_impact(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("functional_impact", "HIGH") == "tier_1"

    def test_tier2_moderate_impact(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("functional_impact", "MODERATE") == "tier_2"

    def test_tier2_low_impact(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("functional_impact", "LOW") == "tier_2"

    def test_tier3_modifier_impact(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("functional_impact", "MODIFIER") == "tier_3"

    def test_tier3_default(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("variant_class", "SNP") == "tier_3"

    def test_tier3_off_panel(self):
        engine = AcceptanceEngine(_default_config())
        assert engine.assign_tier("gene_panel", "off_panel") == "tier_3"


class TestGetThresholds:
    def test_per_panel_override(self):
        engine = AcceptanceEngine(_default_config())
        t = engine.get_thresholds("gene_panel", "ACMG59")
        assert t["min_sensitivity"] == 0.999
        assert t["min_precision"] == 0.9999

    def test_per_stratum_override(self):
        engine = AcceptanceEngine(_default_config())
        t = engine.get_thresholds("variant_class", "SNP")
        assert t["min_sensitivity"] == 0.995

    def test_global_fallback(self):
        engine = AcceptanceEngine(_default_config())
        t = engine.get_thresholds("vaf_bins", "<0.05")
        assert t["min_sensitivity"] == 0.99
        assert t["min_precision"] == 0.999

    def test_indel_thresholds(self):
        engine = AcceptanceEngine(_default_config())
        t = engine.get_thresholds("variant_class", "INDEL")
        assert t["min_sensitivity"] == 0.98
        assert t["min_precision"] == 0.995


class TestEvaluateStratum:
    def test_pass(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(precision=0.9995, recall=0.998)
        )
        assert result["result"] == "PASS"
        assert result["violations"] == []

    def test_fail_sensitivity(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(precision=0.9995, recall=0.90)
        )
        assert result["result"] == "FAIL"
        assert any("sensitivity" in v for v in result["violations"])

    def test_fail_precision(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(precision=0.90, recall=0.998)
        )
        assert result["result"] == "FAIL"
        assert any("precision" in v for v in result["violations"])

    def test_fail_both(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(precision=0.90, recall=0.90)
        )
        assert result["result"] == "FAIL"
        assert len(result["violations"]) == 2

    def test_low_count_warn_policy(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(total_variants=5, low_confidence=True)
        )
        assert result["result"] == "WARN"

    def test_low_count_exclude_policy(self):
        config = _default_config()
        config["low_count_policy"] = "exclude"
        engine = AcceptanceEngine(config)
        result = engine.evaluate_stratum(
            _make_metric(total_variants=5, low_confidence=True)
        )
        assert result["result"] == "EXCLUDED"

    def test_low_count_fail_policy(self):
        config = _default_config()
        config["low_count_policy"] = "fail"
        engine = AcceptanceEngine(config)
        result = engine.evaluate_stratum(
            _make_metric(total_variants=5, low_confidence=True)
        )
        assert result["result"] == "FAIL"

    def test_none_metrics_pass(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(precision=None, recall=None)
        )
        assert result["result"] == "PASS"

    def test_tier_assigned(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_stratum(
            _make_metric(dimension="functional_impact", stratum="HIGH")
        )
        assert result["tier"] == "tier_1"


class TestEvaluateAll:
    def test_all_pass(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(stratum="SNP", precision=0.9995, recall=0.998),
            _make_metric(stratum="INDEL", precision=0.998, recall=0.99),
        ]
        result = engine.evaluate_all(metrics)
        assert result["verdict"] == "PASS"

    def test_tier1_fail_is_fail(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(
                dimension="functional_impact", stratum="HIGH",
                precision=0.80, recall=0.80,
            ),
        ]
        result = engine.evaluate_all(metrics)
        assert result["verdict"] == "FAIL"

    def test_tier2_fail_is_review(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(
                dimension="functional_impact", stratum="MODERATE",
                precision=0.80, recall=0.80,
            ),
        ]
        result = engine.evaluate_all(metrics)
        assert result["verdict"] == "REVIEW_REQUIRED"

    def test_tier2_warn_is_conditional(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(
                dimension="functional_impact", stratum="MODERATE",
                precision=0.9995, recall=0.998,
                total_variants=5, low_confidence=True,
            ),
        ]
        result = engine.evaluate_all(metrics)
        assert result["verdict"] == "CONDITIONAL_PASS"

    def test_tier3_fail_still_pass(self):
        """Tier 3 failures don't affect overall verdict."""
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(
                dimension="functional_impact", stratum="MODIFIER",
                precision=0.50, recall=0.50,
            ),
        ]
        result = engine.evaluate_all(metrics)
        assert result["verdict"] == "PASS"

    def test_mixed_tiers(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [
            _make_metric(
                dimension="functional_impact", stratum="HIGH",
                precision=0.9999, recall=0.9999,
            ),
            _make_metric(
                dimension="functional_impact", stratum="MODERATE",
                precision=0.80, recall=0.80,
            ),
            _make_metric(
                dimension="functional_impact", stratum="MODIFIER",
                precision=0.50, recall=0.50,
            ),
        ]
        result = engine.evaluate_all(metrics)
        # Tier 1 passes, tier 2 fails → REVIEW_REQUIRED
        assert result["verdict"] == "REVIEW_REQUIRED"

    def test_empty_metrics(self):
        engine = AcceptanceEngine(_default_config())
        result = engine.evaluate_all([])
        assert result["verdict"] == "PASS"

    def test_summary_included(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [_make_metric()]
        result = engine.evaluate_all(metrics)
        assert "summary" in result
        assert "PASS" in result["summary"]

    def test_strata_results_included(self):
        engine = AcceptanceEngine(_default_config())
        metrics = [_make_metric()]
        result = engine.evaluate_all(metrics)
        assert len(result["strata_results"]) == 1


class TestDefaultConfig:
    def test_no_config(self):
        engine = AcceptanceEngine()
        result = engine.evaluate_stratum(
            _make_metric(precision=0.9995, recall=0.995)
        )
        assert result["result"] == "PASS"

    def test_defaults_applied(self):
        engine = AcceptanceEngine()
        t = engine.get_thresholds("unknown", "unknown")
        assert t["min_sensitivity"] == 0.99
        assert t["min_precision"] == 0.999
