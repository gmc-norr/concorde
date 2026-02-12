"""Tests for AcceptanceEngine + BRAM integration (S17.5)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.acceptance import AcceptanceEngine  # noqa: E402
from analysis.bayesian_risk import BRAM_FLAG, BRAM_NOT_RUN, BRAM_PASS  # noqa: E402


def _acceptance_config():
    return {
        "global": {"min_sensitivity": 0.99, "min_precision": 0.999},
        "per_stratum": {},
        "per_panel": {},
        "tiers": {
            "tier_1": {"panels": ["ACMG59"], "impacts": ["HIGH"]},
            "tier_2": {"impacts": ["MODERATE", "LOW"]},
            "tier_3": {"impacts": ["MODIFIER"]},
        },
        "low_count_threshold": 20,
        "low_count_policy": "warn",
    }


def _passing_metrics():
    return [
        {
            "dimension": "variant_class",
            "stratum": "SNP",
            "precision": 0.999,
            "recall": 0.998,
            "total_variants": 2000,
            "low_confidence": False,
        },
    ]


class TestEvaluateWithBram:
    """Tests for AcceptanceEngine.evaluate_with_bram()."""

    def test_pass_bram_pass(self):
        engine = AcceptanceEngine(_acceptance_config())
        bram_result = {"verdict": BRAM_PASS, "flagged_stratum_count": 0}
        result = engine.evaluate_with_bram(_passing_metrics(), bram_result)
        assert result["verdict"] == "PASS"
        assert result["deterministic_verdict"] == "PASS"
        assert result["bram_verdict"] == BRAM_PASS

    def test_pass_bram_flag_escalates(self):
        engine = AcceptanceEngine(_acceptance_config())
        bram_result = {"verdict": BRAM_FLAG, "flagged_stratum_count": 1}
        result = engine.evaluate_with_bram(_passing_metrics(), bram_result)
        assert result["verdict"] == "REVIEW_REQUIRED"
        assert result["deterministic_verdict"] == "PASS"
        assert result["bram_verdict"] == BRAM_FLAG

    def test_fail_bram_pass_stays_fail(self):
        engine = AcceptanceEngine(_acceptance_config())
        # Failing metric: HIGH impact with low sensitivity
        metrics = [
            {
                "dimension": "functional_impact",
                "stratum": "HIGH",
                "precision": 0.999,
                "recall": 0.950,  # Below 0.99
                "total_variants": 200,
                "low_confidence": False,
            },
        ]
        bram_result = {"verdict": BRAM_PASS, "flagged_stratum_count": 0}
        result = engine.evaluate_with_bram(metrics, bram_result)
        assert result["verdict"] == "FAIL"

    def test_fail_bram_flag_stays_fail(self):
        engine = AcceptanceEngine(_acceptance_config())
        metrics = [
            {
                "dimension": "functional_impact",
                "stratum": "HIGH",
                "precision": 0.999,
                "recall": 0.950,
                "total_variants": 200,
                "low_confidence": False,
            },
        ]
        bram_result = {"verdict": BRAM_FLAG, "flagged_stratum_count": 1}
        result = engine.evaluate_with_bram(metrics, bram_result)
        assert result["verdict"] == "FAIL"

    def test_no_bram_result(self):
        engine = AcceptanceEngine(_acceptance_config())
        result = engine.evaluate_with_bram(_passing_metrics(), None)
        assert result["verdict"] == "PASS"
        assert result["bram_verdict"] == BRAM_NOT_RUN

    def test_bram_result_included_in_output(self):
        engine = AcceptanceEngine(_acceptance_config())
        bram_result = {
            "verdict": BRAM_PASS,
            "aggregate_risk_score": 0.15,
            "flagged_stratum_count": 0,
        }
        result = engine.evaluate_with_bram(_passing_metrics(), bram_result)
        assert "bram_result" in result
        assert result["bram_result"]["aggregate_risk_score"] == 0.15

    def test_strata_results_preserved(self):
        engine = AcceptanceEngine(_acceptance_config())
        bram_result = {"verdict": BRAM_PASS}
        result = engine.evaluate_with_bram(_passing_metrics(), bram_result)
        assert len(result["strata_results"]) == 1
        assert result["strata_results"][0]["stratum"] == "SNP"
