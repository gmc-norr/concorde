"""Tests for BRAM reporting integration (S17.8)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from reporting.html_report import _status_color, generate_html_report  # noqa: E402
from reporting.json_report import generate_json_report  # noqa: E402


def _sample_bram_result():
    return {
        "verdict": "BRAM_PASS",
        "aggregate_risk_score": 0.142,
        "mean_risk_score": 0.142,
        "flagged_stratum_count": 0,
        "alert_threshold": 0.80,
        "degradation_threshold": 0.01,
        "aggregation_method": "max",
        "per_stratum": [
            {
                "dimension": "gene_panel",
                "stratum": "ACMG59",
                "metric_name": "sensitivity",
                "delta_observed": -0.02,
                "prior_mu": 0.999,
                "prior_sigma": 0.007,
                "sigma_obs": None,
                "posterior_mu": 0.982,
                "posterior_sigma": 0.004,
                "tail_probability": 0.142,
                "risk_weight": 1.0,
                "weighted_risk": 0.142,
                "flagged": False,
                "tier": "tier_1",
                "degradation_threshold": 0.005,
            },
        ],
    }


class TestJsonReportBRAM:
    """Tests for BRAM in JSON reports."""

    def test_report_without_bram(self):
        report = generate_json_report()
        assert "bayesian_risk" not in report

    def test_report_with_bram(self):
        report = generate_json_report(bayesian_risk=_sample_bram_result())
        assert "bayesian_risk" in report
        bram = report["bayesian_risk"]
        assert bram["enabled"] is True
        assert bram["verdict"] == "BRAM_PASS"
        assert bram["aggregate_risk_score"] == 0.142

    def test_bram_per_stratum_formatted(self):
        report = generate_json_report(bayesian_risk=_sample_bram_result())
        per_stratum = report["bayesian_risk"]["per_stratum"]
        assert len(per_stratum) == 1
        s = per_stratum[0]
        assert s["dimension"] == "gene_panel"
        assert s["stratum"] == "ACMG59"
        assert s["metric"] == "sensitivity"
        assert s["delta_observed"] == -0.02
        assert s["tail_probability"] == 0.142
        assert s["flagged"] is False

    def test_bram_flagged_report(self):
        bram = _sample_bram_result()
        bram["verdict"] = "BRAM_FLAG"
        bram["flagged_stratum_count"] = 1
        bram["per_stratum"][0]["flagged"] = True
        report = generate_json_report(bayesian_risk=bram)
        assert report["bayesian_risk"]["verdict"] == "BRAM_FLAG"
        assert report["bayesian_risk"]["flagged_count"] == 1

    def test_bram_thresholds_in_report(self):
        report = generate_json_report(bayesian_risk=_sample_bram_result())
        bram = report["bayesian_risk"]
        assert bram["alert_threshold"] == 0.80
        assert bram["degradation_threshold"] == 0.01


class TestHtmlReportBRAM:
    """Tests for BRAM in HTML reports."""

    def test_html_contains_bram_section(self):
        report = generate_json_report(bayesian_risk=_sample_bram_result())
        html = generate_html_report(report)
        assert "Bayesian Risk Assessment" in html
        assert "BRAM Verdict" in html

    def test_html_bram_shows_stratum(self):
        report = generate_json_report(bayesian_risk=_sample_bram_result())
        html = generate_html_report(report)
        assert "ACMG59" in html
        assert "sensitivity" in html

    def test_html_no_bram_section_without_data(self):
        report = generate_json_report()
        html = generate_html_report(report)
        assert "Bayesian Risk Assessment" not in html

    def test_html_bram_flagged_shows_red(self):
        bram = _sample_bram_result()
        bram["verdict"] = "BRAM_FLAG"
        bram["per_stratum"][0]["flagged"] = True
        report = generate_json_report(bayesian_risk=bram)
        html = generate_html_report(report)
        assert "FLAGGED" in html


class TestBRAMStatusColors:
    """Tests for BRAM verdict color mapping."""

    def test_bram_pass_green(self):
        assert _status_color("BRAM_PASS") == "green"

    def test_bram_flag_red(self):
        assert _status_color("BRAM_FLAG") == "red"

    def test_bram_not_run_gray(self):
        assert _status_color("BRAM_NOT_RUN") == "gray"
