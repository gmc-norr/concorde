"""Tests for HTML report generation (Spec S16.2)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from reporting.html_report import (  # noqa: E402
    _format_float,
    _format_pct,
    _status_color,
    generate_html_report,
)
from reporting.json_report import generate_json_report  # noqa: E402


def _sample_report():
    return generate_json_report(
        mode="germline",
        run_metadata={"sample": "NA12878", "pipeline_version": "v4.5.0"},
        global_metrics={"sensitivity": 0.998, "precision": 0.999, "f1": 0.9985},
        stratified_metrics=[
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "precision": 0.999,
                "recall": 0.998,
                "f1": 0.9985,
                "tp_count": 2950,
                "fp_count": 3,
                "fn_count": 6,
                "total_variants": 2959,
                "low_confidence": False,
            }
        ],
        acceptance_result={
            "verdict": "PASS",
            "strata_results": [
                {
                    "dimension": "variant_class",
                    "stratum": "SNP",
                    "tier": "tier_3",
                    "result": "PASS",
                    "violations": [],
                }
            ],
        },
    )


class TestGenerateHtmlReport:
    def test_generates_html(self):
        report = _sample_report()
        html = generate_html_report(report)
        assert "<!DOCTYPE html>" in html
        assert "Concorde Validation Report" in html

    def test_contains_metadata(self):
        report = _sample_report()
        html = generate_html_report(report)
        assert "NA12878" in html
        assert "v4.5.0" in html

    def test_contains_metrics(self):
        report = _sample_report()
        html = generate_html_report(report)
        assert "SNP" in html
        assert "2950" in html  # tp_count

    def test_contains_acceptance_verdict(self):
        report = _sample_report()
        html = generate_html_report(report)
        assert "PASS" in html

    def test_fail_report_red_banner(self):
        report = generate_json_report(
            acceptance_result={
                "verdict": "FAIL",
                "strata_results": [
                    {
                        "dimension": "functional_impact",
                        "stratum": "HIGH",
                        "tier": "tier_1",
                        "result": "FAIL",
                        "violations": ["sensitivity below threshold"],
                    }
                ],
            },
        )
        html = generate_html_report(report)
        assert "red" in html

    def test_verification_report_with_diff(self):
        report = generate_json_report(
            report_type="continuous_verification",
            variant_diff=[
                {
                    "chrom": "chr1",
                    "pos": 100,
                    "ref": "A",
                    "alt": "T",
                    "baseline_classification": "TP",
                    "verification_classification": "FN",
                    "transition": "TP→FN",
                }
            ],
        )
        html = generate_html_report(report)
        assert "Variant Differences" in html
        assert "TP→FN" in html

    def test_empty_report(self):
        report = generate_json_report()
        html = generate_html_report(report)
        assert "<!DOCTYPE html>" in html

    def test_custom_template_dir(self, tmp_path):
        tpl = tmp_path / "custom.html.j2"
        tpl.write_text("<html>{{ report.mode }}</html>")

        report = generate_json_report(mode="somatic")
        html = generate_html_report(report, "custom.html.j2", str(tmp_path))
        assert "somatic" in html


class TestFilters:
    def test_format_pct(self):
        assert _format_pct(0.998) == "99.80%"
        assert _format_pct(1.0) == "100.00%"
        assert _format_pct(None) == "N/A"

    def test_format_float(self):
        assert _format_float(0.9985) == "0.9985"
        assert _format_float(0.9985, 2) == "1.00"
        assert _format_float(None) == "N/A"

    def test_status_color(self):
        assert _status_color("PASS") == "green"
        assert _status_color("FAIL") == "red"
        assert _status_color("REVIEW_REQUIRED") == "amber"
        assert _status_color("UNKNOWN") == "gray"
