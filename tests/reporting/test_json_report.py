"""Tests for JSON report generation (Spec S16.1)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from reporting.json_report import (  # noqa: E402
    compute_report_checksum,
    generate_json_report,
    serialize_report,
)


class TestGenerateJsonReport:
    def test_minimal_report(self):
        report = generate_json_report()
        assert report["report_type"] == "initial_validation"
        assert report["report_version"] == "1.0"
        assert report["mode"] == "germline"
        assert "generated_at" in report
        assert "metrics" in report
        assert "acceptance" in report

    def test_report_with_metadata(self):
        report = generate_json_report(
            mode="somatic",
            run_metadata={"sample": "TUMOR", "caller": "Mutect2"},
        )
        assert report["mode"] == "somatic"
        assert report["run_metadata"]["sample"] == "TUMOR"

    def test_verification_report_includes_diff(self):
        diff = [{"chrom": "chr1", "pos": 100, "transition": "TP→FN"}]
        evidence = [{"category": "coverage", "score": "STRONG"}]
        report = generate_json_report(
            report_type="continuous_verification",
            variant_diff=diff,
            root_cause_evidence=evidence,
        )
        assert "variant_diff" in report
        assert len(report["variant_diff"]) == 1
        assert "root_cause_evidence" in report

    def test_initial_report_no_diff(self):
        report = generate_json_report(report_type="initial_validation")
        assert "variant_diff" not in report

    def test_stratified_metrics_formatted(self):
        metrics = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "precision": 0.99,
                "recall": 0.98,
                "f1": 0.985,
                "tp_count": 100,
                "fp_count": 1,
                "fn_count": 2,
                "total_variants": 103,
                "low_confidence": False,
            }
        ]
        report = generate_json_report(stratified_metrics=metrics)
        per_stratum = report["metrics"]["per_stratum"]
        assert len(per_stratum) == 1
        assert per_stratum[0]["dimension"] == "variant_class"

    def test_acceptance_formatted(self):
        acceptance = {
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
        }
        report = generate_json_report(acceptance_result=acceptance)
        assert report["acceptance"]["overall_result"] == "PASS"
        assert len(report["acceptance"]["per_stratum_results"]) == 1

    def test_acceptance_violations(self):
        acceptance = {
            "verdict": "FAIL",
            "strata_results": [
                {
                    "dimension": "functional_impact",
                    "stratum": "HIGH",
                    "tier": "tier_1",
                    "result": "FAIL",
                    "violations": ["sensitivity 0.9500 < 0.9990"],
                }
            ],
        }
        report = generate_json_report(acceptance_result=acceptance)
        assert len(report["acceptance"]["violations"]) == 1

    def test_no_acceptance(self):
        report = generate_json_report()
        assert report["acceptance"]["overall_result"] == "NOT_EVALUATED"

    def test_baseline_name_in_metadata(self):
        report = generate_json_report(baseline_name="v1.0_germline")
        assert report["run_metadata"]["baseline_name"] == "v1.0_germline"


class TestSerializeReport:
    def test_serialize_valid_json(self):
        report = generate_json_report()
        json_str = serialize_report(report)
        parsed = json.loads(json_str)
        assert parsed["report_version"] == "1.0"

    def test_serialize_pretty_printed(self):
        report = generate_json_report()
        json_str = serialize_report(report)
        assert "\n" in json_str  # Pretty-printed


class TestComputeChecksum:
    def test_deterministic(self):
        json_str = '{"test": "value"}'
        c1 = compute_report_checksum(json_str)
        c2 = compute_report_checksum(json_str)
        assert c1 == c2
        assert len(c1) == 64

    def test_different_content(self):
        c1 = compute_report_checksum('{"a": 1}')
        c2 = compute_report_checksum('{"b": 2}')
        assert c1 != c2
