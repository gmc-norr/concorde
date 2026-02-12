"""Tests for diff classifier (Spec S13.3)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.diff_classifier import classify_transitions, determine_verdict  # noqa: E402


def _make_transition(chrom="chr1", pos=100, ref="A", alt="T",
                     transition="TP→FN", strata=None):
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "transition": transition,
        "strata": strata or [],
    }


class TestClassifyTransitions:
    def test_drift_when_no_violations(self):
        transitions = [_make_transition()]
        result = classify_transitions(transitions, envelope_violations=[])
        assert result[0]["diff_class"] == "numerical_drift"

    def test_biological_when_envelope_breached(self):
        transitions = [
            _make_transition(strata=[("variant_class", "SNP")])
        ]
        violations = [
            {"dimension": "variant_class", "stratum": "SNP",
             "metric_name": "recall", "current_value": 0.95}
        ]
        result = classify_transitions(transitions, violations)
        assert result[0]["diff_class"] == "biological_difference"

    def test_biological_when_high_impact(self):
        transitions = [_make_transition()]
        high_impact = {("chr1", 100, "A", "T")}
        result = classify_transitions(
            transitions, [], high_impact_variants=high_impact
        )
        assert result[0]["diff_class"] == "biological_difference"

    def test_biological_when_regulated_panel(self):
        transitions = [_make_transition()]
        regulated = {("chr1", 100, "A", "T")}
        result = classify_transitions(
            transitions, [], regulated_panel_variants=regulated
        )
        assert result[0]["diff_class"] == "biological_difference"

    def test_mixed_classifications(self):
        transitions = [
            _make_transition(pos=100, strata=[("variant_class", "SNP")]),
            _make_transition(pos=200),
        ]
        violations = [
            {"dimension": "variant_class", "stratum": "SNP",
             "metric_name": "precision", "current_value": 0.90}
        ]
        result = classify_transitions(transitions, violations)
        assert result[0]["diff_class"] == "biological_difference"
        assert result[1]["diff_class"] == "numerical_drift"

    def test_empty_transitions(self):
        result = classify_transitions([], [])
        assert result == []


class TestDetermineVerdict:
    def test_pass_no_changes(self):
        assert determine_verdict([], []) == "PASS"

    def test_fail_biological_difference(self):
        transitions = [{"diff_class": "biological_difference"}]
        assert determine_verdict(transitions, []) == "FAIL"

    def test_fail_envelope_violations(self):
        violations = [{"dimension": "variant_class", "stratum": "SNP"}]
        assert determine_verdict([], violations) == "FAIL"

    def test_review_drift_only(self):
        transitions = [{"diff_class": "numerical_drift"}]
        assert determine_verdict(transitions, []) == "REVIEW_REQUIRED"

    def test_fail_both_bio_and_drift(self):
        transitions = [
            {"diff_class": "biological_difference"},
            {"diff_class": "numerical_drift"},
        ]
        assert determine_verdict(transitions, []) == "FAIL"
