"""Tests for variant-level diff engine (Spec S13.2)."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.variant_diff import compute_variant_diff, assign_strata_to_transitions  # noqa: E402


def _make_variant(chrom, pos, ref, alt, classification):
    v = MagicMock()
    v.chrom = chrom
    v.pos = pos
    v.ref = ref
    v.alt = alt
    v.classification = classification
    return v


class TestComputeVariantDiff:
    def test_no_changes(self):
        base = [_make_variant("chr1", 100, "A", "T", "TP")]
        verify = [_make_variant("chr1", 100, "A", "T", "TP")]
        diff = compute_variant_diff(base, verify)
        assert len(diff) == 0

    def test_gained_variant(self):
        base = []
        verify = [_make_variant("chr1", 100, "A", "T", "TP")]
        diff = compute_variant_diff(base, verify)
        assert len(diff) == 1
        assert diff[0]["transition"] == "None→TP"
        assert diff[0]["baseline_classification"] is None
        assert diff[0]["verification_classification"] == "TP"

    def test_lost_variant(self):
        base = [_make_variant("chr1", 100, "A", "T", "TP")]
        verify = []
        diff = compute_variant_diff(base, verify)
        assert len(diff) == 1
        assert diff[0]["transition"] == "TP→None"

    def test_changed_classification(self):
        base = [_make_variant("chr1", 100, "A", "T", "TP")]
        verify = [_make_variant("chr1", 100, "A", "T", "FN")]
        diff = compute_variant_diff(base, verify)
        assert len(diff) == 1
        assert diff[0]["transition"] == "TP→FN"

    def test_multiple_transitions(self):
        base = [
            _make_variant("chr1", 100, "A", "T", "TP"),
            _make_variant("chr1", 200, "G", "C", "FP"),
            _make_variant("chr1", 300, "C", "A", "FN"),
        ]
        verify = [
            _make_variant("chr1", 100, "A", "T", "TP"),  # Same
            _make_variant("chr1", 200, "G", "C", "TP"),  # Changed
            _make_variant("chr1", 400, "T", "G", "TP"),  # New
        ]
        diff = compute_variant_diff(base, verify)
        # chr1:200 changed, chr1:300 lost, chr1:400 gained
        assert len(diff) == 3

    def test_dict_input(self):
        base = [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "classification": "TP"}]
        verify = [{"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "classification": "FP"}]
        diff = compute_variant_diff(base, verify)
        assert len(diff) == 1
        assert diff[0]["transition"] == "TP→FP"

    def test_sorted_output(self):
        base = [
            _make_variant("chr2", 100, "A", "T", "TP"),
            _make_variant("chr1", 100, "A", "T", "TP"),
        ]
        verify = [
            _make_variant("chr2", 100, "A", "T", "FN"),
            _make_variant("chr1", 100, "A", "T", "FN"),
        ]
        diff = compute_variant_diff(base, verify)
        assert diff[0]["chrom"] == "chr1"
        assert diff[1]["chrom"] == "chr2"

    def test_empty_both(self):
        diff = compute_variant_diff([], [])
        assert diff == []


class TestAssignStrataToTransitions:
    def test_assigns_strata(self):
        v = _make_variant("chr1", 100, "A", "T", "TP")
        strata = {
            ("variant_class", "SNP"): [v],
            ("functional_impact", "HIGH"): [v],
        }
        transitions = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"}
        ]

        result = assign_strata_to_transitions(transitions, [v], strata)
        assert len(result[0]["strata"]) == 2
        assert ("variant_class", "SNP") in result[0]["strata"]
        assert ("functional_impact", "HIGH") in result[0]["strata"]

    def test_no_strata(self):
        transitions = [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"}
        ]
        result = assign_strata_to_transitions(transitions, [], {})
        assert result[0]["strata"] == []
