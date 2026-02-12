"""Tests for root cause analysis (Spec S14)."""

from __future__ import annotations

import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.root_cause import (  # noqa: E402
    collect_evidence,
    determine_likely_cause,
    requires_review,
)


class TestCoverageEvidence:
    def test_strong_coverage_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"dp": 30},
            verification_variant={"dp": 5},
        )
        cov = [e for e in evidence if e["category"] == "coverage"][0]
        assert cov["score"] == "STRONG"

    def test_moderate_coverage_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"dp": 30},
            verification_variant={"dp": 20},
        )
        cov = [e for e in evidence if e["category"] == "coverage"][0]
        assert cov["score"] == "MODERATE"

    def test_weak_coverage_change(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"dp": 30},
            verification_variant={"dp": 28},
        )
        cov = [e for e in evidence if e["category"] == "coverage"][0]
        assert cov["score"] == "WEAK"

    def test_no_coverage_change(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"dp": 30},
            verification_variant={"dp": 30},
        )
        cov = [e for e in evidence if e["category"] == "coverage"][0]
        assert cov["score"] == "NONE"

    def test_bam_not_available(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            bam_available=False,
        )
        cov = [e for e in evidence if e["category"] == "coverage"][0]
        assert cov["score"] == "NONE"
        assert "BAM not available" in cov["change_description"]


class TestMappingQualityEvidence:
    def test_strong_mq_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"mq": 60.0},
            verification_variant={"mq": 40.0},
        )
        mq = [e for e in evidence if e["category"] == "mapping_quality"][0]
        assert mq["score"] == "STRONG"

    def test_moderate_mq_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"mq": 60.0},
            verification_variant={"mq": 53.0},
        )
        mq = [e for e in evidence if e["category"] == "mapping_quality"][0]
        assert mq["score"] == "MODERATE"

    def test_no_mq_change(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"mq": 60.0},
            verification_variant={"mq": 60.0},
        )
        mq = [e for e in evidence if e["category"] == "mapping_quality"][0]
        assert mq["score"] == "NONE"


class TestFilterEvidence:
    def test_pass_to_lowqual(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"filter_status": "PASS"},
            verification_variant={"filter_status": "LowQual"},
        )
        filt = [e for e in evidence if e["category"] == "filter"][0]
        assert filt["score"] == "STRONG"

    def test_lowqual_to_pass(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "FN→TP"},
            baseline_variant={"filter_status": "LowDP"},
            verification_variant={"filter_status": "PASS"},
        )
        filt = [e for e in evidence if e["category"] == "filter"][0]
        assert filt["score"] == "STRONG"

    def test_no_filter_change(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"filter_status": "PASS"},
            verification_variant={"filter_status": "PASS"},
        )
        filt = [e for e in evidence if e["category"] == "filter"][0]
        assert filt["score"] == "NONE"

    def test_non_pass_filter_change(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"filter_status": "LowDP"},
            verification_variant={"filter_status": "LowQual"},
        )
        filt = [e for e in evidence if e["category"] == "filter"][0]
        assert filt["score"] == "MODERATE"


class TestQualityScoreEvidence:
    def test_strong_qual_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"qual": 100.0},
            verification_variant={"qual": 30.0},
        )
        qual = [e for e in evidence if e["category"] == "quality_score"][0]
        assert qual["score"] == "STRONG"

    def test_moderate_qual_drop(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            baseline_variant={"qual": 100.0},
            verification_variant={"qual": 70.0},
        )
        qual = [e for e in evidence if e["category"] == "quality_score"][0]
        assert qual["score"] == "MODERATE"

    def test_no_qual_data(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
        )
        qual = [e for e in evidence if e["category"] == "quality_score"][0]
        assert qual["score"] == "NONE"


class TestToolVersionEvidence:
    def test_tool_changes(self):
        changes = [
            {"tool": "GATK", "old_version": "4.4.0", "new_version": "4.5.0"},
        ]
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            tool_changes=changes,
        )
        tool = [e for e in evidence if e["category"] == "tool_version"][0]
        assert tool["score"] == "WEAK"
        assert "GATK" in tool["change_description"]

    def test_no_tool_changes(self):
        evidence = collect_evidence(
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "T", "transition": "TP→FN"},
            tool_changes=[],
        )
        tool = [e for e in evidence if e["category"] == "tool_version"][0]
        assert tool["score"] == "NONE"


class TestDetermineLikelyCause:
    def test_coverage_drop(self):
        evidence = [
            {"category": "coverage", "score": "STRONG"},
            {"category": "filter", "score": "NONE"},
        ]
        assert determine_likely_cause(evidence) == "coverage_drop"

    def test_filter_change(self):
        evidence = [
            {"category": "coverage", "score": "NONE"},
            {"category": "filter", "score": "STRONG"},
        ]
        assert determine_likely_cause(evidence) == "filter_change"

    def test_unknown_when_all_none(self):
        evidence = [
            {"category": "coverage", "score": "NONE"},
            {"category": "filter", "score": "NONE"},
        ]
        assert determine_likely_cause(evidence) == "unknown"

    def test_strongest_wins(self):
        evidence = [
            {"category": "coverage", "score": "MODERATE"},
            {"category": "filter", "score": "STRONG"},
            {"category": "tool_version", "score": "WEAK"},
        ]
        assert determine_likely_cause(evidence) == "filter_change"

    def test_empty_evidence(self):
        assert determine_likely_cause([]) == "unknown"


class TestRequiresReview:
    def test_requires_review_when_no_strong(self):
        evidence = [
            {"category": "coverage", "score": "MODERATE"},
            {"category": "filter", "score": "WEAK"},
        ]
        assert requires_review(evidence) is True

    def test_no_review_when_strong(self):
        evidence = [
            {"category": "coverage", "score": "STRONG"},
        ]
        assert requires_review(evidence) is False

    def test_requires_review_empty(self):
        assert requires_review([]) is True
