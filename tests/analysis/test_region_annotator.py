"""Tests for region annotator and new stratification dimensions."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.region_annotator import RegionAnnotator  # noqa: E402
from analysis.stratification import (  # noqa: E402
    _assign_coverage_depth,
    _assign_gc_content,
    _assign_low_complexity,
    _assign_mappability,
    _assign_segdup,
)


class _MockVariant:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __getattr__(self, name):
        return None


# ── New Assign Functions ──


class TestAssignLowComplexity:
    def test_in_lcr(self):
        v = _MockVariant(in_low_complexity=True)
        assert _assign_low_complexity(v) == "low_complexity"

    def test_not_in_lcr(self):
        v = _MockVariant(in_low_complexity=False)
        assert _assign_low_complexity(v) == "non_lcr"

    def test_none(self):
        v = _MockVariant()
        assert _assign_low_complexity(v) is None


class TestAssignGcContent:
    BINS = [
        {"name": "0-25%", "min": 0.0, "max": 0.25},
        {"name": "25-55%", "min": 0.25, "max": 0.55},
        {"name": "55-100%", "min": 0.55, "max": 1.00},
    ]

    def test_low_gc(self):
        v = _MockVariant(gc_content=0.15)
        assert _assign_gc_content(v, self.BINS) == "0-25%"

    def test_mid_gc(self):
        v = _MockVariant(gc_content=0.40)
        assert _assign_gc_content(v, self.BINS) == "25-55%"

    def test_high_gc(self):
        v = _MockVariant(gc_content=0.70)
        assert _assign_gc_content(v, self.BINS) == "55-100%"

    def test_exact_boundary(self):
        v = _MockVariant(gc_content=0.25)
        assert _assign_gc_content(v, self.BINS) == "25-55%"

    def test_gc_1_0(self):
        v = _MockVariant(gc_content=1.0)
        assert _assign_gc_content(v, self.BINS) == "55-100%"

    def test_none_gc(self):
        v = _MockVariant()
        assert _assign_gc_content(v, self.BINS) is None

    def test_no_bins(self):
        v = _MockVariant(gc_content=0.5)
        assert _assign_gc_content(v, None) is None


class TestAssignSegdup:
    def test_in_segdup(self):
        v = _MockVariant(in_segdup=True)
        assert _assign_segdup(v) == "in_segdup"

    def test_not_in_segdup(self):
        v = _MockVariant(in_segdup=False)
        assert _assign_segdup(v) == "non_segdup"

    def test_none(self):
        v = _MockVariant()
        assert _assign_segdup(v) is None


class TestAssignMappability:
    BINS = [
        {"name": "low", "min": 0.0, "max": 0.5},
        {"name": "medium", "min": 0.5, "max": 0.9},
        {"name": "high", "min": 0.9, "max": 1.01},
    ]

    def test_low_mappability(self):
        v = _MockVariant(mappability_score=0.3)
        assert _assign_mappability(v, self.BINS) == "low"

    def test_medium_mappability(self):
        v = _MockVariant(mappability_score=0.7)
        assert _assign_mappability(v, self.BINS) == "medium"

    def test_high_mappability(self):
        v = _MockVariant(mappability_score=1.0)
        assert _assign_mappability(v, self.BINS) == "high"

    def test_none(self):
        v = _MockVariant()
        assert _assign_mappability(v, self.BINS) is None

    def test_no_bins(self):
        v = _MockVariant(mappability_score=0.5)
        assert _assign_mappability(v, None) is None


class TestAssignCoverageDepth:
    BINS = [
        {"name": "<10x", "min": 0, "max": 10},
        {"name": "10-30x", "min": 10, "max": 30},
        {"name": "30-50x", "min": 30, "max": 50},
        {"name": ">50x", "min": 50, "max": None},
    ]

    def test_low_coverage(self):
        v = _MockVariant(dp=5)
        assert _assign_coverage_depth(v, self.BINS) == "<10x"

    def test_medium_coverage(self):
        v = _MockVariant(dp=25)
        assert _assign_coverage_depth(v, self.BINS) == "10-30x"

    def test_high_coverage(self):
        v = _MockVariant(dp=45)
        assert _assign_coverage_depth(v, self.BINS) == "30-50x"

    def test_very_high_coverage(self):
        v = _MockVariant(dp=100)
        assert _assign_coverage_depth(v, self.BINS) == ">50x"

    def test_boundary(self):
        v = _MockVariant(dp=10)
        assert _assign_coverage_depth(v, self.BINS) == "10-30x"

    def test_none_dp(self):
        v = _MockVariant()
        assert _assign_coverage_depth(v, self.BINS) is None

    def test_no_bins(self):
        v = _MockVariant(dp=30)
        assert _assign_coverage_depth(v, None) is None


# ── RegionAnnotator ──


class TestRegionAnnotator:
    def test_no_config_skips(self):
        ann = RegionAnnotator()
        variants = [_MockVariant(chrom="chr1", pos=100)]
        result = ann.annotate_variants(variants)
        assert len(result) == 1
        # No annotations should be set since no sources configured

    def test_gc_content_computation(self):
        """Test GC content with mocked pysam FastaFile."""
        mock_fasta = MagicMock()
        mock_fasta.fetch.return_value = "GCGCGCGCGC"  # 100% GC

        ann = RegionAnnotator(reference_path="/fake/ref.fa", gc_window=10)
        ann._fasta = mock_fasta

        gc = ann.annotate_gc_content("chr1", 100)
        assert gc == pytest.approx(1.0)

    def test_gc_content_mixed(self):
        mock_fasta = MagicMock()
        mock_fasta.fetch.return_value = "GCATATAT"  # 2 GC out of 8 = 0.25

        ann = RegionAnnotator(reference_path="/fake/ref.fa", gc_window=8)
        ann._fasta = mock_fasta

        gc = ann.annotate_gc_content("chr1", 100)
        assert gc == pytest.approx(0.25)

    def test_tabix_overlaps_found(self):
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = ["chr1\t99\t200"]

        ann = RegionAnnotator(config={"low_complexity": "/fake/lcr.bed.gz"})
        ann._tabix_cache["/fake/lcr.bed.gz"] = mock_tabix

        result = ann.annotate_low_complexity("chr1", 100)
        assert result is True

    def test_tabix_overlaps_not_found(self):
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = []

        ann = RegionAnnotator(config={"low_complexity": "/fake/lcr.bed.gz"})
        ann._tabix_cache["/fake/lcr.bed.gz"] = mock_tabix

        result = ann.annotate_low_complexity("chr1", 100)
        assert result is False

    def test_missing_bed_returns_none(self):
        ann = RegionAnnotator(config={"low_complexity": ""})
        result = ann.annotate_low_complexity("chr1", 100)
        assert result is None

    def test_tabix_score(self):
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = ["chr1\t99\t200\t.\t0.85"]

        ann = RegionAnnotator(config={"mappability": "/fake/mapp.bed.gz"})
        ann._tabix_cache["/fake/mapp.bed.gz"] = mock_tabix

        score = ann.annotate_mappability("chr1", 100)
        assert score == pytest.approx(0.85)

    def test_annotate_variants_sets_attributes(self):
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = ["chr1\t99\t200"]

        ann = RegionAnnotator(config={"low_complexity": "/fake/lcr.bed.gz"})
        ann._tabix_cache["/fake/lcr.bed.gz"] = mock_tabix

        v = _MockVariant(chrom="chr1", pos=100)
        ann.annotate_variants([v])
        assert v.in_low_complexity is True
