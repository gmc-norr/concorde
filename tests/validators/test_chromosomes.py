"""Tests for chromosome normalization utilities."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from validators import (
    ValidationError,
    normalize_chromosome,
    get_chromosome_format,
    validate_chromosome_format,
)


class TestNormalizeChromosome:
    def test_adds_chr_prefix(self):
        assert normalize_chromosome("1") == "chr1"
        assert normalize_chromosome("X") == "chrX"

    def test_preserves_chr_prefix(self):
        assert normalize_chromosome("chr1") == "chr1"
        assert normalize_chromosome("chrX") == "chrX"

    def test_normalizes_mt(self):
        assert normalize_chromosome("MT") == "chrM"
        assert normalize_chromosome("chrMT") == "chrM"
        assert normalize_chromosome("M") == "chrM"


class TestGetChromosomeFormat:
    def test_detects_with_chr(self):
        chroms = {"chr1", "chr2", "chrX"}
        assert get_chromosome_format(chroms) == "with_chr"

    def test_detects_without_chr(self):
        chroms = {"1", "2", "X"}
        assert get_chromosome_format(chroms) == "without_chr"

    def test_detects_mixed(self):
        chroms = {"chr1", "2", "X"}
        assert get_chromosome_format(chroms) == "mixed"


class TestNormalizeChromosomeEdgeCases:
    """Edge case tests for chromosome normalization."""

    def test_mitochondrial_variants(self):
        """Test normalization of mitochondrial chromosome variants."""
        # All MT variations normalize to chrM (standard behavior)
        assert normalize_chromosome("M") == "chrM"
        assert normalize_chromosome("MT") == "chrM"  # MT -> chrM
        assert normalize_chromosome("chrM") == "chrM"
        assert normalize_chromosome("chrMT") == "chrM"  # chrMT -> chrM

    def test_leading_zeros(self):
        """Test normalization behavior with leading zeros."""
        # Leading zeros are currently preserved (documents actual behavior)
        assert normalize_chromosome("01") == "chr01"
        assert normalize_chromosome("chr01") == "chr01"
        assert normalize_chromosome("09") == "chr09"
        # Note: Future enhancement could strip leading zeros

    def test_nonstandard_contigs(self):
        """Test normalization of non-standard contigs."""
        # HLA contigs should pass through unchanged
        hla = "HLA-A*01:01"
        assert normalize_chromosome(hla) == f"chr{hla}"

        # Alt scaffolds
        assert normalize_chromosome("chr1_KI270706v1_random") == "chr1_KI270706v1_random"


class TestValidateChromosomeFormat:
    def test_valid_with_chr_passes(self):
        df = pd.DataFrame({"chrom": ["chr1", "chr2", "chrX"]})
        validate_chromosome_format(df, "with_chr", "test.tsv")

    def test_valid_without_chr_passes(self):
        df = pd.DataFrame({"chrom": ["1", "2", "X"]})
        validate_chromosome_format(df, "without_chr", "test.tsv")

    def test_mixed_format_raises(self):
        df = pd.DataFrame({"chrom": ["chr1", "2", "X"]})
        with pytest.raises(ValidationError, match="mixed chromosome formats"):
            validate_chromosome_format(df, "with_chr", "test.tsv")

    def test_wrong_format_raises(self):
        df = pd.DataFrame({"chrom": ["chr1", "chr2"]})
        with pytest.raises(ValidationError, match="expected without_chr"):
            validate_chromosome_format(df, "without_chr", "test.tsv")
