"""Tests for data quality validators."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from validators import (
    ValidationError,
    validate_dna_sequence,
    validate_probability,
    validate_percentage,
    validate_coverage_metrics,
)


class TestValidateDNASequence:
    def test_valid_dna_passes(self):
        validate_dna_sequence("ACGT")
        validate_dna_sequence("ACGTN")
        validate_dna_sequence("acgtn")  # Case insensitive

    def test_invalid_character_raises(self):
        with pytest.raises(ValidationError, match="invalid DNA characters"):
            validate_dna_sequence("ACGX")

    def test_empty_sequence_raises(self):
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_dna_sequence("")


class TestValidateProbability:
    def test_valid_probability_passes(self):
        validate_probability(0.0, "test")
        validate_probability(0.5, "test")
        validate_probability(1.0, "test")

    def test_too_large_raises(self):
        with pytest.raises(ValidationError, match="between 0.0 and 1.0"):
            validate_probability(1.5, "test")

    def test_negative_raises(self):
        with pytest.raises(ValidationError, match="between 0.0 and 1.0"):
            validate_probability(-0.1, "test")


class TestValidatePercentage:
    def test_valid_percentage_passes(self):
        validate_percentage(0, "test")
        validate_percentage(50.5, "test")
        validate_percentage(100, "test")

    def test_too_large_raises(self):
        with pytest.raises(ValidationError, match="between 0 and 100"):
            validate_percentage(150, "test")

    def test_negative_raises(self):
        with pytest.raises(ValidationError, match="between 0 and 100"):
            validate_percentage(-10, "test")


class TestValidateCoverageMetrics:
    def test_valid_metrics_pass(self):
        metrics = {
            "mean_coverage": 32.0,
            "median_coverage": 30.5,
            "pct_bases_10x": 99.1,
            "mapping_rate": 99.8
        }
        validate_coverage_metrics(metrics)

    def test_negative_coverage_raises(self):
        metrics = {"mean_coverage": -5}
        with pytest.raises(ValidationError, match="non-negative"):
            validate_coverage_metrics(metrics)

    def test_invalid_percentage_raises(self):
        metrics = {"mapping_rate": 150}
        with pytest.raises(ValidationError, match="between 0 and 100"):
            validate_coverage_metrics(metrics)

    def test_mean_less_than_median_warns(self):
        metrics = {"mean_coverage": 20, "median_coverage": 30}
        with pytest.raises(ValidationError, match="less than"):
            validate_coverage_metrics(metrics)


class TestEdgeCases:
    """Additional edge case tests for data quality validators."""

    def test_validate_dna_iupac_ambiguity_codes(self):
        """Test DNA validation with IUPAC ambiguity codes."""
        # Standard codes plus N should work
        validate_dna_sequence("ACGTN")
        # Full IUPAC ambiguity codes - currently not supported
        # If we want to support them, add: RYWSMKHBVDN
        with pytest.raises(ValidationError, match="invalid DNA characters"):
            validate_dna_sequence("ACGTR")  # R = puRine (A or G)

    def test_validate_probability_infinity(self):
        """Test probability validation rejects infinity."""
        with pytest.raises(ValidationError, match="between 0.0 and 1.0"):
            validate_probability(float('inf'), "test")
        with pytest.raises(ValidationError, match="between 0.0 and 1.0"):
            validate_probability(float('-inf'), "test")

    def test_validate_percentage_extreme_boundary(self):
        """Test percentage validation at extreme boundaries."""
        validate_percentage(0.0, "test")   # Exactly 0
        validate_percentage(100.0, "test")  # Exactly 100
        validate_percentage(0.001, "test")  # Very small positive
        validate_percentage(99.999, "test") # Very close to 100

    def test_validate_coverage_metrics_empty_dict(self):
        """Test coverage validation with empty metrics dict."""
        # Empty dict should be handled gracefully
        validate_coverage_metrics({})

    def test_validate_coverage_metrics_partial_data(self):
        """Test coverage validation with only some metrics."""
        # Should handle partial data
        validate_coverage_metrics({"mean_coverage": 30.0})
        validate_coverage_metrics({"mapping_rate": 95.0})
