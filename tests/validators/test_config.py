"""Tests for configuration validators."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from validators import (
    ValidationError,
    validate_comparison_tool,
    validate_decomposition_modes,
    validate_gene_sets_config,
    validate_config,
)


class TestValidateComparisonTool:
    """Tests for validate_comparison_tool function."""

    def test_happy_passes(self):
        """Test 'happy' is valid."""
        validate_comparison_tool("happy")  # Should not raise

    def test_rtg_passes(self):
        """Test 'rtg' is valid."""
        validate_comparison_tool("rtg")  # Should not raise

    def test_invalid_tool_raises(self):
        """Test invalid tool raises error."""
        with pytest.raises(ValidationError, match="Invalid comparison_tool"):
            validate_comparison_tool("invalid")


class TestValidateDecompositionModes:
    """Tests for validate_decomposition_modes function."""

    def test_valid_modes_pass(self):
        """Test valid modes pass."""
        validate_decomposition_modes(["decomposed"])
        validate_decomposition_modes(["non-decomposed"])
        validate_decomposition_modes(["decomposed", "non-decomposed"])

    def test_invalid_mode_raises(self):
        """Test invalid mode raises error."""
        with pytest.raises(ValidationError, match="Invalid decomposition_modes"):
            validate_decomposition_modes(["invalid"])

    def test_empty_list_raises(self):
        """Test empty list raises error."""
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_decomposition_modes([])


class TestValidateGeneSetsConfig:
    """Tests for validate_gene_sets_config function."""

    def test_valid_gene_set_passes(self):
        """Test valid gene set configuration passes."""
        gene_sets = [{
            "name": "ACMG59",
            "description": "ACMG genes",
            "version": "v3.1",
            "bed": "/path/to/acmg.bed"
        }]
        validate_gene_sets_config(gene_sets)  # Should not raise

    def test_missing_required_field_raises(self):
        """Test missing required field raises error."""
        gene_sets = [{
            "name": "ACMG59",
            "description": "ACMG genes",
            # Missing version and bed
        }]
        with pytest.raises(ValidationError, match="missing required fields"):
            validate_gene_sets_config(gene_sets)

    def test_empty_name_raises(self):
        """Test empty name raises error."""
        gene_sets = [{
            "name": "",
            "description": "ACMG genes",
            "version": "v3.1",
            "bed": "/path/to/bed"
        }]
        with pytest.raises(ValidationError, match="invalid name"):
            validate_gene_sets_config(gene_sets)

    def test_non_dict_gene_set_raises(self):
        """Test non-dictionary gene set raises error."""
        with pytest.raises(ValidationError, match="not a dictionary"):
            validate_gene_sets_config(["not a dict"])


class TestValidateConfig:
    """Tests for validate_config function."""

    def test_valid_config_passes(self, sample_config):
        """Test valid configuration passes."""
        config = sample_config()
        validate_config(config)  # Should not raise

    def test_missing_required_field_raises(self, sample_config):
        """Test missing required field raises error."""
        config = sample_config()
        del config["sample"]

        with pytest.raises(ValidationError, match="missing required fields"):
            validate_config(config)

    def test_empty_sample_raises(self, sample_config):
        """Test empty sample name raises error."""
        config = sample_config(sample="")

        with pytest.raises(ValidationError, match="must be a non-empty string"):
            validate_config(config)

    def test_invalid_comparison_tool_raises(self, sample_config):
        """Test invalid comparison tool raises error."""
        config = sample_config(comparison_tool="invalid")

        with pytest.raises(ValidationError, match="Invalid comparison_tool"):
            validate_config(config)

    def test_invalid_decomposition_modes_raises(self, sample_config):
        """Test invalid decomposition modes raise error."""
        config = sample_config(decomposition_modes=["invalid"])

        with pytest.raises(ValidationError, match="Invalid decomposition_modes"):
            validate_config(config)

    def test_empty_query_vcfs_raises(self, sample_config):
        """Test empty query_vcfs raises error."""
        config = sample_config(query_vcfs=[])

        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_config(config)

    def test_non_list_decomposition_modes_raises(self, sample_config):
        """Test non-list decomposition_modes raises error."""
        config = sample_config(decomposition_modes="decomposed")

        with pytest.raises(ValidationError, match="must be a list"):
            validate_config(config)

    def test_valid_gene_sets_passes(self, sample_config):
        """Test valid gene sets pass."""
        config = sample_config(gene_sets=[{
            "name": "ACMG59",
            "description": "ACMG genes",
            "version": "v3.1",
            "bed": "/path/to/bed"
        }])
        validate_config(config)  # Should not raise

    def test_invalid_gene_sets_raises(self, sample_config):
        """Test invalid gene sets raise error."""
        config = sample_config(gene_sets=[{
            "name": "ACMG59",
            # Missing required fields
        }])

        with pytest.raises(ValidationError, match="missing required fields"):
            validate_config(config)
