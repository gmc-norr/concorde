"""Tests for somatic mode configuration validation."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from validators.config import (  # noqa: E402
    validate_mode,
    validate_somatic_config,
    validate_truth_matching_config,
    validate_comparison_tool,
    validate_config,
)
from validators.base import ValidationError  # noqa: E402


class TestValidateMode:
    """Tests for mode validation."""

    def test_germline_passes(self):
        validate_mode("germline")

    def test_somatic_passes(self):
        validate_mode("somatic")

    def test_invalid_mode_raises(self):
        with pytest.raises(ValidationError, match="Invalid mode"):
            validate_mode("invalid")

    def test_empty_mode_raises(self):
        with pytest.raises(ValidationError, match="Invalid mode"):
            validate_mode("")


class TestValidateComparisonToolWithMode:
    """Tests for mode-aware comparison tool validation."""

    def test_happy_germline_passes(self):
        validate_comparison_tool("happy", "germline")

    def test_rtg_germline_passes(self):
        validate_comparison_tool("rtg", "germline")

    def test_sompy_somatic_passes(self):
        validate_comparison_tool("sompy", "somatic")

    def test_rtg_somatic_passes(self):
        validate_comparison_tool("rtg", "somatic")

    def test_internal_somatic_passes(self):
        validate_comparison_tool("internal", "somatic")

    def test_happy_somatic_raises(self):
        with pytest.raises(ValidationError, match="Invalid comparison_tool.*somatic"):
            validate_comparison_tool("happy", "somatic")

    def test_sompy_germline_raises(self):
        with pytest.raises(ValidationError, match="Invalid comparison_tool.*germline"):
            validate_comparison_tool("sompy", "germline")


class TestValidateSomaticConfig:
    """Tests for somatic config block validation."""

    def test_valid_tumor_normal_passes(self):
        validate_somatic_config({
            "tumor_sample": "TUMOR",
            "normal_sample": "NORMAL",
            "calling_mode": "tumor_normal",
        })

    def test_valid_tumor_only_passes(self):
        validate_somatic_config({
            "tumor_sample": "TUMOR",
            "calling_mode": "tumor_only",
        })

    def test_missing_tumor_sample_raises(self):
        with pytest.raises(ValidationError, match="tumor_sample"):
            validate_somatic_config({
                "calling_mode": "tumor_only",
            })

    def test_empty_tumor_sample_raises(self):
        with pytest.raises(ValidationError, match="tumor_sample"):
            validate_somatic_config({
                "tumor_sample": "",
                "calling_mode": "tumor_only",
            })

    def test_missing_calling_mode_raises(self):
        with pytest.raises(ValidationError, match="calling_mode"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
            })

    def test_invalid_calling_mode_raises(self):
        with pytest.raises(ValidationError, match="Invalid somatic calling_mode"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
                "calling_mode": "invalid",
            })

    def test_tumor_normal_without_normal_raises(self):
        with pytest.raises(ValidationError, match="normal_sample"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_normal",
            })

    def test_vaf_tolerance_valid(self):
        validate_somatic_config({
            "tumor_sample": "TUMOR",
            "calling_mode": "tumor_only",
            "vaf_tolerance": 0.10,
        })

    def test_vaf_tolerance_invalid_raises(self):
        with pytest.raises(ValidationError, match="vaf_tolerance"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_only",
                "vaf_tolerance": 1.5,
            })

    def test_vaf_tolerance_negative_raises(self):
        with pytest.raises(ValidationError, match="vaf_tolerance"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_only",
                "vaf_tolerance": -0.1,
            })

    def test_min_vaf_valid(self):
        validate_somatic_config({
            "tumor_sample": "TUMOR",
            "calling_mode": "tumor_only",
            "min_vaf": 0.05,
        })

    def test_min_vaf_invalid_raises(self):
        with pytest.raises(ValidationError, match="min_vaf"):
            validate_somatic_config({
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_only",
                "min_vaf": 2.0,
            })

    def test_non_dict_raises(self):
        with pytest.raises(ValidationError, match="must be a dictionary"):
            validate_somatic_config("not a dict")


class TestValidateTruthMatchingConfig:
    """Tests for truth matching config validation."""

    def test_empty_config_passes(self):
        validate_truth_matching_config({}, "germline")

    def test_valid_germline_config_passes(self):
        validate_truth_matching_config({
            "germline": {
                "tool": "happy",
                "genotype_match": True,
            }
        }, "germline")

    def test_valid_somatic_config_passes(self):
        validate_truth_matching_config({
            "somatic": {
                "tool": "sompy",
                "vaf_tolerance": 0.10,
            }
        }, "somatic")

    def test_invalid_germline_tool_raises(self):
        with pytest.raises(ValidationError, match="Invalid truth_matching"):
            validate_truth_matching_config({
                "germline": {
                    "tool": "sompy",
                }
            }, "germline")

    def test_invalid_somatic_tool_raises(self):
        with pytest.raises(ValidationError, match="Invalid truth_matching"):
            validate_truth_matching_config({
                "somatic": {
                    "tool": "happy",
                }
            }, "somatic")

    def test_non_dict_raises(self):
        with pytest.raises(ValidationError, match="must be a dictionary"):
            validate_truth_matching_config("not a dict", "germline")

    def test_missing_mode_key_passes(self):
        """Config without a matching mode key should pass (use defaults)."""
        validate_truth_matching_config({
            "somatic": {"tool": "sompy"}
        }, "germline")


class TestValidateConfigWithMode:
    """Tests for validate_config with mode support."""

    def test_germline_config_backward_compat(self, sample_config):
        """Config without mode key should default to germline."""
        config = sample_config()
        # No mode key - should default to germline
        validate_config(config)

    def test_explicit_germline_passes(self, sample_config):
        config = sample_config(mode="germline")
        validate_config(config)

    def test_somatic_with_valid_config_passes(self, sample_config):
        config = sample_config(
            mode="somatic",
            comparison_tool="sompy",
            somatic={
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_only",
            },
        )
        validate_config(config)

    def test_somatic_without_somatic_block_raises(self, sample_config):
        config = sample_config(
            mode="somatic",
            comparison_tool="sompy",
        )
        with pytest.raises(ValidationError, match="somatic"):
            validate_config(config)

    def test_somatic_with_germline_tool_raises(self, sample_config):
        config = sample_config(
            mode="somatic",
            comparison_tool="happy",
            somatic={
                "tumor_sample": "TUMOR",
                "calling_mode": "tumor_only",
            },
        )
        with pytest.raises(ValidationError, match="Invalid comparison_tool.*somatic"):
            validate_config(config)
