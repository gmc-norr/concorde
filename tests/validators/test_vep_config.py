"""Tests for VEP configuration validation."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from validators.config import validate_vep_config  # noqa: E402
from validators.base import ValidationError  # noqa: E402


class TestValidateVepConfig:
    """Tests for VEP config validation."""

    def test_disabled_passes_without_fields(self):
        validate_vep_config({"enabled": False})

    def test_enabled_minimal_passes(self):
        validate_vep_config({"enabled": True})

    def test_enabled_full_config_passes(self):
        validate_vep_config({
            "enabled": True,
            "vep_version": "110",
            "cache_version": "110",
            "genome_build": "GRCh38",
            "transcript_source": "ensembl",
            "transcript_selection": "canonical",
            "cache_dir": "/data/vep/cache",
            "plugins": ["CADD", "SpliceAI"],
        })

    def test_invalid_genome_build_raises(self):
        with pytest.raises(ValidationError, match="genome_build"):
            validate_vep_config({
                "enabled": True,
                "genome_build": "hg38",
            })

    def test_invalid_transcript_source_raises(self):
        with pytest.raises(ValidationError, match="transcript_source"):
            validate_vep_config({
                "enabled": True,
                "transcript_source": "gencode",
            })

    def test_invalid_transcript_selection_raises(self):
        with pytest.raises(ValidationError, match="transcript_selection"):
            validate_vep_config({
                "enabled": True,
                "transcript_selection": "longest",
            })

    def test_plugins_not_list_raises(self):
        with pytest.raises(ValidationError, match="plugins.*list"):
            validate_vep_config({
                "enabled": True,
                "plugins": "CADD",
            })

    def test_non_dict_raises(self):
        with pytest.raises(ValidationError, match="must be a dictionary"):
            validate_vep_config("not a dict")

    def test_grch37_passes(self):
        validate_vep_config({
            "enabled": True,
            "genome_build": "GRCh37",
        })

    def test_refseq_source_passes(self):
        validate_vep_config({
            "enabled": True,
            "transcript_source": "refseq",
        })

    def test_mane_select_passes(self):
        validate_vep_config({
            "enabled": True,
            "transcript_selection": "mane_select",
        })

    def test_most_severe_passes(self):
        validate_vep_config({
            "enabled": True,
            "transcript_selection": "most_severe",
        })
