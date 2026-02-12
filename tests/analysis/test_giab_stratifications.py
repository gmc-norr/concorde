"""Tests for GIAB stratification BED file loader."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from analysis.giab_stratifications import GIABStratificationLoader  # noqa: E402


@pytest.fixture
def giab_dir(tmp_path):
    """Create a mock GIAB stratification directory structure."""
    # LowComplexity
    lcr_dir = tmp_path / "LowComplexity"
    lcr_dir.mkdir()
    (lcr_dir / "GRCh38_AllHomopolymers.bed.gz").touch()
    (lcr_dir / "GRCh38_SimpleRepeats.bed.gz").touch()

    # SegmentalDuplications
    seg_dir = tmp_path / "SegmentalDuplications"
    seg_dir.mkdir()
    (seg_dir / "GRCh38_segdups.bed.gz").touch()

    # GCcontent
    gc_dir = tmp_path / "GCcontent"
    gc_dir.mkdir()
    (gc_dir / "GRCh38_gc15_slop50.bed.gz").touch()

    # Mappability
    mapp_dir = tmp_path / "Mappability"
    mapp_dir.mkdir()
    (mapp_dir / "GRCh38_lowmappabilityall.bed.gz").touch()

    # Irrelevant file (should be ignored)
    (tmp_path / "README.md").touch()

    return tmp_path


class TestGIABStratificationLoader:
    def test_discover_all_categories(self, giab_dir):
        loader = GIABStratificationLoader(giab_dir)
        discovered = loader.discover_bed_files()

        assert "LowComplexity" in discovered
        assert len(discovered["LowComplexity"]) == 2
        assert "SegmentalDuplications" in discovered
        assert "GCcontent" in discovered
        assert "Mappability" in discovered

    def test_filter_categories(self, giab_dir):
        loader = GIABStratificationLoader(giab_dir, categories=["LowComplexity"])
        discovered = loader.discover_bed_files()

        assert "LowComplexity" in discovered
        assert "SegmentalDuplications" not in discovered

    def test_get_region_beds_config(self, giab_dir):
        loader = GIABStratificationLoader(giab_dir)
        config = loader.get_region_beds_config()

        assert "low_complexity" in config
        assert "segmental_duplications" in config
        assert "gc_content" in config
        assert "mappability" in config
        # Each should be a path string
        assert config["low_complexity"].endswith(".bed.gz")

    def test_nonexistent_dir(self, tmp_path):
        loader = GIABStratificationLoader(tmp_path / "nonexistent")
        discovered = loader.discover_bed_files()
        assert discovered == {}

    def test_empty_dir(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        loader = GIABStratificationLoader(empty)
        discovered = loader.discover_bed_files()
        assert discovered == {}

    def test_none_bed_dir(self):
        loader = GIABStratificationLoader(None)
        discovered = loader.discover_bed_files()
        assert discovered == {}

    def test_get_all_bed_files(self, giab_dir):
        loader = GIABStratificationLoader(giab_dir)
        all_files = loader.get_all_bed_files()

        # Should have 5 total files across all categories
        assert len(all_files) == 5
        assert all(isinstance(k, tuple) for k in all_files)
