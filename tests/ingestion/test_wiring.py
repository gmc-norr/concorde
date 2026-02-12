"""Integration tests for data wiring: parser → ingest_helpers → models.

These tests verify that fields flow correctly through the full pipeline
rather than testing individual functions in isolation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import pandas as pd  # noqa: E402
from sqlalchemy import create_engine  # noqa: E402
from sqlalchemy.orm import sessionmaker  # noqa: E402

from models.base import Base  # noqa: E402
from models.variant import Variant  # noqa: E402

from ingestion.ingest_helpers import (  # noqa: E402
    _derive_zygosity,
    load_and_insert_variants,
)


# ── Zygosity derivation ──


class TestDeriveZygosity:
    """Verify zygosity is correctly derived from genotype strings."""

    def test_het(self):
        assert _derive_zygosity("0/1") == "HET"

    def test_het_phased(self):
        assert _derive_zygosity("0|1") == "HET"

    def test_het_multiallelic(self):
        assert _derive_zygosity("1/2") == "HET"

    def test_hom_alt(self):
        assert _derive_zygosity("1/1") == "HOM_ALT"

    def test_hom_ref(self):
        assert _derive_zygosity("0/0") == "HOM_REF"

    def test_none(self):
        assert _derive_zygosity(None) is None

    def test_missing_alleles(self):
        assert _derive_zygosity("./.") is None

    def test_partial_missing(self):
        assert _derive_zygosity("0/.") is None


# ── End-to-end variant loading with zygosity and indel_size ──


@pytest.fixture
def db_session(tmp_path):
    """Create an in-memory database session."""
    db_path = tmp_path / "wiring.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    yield session
    session.close()


class TestVariantLoadingWiring:
    """Verify that load_and_insert_variants correctly wires all fields."""

    def _write_variants_tsv(self, tmp_path, rows):
        """Write a variants TSV from a list of dicts."""
        df = pd.DataFrame(rows)
        path = tmp_path / "variants.tsv"
        df.to_csv(path, sep="\t", index=False)
        return str(path)

    def test_zygosity_derived_from_query_gt(self, tmp_path, db_session):
        """Zygosity should be derived from query_gt when not explicitly set."""
        from models.run import Run

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "0/1"},
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T",
             "type": "SNP", "classification": "TP",
             "truth_gt": "1/1", "query_gt": "1/1"},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        assert variants[0].zygosity == "HET"
        assert variants[1].zygosity == "HOM_ALT"

    def test_zygosity_falls_back_to_truth_gt(self, tmp_path, db_session):
        """For FN variants with only truth_gt, derive zygosity from truth."""
        from models.run import Run

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
             "type": "SNP", "classification": "FN",
             "truth_gt": "0/1"},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        assert variants[0].zygosity == "HET"
        assert variants[0].query_gt is None

    def test_indel_size_loaded(self, tmp_path, db_session):
        """Verify indel_size is loaded from TSV when present."""
        from models.run import Run

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "AT", "alt": "A",
             "type": "INDEL", "classification": "TP",
             "indel_size": 1},
            {"chrom": "chr1", "pos": 200, "ref": "A", "alt": "ATG",
             "type": "INDEL", "classification": "TP",
             "indel_size": 2},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        assert variants[0].indel_size == 1
        assert variants[1].indel_size == 2

    def test_gt_concordant_computed(self, tmp_path, db_session):
        """Verify gt_concordant is correctly computed from truth/query GT."""
        from models.run import Run

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "0/1"},
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "1/1"},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        assert variants[0].gt_concordant is True
        assert variants[1].gt_concordant is False

    def test_explicit_zygosity_preserved(self, tmp_path, db_session):
        """If zygosity is in TSV, use it rather than deriving from GT."""
        from models.run import Run

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "0/1",
             "zygosity": "HET"},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        assert variants[0].zygosity == "HET"


# ── Stratification uses zygosity correctly ──


class TestStratificationWiring:
    """Verify stratification sees derived zygosity for Het/Hom ratio."""

    def test_het_hom_ratio_with_derived_zygosity(self, tmp_path, db_session):
        """Extended metrics should compute correct Het/Hom from derived zygosity."""
        from models.run import Run
        from analysis.extended_metrics import compute_het_hom_ratio

        run = Run(sample="S", caller="C", pipeline_version="v1",
                  comparison_tool="h", mode="germline")
        db_session.add(run)
        db_session.flush()

        path = self._write_variants_tsv(tmp_path, [
            {"chrom": "chr1", "pos": 100, "ref": "A", "alt": "G",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "0/1"},
            {"chrom": "chr1", "pos": 200, "ref": "C", "alt": "T",
             "type": "SNP", "classification": "TP",
             "truth_gt": "0/1", "query_gt": "0/1"},
            {"chrom": "chr1", "pos": 300, "ref": "G", "alt": "A",
             "type": "SNP", "classification": "TP",
             "truth_gt": "1/1", "query_gt": "1/1"},
        ])

        variants = load_and_insert_variants(path, run.id, db_session)
        db_session.flush()

        # Verify zygosity was derived
        assert variants[0].zygosity == "HET"
        assert variants[1].zygosity == "HET"
        assert variants[2].zygosity == "HOM_ALT"

        # Verify Het/Hom ratio computation
        result = compute_het_hom_ratio(variants)
        assert result["het_count"] == 2
        assert result["hom_count"] == 1
        assert result["het_hom_ratio"] == pytest.approx(2.0)

    def _write_variants_tsv(self, tmp_path, rows):
        df = pd.DataFrame(rows)
        path = tmp_path / "variants.tsv"
        df.to_csv(path, sep="\t", index=False)
        return str(path)
