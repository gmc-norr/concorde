"""Tests for multi-caller ensemble engine."""

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

from sqlalchemy import create_engine  # noqa: E402
from sqlalchemy.orm import sessionmaker  # noqa: E402

from models.base import Base  # noqa: E402
from models.run import Run  # noqa: E402
from models.variant import Variant  # noqa: E402

from analysis.ensemble import EnsembleEngine  # noqa: E402


@pytest.fixture
def db_session(tmp_path):
    """Create database with multi-caller runs."""
    db_path = tmp_path / "ensemble.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Run 1: GATK
    run1 = Run(
        sample="NA12878", caller="GATK", pipeline_version="v1",
        comparison_tool="happy", mode="germline",
    )
    session.add(run1)
    session.flush()

    # Variants: chr1:100 (shared), chr1:200 (shared), chr1:300 (GATK only)
    session.add_all([
        Variant(run_id=run1.id, chrom="chr1", pos=100, ref="A", alt="G",
                type="SNP", classification="TP"),
        Variant(run_id=run1.id, chrom="chr1", pos=200, ref="C", alt="T",
                type="SNP", classification="TP"),
        Variant(run_id=run1.id, chrom="chr1", pos=300, ref="G", alt="A",
                type="SNP", classification="FP"),
    ])

    # Run 2: DeepVariant
    run2 = Run(
        sample="NA12878", caller="DeepVariant", pipeline_version="v1",
        comparison_tool="happy", mode="germline",
    )
    session.add(run2)
    session.flush()

    # Variants: chr1:100 (shared), chr1:200 (shared), chr1:400 (DV only)
    session.add_all([
        Variant(run_id=run2.id, chrom="chr1", pos=100, ref="A", alt="G",
                type="SNP", classification="TP"),
        Variant(run_id=run2.id, chrom="chr1", pos=200, ref="C", alt="T",
                type="SNP", classification="TP"),
        Variant(run_id=run2.id, chrom="chr1", pos=400, ref="T", alt="C",
                type="SNP", classification="TP"),
    ])

    # Run 3: Strelka
    run3 = Run(
        sample="NA12878", caller="Strelka", pipeline_version="v1",
        comparison_tool="happy", mode="germline",
    )
    session.add(run3)
    session.flush()

    # Variants: chr1:100 (shared by all), chr1:400 (shared with DV)
    session.add_all([
        Variant(run_id=run3.id, chrom="chr1", pos=100, ref="A", alt="G",
                type="SNP", classification="TP"),
        Variant(run_id=run3.id, chrom="chr1", pos=400, ref="T", alt="C",
                type="SNP", classification="TP"),
    ])

    session.commit()
    yield session, [run1.id, run2.id, run3.id]
    session.close()


class TestEnsembleEngine:
    def test_find_runs_for_sample(self, db_session):
        session, run_ids = db_session
        engine = EnsembleEngine(session)
        runs = engine.find_runs_for_sample("NA12878")
        assert len(runs) == 3

    def test_cross_caller_concordance(self, db_session):
        session, run_ids = db_session
        engine = EnsembleEngine(session)
        result = engine.compute_cross_caller_concordance(run_ids)

        assert len(result["callers"]) == 3
        assert result["total_unique_variants"] == 4  # chr1:100, 200, 300, 400
        assert 1 in result["variants_by_caller_count"]  # Some unique to one caller
        assert 3 in result["variants_by_caller_count"]  # chr1:100 in all three

    def test_ensemble_intersection(self, db_session):
        session, run_ids = db_session
        engine = EnsembleEngine(session)
        result = engine.apply_ensemble_filter(run_ids, method="intersection")

        # Only chr1:100 is in all 3 callers
        assert len(result) == 1
        assert result[0]["pos"] == 100
        assert result[0]["caller_support"] == 3

    def test_ensemble_union(self, db_session):
        session, run_ids = db_session
        engine = EnsembleEngine(session)
        result = engine.apply_ensemble_filter(run_ids, method="union")

        # All 4 unique variants
        assert len(result) == 4

    def test_ensemble_majority_vote(self, db_session):
        session, run_ids = db_session
        engine = EnsembleEngine(session)
        result = engine.apply_ensemble_filter(run_ids, method="majority_vote")

        # Majority = 2+ out of 3
        # chr1:100 (3/3), chr1:200 (2/3), chr1:400 (2/3)
        # chr1:300 is only in 1/3 caller = excluded
        assert len(result) == 3
        positions = {r["pos"] for r in result}
        assert 100 in positions
        assert 200 in positions
        assert 400 in positions
        assert 300 not in positions
