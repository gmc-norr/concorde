"""Tests for the verification workflow orchestrator."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Add scripts and project root to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from sqlalchemy import create_engine  # noqa: E402
from sqlalchemy.orm import sessionmaker  # noqa: E402

from models.base import Base  # noqa: E402
from models.baseline import Baseline  # noqa: E402
from models.run import Run  # noqa: E402
from models.variant import Variant  # noqa: E402
from models.verification import VariantTransition, VerificationResult  # noqa: E402


@pytest.fixture
def db_session(tmp_path):
    """Create an in-memory database with test data."""
    db_path = tmp_path / "test.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Create baseline run
    base_run = Run(
        sample="NA12878", caller="GATK", pipeline_version="v1.0",
        comparison_tool="happy", mode="germline",
    )
    session.add(base_run)
    session.flush()

    # Create baseline variants
    base_variants = [
        Variant(run_id=base_run.id, chrom="chr1", pos=100, ref="A", alt="G",
                type="SNP", classification="TP", dp=30, mq=60.0, qual=100.0),
        Variant(run_id=base_run.id, chrom="chr1", pos=200, ref="C", alt="T",
                type="SNP", classification="TP", dp=25, mq=55.0, qual=80.0),
        Variant(run_id=base_run.id, chrom="chr1", pos=300, ref="G", alt="A",
                type="SNP", classification="FP", dp=10, mq=30.0, qual=20.0),
    ]
    session.add_all(base_variants)

    # Create locked baseline
    baseline = Baseline(
        name="test_baseline", run_id=base_run.id, locked=True,
    )
    session.add(baseline)

    # Create verification run
    verify_run = Run(
        sample="NA12878", caller="GATK", pipeline_version="v1.1",
        comparison_tool="happy", mode="germline",
    )
    session.add(verify_run)
    session.flush()

    # Create verification variants (with some differences)
    verify_variants = [
        Variant(run_id=verify_run.id, chrom="chr1", pos=100, ref="A", alt="G",
                type="SNP", classification="TP", dp=30, mq=60.0, qual=100.0),
        Variant(run_id=verify_run.id, chrom="chr1", pos=200, ref="C", alt="T",
                type="SNP", classification="FN", dp=5, mq=20.0, qual=10.0),  # Changed!
        Variant(run_id=verify_run.id, chrom="chr1", pos=300, ref="G", alt="A",
                type="SNP", classification="FP", dp=10, mq=30.0, qual=20.0),
    ]
    session.add_all(verify_variants)

    session.commit()
    yield session, str(db_path)
    session.close()


class TestRunVerification:
    def test_basic_verification(self, db_session):
        from analysis.run_verification import run_verification

        session, db_path = db_session
        result = run_verification(
            database=db_path,
            baseline_name="test_baseline",
            mode="germline",
        )

        assert result["verdict"] in ("PASS", "FAIL", "REVIEW_REQUIRED")
        assert result["total_transitions"] >= 1  # At least the TP→FN change
        assert isinstance(result["drift_count"], int)
        assert isinstance(result["biological_count"], int)

    def test_verification_writes_to_db(self, db_session):
        from analysis.run_verification import run_verification

        session, db_path = db_session
        run_verification(
            database=db_path,
            baseline_name="test_baseline",
        )

        # Check that VerificationResult was created
        vr = session.query(VerificationResult).first()
        assert vr is not None
        assert vr.verdict in ("PASS", "FAIL", "REVIEW_REQUIRED")
        assert vr.total_transitions >= 1

        # Check that VariantTransitions were created
        transitions = session.query(VariantTransition).all()
        assert len(transitions) >= 1

    def test_missing_baseline_raises(self, db_session):
        from analysis.run_verification import run_verification

        _, db_path = db_session
        with pytest.raises(ValueError, match="No locked baseline found"):
            run_verification(
                database=db_path,
                baseline_name="nonexistent_baseline",
            )

    def test_no_changes_pass(self, db_session, tmp_path):
        """When baseline and verification have identical variants, expect PASS."""
        from analysis.run_verification import run_verification

        session, _ = db_session

        # Create a new database with identical runs
        db_path = tmp_path / "identical.db"
        engine = create_engine(f"sqlite:///{db_path}")
        Base.metadata.create_all(engine)
        Session2 = sessionmaker(bind=engine)
        session2 = Session2()

        run1 = Run(sample="S1", caller="C1", pipeline_version="v1",
                    comparison_tool="happy", mode="germline")
        session2.add(run1)
        session2.flush()
        session2.add(Variant(run_id=run1.id, chrom="chr1", pos=100, ref="A",
                             alt="G", type="SNP", classification="TP"))
        baseline = Baseline(name="bl", run_id=run1.id, locked=True)
        session2.add(baseline)

        run2 = Run(sample="S1", caller="C1", pipeline_version="v1",
                    comparison_tool="happy", mode="germline")
        session2.add(run2)
        session2.flush()
        session2.add(Variant(run_id=run2.id, chrom="chr1", pos=100, ref="A",
                             alt="G", type="SNP", classification="TP"))
        session2.commit()

        result = run_verification(database=str(db_path), baseline_name="bl")
        assert result["verdict"] == "PASS"
        assert result["total_transitions"] == 0
        session2.close()
