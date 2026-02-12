"""Tests for QC metric trending z-score analysis."""

from __future__ import annotations

import sys
from pathlib import Path

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
from models.qc_metric import QCMetric  # noqa: E402
from models.run import Run  # noqa: E402

from analysis.qc_trending import compute_historical_stats, compute_trending_report  # noqa: E402


@pytest.fixture
def db_session(tmp_path):
    """Create an in-memory database with historical QC data."""
    db_path = tmp_path / "trending.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Create historical runs with consistent metrics
    for i in range(10):
        run = Run(
            sample="NA12878", caller="GATK", pipeline_version=f"v{i}",
            comparison_tool="happy", mode="germline",
        )
        session.add(run)
        session.flush()

        # Ti/Tv ratio hovering around 2.1 with small variance
        session.add(QCMetric(
            run_id=run.id, metric_name="titv_ratio",
            metric_value_float=2.1 + (i - 5) * 0.01,
            metric_source="concorde", metric_category="variant_quality",
        ))
        # Het/Hom ratio around 1.5
        session.add(QCMetric(
            run_id=run.id, metric_name="het_hom_ratio",
            metric_value_float=1.5 + (i - 5) * 0.02,
            metric_source="concorde", metric_category="variant_quality",
        ))

    # Create "current" run with a normal Ti/Tv but outlier Het/Hom
    current_run = Run(
        sample="NA12878", caller="GATK", pipeline_version="v_current",
        comparison_tool="happy", mode="germline",
    )
    session.add(current_run)
    session.flush()
    session.add(QCMetric(
        run_id=current_run.id, metric_name="titv_ratio",
        metric_value_float=2.1,
        metric_source="concorde", metric_category="variant_quality",
    ))
    session.add(QCMetric(
        run_id=current_run.id, metric_name="het_hom_ratio",
        metric_value_float=3.0,  # Far outlier!
        metric_source="concorde", metric_category="variant_quality",
    ))

    session.commit()
    yield session, current_run.id
    session.close()


class TestComputeHistoricalStats:
    def test_sufficient_history(self, db_session):
        session, current_run_id = db_session
        stats = compute_historical_stats(session, "titv_ratio", current_run_id, min_runs=5)

        assert stats is not None
        assert stats["count"] == 10
        assert abs(stats["mean"] - 2.1) < 0.05
        assert stats["stdev"] > 0

    def test_insufficient_history(self, db_session):
        session, current_run_id = db_session
        stats = compute_historical_stats(
            session, "titv_ratio", current_run_id, min_runs=50
        )
        assert stats is None

    def test_nonexistent_metric(self, db_session):
        session, current_run_id = db_session
        stats = compute_historical_stats(
            session, "nonexistent_metric", current_run_id, min_runs=1
        )
        assert stats is None


class TestComputeTrendingReport:
    def test_basic_trending(self, db_session):
        session, current_run_id = db_session
        report = compute_trending_report(session, current_run_id)

        assert len(report) == 2  # titv_ratio and het_hom_ratio
        names = {r["metric_name"] for r in report}
        assert "titv_ratio" in names
        assert "het_hom_ratio" in names

    def test_outlier_flagged(self, db_session):
        session, current_run_id = db_session
        report = compute_trending_report(
            session, current_run_id,
            config={"z_score_threshold": 2.0},
        )

        hethom = next(r for r in report if r["metric_name"] == "het_hom_ratio")
        assert hethom["flagged"] is True
        assert hethom["z_score"] is not None
        assert abs(hethom["z_score"]) > 2.0

    def test_normal_not_flagged(self, db_session):
        session, current_run_id = db_session
        report = compute_trending_report(
            session, current_run_id,
            config={"z_score_threshold": 2.0},
        )

        titv = next(r for r in report if r["metric_name"] == "titv_ratio")
        assert titv["flagged"] is False

    def test_custom_threshold(self, db_session):
        session, current_run_id = db_session
        # Very strict threshold should flag both
        report = compute_trending_report(
            session, current_run_id,
            config={"z_score_threshold": 0.01},
        )

        flagged = [r for r in report if r["flagged"]]
        assert len(flagged) >= 1
