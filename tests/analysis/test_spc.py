"""Tests for Statistical Process Control engine."""

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
from models.stratified_metric import StratifiedMetric  # noqa: E402

from analysis.spc import SPCEngine  # noqa: E402


@pytest.fixture
def db_session(tmp_path):
    """Create database with historical stratified metrics."""
    db_path = tmp_path / "spc.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Create 15 runs with consistent precision ~0.95 for SNP stratum
    for i in range(15):
        run = Run(
            sample="NA12878", caller="GATK", pipeline_version=f"v{i}",
            comparison_tool="happy", mode="germline",
        )
        session.add(run)
        session.flush()

        session.add(StratifiedMetric(
            run_id=run.id,
            dimension="variant_class",
            stratum="SNP",
            variant_type="ALL",
            precision=0.95 + (i - 7) * 0.002,
            recall=0.90 + (i - 7) * 0.001,
            f1=0.92 + (i - 7) * 0.001,
            tp_count=100,
            fp_count=5,
            fn_count=10,
            total_variants=115,
            low_confidence=False,
        ))

    session.commit()
    yield session
    session.close()


class TestSPCEngine:
    def test_compute_control_limits(self, db_session):
        engine = SPCEngine(db_session, config={"min_data_points": 5})
        limits = engine.compute_control_limits(
            "precision", dimension="variant_class", stratum="SNP"
        )

        assert limits is not None
        assert limits["data_points"] == 15
        assert limits["mean"] == pytest.approx(0.95, abs=0.01)
        assert limits["stdev"] > 0
        assert limits["ucl"] > limits["mean"]
        assert limits["lcl"] < limits["mean"]

    def test_insufficient_data(self, db_session):
        engine = SPCEngine(db_session, config={"min_data_points": 100})
        limits = engine.compute_control_limits("precision")
        assert limits is None

    def test_evaluate_run_in_control(self, db_session):
        engine = SPCEngine(db_session, config={"min_data_points": 5})
        # The latest run (id=15) should be in control
        results = engine.evaluate_run(15)
        assert len(results) > 0
        # All should be in control since data is within normal range
        for r in results:
            if not r.get("insufficient_data"):
                assert r["in_control"] is True

    def test_generate_chart_data(self, db_session):
        engine = SPCEngine(db_session, config={"min_data_points": 5})
        chart = engine.generate_chart_data(
            "precision", dimension="variant_class", stratum="SNP"
        )

        assert len(chart["run_ids"]) == 15
        assert len(chart["values"]) == 15
        assert chart["mean"] is not None
        assert chart["ucl"] is not None
        assert chart["lcl"] is not None

    def test_custom_sigma_multiplier(self, db_session):
        engine_tight = SPCEngine(db_session, config={"sigma_multiplier": 1.0, "min_data_points": 5})
        engine_loose = SPCEngine(db_session, config={"sigma_multiplier": 5.0, "min_data_points": 5})

        limits_tight = engine_tight.compute_control_limits("precision")
        limits_loose = engine_loose.compute_control_limits("precision")

        assert limits_tight["ucl"] - limits_tight["lcl"] < limits_loose["ucl"] - limits_loose["lcl"]
