"""Tests for baseline management (Spec S12)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Add project root and scripts directory to path
PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

from models.base import Base  # noqa: E402
from models.baseline import Baseline, BaselineEnvelope, FixtureRecord  # noqa: E402
from models.run import Run  # noqa: E402
from models.stratified_metric import StratifiedMetric  # noqa: E402
from analysis.baseline_manager import BaselineManager  # noqa: E402


@pytest.fixture
def db_session():
    """Create in-memory test database with schema."""
    engine = create_engine(
        "sqlite:///:memory:",
        connect_args={"check_same_thread": False},
    )

    @sqlalchemy.event.listens_for(engine, "connect")
    def set_sqlite_pragma(dbapi_conn, connection_record):
        cursor = dbapi_conn.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    yield session
    session.close()
    Base.metadata.drop_all(engine)
    engine.dispose()


@pytest.fixture
def run_with_metrics(db_session):
    """Create a run with stratified metrics for baseline testing."""
    run = Run(
        sample="NA12878",
        caller="GATK-HC",
        pipeline_version="v4.5.0",
        decomposition_mode="decomposed",
        comparison_tool="happy",
        mode="germline",
    )
    db_session.add(run)
    db_session.flush()

    # Add stratified metrics
    metrics = [
        StratifiedMetric(
            run_id=run.id,
            dimension="variant_class",
            stratum="SNP",
            variant_type="ALL",
            precision=0.9975,
            recall=0.9960,
            f1=0.9968,
            tp_count=2950,
            fp_count=7,
            fn_count=12,
            total_variants=2969,
            low_confidence=False,
        ),
        StratifiedMetric(
            run_id=run.id,
            dimension="variant_class",
            stratum="INDEL",
            variant_type="ALL",
            precision=0.9900,
            recall=0.9850,
            f1=0.9875,
            tp_count=490,
            fp_count=5,
            fn_count=7,
            total_variants=502,
            low_confidence=False,
        ),
        StratifiedMetric(
            run_id=run.id,
            dimension="functional_impact",
            stratum="HIGH",
            variant_type="ALL",
            precision=1.0,
            recall=0.9500,
            f1=0.9744,
            tp_count=19,
            fp_count=0,
            fn_count=1,
            total_variants=20,
            low_confidence=False,
        ),
    ]
    db_session.add_all(metrics)
    db_session.commit()
    db_session.refresh(run)
    return run


class TestBaselineCreation:
    def test_create_baseline(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="v1.0_germline_GRCh38",
            run_id=run_with_metrics.id,
            mode="germline",
        )
        db_session.commit()

        assert baseline.id is not None
        assert baseline.name == "v1.0_germline_GRCh38"
        assert baseline.run_id == run_with_metrics.id
        assert baseline.mode == "germline"
        assert baseline.locked is False
        assert baseline.created_at is not None

    def test_create_baseline_with_config(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="v1.0_test",
            run_id=run_with_metrics.id,
            config_snapshot="mode: germline\nreference: GRCh38",
            stratification_snapshot="dimensions:\n  variant_class: true",
            pipeline_version="v4.5.0",
        )
        db_session.commit()

        assert baseline.config_snapshot == "mode: germline\nreference: GRCh38"
        assert baseline.stratification_snapshot is not None
        assert baseline.pipeline_version == "v4.5.0"

    def test_create_baseline_duplicate_name_raises(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        manager.create_baseline(name="dup_test", run_id=run_with_metrics.id)
        db_session.commit()

        with pytest.raises(ValueError, match="already exists"):
            manager.create_baseline(name="dup_test", run_id=run_with_metrics.id)

    def test_create_baseline_invalid_run_raises(self, db_session):
        manager = BaselineManager(db_session)
        with pytest.raises(ValueError, match="not found"):
            manager.create_baseline(name="bad_run", run_id=99999)

    def test_create_baseline_inherits_pipeline_version(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="v_inherit",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        assert baseline.pipeline_version == "v4.5.0"


class TestAutoDerivation:
    def test_envelopes_auto_derived(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="auto_test",
            run_id=run_with_metrics.id,
            tolerance=0.005,
        )
        db_session.commit()

        # 3 strat metrics × 3 envelope metrics (precision, recall, f1) = 9
        assert len(baseline.envelopes) == 9

    def test_envelope_values_with_tolerance(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="tol_test",
            run_id=run_with_metrics.id,
            tolerance=0.01,
        )
        db_session.commit()

        # Find the SNP precision envelope
        snp_prec = [
            e
            for e in baseline.envelopes
            if e.dimension == "variant_class"
            and e.stratum == "SNP"
            and e.metric_name == "precision"
        ]
        assert len(snp_prec) == 1
        env = snp_prec[0]

        assert env.expected_value == pytest.approx(0.9975)
        assert env.lower_bound == pytest.approx(0.9875)
        assert env.upper_bound == pytest.approx(1.0)  # Clamped to 1.0
        assert env.manually_set is False

    def test_envelope_lower_bound_clamped_to_zero(self, db_session):
        """Test that lower bounds are clamped to 0.0."""
        run = Run(
            sample="LOW",
            caller="GATK-HC",
            pipeline_version="v1.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
        )
        db_session.add(run)
        db_session.flush()

        sm = StratifiedMetric(
            run_id=run.id,
            dimension="variant_class",
            stratum="SNP",
            variant_type="ALL",
            precision=0.002,
            recall=0.001,
            f1=0.001,
            total_variants=5,
            low_confidence=True,
        )
        db_session.add(sm)
        db_session.commit()

        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="low_test",
            run_id=run.id,
            tolerance=0.005,
        )
        db_session.commit()

        prec_env = [
            e for e in baseline.envelopes if e.metric_name == "recall"
        ][0]
        assert prec_env.lower_bound == 0.0  # Clamped, not negative

    def test_null_metrics_skip_envelope(self, db_session):
        """Test that null metric values don't create envelopes."""
        run = Run(
            sample="NULL",
            caller="GATK-HC",
            pipeline_version="v1.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
        )
        db_session.add(run)
        db_session.flush()

        sm = StratifiedMetric(
            run_id=run.id,
            dimension="variant_class",
            stratum="SNP",
            variant_type="ALL",
            precision=None,
            recall=None,
            f1=None,
            total_variants=0,
            low_confidence=True,
        )
        db_session.add(sm)
        db_session.commit()

        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="null_test",
            run_id=run.id,
        )
        db_session.commit()

        assert len(baseline.envelopes) == 0


class TestManualEnvelopes:
    def test_set_envelope_creates_new(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="manual_test",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        env = manager.set_envelope(
            baseline_id=baseline.id,
            dimension="custom",
            stratum="test",
            metric_name="precision",
            expected_value=0.99,
            lower_bound=0.95,
            upper_bound=1.0,
        )
        db_session.commit()

        assert env.manually_set is True
        assert env.expected_value == 0.99
        assert env.lower_bound == 0.95

    def test_set_envelope_updates_existing(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="update_test",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        # Override an auto-derived envelope
        env = manager.set_envelope(
            baseline_id=baseline.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            expected_value=0.995,
            lower_bound=0.990,
        )
        db_session.commit()

        assert env.manually_set is True
        assert env.expected_value == 0.995

    def test_set_envelope_on_locked_raises(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="locked_set_test",
            run_id=run_with_metrics.id,
        )
        manager.lock_baseline(baseline.id, approver="test")
        db_session.commit()

        with pytest.raises(ValueError, match="locked"):
            manager.set_envelope(
                baseline_id=baseline.id,
                dimension="test",
                stratum="test",
                metric_name="precision",
                expected_value=0.99,
            )

    def test_set_envelope_not_found_raises(self, db_session):
        manager = BaselineManager(db_session)
        with pytest.raises(ValueError, match="not found"):
            manager.set_envelope(
                baseline_id=99999,
                dimension="test",
                stratum="test",
                metric_name="precision",
                expected_value=0.99,
            )


class TestBaselineLocking:
    def test_lock_baseline(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="lock_test",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        locked = manager.lock_baseline(
            baseline.id,
            approver="Dr. Smith",
            comment="Approved for production use",
        )
        db_session.commit()

        assert locked.locked is True
        assert locked.locked_at is not None
        assert locked.approver == "Dr. Smith"
        assert locked.approval_comment == "Approved for production use"
        assert locked.artifact_hash is not None
        assert len(locked.artifact_hash) == 64  # SHA-256

    def test_lock_already_locked_raises(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="double_lock",
            run_id=run_with_metrics.id,
        )
        manager.lock_baseline(baseline.id, approver="test")
        db_session.commit()

        with pytest.raises(ValueError, match="already locked"):
            manager.lock_baseline(baseline.id, approver="test2")

    def test_lock_not_found_raises(self, db_session):
        manager = BaselineManager(db_session)
        with pytest.raises(ValueError, match="not found"):
            manager.lock_baseline(99999, approver="test")

    def test_artifact_hash_deterministic(self, db_session, run_with_metrics):
        """The same baseline should always produce the same hash."""
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="hash_test",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        hash1 = manager._compute_artifact_hash(baseline)
        hash2 = manager._compute_artifact_hash(baseline)
        assert hash1 == hash2


class TestBaselineLookup:
    def test_get_baseline(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        manager.create_baseline(name="lookup_test", run_id=run_with_metrics.id)
        db_session.commit()

        found = manager.get_baseline("lookup_test")
        assert found is not None
        assert found.name == "lookup_test"

    def test_get_baseline_not_found(self, db_session):
        manager = BaselineManager(db_session)
        assert manager.get_baseline("nonexistent") is None


class TestEnvelopeChecking:
    def test_check_no_violations(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="check_pass",
            run_id=run_with_metrics.id,
            tolerance=0.01,
        )
        db_session.commit()

        # Pass current metrics that are within envelope
        current = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "precision": 0.9970,
                "recall": 0.9955,
                "f1": 0.9963,
            }
        ]
        violations = manager.check_envelopes(baseline.id, current)
        assert len(violations) == 0

    def test_check_below_lower_bound(self, db_session, run_with_metrics):
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="check_fail",
            run_id=run_with_metrics.id,
            tolerance=0.005,
        )
        db_session.commit()

        current = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "precision": 0.98,  # Well below expected 0.9975
                "recall": 0.9955,
                "f1": 0.99,
            }
        ]
        violations = manager.check_envelopes(baseline.id, current)

        # precision should be violated
        prec_violations = [v for v in violations if v["metric_name"] == "precision"]
        assert len(prec_violations) == 1
        assert prec_violations[0]["violation"] == "below_lower_bound"
        assert prec_violations[0]["current_value"] == 0.98

    def test_check_above_upper_bound(self, db_session, run_with_metrics):
        """Not common for precision/recall, but tests the logic."""
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="check_above",
            run_id=run_with_metrics.id,
            tolerance=0.001,  # Very tight
        )
        db_session.commit()

        # Manually set a very tight upper bound
        manager.set_envelope(
            baseline_id=baseline.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            expected_value=0.99,
            lower_bound=0.98,
            upper_bound=0.995,
        )
        db_session.commit()

        current = [
            {
                "dimension": "variant_class",
                "stratum": "SNP",
                "precision": 0.999,
            }
        ]
        violations = manager.check_envelopes(baseline.id, current)
        assert any(v["violation"] == "above_upper_bound" for v in violations)

    def test_check_unmatched_stratum_ignored(self, db_session, run_with_metrics):
        """Current metrics for strata not in the baseline are ignored."""
        manager = BaselineManager(db_session)
        baseline = manager.create_baseline(
            name="check_unmatched",
            run_id=run_with_metrics.id,
        )
        db_session.commit()

        current = [
            {
                "dimension": "unknown_dim",
                "stratum": "unknown",
                "precision": 0.5,
            }
        ]
        violations = manager.check_envelopes(baseline.id, current)
        assert len(violations) == 0
