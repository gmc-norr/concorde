"""Tests for audit trail recording."""

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
from models.audit_event import AuditEvent  # noqa: E402
from models.run import Run  # noqa: E402

from analysis.audit import get_audit_trail, record_event  # noqa: E402


@pytest.fixture
def db_session(tmp_path):
    db_path = tmp_path / "audit.db"
    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    run = Run(
        sample="NA12878", caller="GATK", pipeline_version="v1",
        comparison_tool="happy", mode="germline",
    )
    session.add(run)
    session.commit()

    yield session
    session.close()


class TestRecordEvent:
    def test_basic_event(self, db_session):
        event = record_event(
            db_session,
            event_type="ingestion_complete",
            run_id=1,
            actor="test_user",
            detail="Ingested 100 variants",
        )
        db_session.flush()

        assert event.event_type == "ingestion_complete"
        assert event.run_id == 1
        assert event.actor == "test_user"
        assert event.detail == "Ingested 100 variants"
        assert event.created_at is not None

    def test_default_actor(self, db_session):
        event = record_event(db_session, event_type="test")
        db_session.flush()
        assert event.actor == "system"

    def test_baseline_event(self, db_session):
        event = record_event(
            db_session,
            event_type="baseline_locked",
            baseline_id=42,
            actor="admin",
        )
        db_session.flush()
        assert event.baseline_id == 42


class TestGetAuditTrail:
    def test_filter_by_run(self, db_session):
        record_event(db_session, "event_a", run_id=1)
        record_event(db_session, "event_b", run_id=2)
        db_session.flush()

        trail = get_audit_trail(db_session, run_id=1)
        assert len(trail) == 1
        assert trail[0].event_type == "event_a"

    def test_filter_by_event_type(self, db_session):
        record_event(db_session, "ingestion_complete", run_id=1)
        record_event(db_session, "verification_complete", run_id=1)
        db_session.flush()

        trail = get_audit_trail(db_session, event_type="ingestion_complete")
        assert len(trail) == 1

    def test_limit(self, db_session):
        for i in range(20):
            record_event(db_session, f"event_{i}")
        db_session.flush()

        trail = get_audit_trail(db_session, limit=5)
        assert len(trail) == 5

    def test_empty_trail(self, db_session):
        trail = get_audit_trail(db_session, run_id=999)
        assert len(trail) == 0
