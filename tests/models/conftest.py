"""Shared fixtures for model tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Add project root to path
PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from models.base import Base


@pytest.fixture
def test_db():
    """Create in-memory test database with schema.

    Each test gets a fresh database to avoid conflicts.
    """
    # Create unique in-memory database for this test
    engine = create_engine(
        "sqlite:///:memory:",
        connect_args={"check_same_thread": False},
        # Enable foreign key constraints in SQLite
        poolclass=None
    )

    # Enable foreign key constraints
    @sqlalchemy.event.listens_for(engine, "connect")
    def set_sqlite_pragma(dbapi_conn, connection_record):
        cursor = dbapi_conn.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()

    # Create all tables
    Base.metadata.create_all(engine)

    # Create session
    Session = sessionmaker(bind=engine)
    session = Session()

    yield session

    # Cleanup
    session.close()
    Base.metadata.drop_all(engine)
    engine.dispose()


@pytest.fixture
def sample_run(test_db):
    """Create a sample Run for testing."""
    # Import here to avoid circular dependency
    sys.path.insert(0, PROJECT_ROOT)
    from models import Run

    run = Run(
        sample="NA12878",
        caller="GATK-HC",
        pipeline_version="v4.5.0",
        decomposition_mode="decomposed",
        comparison_tool="happy",
        mean_coverage=35.2,
        median_coverage=34.8,
        pct_bases_10x=99.1,
        pct_bases_30x=95.3
    )
    test_db.add(run)
    test_db.commit()
    test_db.refresh(run)
    return run
