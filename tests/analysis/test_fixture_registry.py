"""Tests for fixture registry (Spec S12.4)."""

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
from analysis.fixture_registry import FixtureRegistry, compute_file_checksum  # noqa: E402


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


class TestComputeFileChecksum:
    def test_sha256(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        checksum = compute_file_checksum(str(f), "sha256")
        assert len(checksum) == 64  # SHA-256 hex digest length
        assert checksum.isalnum()

    def test_md5(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("hello world")
        checksum = compute_file_checksum(str(f), "md5")
        assert len(checksum) == 32  # MD5 hex digest length

    def test_deterministic(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("reproducible")
        c1 = compute_file_checksum(str(f))
        c2 = compute_file_checksum(str(f))
        assert c1 == c2

    def test_different_content_different_hash(self, tmp_path):
        f1 = tmp_path / "a.txt"
        f2 = tmp_path / "b.txt"
        f1.write_text("content A")
        f2.write_text("content B")
        assert compute_file_checksum(str(f1)) != compute_file_checksum(str(f2))

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            compute_file_checksum(str(tmp_path / "missing.txt"))

    def test_unsupported_algorithm(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("test")
        with pytest.raises(ValueError, match="Unsupported"):
            compute_file_checksum(str(f), "sha512")


class TestFixtureRegistration:
    def test_register_fixture(self, db_session, tmp_path):
        f = tmp_path / "truth.vcf.gz"
        f.write_bytes(b"fake VCF content")

        registry = FixtureRegistry(db_session)
        record = registry.register(
            file_path=str(f),
            source="GIAB",
            genome_build="GRCh38",
            sample_identifier="NA12878",
            truth_set_version="4.2.1",
        )
        db_session.commit()

        assert record.id is not None
        assert record.file_path == str(f)
        assert record.source == "GIAB"
        assert record.genome_build == "GRCh38"
        assert record.sample_identifier == "NA12878"
        assert record.truth_set_version == "4.2.1"
        assert record.checksum_type == "sha256"
        assert len(record.checksum) == 64
        assert record.date_registered is not None

    def test_register_with_md5(self, db_session, tmp_path):
        f = tmp_path / "ref.fa"
        f.write_bytes(b"ACGTACGT")

        registry = FixtureRegistry(db_session)
        record = registry.register(str(f), algorithm="md5")
        db_session.commit()

        assert record.checksum_type == "md5"
        assert len(record.checksum) == 32

    def test_register_duplicate_same_content(self, db_session, tmp_path):
        f = tmp_path / "test.bed"
        f.write_text("chr1\t100\t200")

        registry = FixtureRegistry(db_session)
        r1 = registry.register(str(f))
        db_session.commit()

        r2 = registry.register(str(f))
        db_session.commit()

        assert r1.id == r2.id  # Same record returned

    def test_register_duplicate_changed_content(self, db_session, tmp_path):
        f = tmp_path / "test.bed"
        f.write_text("chr1\t100\t200")

        registry = FixtureRegistry(db_session)
        r1 = registry.register(str(f))
        db_session.commit()
        old_checksum = r1.checksum

        # Modify file
        f.write_text("chr1\t100\t300")
        r2 = registry.register(str(f))
        db_session.commit()

        assert r2.id == r1.id
        assert r2.checksum != old_checksum

    def test_register_missing_file_raises(self, db_session, tmp_path):
        registry = FixtureRegistry(db_session)
        with pytest.raises(FileNotFoundError):
            registry.register(str(tmp_path / "missing.vcf"))


class TestFixtureVerification:
    def test_verify_passes(self, db_session, tmp_path):
        f = tmp_path / "truth.vcf.gz"
        f.write_bytes(b"consistent content")

        registry = FixtureRegistry(db_session)
        registry.register(str(f))
        db_session.commit()

        assert registry.verify(str(f)) is True

    def test_verify_fails_after_modification(self, db_session, tmp_path):
        f = tmp_path / "truth.vcf.gz"
        f.write_bytes(b"original content")

        registry = FixtureRegistry(db_session)
        registry.register(str(f))
        db_session.commit()

        # Modify the file
        f.write_bytes(b"tampered content")
        assert registry.verify(str(f)) is False

    def test_verify_unregistered_raises(self, db_session, tmp_path):
        f = tmp_path / "unknown.vcf"
        f.write_text("test")

        registry = FixtureRegistry(db_session)
        with pytest.raises(ValueError, match="not registered"):
            registry.verify(str(f))

    def test_verify_missing_file_raises(self, db_session, tmp_path):
        f = tmp_path / "truth.vcf.gz"
        f.write_bytes(b"will be deleted")

        registry = FixtureRegistry(db_session)
        registry.register(str(f))
        db_session.commit()

        # Delete the file
        f.unlink()
        with pytest.raises(FileNotFoundError):
            registry.verify(str(f))


class TestVerifyAll:
    def test_verify_all_passes(self, db_session, tmp_path):
        f1 = tmp_path / "a.vcf"
        f2 = tmp_path / "b.bed"
        f1.write_text("content a")
        f2.write_text("content b")

        registry = FixtureRegistry(db_session)
        registry.register(str(f1))
        registry.register(str(f2))
        db_session.commit()

        results = registry.verify_all()
        assert len(results) == 2
        assert all(r["passed"] for r in results)

    def test_verify_all_mixed(self, db_session, tmp_path):
        f1 = tmp_path / "good.vcf"
        f2 = tmp_path / "bad.vcf"
        f1.write_text("good")
        f2.write_text("original")

        registry = FixtureRegistry(db_session)
        registry.register(str(f1))
        registry.register(str(f2))
        db_session.commit()

        # Tamper with one file
        f2.write_text("tampered")

        results = registry.verify_all()
        assert len(results) == 2
        good = [r for r in results if r["file_path"] == str(f1)][0]
        bad = [r for r in results if r["file_path"] == str(f2)][0]
        assert good["passed"] is True
        assert bad["passed"] is False

    def test_verify_all_missing_file(self, db_session, tmp_path):
        f = tmp_path / "missing.vcf"
        f.write_text("will delete")

        registry = FixtureRegistry(db_session)
        registry.register(str(f))
        db_session.commit()

        f.unlink()

        results = registry.verify_all()
        assert len(results) == 1
        assert results[0]["passed"] is False
        assert results[0]["actual_checksum"] is None

    def test_verify_all_empty(self, db_session):
        registry = FixtureRegistry(db_session)
        results = registry.verify_all()
        assert results == []


class TestGetRecord:
    def test_get_existing_record(self, db_session, tmp_path):
        f = tmp_path / "test.vcf"
        f.write_text("content")

        registry = FixtureRegistry(db_session)
        registry.register(str(f), source="GIAB")
        db_session.commit()

        record = registry.get_record(str(f))
        assert record is not None
        assert record.source == "GIAB"

    def test_get_nonexistent_record(self, db_session):
        registry = FixtureRegistry(db_session)
        assert registry.get_record("/nonexistent") is None
