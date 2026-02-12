"""Fixture registry for test fixture management (Spec S12.4).

Manages checksummed registration and verification of validation fixtures
(truth VCFs, references, BED files, gene panels).
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

from sqlalchemy.orm import Session

log = logging.getLogger(__name__)

CHUNK_SIZE = 8192


def compute_file_checksum(file_path: str, algorithm: str = "sha256") -> str:
    """Compute checksum of a file.

    Args:
        file_path: Path to the file
        algorithm: Hash algorithm (sha256 or md5)

    Returns:
        Hex digest string

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If algorithm is unsupported
    """
    if algorithm not in ("sha256", "md5"):
        raise ValueError(f"Unsupported checksum algorithm: {algorithm}")

    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        while True:
            chunk = f.read(CHUNK_SIZE)
            if not chunk:
                break
            h.update(chunk)

    return h.hexdigest()


class FixtureRegistry:
    """Registry for managing validation fixture files.

    Args:
        session: SQLAlchemy database session
    """

    def __init__(self, session: Session):
        self.session = session

    def register(
        self,
        file_path: str,
        source: str | None = None,
        genome_build: str | None = None,
        sample_identifier: str | None = None,
        truth_set_version: str | None = None,
        algorithm: str = "sha256",
    ) -> "FixtureRecord":
        """Register a fixture file with its checksum and metadata.

        Args:
            file_path: Path to the fixture file
            source: Source of the fixture (e.g. GIAB, Horizon, in-house)
            genome_build: Reference genome build (e.g. GRCh38)
            sample_identifier: Sample name/ID
            truth_set_version: Version of the truth set
            algorithm: Checksum algorithm (sha256 or md5)

        Returns:
            Created FixtureRecord

        Raises:
            FileNotFoundError: If file does not exist
        """
        from models.baseline import FixtureRecord

        checksum = compute_file_checksum(file_path, algorithm)

        # Check if already registered with same path
        existing = (
            self.session.query(FixtureRecord)
            .filter(FixtureRecord.file_path == file_path)
            .first()
        )
        if existing:
            if existing.checksum == checksum:
                log.info("Fixture already registered: %s", file_path)
                return existing
            # File changed - update checksum
            existing.checksum = checksum
            existing.checksum_type = algorithm
            existing.source = source or existing.source
            existing.genome_build = genome_build or existing.genome_build
            existing.sample_identifier = sample_identifier or existing.sample_identifier
            existing.truth_set_version = truth_set_version or existing.truth_set_version
            log.warning(
                "Fixture checksum changed, updated: %s (new=%s)",
                file_path,
                checksum[:12],
            )
            return existing

        record = FixtureRecord(
            file_path=file_path,
            checksum_type=algorithm,
            checksum=checksum,
            source=source,
            genome_build=genome_build,
            sample_identifier=sample_identifier,
            truth_set_version=truth_set_version,
        )
        self.session.add(record)
        log.info(
            "Registered fixture: %s (%s=%s)",
            file_path,
            algorithm,
            checksum[:12],
        )
        return record

    def verify(self, file_path: str) -> bool:
        """Verify a fixture file's integrity against its registered checksum.

        Args:
            file_path: Path to the fixture file

        Returns:
            True if checksum matches, False otherwise

        Raises:
            ValueError: If file is not registered
            FileNotFoundError: If file does not exist
        """
        from models.baseline import FixtureRecord

        record = (
            self.session.query(FixtureRecord)
            .filter(FixtureRecord.file_path == file_path)
            .first()
        )
        if not record:
            raise ValueError(f"Fixture not registered: {file_path}")

        current_checksum = compute_file_checksum(file_path, record.checksum_type)
        matches = current_checksum == record.checksum

        if not matches:
            log.warning(
                "Fixture integrity check FAILED: %s (expected=%s, got=%s)",
                file_path,
                record.checksum[:12],
                current_checksum[:12],
            )
        else:
            log.info("Fixture integrity check passed: %s", file_path)

        return matches

    def verify_all(self) -> list[dict]:
        """Verify all registered fixtures.

        Returns:
            List of dicts with file_path, expected_checksum, actual_checksum, passed
        """
        from models.baseline import FixtureRecord

        records = self.session.query(FixtureRecord).all()
        results = []

        for record in records:
            try:
                current = compute_file_checksum(record.file_path, record.checksum_type)
                passed = current == record.checksum
            except FileNotFoundError:
                current = None
                passed = False

            results.append(
                {
                    "file_path": record.file_path,
                    "expected_checksum": record.checksum,
                    "actual_checksum": current,
                    "passed": passed,
                }
            )

        return results

    def get_record(self, file_path: str) -> "FixtureRecord | None":
        """Look up a fixture record by file path.

        Args:
            file_path: Path to the fixture file

        Returns:
            FixtureRecord or None
        """
        from models.baseline import FixtureRecord

        return (
            self.session.query(FixtureRecord)
            .filter(FixtureRecord.file_path == file_path)
            .first()
        )
