"""Baseline management models for initial validation.

Baselines capture a snapshot of validation performance at a point in time,
including configuration, input checksums, performance envelopes, and
approval metadata. See spec S12.
"""

from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import Boolean, DateTime, Float, ForeignKey, Index, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class Baseline(Base):
    """An immutable, versioned validation baseline.

    Contains the configuration snapshot, computed metrics for every stratum,
    and approval metadata. Once locked, a baseline cannot be modified.
    """

    __tablename__ = "baselines"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, unique=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=False)

    # Configuration snapshot
    config_snapshot: Mapped[str | None] = mapped_column(Text, nullable=True)
    stratification_snapshot: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Pipeline metadata
    pipeline_version: Mapped[str | None] = mapped_column(String(100), nullable=True)
    mode: Mapped[str] = mapped_column(String(20), nullable=False, default="germline")

    # Signing and locking (S12.3)
    locked: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    locked_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)
    approver: Mapped[str | None] = mapped_column(String(255), nullable=True)
    approval_comment: Mapped[str | None] = mapped_column(Text, nullable=True)
    artifact_hash: Mapped[str | None] = mapped_column(String(64), nullable=True)

    created_at: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(UTC), nullable=False
    )

    # Relationships
    run: Mapped["Run"] = relationship(back_populates="baselines")
    envelopes: Mapped[list["BaselineEnvelope"]] = relationship(
        back_populates="baseline", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_baselines_name", "name"),
        Index("ix_baselines_run_id", "run_id"),
    )

    def __repr__(self) -> str:
        return f"<Baseline(id={self.id}, name={self.name}, locked={self.locked})>"


class BaselineEnvelope(Base):
    """Performance envelope for a single metric in a single stratum.

    Defines expected value, lower bound, and optional upper bound (S12.2).
    """

    __tablename__ = "baseline_envelopes"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    baseline_id: Mapped[int] = mapped_column(ForeignKey("baselines.id"), nullable=False)

    # Stratum identity
    dimension: Mapped[str] = mapped_column(String(50), nullable=False)
    stratum: Mapped[str] = mapped_column(String(100), nullable=False)
    metric_name: Mapped[str] = mapped_column(String(50), nullable=False)

    # Envelope values
    expected_value: Mapped[float] = mapped_column(Float, nullable=False)
    lower_bound: Mapped[float | None] = mapped_column(Float, nullable=True)
    upper_bound: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Whether this was manually overridden vs auto-derived
    manually_set: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    baseline: Mapped["Baseline"] = relationship(back_populates="envelopes")

    __table_args__ = (
        Index("ix_envelope_baseline_id", "baseline_id"),
        Index("ix_envelope_dimension_stratum", "dimension", "stratum"),
    )

    def __repr__(self) -> str:
        return (
            f"<BaselineEnvelope({self.dimension}/{self.stratum}: "
            f"{self.metric_name}={self.expected_value})>"
        )


class FixtureRecord(Base):
    """Registry entry for a test fixture used in validation (S12.4).

    Tracks checksums and metadata for truth sets, references, and BED files.
    """

    __tablename__ = "fixture_records"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    file_path: Mapped[str] = mapped_column(String(500), nullable=False)
    checksum_type: Mapped[str] = mapped_column(String(10), nullable=False, default="sha256")
    checksum: Mapped[str] = mapped_column(String(64), nullable=False)

    # Metadata
    source: Mapped[str | None] = mapped_column(String(100), nullable=True)
    genome_build: Mapped[str | None] = mapped_column(String(20), nullable=True)
    sample_identifier: Mapped[str | None] = mapped_column(String(255), nullable=True)
    truth_set_version: Mapped[str | None] = mapped_column(String(50), nullable=True)
    date_registered: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(UTC), nullable=False
    )

    __table_args__ = (
        Index("ix_fixture_file_path", "file_path"),
        Index("ix_fixture_checksum", "checksum"),
    )

    def __repr__(self) -> str:
        return f"<FixtureRecord(path={self.file_path}, checksum={self.checksum[:12]}...)>"
