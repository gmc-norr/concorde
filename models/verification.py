"""Verification models for continuous verification (Spec S13).

Stores verification run results, variant-level diffs, and per-stratum
comparisons against a locked baseline.
"""

from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import DateTime, ForeignKey, Index, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class VerificationResult(Base):
    """Result of a verification run compared against a baseline.

    Links a new run to a baseline and stores the overall verdict.
    """

    __tablename__ = "verification_results"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=False)
    baseline_id: Mapped[int] = mapped_column(ForeignKey("baselines.id"), nullable=False)

    # Overall verdict: PASS, FAIL, REVIEW_REQUIRED
    verdict: Mapped[str] = mapped_column(String(20), nullable=False)
    summary: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Counts
    total_transitions: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    equivalent_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    drift_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    biological_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    envelope_breaches: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    created_at: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(UTC), nullable=False
    )

    transitions: Mapped[list["VariantTransition"]] = relationship(
        back_populates="verification_result", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_verification_run_id", "run_id"),
        Index("ix_verification_baseline_id", "baseline_id"),
    )

    def __repr__(self) -> str:
        return (
            f"<VerificationResult(id={self.id}, verdict={self.verdict}, "
            f"transitions={self.total_transitions})>"
        )


class VariantTransition(Base):
    """A single variant-level change between baseline and verification run.

    Records what happened to a specific variant (gained, lost, changed).
    """

    __tablename__ = "variant_transitions"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    verification_id: Mapped[int] = mapped_column(
        ForeignKey("verification_results.id"), nullable=False
    )

    # Variant identity
    chrom: Mapped[str] = mapped_column(String(10), nullable=False)
    pos: Mapped[int] = mapped_column(Integer, nullable=False)
    ref: Mapped[str] = mapped_column(String(500), nullable=False)
    alt: Mapped[str] = mapped_column(String(500), nullable=False)

    # Classifications
    baseline_classification: Mapped[str | None] = mapped_column(String(20), nullable=True)
    verification_classification: Mapped[str | None] = mapped_column(String(20), nullable=True)
    transition: Mapped[str] = mapped_column(String(20), nullable=False)

    # Classification of the difference (S13.3)
    diff_class: Mapped[str] = mapped_column(String(20), nullable=False, default="equivalent")

    # Affected strata (JSON serialized list)
    strata_json: Mapped[str | None] = mapped_column(Text, nullable=True)

    verification_result: Mapped["VerificationResult"] = relationship(
        back_populates="transitions"
    )
    evidence: Mapped[list["RootCauseEvidence"]] = relationship(
        back_populates="transition", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_transition_verification_id", "verification_id"),
        Index("ix_transition_chrom_pos", "chrom", "pos"),
    )

    def __repr__(self) -> str:
        return (
            f"<VariantTransition({self.chrom}:{self.pos}:{self.ref}>{self.alt} "
            f"{self.transition}, class={self.diff_class})>"
        )
