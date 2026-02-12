"""Root cause evidence model (Spec S14).

Stores evidence items linked to variant transitions, with scoring
for explanatory power.
"""

from __future__ import annotations

from sqlalchemy import Float, ForeignKey, Index, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class RootCauseEvidence(Base):
    """A single piece of root cause evidence for a variant transition.

    Each evidence item has a category (coverage, mapping_quality, filter,
    quality_score, tool_version), values, and an explanatory score.
    """

    __tablename__ = "root_cause_evidence"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    transition_id: Mapped[int] = mapped_column(
        ForeignKey("variant_transitions.id"), nullable=False
    )

    # Evidence category
    category: Mapped[str] = mapped_column(String(50), nullable=False)

    # Evidence values (JSON or key-value pairs)
    baseline_value: Mapped[str | None] = mapped_column(String(100), nullable=True)
    verification_value: Mapped[str | None] = mapped_column(String(100), nullable=True)
    change_description: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Numeric change for thresholding
    change_magnitude: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Score: STRONG, MODERATE, WEAK, NONE
    score: Mapped[str] = mapped_column(String(10), nullable=False, default="NONE")

    transition: Mapped["VariantTransition"] = relationship(back_populates="evidence")

    __table_args__ = (
        Index("ix_evidence_transition_id", "transition_id"),
        Index("ix_evidence_category", "category"),
    )

    def __repr__(self) -> str:
        return f"<RootCauseEvidence(category={self.category}, score={self.score})>"
