"""Bayesian Risk Assessment models (Spec S17).

Stores per-assessment and per-stratum posterior risk estimates from
the Bayesian Risk Assessment Module (BRAM).
"""

from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import Boolean, DateTime, Float, ForeignKey, Index, Integer, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class BayesianRiskAssessment(Base):
    """Result of a BRAM assessment for a verification run.

    Links a run to its baseline and stores aggregate risk scores
    and the BRAM verdict (S17.7.1).
    """

    __tablename__ = "bayesian_risk_assessments"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=False)
    baseline_id: Mapped[int] = mapped_column(ForeignKey("baselines.id"), nullable=False)

    # Aggregate scores
    aggregate_risk_score: Mapped[float] = mapped_column(Float, nullable=False)
    mean_risk_score: Mapped[float] = mapped_column(Float, nullable=False)
    flagged_stratum_count: Mapped[int] = mapped_column(Integer, nullable=False, default=0)

    # Configuration used
    alert_threshold: Mapped[float] = mapped_column(Float, nullable=False)
    degradation_threshold: Mapped[float] = mapped_column(Float, nullable=False)
    aggregation_method: Mapped[str] = mapped_column(String(20), nullable=False, default="max")

    # Verdict: BRAM_PASS, BRAM_FLAG, BRAM_NOT_RUN
    verdict: Mapped[str] = mapped_column(String(20), nullable=False)

    created_at: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(UTC), nullable=False
    )

    # Relationships
    stratum_posteriors: Mapped[list["StratumPosterior"]] = relationship(
        back_populates="assessment", cascade="all, delete-orphan"
    )

    __table_args__ = (
        Index("ix_bram_run_id", "run_id"),
        Index("ix_bram_baseline_id", "baseline_id"),
    )

    def __repr__(self) -> str:
        return (
            f"<BayesianRiskAssessment(id={self.id}, verdict={self.verdict}, "
            f"aggregate={self.aggregate_risk_score:.3f})>"
        )


class StratumPosterior(Base):
    """Per-stratum posterior risk estimate from BRAM (S17.7.2).

    Stores the full posterior computation for a single metric
    in a single stratum, including prior, likelihood, and posterior
    parameters plus the weighted risk quantity.
    """

    __tablename__ = "stratum_posteriors"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    assessment_id: Mapped[int] = mapped_column(
        ForeignKey("bayesian_risk_assessments.id"), nullable=False
    )

    # Stratum identity
    dimension: Mapped[str] = mapped_column(String(50), nullable=False)
    stratum: Mapped[str] = mapped_column(String(100), nullable=False)
    metric_name: Mapped[str] = mapped_column(String(50), nullable=False)

    # Observed data
    delta_observed: Mapped[float] = mapped_column(Float, nullable=False)

    # Prior parameters
    prior_mu: Mapped[float] = mapped_column(Float, nullable=False)
    prior_sigma: Mapped[float] = mapped_column(Float, nullable=False)

    # Observation noise
    sigma_obs: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Posterior parameters
    posterior_mu: Mapped[float] = mapped_column(Float, nullable=False)
    posterior_sigma: Mapped[float] = mapped_column(Float, nullable=False)

    # Risk quantities
    tail_probability: Mapped[float] = mapped_column(Float, nullable=False)
    risk_weight: Mapped[float] = mapped_column(Float, nullable=False)
    weighted_risk: Mapped[float] = mapped_column(Float, nullable=False)
    flagged: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    # Tier
    tier: Mapped[str] = mapped_column(String(10), nullable=False, default="tier_3")

    assessment: Mapped["BayesianRiskAssessment"] = relationship(
        back_populates="stratum_posteriors"
    )

    __table_args__ = (
        Index("ix_stratum_post_assessment_id", "assessment_id"),
        Index("ix_stratum_post_dimension_stratum", "dimension", "stratum"),
    )

    def __repr__(self) -> str:
        return (
            f"<StratumPosterior({self.dimension}/{self.stratum}:{self.metric_name} "
            f"risk={self.weighted_risk:.3f}, flagged={self.flagged})>"
        )
