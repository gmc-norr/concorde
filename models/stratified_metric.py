from sqlalchemy import Boolean, Float, ForeignKey, Index, Integer, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class StratifiedMetric(Base):
    __tablename__ = "stratified_metrics"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=False)

    # Stratification dimension and value
    dimension: Mapped[str] = mapped_column(String(50), nullable=False)
    stratum: Mapped[str] = mapped_column(String(100), nullable=False)

    # Variant type (ALL, SNP, INDEL, etc.)
    variant_type: Mapped[str] = mapped_column(String(20), nullable=False, default="ALL")

    # Metrics
    precision: Mapped[float] = mapped_column(Float, nullable=True)
    recall: Mapped[float] = mapped_column(Float, nullable=True)
    f1: Mapped[float] = mapped_column(Float, nullable=True)
    tp_count: Mapped[int] = mapped_column(Integer, nullable=True)
    fp_count: Mapped[int] = mapped_column(Integer, nullable=True)
    fn_count: Mapped[int] = mapped_column(Integer, nullable=True)

    # Total variants in this stratum
    total_variants: Mapped[int] = mapped_column(Integer, nullable=True)

    # Confidence intervals (Wilson score)
    precision_ci_lower: Mapped[float | None] = mapped_column(Float, nullable=True)
    precision_ci_upper: Mapped[float | None] = mapped_column(Float, nullable=True)
    recall_ci_lower: Mapped[float | None] = mapped_column(Float, nullable=True)
    recall_ci_upper: Mapped[float | None] = mapped_column(Float, nullable=True)
    f1_ci_lower: Mapped[float | None] = mapped_column(Float, nullable=True)
    f1_ci_upper: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Matthews Correlation Coefficient
    mcc: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Genotype concordance rate
    genotype_concordance: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Low-confidence flag (when total_variants < min_variants_per_stratum)
    low_confidence: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

    run: Mapped["Run"] = relationship(back_populates="stratified_metrics")

    __table_args__ = (
        Index("ix_strat_metrics_run_id", "run_id"),
        Index("ix_strat_metrics_dimension", "dimension"),
        Index("ix_strat_metrics_stratum", "stratum"),
    )

    def __repr__(self) -> str:
        return f"<StratifiedMetric(id={self.id}, {self.dimension}={self.stratum}, F1={self.f1})>"
