from sqlalchemy import Float, ForeignKey, Integer, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class Metric(Base):
    __tablename__ = "metrics"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        ForeignKey("runs.id"), nullable=False, index=True
    )
    variant_type: Mapped[str] = mapped_column(String(20), nullable=False)
    stratification: Mapped[str | None] = mapped_column(String(100), nullable=True)
    precision: Mapped[float] = mapped_column(Float, nullable=False)
    recall: Mapped[float] = mapped_column(Float, nullable=False)
    f1: Mapped[float] = mapped_column(Float, nullable=False)
    tp_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    fp_count: Mapped[int | None] = mapped_column(Integer, nullable=True)
    fn_count: Mapped[int | None] = mapped_column(Integer, nullable=True)

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

    run: Mapped["Run"] = relationship(back_populates="metrics")

    def __repr__(self) -> str:
        return f"<Metric(id={self.id}, type={self.variant_type}, F1={self.f1})>"
