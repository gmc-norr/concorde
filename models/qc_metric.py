from __future__ import annotations

from sqlalchemy import Float, ForeignKey, Index, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class QCMetric(Base):
    """Flexible storage for QC metrics from various tools.

    Stores QC metrics in key-value format, allowing different tools to store
    different metrics without schema changes. Common metrics include coverage
    statistics, alignment metrics, contamination estimates, etc.
    """

    __tablename__ = "qc_metrics"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        ForeignKey("runs.id"), nullable=False
    )
    metric_name: Mapped[str] = mapped_column(String(100), nullable=False)
    metric_value_float: Mapped[float | None] = mapped_column(Float, nullable=True)
    metric_value_str: Mapped[str | None] = mapped_column(String(500), nullable=True)
    metric_source: Mapped[str] = mapped_column(String(50), nullable=False)
    metric_category: Mapped[str | None] = mapped_column(String(50), nullable=True)
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)

    run: Mapped["Run"] = relationship(back_populates="qc_metrics")

    __table_args__ = (
        Index("ix_qc_metrics_run_id", "run_id"),
        Index("ix_qc_metrics_name", "metric_name"),
        Index("ix_qc_metrics_category", "metric_category"),
        Index("ix_qc_metrics_source", "metric_source"),
    )

    def __repr__(self) -> str:
        return f"<QCMetric(id={self.id}, name={self.metric_name}, source={self.metric_source})>"
