from __future__ import annotations

from datetime import UTC, datetime

from sqlalchemy import DateTime, Float, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class Run(Base):
    __tablename__ = "runs"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    sample: Mapped[str] = mapped_column(String(255), nullable=False)
    pipeline_version: Mapped[str] = mapped_column(String(100), nullable=False)
    caller: Mapped[str] = mapped_column(String(100), nullable=False)
    parameters: Mapped[str | None] = mapped_column(Text, nullable=True)
    decomposition_mode: Mapped[str | None] = mapped_column(String(50), nullable=True)
    comparison_tool: Mapped[str | None] = mapped_column(String(50), nullable=True)
    mode: Mapped[str] = mapped_column(String(20), nullable=False, default="germline")
    date: Mapped[datetime] = mapped_column(
        DateTime, default=lambda: datetime.now(UTC), nullable=False
    )

    # Somatic mode fields
    tumor_sample: Mapped[str | None] = mapped_column(String(255), nullable=True)
    normal_sample: Mapped[str | None] = mapped_column(String(255), nullable=True)
    calling_mode: Mapped[str | None] = mapped_column(String(20), nullable=True)

    # QC Summary Fields (from nf-core/raredisease or other pipelines)
    mean_coverage: Mapped[float | None] = mapped_column(Float, nullable=True)
    median_coverage: Mapped[float | None] = mapped_column(Float, nullable=True)
    pct_bases_10x: Mapped[float | None] = mapped_column(Float, nullable=True)
    pct_bases_20x: Mapped[float | None] = mapped_column(Float, nullable=True)
    pct_bases_30x: Mapped[float | None] = mapped_column(Float, nullable=True)
    mapping_rate: Mapped[float | None] = mapped_column(Float, nullable=True)
    duplicate_rate: Mapped[float | None] = mapped_column(Float, nullable=True)
    median_insert_size: Mapped[float | None] = mapped_column(Float, nullable=True)
    contamination_pct: Mapped[float | None] = mapped_column(Float, nullable=True)
    sex_inferred: Mapped[str | None] = mapped_column(String(10), nullable=True)

    # Pipeline metadata (from nf-core/raredisease pipeline_info)
    nfcore_pipeline_version: Mapped[str | None] = mapped_column(
        String(50), nullable=True
    )
    pipeline_params_json: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Ensemble caller support
    callers_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    ensemble_method: Mapped[str | None] = mapped_column(String(50), nullable=True)

    variants: Mapped[list["Variant"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )
    metrics: Mapped[list["Metric"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )
    qc_metrics: Mapped[list["QCMetric"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )
    software_versions: Mapped[list["SoftwareVersion"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )
    stratified_metrics: Mapped[list["StratifiedMetric"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )
    baselines: Mapped[list["Baseline"]] = relationship(
        back_populates="run", cascade="all, delete-orphan"
    )

    def __repr__(self) -> str:
        return f"<Run(id={self.id}, sample={self.sample}, caller={self.caller})>"
