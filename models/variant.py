from sqlalchemy import Boolean, Float, ForeignKey, Index, Integer, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class Variant(Base):
    __tablename__ = "variants"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(ForeignKey("runs.id"), nullable=False)
    chrom: Mapped[str] = mapped_column(String(10), nullable=False)
    pos: Mapped[int] = mapped_column(Integer, nullable=False)
    ref: Mapped[str] = mapped_column(String(500), nullable=False)
    alt: Mapped[str] = mapped_column(String(500), nullable=False)
    type: Mapped[str] = mapped_column(String(20), nullable=False)
    classification: Mapped[str] = mapped_column(String(10), nullable=False)
    dp: Mapped[int | None] = mapped_column(Integer, nullable=True)
    mq: Mapped[float | None] = mapped_column(Float, nullable=True)
    af: Mapped[float | None] = mapped_column(Float, nullable=True)
    gq: Mapped[int | None] = mapped_column(Integer, nullable=True)
    qual: Mapped[float | None] = mapped_column(Float, nullable=True)
    filter_status: Mapped[str | None] = mapped_column(String(100), nullable=True)
    gene: Mapped[str | None] = mapped_column(String(100), nullable=True)
    impact: Mapped[str | None] = mapped_column(String(50), nullable=True)

    # Somatic-specific fields
    tumor_dp: Mapped[int | None] = mapped_column(Integer, nullable=True)
    tumor_af: Mapped[float | None] = mapped_column(Float, nullable=True)
    normal_dp: Mapped[int | None] = mapped_column(Integer, nullable=True)
    normal_af: Mapped[float | None] = mapped_column(Float, nullable=True)
    somatic_quality: Mapped[float | None] = mapped_column(Float, nullable=True)
    caller_filter: Mapped[str | None] = mapped_column(String(100), nullable=True)

    # Genotype concordance fields
    truth_gt: Mapped[str | None] = mapped_column(String(10), nullable=True)
    query_gt: Mapped[str | None] = mapped_column(String(10), nullable=True)
    gt_concordant: Mapped[bool | None] = mapped_column(Boolean, nullable=True)

    # Stratification support fields
    zygosity: Mapped[str | None] = mapped_column(String(10), nullable=True)
    indel_size: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # Region annotation fields
    gc_content: Mapped[float | None] = mapped_column(Float, nullable=True)
    in_low_complexity: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    in_segdup: Mapped[bool | None] = mapped_column(Boolean, nullable=True)
    mappability_score: Mapped[float | None] = mapped_column(Float, nullable=True)

    # VEP annotation fields
    consequence: Mapped[str | None] = mapped_column(String(100), nullable=True)
    gene_symbol: Mapped[str | None] = mapped_column(String(50), nullable=True)
    gene_id: Mapped[str | None] = mapped_column(String(30), nullable=True)
    transcript_id: Mapped[str | None] = mapped_column(String(30), nullable=True)
    hgvsc: Mapped[str | None] = mapped_column(String(200), nullable=True)
    hgvsp: Mapped[str | None] = mapped_column(String(200), nullable=True)
    exon: Mapped[str | None] = mapped_column(String(20), nullable=True)
    protein_position: Mapped[int | None] = mapped_column(Integer, nullable=True)
    sift_prediction: Mapped[str | None] = mapped_column(String(50), nullable=True)
    polyphen_prediction: Mapped[str | None] = mapped_column(String(50), nullable=True)

    run: Mapped["Run"] = relationship(back_populates="variants")
    gene_sets: Mapped[list["GeneSet"]] = relationship(
        secondary="variant_gene_sets", back_populates="variants"
    )

    __table_args__ = (
        Index("ix_variants_chrom_pos", "chrom", "pos"),
        Index("ix_variants_run_id", "run_id"),
        Index("ix_variants_classification", "classification"),
        Index("ix_variants_tumor_af", "tumor_af"),
        Index("ix_variants_consequence", "consequence"),
        Index("ix_variants_impact", "impact"),
    )

    def __repr__(self) -> str:
        return f"<Variant(id={self.id}, {self.chrom}:{self.pos} {self.ref}>{self.alt}, {self.classification})>"
