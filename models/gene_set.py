from sqlalchemy import Column, ForeignKey, Index, Integer, String, Table, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base

variant_gene_sets = Table(
    "variant_gene_sets",
    Base.metadata,
    Column(
        "variant_id",
        Integer,
        ForeignKey("variants.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    Column(
        "gene_set_id",
        Integer,
        ForeignKey("gene_sets.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    Index("ix_vgs_gene_set_id", "gene_set_id"),
)


class GeneSet(Base):
    __tablename__ = "gene_sets"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), unique=True, nullable=False)
    description: Mapped[str | None] = mapped_column(Text, nullable=True)
    version: Mapped[str | None] = mapped_column(String(50), nullable=True)

    variants: Mapped[list["Variant"]] = relationship(
        secondary=variant_gene_sets, back_populates="gene_sets"
    )

    def __repr__(self) -> str:
        return f"<GeneSet(id={self.id}, name={self.name})>"
