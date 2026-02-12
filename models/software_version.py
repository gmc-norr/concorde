from __future__ import annotations

from sqlalchemy import ForeignKey, Index, String, Text
from sqlalchemy.orm import Mapped, mapped_column, relationship

from models.base import Base


class SoftwareVersion(Base):
    """Tracks software tool versions used in pipeline runs.

    Stores version information for each tool used in the pipeline execution.
    Enables tracking version changes across runs and identifying version-specific
    issues in variant calling results.
    """

    __tablename__ = "software_versions"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    run_id: Mapped[int] = mapped_column(
        ForeignKey("runs.id"), nullable=False
    )
    tool_name: Mapped[str] = mapped_column(String(100), nullable=False)
    version: Mapped[str] = mapped_column(String(100), nullable=False)
    category: Mapped[str | None] = mapped_column(String(50), nullable=True)
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)

    run: Mapped["Run"] = relationship(back_populates="software_versions")

    __table_args__ = (
        Index("ix_software_versions_run_id", "run_id"),
        Index("ix_software_versions_tool_name", "tool_name"),
        Index("ix_software_versions_version", "version"),
    )

    def __repr__(self) -> str:
        return f"<SoftwareVersion(id={self.id}, tool={self.tool_name}, version={self.version})>"
