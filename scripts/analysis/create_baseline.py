"""Snakemake script for creating a validation baseline.

Creates a baseline from the most recent run, auto-deriving performance
envelopes from stratified metrics.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Add scripts directory to path for imports
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)


def create_baseline_from_run(
    database: str,
    baseline_name: str,
    tolerance: float = 0.005,
    mode: str = "germline",
    pipeline_version: str | None = None,
) -> None:
    """Create a baseline from the most recent run in the database.

    Args:
        database: Path to SQLite database
        baseline_name: Name for the new baseline
        tolerance: Tolerance for auto-derived envelopes
        mode: Execution mode
        pipeline_version: Pipeline version string
    """
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker

    from models.base import Base
    from models.run import Run
    from analysis.baseline_manager import BaselineManager

    engine = create_engine(f"sqlite:///{database}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        # Get the most recent run
        run = session.query(Run).order_by(Run.id.desc()).first()
        if not run:
            raise ValueError("No runs found in database")

        manager = BaselineManager(session)
        baseline = manager.create_baseline(
            name=baseline_name,
            run_id=run.id,
            pipeline_version=pipeline_version or run.pipeline_version,
            mode=mode,
            tolerance=tolerance,
        )

        session.commit()
        logging.info(
            "Baseline '%s' created successfully (id=%d, envelopes=%d)",
            baseline.name,
            baseline.id,
            len(baseline.envelopes),
        )
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


try:
    from typing import TYPE_CHECKING

    if TYPE_CHECKING:
        from snakemake.script import Snakemake

        snakemake: Snakemake
    else:
        snakemake = snakemake  # type: ignore  # noqa: F821

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    params = snakemake.params
    create_baseline_from_run(
        database=params.database,
        baseline_name=params.baseline_name,
        tolerance=float(params.tolerance),
        mode=params.execution_mode,
        pipeline_version=params.pipeline_version,
    )

    # Touch output
    Path(snakemake.output.done).touch()

except NameError:
    pass  # Not running via Snakemake
