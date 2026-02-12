"""Ingest parsed pipeline outputs into the Concorde SQLite database.

Called via Snakemake `script:` directive.
Uses SQLAlchemy models directly to ensure schema consistency.
"""

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # snakemake is injected by Snakemake at runtime when using script: directive
    from snakemake.script import Snakemake
    snakemake: Snakemake
else:
    # Declare snakemake as a runtime global to avoid F821 linting errors
    snakemake = snakemake  # type: ignore  # noqa: F821

from sqlalchemy import create_engine  # noqa: E402
from sqlalchemy.exc import IntegrityError, OperationalError  # noqa: E402
from sqlalchemy.orm import sessionmaker  # noqa: E402

from models.base import Base  # noqa: E402

from ingest_helpers import (  # noqa: E402
    create_run_from_qc,
    insert_metrics,
    insert_qc_metrics,
    insert_software_versions,
    insert_stratified_metrics,
    load_and_insert_variants,
    load_pipeline_metadata,
    load_qc_summary,
    merge_vep_annotations,
    process_gene_sets,
)
from validators import ValidationError  # noqa: E402

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def resolve_db_path(db_path_str: str) -> Path:
    """Resolve the database path relative to the pipeline directory."""
    p = Path(db_path_str)
    if not p.is_absolute():
        # Resolve relative to the Snakefile's working directory (pipeline/)
        p = Path(__file__).resolve().parent.parent / p
    return p.resolve()


# --- Main execution ---
db_path = resolve_db_path(snakemake.params.database)
db_path.parent.mkdir(parents=True, exist_ok=True)
db_url = f"sqlite:///{db_path}"

log.info("Connecting to database: %s", db_url)

engine = create_engine(
    db_url,
    connect_args={"check_same_thread": False},
)

# Ensure tables exist
Base.metadata.create_all(bind=engine)

Session = sessionmaker(bind=engine)
session = Session()

try:
    # Load input data
    qc_summary = load_qc_summary(snakemake.input.qc_summary_tsv)
    pipeline_metadata = load_pipeline_metadata(snakemake.input.pipeline_metadata_tsv)

    # Create Run record from QC data and Snakemake params
    execution_mode = snakemake.params.get("execution_mode", "germline")
    somatic_config = snakemake.params.get("somatic_config", {})

    params = {
        "sample": snakemake.params.sample,
        "pipeline_version": snakemake.params.pipeline_version,
        "caller": snakemake.params.caller,
        "parameters": snakemake.params.parameters,
        "decomposition_mode": snakemake.params.decomposition_mode,
        "comparison_tool": snakemake.params.get("comparison_tool", "happy"),
        "mode": execution_mode,
        "tumor_sample": somatic_config.get("tumor_sample") if execution_mode == "somatic" else None,
        "normal_sample": somatic_config.get("normal_sample") if execution_mode == "somatic" else None,
        "calling_mode": somatic_config.get("calling_mode") if execution_mode == "somatic" else None,
    }
    run = create_run_from_qc(qc_summary, pipeline_metadata, params)
    session.add(run)
    session.flush()
    log.info(
        "Created Run id=%d (sample=%s, caller=%s, comparison_tool=%s)",
        run.id,
        run.sample,
        run.caller,
        run.comparison_tool,
    )

    # Insert variants and get variant objects for gene set processing
    variant_objects = load_and_insert_variants(
        snakemake.input.variants_tsv, run.id, session
    )

    # Insert metrics and related data
    insert_metrics(snakemake.input.metrics_tsv, run.id, session)
    insert_qc_metrics(snakemake.input.qc_metrics_tsv, run.id, session)
    insert_software_versions(snakemake.input.software_versions_tsv, run.id, session)

    # Merge VEP annotations onto variant objects
    merge_vep_annotations(
        snakemake.input.vep_annotations_tsv,
        variant_objects,
    )

    # Process gene set assignments
    process_gene_sets(
        snakemake.input.gene_set_tsv,
        run.id,
        variant_objects,
        snakemake.params.gene_sets_config,
        session,
    )

    # Annotate variants with region information (GC content, LCR, segdup, mappability)
    region_beds_config = snakemake.params.get("region_beds", {})
    reference_path = snakemake.params.get("reference", None)
    if region_beds_config and any(region_beds_config.values()):
        from analysis.region_annotator import RegionAnnotator  # noqa: E402

        annotator = RegionAnnotator(
            config=region_beds_config, reference_path=reference_path
        )
        annotator.annotate_variants(variant_objects)
        log.info("Applied region annotations to %d variants", len(variant_objects))

    # Compute and insert stratified metrics
    strat_config = snakemake.params.get("stratification_config", {})
    gene_set_names = [gs["name"] for gs in (snakemake.params.gene_sets_config or []) if gs.get("name")]
    insert_stratified_metrics(
        variant_objects,
        run.id,
        session,
        stratification_config=strat_config,
        mode=execution_mode,
        gene_set_names=gene_set_names,
    )

    # Compute and insert extended QC metrics (Ti/Tv and Het/Hom ratios)
    from analysis.extended_metrics import compute_het_hom_ratio, compute_titv_ratio  # noqa: E402
    from models.qc_metric import QCMetric  # noqa: E402

    titv = compute_titv_ratio(variant_objects)
    if titv["titv_ratio"] is not None:
        session.add(
            QCMetric(
                run_id=run.id,
                metric_name="titv_ratio",
                metric_value_float=titv["titv_ratio"],
                metric_source="concorde",
                metric_category="variant_quality",
            )
        )
        log.info("Inserted Ti/Tv ratio: %.4f", titv["titv_ratio"])

    hethom = compute_het_hom_ratio(variant_objects)
    if hethom["het_hom_ratio"] is not None:
        session.add(
            QCMetric(
                run_id=run.id,
                metric_name="het_hom_ratio",
                metric_value_float=hethom["het_hom_ratio"],
                metric_source="concorde",
                metric_category="variant_quality",
            )
        )
        log.info("Inserted Het/Hom ratio: %.4f", hethom["het_hom_ratio"])

    # Commit transaction
    session.commit()
    log.info("Database ingestion complete for run id=%d", run.id)

except ValidationError as e:
    session.rollback()
    log.error("Data validation failed: %s", e)
    raise
except IntegrityError as e:
    session.rollback()
    log.error("Database constraint violation: %s", e)
    raise
except OperationalError as e:
    session.rollback()
    log.error("Database operational error: %s", e)
    raise
except Exception:
    session.rollback()
    log.exception("Ingestion failed, rolling back transaction")
    raise
finally:
    session.close()

# Touch sentinel file
Path(snakemake.output.done).touch()
log.info("Wrote sentinel: %s", snakemake.output.done)
