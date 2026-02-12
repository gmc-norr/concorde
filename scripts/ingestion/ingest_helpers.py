"""Helper functions for database ingestion script.

This module contains modular functions for loading, transforming,
and inserting data into the database. Extracted from ingest.py
to improve maintainability and testability.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd
from sqlalchemy.orm import Session

from utils import safe_float, safe_int, safe_str
from validators import (
    validate_required_columns,
    validate_tsv_schema,
    validate_no_duplicate_variants,
    validate_no_nulls,
)

log = logging.getLogger(__name__)

# Required columns for each data type
VARIANTS_REQUIRED = ["chrom", "pos", "ref", "alt", "type", "classification"]
METRICS_REQUIRED = ["variant_type", "precision", "recall", "f1"]
QC_METRICS_REQUIRED = ["metric_name", "metric_source"]
SOFTWARE_VERSIONS_REQUIRED = ["tool_name", "version"]


def _derive_zygosity(gt_str: str | None) -> str | None:
    """Derive zygosity from a genotype string (e.g. '0/1', '1/1').

    Args:
        gt_str: Genotype string like '0/1', '1/1', '0|1', etc.

    Returns:
        'HET', 'HOM_ALT', 'HOM_REF', or None if GT is missing/unparseable
    """
    if gt_str is None:
        return None

    # Normalize phased separator to unphased
    gt = gt_str.replace("|", "/")
    alleles = gt.split("/")

    if len(alleles) != 2:
        return None

    try:
        a1, a2 = alleles[0].strip(), alleles[1].strip()
        if a1 == "." or a2 == ".":
            return None
        if a1 == a2:
            return "HOM_REF" if a1 == "0" else "HOM_ALT"
        return "HET"
    except (ValueError, IndexError):
        return None


class TSVBulkInserter:
    """Generic TSV loader and bulk inserter for database models.

    This class encapsulates the common pattern of:
    1. Loading and validating TSV files
    2. Converting rows to model objects
    3. Bulk inserting into database

    Attributes:
        session: SQLAlchemy session
        run_id: Run ID to associate records with
    """

    def __init__(self, session: Session, run_id: int):
        """Initialize inserter.

        Args:
            session: Database session
            run_id: Run ID for associating records
        """
        self.session = session
        self.run_id = run_id

    def insert_from_tsv(
        self,
        tsv_path: str,
        model_class: type,
        required_columns: list[str],
        row_mapper: callable,
        entity_name: str = "records",
    ) -> int:
        """Load TSV and bulk insert into database.

        Args:
            tsv_path: Path to TSV file
            model_class: SQLAlchemy model class
            required_columns: List of required column names
            row_mapper: Function to map DataFrame row to model kwargs
            entity_name: Name of entities for logging (e.g., "metrics")

        Returns:
            Number of records inserted

        Raises:
            ValidationError: If data validation fails
        """
        # Validate and load TSV
        df = validate_tsv_schema(tsv_path, required_columns=required_columns)

        log.info("Loading %d %s...", len(df), entity_name)

        # Convert rows to model objects
        objects = []
        for _, row in df.iterrows():
            kwargs = row_mapper(row, self.run_id)
            obj = model_class(**kwargs)
            objects.append(obj)

        # Bulk insert
        self.session.add_all(objects)
        log.info("Inserted %d %s", len(objects), entity_name)

        return len(objects)


def load_qc_summary(qc_summary_path: str) -> dict[str, Any]:
    """Load QC summary data from TSV file.

    Args:
        qc_summary_path: Path to qc_summary.tsv

    Returns:
        Dictionary of QC summary fields

    Raises:
        ValidationError: If file is invalid
    """
    try:
        df = pd.read_csv(qc_summary_path, sep="\t")
        if len(df) > 0:
            summary = df.iloc[0].to_dict()
            log.info("Loaded QC summary with %d fields", len(summary))
            return summary
        return {}
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        log.warning("Failed to load QC summary: %s", e)
        return {}


def load_pipeline_metadata(metadata_path: str) -> dict[str, Any]:
    """Load pipeline metadata from TSV file.

    Args:
        metadata_path: Path to pipeline_metadata.tsv

    Returns:
        Dictionary of pipeline metadata fields
    """
    try:
        df = pd.read_csv(metadata_path, sep="\t")
        if len(df) > 0:
            metadata = df.iloc[0].to_dict()
            log.info("Loaded pipeline metadata with %d fields", len(metadata))
            return metadata
        return {}
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        log.warning("Failed to load pipeline metadata: %s", e)
        return {}


def create_run_from_qc(
    qc_summary: dict,
    pipeline_metadata: dict,
    params: dict,
) -> Any:  # Returns Run object
    """Create Run object from QC data and parameters.

    Args:
        qc_summary: QC summary dictionary
        pipeline_metadata: Pipeline metadata dictionary
        params: Snakemake parameters

    Returns:
        Run ORM object (not yet added to session)
    """
    from models.run import Run
    from datetime import UTC, datetime

    # Helper to safely extract and convert QC values
    def get_float(key: str) -> float | None:
        val = qc_summary.get(key)
        return safe_float(val) if pd.notna(val) else None

    def get_str(key: str) -> str | None:
        val = qc_summary.get(key)
        return safe_str(val) if pd.notna(val) else None

    def get_meta_str(key: str) -> str | None:
        val = pipeline_metadata.get(key)
        return safe_str(val) if pd.notna(val) else None

    run = Run(
        sample=params["sample"],
        pipeline_version=params["pipeline_version"],
        caller=params["caller"],
        parameters=params.get("parameters") or None,
        decomposition_mode=params.get("decomposition_mode") or None,
        comparison_tool=params.get("comparison_tool", "happy"),
        mode=params.get("mode", "germline"),
        tumor_sample=params.get("tumor_sample"),
        normal_sample=params.get("normal_sample"),
        calling_mode=params.get("calling_mode"),
        mean_coverage=get_float("mean_coverage"),
        median_coverage=get_float("median_coverage"),
        pct_bases_10x=get_float("pct_bases_10x"),
        pct_bases_20x=get_float("pct_bases_20x"),
        pct_bases_30x=get_float("pct_bases_30x"),
        mapping_rate=get_float("mapping_rate"),
        duplicate_rate=get_float("duplicate_rate"),
        median_insert_size=get_float("median_insert_size"),
        contamination_pct=get_float("contamination_pct"),
        sex_inferred=get_str("sex_inferred"),
        nfcore_pipeline_version=get_meta_str("nfcore_pipeline_version"),
        pipeline_params_json=get_meta_str("pipeline_params_json"),
        date=datetime.now(UTC),
    )

    log.info(
        "Created Run object (sample=%s, caller=%s, mean_cov=%.1f)",
        run.sample,
        run.caller,
        run.mean_coverage or 0,
    )

    return run


def load_and_insert_variants(
    variants_path: str, run_id: int, session: Session
) -> list[Any]:  # Returns list of Variant objects
    """Load variants from TSV and insert into database.

    Args:
        variants_path: Path to variants.tsv
        run_id: Run ID to associate variants with
        session: Database session

    Returns:
        List of inserted Variant ORM objects

    Raises:
        ValidationError: If data validation fails
        IntegrityError: If database constraint violated
    """
    from models.variant import Variant

    # Validate and load
    valid_classifications = {"TP", "FP", "FN", "TP_GT_MISMATCH", "TP_VAF_DISCORDANT", "PARTIAL"}
    valid_types = {"SNP", "INDEL", "COMPLEX", "SV"}

    variants_df = validate_tsv_schema(
        variants_path,
        required_columns=VARIANTS_REQUIRED,
        enum_columns={
            "type": valid_types,
            "classification": valid_classifications,
        },
        positive_columns=["pos", "dp", "gq"],
    )

    # Validate no duplicates
    validate_no_duplicate_variants(variants_df, Path(variants_path).name)

    # Validate required fields have no nulls
    validate_no_nulls(variants_df, VARIANTS_REQUIRED, Path(variants_path).name)

    log.info("Loading %d variants...", len(variants_df))

    # Create variant objects
    variant_objects = []
    for _, row in variants_df.iterrows():
        truth_gt = safe_str(row.get("truth_gt"))
        query_gt = safe_str(row.get("query_gt"))

        # Derive zygosity from genotype if not explicitly provided
        zygosity = safe_str(row.get("zygosity"))
        if zygosity is None:
            zygosity = _derive_zygosity(query_gt or truth_gt)

        # Compute gt_concordant from truth/query GTs
        gt_concordant = (
            truth_gt == query_gt
            if truth_gt is not None and query_gt is not None
            else None
        )

        v = Variant(
            run_id=run_id,
            chrom=str(row["chrom"]),
            pos=int(row["pos"]),
            ref=str(row["ref"]),
            alt=str(row["alt"]),
            type=str(row["type"]),
            classification=str(row["classification"]),
            dp=safe_int(row.get("dp")),
            mq=safe_float(row.get("mq")),
            af=safe_float(row.get("af")),
            gq=safe_int(row.get("gq")),
            qual=safe_float(row.get("qual")),
            filter_status=safe_str(row.get("filter_status")),
            # Somatic fields (nullable, populated when present in TSV)
            tumor_dp=safe_int(row.get("tumor_dp")),
            tumor_af=safe_float(row.get("tumor_af")),
            normal_dp=safe_int(row.get("normal_dp")),
            normal_af=safe_float(row.get("normal_af")),
            somatic_quality=safe_float(row.get("somatic_quality")),
            caller_filter=safe_str(row.get("caller_filter")),
            # Stratification fields
            zygosity=zygosity,
            indel_size=safe_int(row.get("indel_size")),
            # Genotype concordance fields
            truth_gt=truth_gt,
            query_gt=query_gt,
            gt_concordant=gt_concordant,
        )
        variant_objects.append(v)

    # Bulk insert
    session.add_all(variant_objects)
    session.flush()  # Assign IDs
    log.info("Inserted %d variants", len(variant_objects))

    return variant_objects


def _map_metric_row(row: pd.Series, run_id: int) -> dict:
    """Map DataFrame row to Metric model kwargs.

    Args:
        row: DataFrame row
        run_id: Run ID

    Returns:
        Dictionary of model fields
    """
    return {
        "run_id": run_id,
        "variant_type": str(row["variant_type"]),
        "stratification": safe_str(row.get("stratification")),
        "precision": float(row["precision"]),
        "recall": float(row["recall"]),
        "f1": float(row["f1"]),
        "tp_count": safe_int(row.get("tp_count")),
        "fp_count": safe_int(row.get("fp_count")),
        "fn_count": safe_int(row.get("fn_count")),
    }


def insert_metrics(metrics_path: str, run_id: int, session: Session) -> None:
    """Insert metrics from TSV into database.

    Args:
        metrics_path: Path to metrics.tsv
        run_id: Run ID to associate metrics with
        session: Database session

    Raises:
        ValidationError: If data validation fails
    """
    from models.metric import Metric

    inserter = TSVBulkInserter(session, run_id)
    inserter.insert_from_tsv(
        metrics_path,
        Metric,
        METRICS_REQUIRED,
        _map_metric_row,
        "metrics"
    )


def _map_qc_metric_row(row: pd.Series, run_id: int) -> dict:
    """Map DataFrame row to QCMetric model kwargs.

    Args:
        row: DataFrame row
        run_id: Run ID

    Returns:
        Dictionary of model fields
    """
    return {
        "run_id": run_id,
        "metric_name": str(row["metric_name"]),
        "metric_value_float": safe_float(row.get("metric_value_float")),
        "metric_value_str": safe_str(row.get("metric_value_str")),
        "metric_source": str(row["metric_source"]),
        "metric_category": safe_str(row.get("metric_category")),
        "notes": safe_str(row.get("notes")),
    }


def insert_qc_metrics(qc_metrics_path: str, run_id: int, session: Session) -> None:
    """Insert QC metrics from TSV into database.

    Args:
        qc_metrics_path: Path to qc_metrics.tsv
        run_id: Run ID to associate QC metrics with
        session: Database session

    Raises:
        ValidationError: If data validation fails
    """
    from models.qc_metric import QCMetric

    inserter = TSVBulkInserter(session, run_id)
    inserter.insert_from_tsv(
        qc_metrics_path,
        QCMetric,
        QC_METRICS_REQUIRED,
        _map_qc_metric_row,
        "QC metrics"
    )


def _map_software_version_row(row: pd.Series, run_id: int) -> dict:
    """Map DataFrame row to SoftwareVersion model kwargs.

    Args:
        row: DataFrame row
        run_id: Run ID

    Returns:
        Dictionary of model fields
    """
    return {
        "run_id": run_id,
        "tool_name": str(row["tool_name"]),
        "version": str(row["version"]),
        "category": safe_str(row.get("category")),
        "notes": safe_str(row.get("notes")),
    }


def insert_software_versions(
    versions_path: str, run_id: int, session: Session
) -> None:
    """Insert software versions from TSV into database.

    Args:
        versions_path: Path to software_versions.tsv
        run_id: Run ID to associate versions with
        session: Database session

    Raises:
        ValidationError: If data validation fails
    """
    from models.software_version import SoftwareVersion

    inserter = TSVBulkInserter(session, run_id)
    inserter.insert_from_tsv(
        versions_path,
        SoftwareVersion,
        SOFTWARE_VERSIONS_REQUIRED,
        _map_software_version_row,
        "software versions"
    )


def merge_vep_annotations(
    vep_annotations_path: str,
    variant_objects: list[Any],
) -> int:
    """Merge VEP annotations onto existing Variant objects by (chrom, pos, ref, alt).

    Args:
        vep_annotations_path: Path to vep_annotations.tsv
        variant_objects: List of Variant objects to annotate

    Returns:
        Number of variants annotated
    """
    try:
        df = pd.read_csv(vep_annotations_path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        log.warning("Failed to load VEP annotations: %s", e)
        return 0

    if len(df) == 0:
        log.info("No VEP annotations to merge")
        return 0

    # Build lookup by (chrom, pos, ref, alt)
    annot_map = {}
    for _, row in df.iterrows():
        key = (str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"]))
        annot_map[key] = row

    annotated = 0
    for variant in variant_objects:
        key = (variant.chrom, variant.pos, variant.ref, variant.alt)
        annot = annot_map.get(key)
        if annot is None:
            continue

        variant.consequence = safe_str(annot.get("consequence"))
        variant.impact = safe_str(annot.get("impact"))
        variant.gene_symbol = safe_str(annot.get("gene_symbol"))
        variant.gene_id = safe_str(annot.get("gene_id"))
        variant.transcript_id = safe_str(annot.get("transcript_id"))
        variant.hgvsc = safe_str(annot.get("hgvsc"))
        variant.hgvsp = safe_str(annot.get("hgvsp"))
        variant.exon = safe_str(annot.get("exon"))
        variant.protein_position = safe_int(annot.get("protein_position"))
        variant.sift_prediction = safe_str(annot.get("sift_prediction"))
        variant.polyphen_prediction = safe_str(annot.get("polyphen_prediction"))
        annotated += 1

    log.info("Merged VEP annotations onto %d / %d variants", annotated, len(variant_objects))
    return annotated


def process_gene_sets(
    gene_set_path: str,
    run_id: int,
    variant_objects: list[Any],
    gene_sets_config: list[dict],
    session: Session,
) -> None:
    """Process gene set assignments and create associations.

    Args:
        gene_set_path: Path to gene_set_assignments.tsv
        run_id: Run ID
        variant_objects: List of Variant objects (in order)
        gene_sets_config: Gene set configuration from Snakemake
        session: Database session

    Raises:
        ValidationError: If data validation fails
    """
    from models.gene_set import GeneSet, variant_gene_sets
    from sqlalchemy import insert

    gs_df = pd.read_csv(gene_set_path, sep="\t")

    if len(gs_df) == 0:
        log.info("No gene set assignments to process")
        return

    log.info("Processing %d gene set assignments...", len(gs_df))

    # Validate required columns
    validate_required_columns(
        gs_df, ["variant_idx", "gene_set_name"], "gene_set_assignments.tsv"
    )

    # Get or create GeneSet records
    gene_set_cache = {}
    gs_config_map = {gs["name"]: gs for gs in (gene_sets_config or [])}

    for gs_name in gs_df["gene_set_name"].unique():
        existing = session.query(GeneSet).filter(GeneSet.name == gs_name).first()
        if existing:
            gene_set_cache[gs_name] = existing
        else:
            gs_cfg = gs_config_map.get(gs_name, {})
            new_gs = GeneSet(
                name=gs_name,
                description=gs_cfg.get("description"),
                version=gs_cfg.get("version"),
            )
            session.add(new_gs)
            session.flush()
            gene_set_cache[gs_name] = new_gs
            log.info("Created GeneSet: %s (id=%d)", gs_name, new_gs.id)

    # Create variant_gene_sets associations
    association_rows = []
    seen = set()

    for _, row in gs_df.iterrows():
        variant_idx = int(row["variant_idx"])
        gs_name = row["gene_set_name"]
        gene = safe_str(row.get("gene"))

        if variant_idx >= len(variant_objects):
            log.warning(
                "Invalid variant_idx %d (max: %d), skipping",
                variant_idx,
                len(variant_objects) - 1,
            )
            continue

        variant = variant_objects[variant_idx]
        gene_set = gene_set_cache[gs_name]

        # Deduplicate
        assoc_key = (variant.id, gene_set.id)
        if assoc_key in seen:
            continue
        seen.add(assoc_key)

        association_rows.append({"variant_id": variant.id, "gene_set_id": gene_set.id})

        # Set gene name on variant if not already set
        if gene and not variant.gene:
            variant.gene = gene

    if association_rows:
        session.execute(insert(variant_gene_sets), association_rows)
        log.info("Created %d variant-gene_set associations", len(association_rows))


def insert_stratified_metrics(
    variant_objects: list[Any],
    run_id: int,
    session: Session,
    stratification_config: dict | None = None,
    mode: str = "germline",
    gene_set_names: list[str] | None = None,
) -> int:
    """Compute stratified metrics and insert into database.

    Args:
        variant_objects: List of Variant objects
        run_id: Run ID to associate metrics with
        session: Database session
        stratification_config: Stratification config dict
        mode: Execution mode (germline or somatic)
        gene_set_names: List of gene set names

    Returns:
        Number of stratified metric records inserted
    """
    from analysis.stratification import StratificationEngine
    from models.stratified_metric import StratifiedMetric

    if not variant_objects:
        log.info("No variants to stratify")
        return 0

    engine = StratificationEngine(
        config=stratification_config or {},
        mode=mode,
        gene_set_names=gene_set_names or [],
    )

    strata = engine.stratify(variant_objects)
    metrics = engine.compute_metrics(strata)

    objects = []
    for m in metrics:
        sm = StratifiedMetric(
            run_id=run_id,
            dimension=m["dimension"],
            stratum=m["stratum"],
            variant_type=m["variant_type"],
            precision=m.get("precision"),
            recall=m.get("recall"),
            f1=m.get("f1"),
            tp_count=m.get("tp_count"),
            fp_count=m.get("fp_count"),
            fn_count=m.get("fn_count"),
            total_variants=m.get("total_variants"),
            precision_ci_lower=m.get("precision_ci_lower"),
            precision_ci_upper=m.get("precision_ci_upper"),
            recall_ci_lower=m.get("recall_ci_lower"),
            recall_ci_upper=m.get("recall_ci_upper"),
            f1_ci_lower=m.get("f1_ci_lower"),
            f1_ci_upper=m.get("f1_ci_upper"),
            mcc=m.get("mcc"),
            genotype_concordance=m.get("genotype_concordance"),
            low_confidence=m.get("low_confidence", False),
        )
        objects.append(sm)

    if objects:
        session.add_all(objects)
        log.info("Inserted %d stratified metric records", len(objects))

    return len(objects)
