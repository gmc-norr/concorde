"""Parse nf-core/raredisease QC outputs into TSV files for database ingestion.

Called via Snakemake `script:` directive.
Parses multiple QC sources:
1. MultiQC JSON - aggregated metrics from all QC tools
2. Mosdepth summary - coverage statistics
3. VerifyBamID2 selfSM - contamination estimates
4. Software versions YAML - tool versions used in pipeline
5. Pipeline params JSON - pipeline parameters

Outputs four TSV files:
- qc_summary.tsv: Key summary metrics for the Run table
- qc_metrics.tsv: Detailed metrics for the QCMetric table
- software_versions.tsv: Tool versions for SoftwareVersion table
- pipeline_metadata.tsv: Pipeline version and parameters
"""

import json
import logging
from pathlib import Path
from typing import Any, TYPE_CHECKING

import pandas as pd
import yaml

from utils import safe_float

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
log = logging.getLogger(__name__)


# Constants for metric sources
class MetricSource:
    """Tool sources for QC metrics."""
    MULTIQC = "multiqc"
    MOSDEPTH = "mosdepth"
    PICARD = "picard"
    SAMTOOLS = "samtools"
    FASTQC = "fastqc"


# Metric name prefixes mapped to sources
METRIC_PREFIXES = {
    "mosdepth_": MetricSource.MOSDEPTH,
    "picard_": MetricSource.PICARD,
    "samtools_": MetricSource.SAMTOOLS,
    "fastqc_": MetricSource.FASTQC,
}


# Constants for metric categories
class MetricCategory:
    """Categories for organizing QC metrics."""
    COVERAGE = "coverage"
    ALIGNMENT = "alignment"
    DUPLICATION = "duplication"
    INSERT_SIZE = "insert_size"
    BASE_QUALITY = "base_quality"
    GC_BIAS = "gc_bias"


# Keywords for categorizing metrics
CATEGORY_KEYWORDS = {
    MetricCategory.COVERAGE: ["coverage", "depth"],
    MetricCategory.ALIGNMENT: ["map", "align"],
    MetricCategory.DUPLICATION: ["dup"],
    MetricCategory.INSERT_SIZE: ["insert"],
    MetricCategory.BASE_QUALITY: ["quality", "qual"],
    MetricCategory.GC_BIAS: ["gc"],
}


def _determine_metric_source(metric_name: str) -> str:
    """Determine the source tool from metric name prefix.

    Args:
        metric_name: Name of the metric

    Returns:
        Source tool name (e.g., "mosdepth", "picard")
    """
    for prefix, source in METRIC_PREFIXES.items():
        if metric_name.startswith(prefix):
            return source
    return MetricSource.MULTIQC


def _determine_metric_category(metric_name: str) -> str | None:
    """Determine the category from metric name.

    Args:
        metric_name: Name of the metric

    Returns:
        Category name or None if no match
    """
    key_lower = metric_name.lower()
    for category, keywords in CATEGORY_KEYWORDS.items():
        if any(keyword in key_lower for keyword in keywords):
            return category
    return None


def _extract_summary_metrics(sample_data: dict) -> dict:
    """Extract summary metrics for Run table from MultiQC data.

    Args:
        sample_data: Sample-specific data from MultiQC

    Returns:
        Dictionary of summary metrics
    """
    return {
        "mean_coverage": safe_float(
            sample_data.get("mosdepth_mean_coverage")
            or sample_data.get("picard_PCT_TARGET_BASES_30X")
        ),
        "duplicate_rate": safe_float(
            sample_data.get("picard_PERCENT_DUPLICATION")
        ),
        "mapping_rate": safe_float(
            sample_data.get("samtools_mapped_passed")
            or sample_data.get("picard_PCT_PF_READS_ALIGNED")
        ),
        "median_insert_size": safe_float(
            sample_data.get("picard_MEDIAN_INSERT_SIZE")
        ),
    }


def _build_detailed_metric(key: str, value: Any) -> dict:
    """Build a detailed metric entry for QCMetric table.

    Args:
        key: Metric name
        value: Metric value

    Returns:
        Dictionary representing a metric entry
    """
    float_val = safe_float(value)
    return {
        "metric_name": key,
        "metric_value_float": float_val,
        "metric_value_str": str(value) if float_val is None else None,
        "metric_source": _determine_metric_source(key),
        "metric_category": _determine_metric_category(key),
        "notes": "Extracted from MultiQC general statistics",
    }


def _extract_detailed_metrics(sample_data: dict) -> list[dict]:
    """Extract detailed metrics for QCMetric table from MultiQC data.

    Args:
        sample_data: Sample-specific data from MultiQC

    Returns:
        List of metric dictionaries
    """
    return [_build_detailed_metric(key, value) for key, value in sample_data.items()]


def _find_sample_data(general_stats: list[dict], sample_id: str) -> dict | None:
    """Find sample-specific data in MultiQC general stats.

    Args:
        general_stats: List of general statistics dictionaries
        sample_id: Sample identifier to find

    Returns:
        Sample data dictionary or None if not found
    """
    # Try exact match first
    for stats_dict in general_stats:
        if sample_id in stats_dict:
            return stats_dict[sample_id]

    # Fall back to first sample if exact match not found
    if general_stats and len(general_stats[0]) > 0:
        first_sample = list(general_stats[0].keys())[0]
        log.warning(
            "Sample %s not found in MultiQC, using first sample: %s",
            sample_id,
            first_sample
        )
        return general_stats[0][first_sample]

    return None


def parse_multiqc_json(json_path, sample_id):
    """Parse MultiQC JSON and extract key metrics.

    Args:
        json_path: Path to MultiQC JSON file
        sample_id: Sample identifier

    Returns:
        tuple: (summary metrics dict, detailed metrics list)
    """
    summary = {}
    detailed_metrics = []

    # Validate file exists
    if not Path(json_path).exists():
        log.warning("MultiQC JSON not found: %s", json_path)
        return summary, detailed_metrics

    log.info("Parsing MultiQC JSON: %s", json_path)

    # Load and parse JSON
    with open(json_path) as f:
        data = json.load(f)

    general_stats = data.get("report_general_stats_data", [])
    if not general_stats:
        log.warning("No general statistics found in MultiQC JSON")
        return summary, detailed_metrics

    # Find sample data
    sample_data = _find_sample_data(general_stats, sample_id)
    if not sample_data:
        log.warning("No sample data found for %s", sample_id)
        return summary, detailed_metrics

    # Extract metrics
    summary = _extract_summary_metrics(sample_data)
    detailed_metrics = _extract_detailed_metrics(sample_data)

    log.info("Extracted %d metrics from MultiQC", len(detailed_metrics))
    return summary, detailed_metrics


def parse_mosdepth_summary(summary_path, sample_id):
    """Parse Mosdepth summary.txt for coverage statistics.

    Returns:
        dict: Summary metrics for Run table
        list: Detailed metrics for QCMetric table
    """
    summary = {}
    detailed_metrics = []

    if not Path(summary_path).exists():
        log.warning("Mosdepth summary not found: %s", summary_path)
        return summary, detailed_metrics

    log.info("Parsing Mosdepth summary: %s", summary_path)

    df = pd.read_csv(summary_path, sep="\t")

    # Get total/autosomal coverage
    total_row = df[df["chrom"] == "total"]
    if len(total_row) > 0:
        mean_cov = safe_float(total_row.iloc[0]["mean"])
        summary["mean_coverage"] = mean_cov

        detailed_metrics.append(
            {
                "metric_name": "mean_coverage_total",
                "metric_value_float": mean_cov,
                "metric_value_str": None,
                "metric_source": "mosdepth",
                "metric_category": "coverage",
                "notes": "Mean coverage across all contigs",
            }
        )

    # Per-chromosome coverage
    for _, row in df.iterrows():
        chrom = row["chrom"]
        if chrom == "total":
            continue

        detailed_metrics.append(
            {
                "metric_name": f"mean_coverage_{chrom}",
                "metric_value_float": safe_float(row["mean"]),
                "metric_value_str": None,
                "metric_source": "mosdepth",
                "metric_category": "coverage",
                "notes": f"Mean coverage for {chrom}",
            }
        )

    log.info("Extracted %d coverage metrics from Mosdepth", len(detailed_metrics))
    return summary, detailed_metrics


def parse_mosdepth_dist(dist_path, sample_id):
    """Parse Mosdepth global.dist.txt for coverage distribution.

    Returns:
        dict: Summary metrics for Run table (% bases at various depths)
        list: Detailed metrics for QCMetric table
    """
    summary = {}
    detailed_metrics = []

    if not Path(dist_path).exists():
        log.warning("Mosdepth distribution file not found: %s", dist_path)
        return summary, detailed_metrics

    log.info("Parsing Mosdepth distribution: %s", dist_path)

    df = pd.read_csv(dist_path, sep="\t", names=["chrom", "depth", "count"])

    # Calculate % bases at key depth thresholds (for total/genome-wide)
    total_data = df[df["chrom"] == "total"]

    if len(total_data) == 0:
        log.warning("No total coverage distribution found in Mosdepth dist file")
        return summary, detailed_metrics

    # Cumulative coverage calculation
    total_data = total_data.sort_values("depth")
    total_bases = total_data["count"].sum()

    for threshold in [10, 20, 30]:
        bases_above = total_data[total_data["depth"] >= threshold]["count"].sum()
        pct_above = (bases_above / total_bases * 100) if total_bases > 0 else 0

        field_name = f"pct_bases_{threshold}x"
        summary[field_name] = pct_above

        detailed_metrics.append(
            {
                "metric_name": field_name,
                "metric_value_float": pct_above,
                "metric_value_str": None,
                "metric_source": "mosdepth",
                "metric_category": "coverage",
                "notes": f"Percentage of bases with coverage >= {threshold}X",
            }
        )

    # Calculate median coverage (50th percentile)
    cumsum = 0
    median_depth = None
    for _, row in total_data.iterrows():
        cumsum += row["count"]
        if cumsum >= total_bases / 2:
            median_depth = row["depth"]
            break

    if median_depth is not None:
        summary["median_coverage"] = float(median_depth)
        detailed_metrics.append(
            {
                "metric_name": "median_coverage",
                "metric_value_float": float(median_depth),
                "metric_value_str": None,
                "metric_source": "mosdepth",
                "metric_category": "coverage",
                "notes": "Median coverage depth across genome",
            }
        )

    log.info(
        "Extracted coverage distribution metrics: %dX, %dX, %dX",
        summary.get("pct_bases_10x", 0),
        summary.get("pct_bases_20x", 0),
        summary.get("pct_bases_30x", 0),
    )
    return summary, detailed_metrics


def parse_verifybamid(selfsm_path, sample_id):
    """Parse VerifyBamID2 selfSM file for contamination estimates.

    Returns:
        dict: Summary metrics for Run table
        list: Detailed metrics for QCMetric table
    """
    summary = {}
    detailed_metrics = []

    if not Path(selfsm_path).exists():
        log.warning("VerifyBamID2 selfSM file not found: %s", selfsm_path)
        return summary, detailed_metrics

    log.info("Parsing VerifyBamID2 selfSM: %s", selfsm_path)

    df = pd.read_csv(selfsm_path, sep="\t")

    if len(df) == 0:
        log.warning("VerifyBamID2 selfSM file is empty")
        return summary, detailed_metrics

    row = df.iloc[0]

    # Extract FREEMIX (contamination estimate)
    contamination = safe_float(row.get("FREEMIX", 0))
    summary["contamination_pct"] = contamination * 100 if contamination else None

    # Store all VerifyBamID2 metrics
    for col in df.columns:
        value = row[col]
        float_val = safe_float(value)

        detailed_metrics.append(
            {
                "metric_name": f"verifybamid_{col.lower()}",
                "metric_value_float": float_val,
                "metric_value_str": str(value) if float_val is None else None,
                "metric_source": "verifybamid2",
                "metric_category": "contamination",
                "notes": f"VerifyBamID2 {col} metric",
            }
        )

    log.info("Contamination estimate: %.4f%%", summary.get("contamination_pct", 0))
    return summary, detailed_metrics


def parse_software_versions(versions_yml_path):
    """Parse nf-core/raredisease software_versions.yml file.

    Returns:
        list: Software version records for SoftwareVersion table
        dict: Metadata including pipeline version
    """
    versions = []
    metadata = {}

    if not Path(versions_yml_path).exists():
        log.warning("Software versions YAML not found: %s", versions_yml_path)
        return versions, metadata

    log.info("Parsing software versions YAML: %s", versions_yml_path)

    with open(versions_yml_path) as f:
        data = yaml.safe_load(f)

    if not data:
        log.warning("Software versions YAML is empty")
        return versions, metadata

    # Extract nf-core/raredisease pipeline version
    if "Workflow" in data:
        workflow_info = data["Workflow"]
        if "nf-core/raredisease" in workflow_info:
            metadata["nfcore_pipeline_version"] = workflow_info["nf-core/raredisease"]
            log.info("Pipeline version: %s", metadata["nfcore_pipeline_version"])

    # Parse tool versions
    for category, tools in data.items():
        if category == "Workflow":
            continue

        if isinstance(tools, dict):
            for tool_name, version in tools.items():
                # Clean up tool name and version
                tool_name_clean = tool_name.strip()
                version_clean = str(version).strip() if version else "unknown"

                versions.append(
                    {
                        "tool_name": tool_name_clean,
                        "version": version_clean,
                        "category": category,
                        "notes": f"From nf-core/raredisease {category} category",
                    }
                )

    log.info("Extracted %d tool versions", len(versions))
    return versions, metadata


def parse_pipeline_params(params_json_path):
    """Parse nf-core/raredisease pipeline params.json file.

    Returns:
        str: JSON string of pipeline parameters
    """
    if not Path(params_json_path).exists():
        log.warning("Pipeline params JSON not found: %s", params_json_path)
        return None

    log.info("Parsing pipeline params JSON: %s", params_json_path)

    with open(params_json_path) as f:
        params = json.load(f)

    # Return as JSON string for storage
    params_json_str = json.dumps(params, indent=2)
    log.info("Loaded pipeline parameters (%d keys)", len(params))

    return params_json_str


# --- Main execution via Snakemake ---
sample_id = snakemake.params.get("sample_id", "unknown")

log.info("Parsing nf-core/raredisease QC outputs for sample: %s", sample_id)

# Parse all QC sources
all_summary = {}
all_detailed = []

# 1. MultiQC JSON
if snakemake.input.get("multiqc_json"):
    mqc_summary, mqc_detailed = parse_multiqc_json(
        snakemake.input.multiqc_json, sample_id
    )
    all_summary.update(mqc_summary)
    all_detailed.extend(mqc_detailed)

# 2. Mosdepth summary
if snakemake.input.get("mosdepth_summary"):
    md_summary, md_detailed = parse_mosdepth_summary(
        snakemake.input.mosdepth_summary, sample_id
    )
    # Prefer Mosdepth values over MultiQC if both exist
    all_summary.update(md_summary)
    all_detailed.extend(md_detailed)

# 3. Mosdepth distribution
if snakemake.input.get("mosdepth_dist"):
    mdd_summary, mdd_detailed = parse_mosdepth_dist(
        snakemake.input.mosdepth_dist, sample_id
    )
    all_summary.update(mdd_summary)
    all_detailed.extend(mdd_detailed)

# 4. VerifyBamID2
if snakemake.input.get("verifybamid_selfsm"):
    vb_summary, vb_detailed = parse_verifybamid(
        snakemake.input.verifybamid_selfsm, sample_id
    )
    all_summary.update(vb_summary)
    all_detailed.extend(vb_detailed)

# 5. Software versions
software_versions = []
pipeline_metadata = {}
if snakemake.input.get("software_versions_yml"):
    software_versions, version_metadata = parse_software_versions(
        snakemake.input.software_versions_yml
    )
    pipeline_metadata.update(version_metadata)

# 6. Pipeline parameters
if snakemake.input.get("pipeline_params_json"):
    params_json = parse_pipeline_params(snakemake.input.pipeline_params_json)
    if params_json:
        pipeline_metadata["pipeline_params_json"] = params_json

# Write summary TSV (for Run table)
summary_df = pd.DataFrame([all_summary])
summary_df.to_csv(snakemake.output.summary_tsv, sep="\t", index=False)
log.info("Wrote QC summary to %s", snakemake.output.summary_tsv)

# Write detailed metrics TSV (for QCMetric table)
detailed_df = pd.DataFrame(all_detailed)
detailed_df.to_csv(snakemake.output.metrics_tsv, sep="\t", index=False)
log.info(
    "Wrote %d detailed QC metrics to %s", len(detailed_df), snakemake.output.metrics_tsv
)

# Write software versions TSV (for SoftwareVersion table)
versions_df = pd.DataFrame(software_versions)
versions_df.to_csv(snakemake.output.versions_tsv, sep="\t", index=False)
log.info(
    "Wrote %d software versions to %s", len(versions_df), snakemake.output.versions_tsv
)

# Write pipeline metadata TSV (for Run table nfcore fields)
metadata_df = pd.DataFrame([pipeline_metadata])
metadata_df.to_csv(snakemake.output.metadata_tsv, sep="\t", index=False)
log.info("Wrote pipeline metadata to %s", snakemake.output.metadata_tsv)

log.info("QC parsing complete")
