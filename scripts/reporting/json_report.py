"""JSON report generation (Spec S16.1).

Produces a structured, versioned JSON report containing run metadata,
metrics, acceptance decisions, variant diffs, and root cause evidence.
"""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import UTC, datetime

log = logging.getLogger(__name__)

REPORT_VERSION = "1.0"


def generate_json_report(
    report_type: str = "initial_validation",
    mode: str = "germline",
    run_metadata: dict | None = None,
    input_checksums: dict | None = None,
    annotation_config: dict | None = None,
    global_metrics: dict | None = None,
    stratified_metrics: list[dict] | None = None,
    acceptance_result: dict | None = None,
    variant_diff: list[dict] | None = None,
    root_cause_evidence: list[dict] | None = None,
    baseline_name: str | None = None,
    bayesian_risk: dict | None = None,
    qc_trending: list[dict] | None = None,
) -> dict:
    """Generate a structured JSON report per S16.1 schema.

    Args:
        report_type: "initial_validation" or "continuous_verification"
        mode: "germline" or "somatic"
        run_metadata: Dict with sample, pipeline_version, caller, etc.
        input_checksums: Dict of file paths to checksums
        annotation_config: VEP/annotation configuration
        global_metrics: Dict with global sensitivity, precision, f1
        stratified_metrics: List of per-stratum metric dicts
        acceptance_result: Acceptance evaluation result dict
        variant_diff: List of variant transition dicts (verification only)
        root_cause_evidence: List of evidence dicts (verification only)
        baseline_name: Baseline name for verification reports
        bayesian_risk: Optional BRAM assessment result dict (S17.8)
        qc_trending: Optional list of QC trending results

    Returns:
        Complete JSON report as a dict
    """
    report = {
        "report_type": report_type,
        "report_version": REPORT_VERSION,
        "generated_at": datetime.now(UTC).isoformat(),
        "mode": mode,
        "run_metadata": run_metadata or {},
        "input_checksums": input_checksums or {},
        "annotation_config": annotation_config or {},
        "metrics": {
            "global": global_metrics or {},
            "per_stratum": _format_stratified_metrics(stratified_metrics or []),
        },
        "acceptance": _format_acceptance(acceptance_result),
    }

    if baseline_name:
        report["run_metadata"]["baseline_name"] = baseline_name

    if report_type == "continuous_verification":
        report["variant_diff"] = variant_diff or []
        report["root_cause_evidence"] = root_cause_evidence or []

    if bayesian_risk:
        report["bayesian_risk"] = _format_bayesian_risk(bayesian_risk)

    if qc_trending:
        report["qc_trending"] = qc_trending

    log.info(
        "Generated %s JSON report (version=%s, mode=%s)",
        report_type,
        REPORT_VERSION,
        mode,
    )
    return report


def _format_stratified_metrics(metrics: list[dict]) -> list[dict]:
    """Format stratified metrics for JSON output."""
    formatted = []
    for m in metrics:
        entry = {
            "dimension": m.get("dimension"),
            "stratum": m.get("stratum"),
            "variant_type": m.get("variant_type", "ALL"),
            "total_variants": m.get("total_variants"),
            "tp_count": m.get("tp_count"),
            "fp_count": m.get("fp_count"),
            "fn_count": m.get("fn_count"),
            "precision": m.get("precision"),
            "recall": m.get("recall"),
            "f1": m.get("f1"),
            "low_confidence": m.get("low_confidence", False),
            "precision_ci_lower": m.get("precision_ci_lower"),
            "precision_ci_upper": m.get("precision_ci_upper"),
            "recall_ci_lower": m.get("recall_ci_lower"),
            "recall_ci_upper": m.get("recall_ci_upper"),
            "f1_ci_lower": m.get("f1_ci_lower"),
            "f1_ci_upper": m.get("f1_ci_upper"),
            "mcc": m.get("mcc"),
            "genotype_concordance": m.get("genotype_concordance"),
        }
        formatted.append(entry)
    return formatted


def _format_acceptance(result: dict | None) -> dict:
    """Format acceptance result for JSON output."""
    if not result:
        return {"overall_result": "NOT_EVALUATED", "per_stratum_results": [], "violations": []}

    violations = []
    for sr in result.get("strata_results", []):
        if sr.get("violations"):
            violations.append(
                {
                    "dimension": sr["dimension"],
                    "stratum": sr["stratum"],
                    "tier": sr["tier"],
                    "result": sr["result"],
                    "violations": sr["violations"],
                }
            )

    return {
        "overall_result": result.get("verdict", "NOT_EVALUATED"),
        "per_stratum_results": [
            {
                "dimension": sr["dimension"],
                "stratum": sr["stratum"],
                "tier": sr["tier"],
                "result": sr["result"],
            }
            for sr in result.get("strata_results", [])
        ],
        "violations": violations,
    }


def _format_bayesian_risk(result: dict) -> dict:
    """Format BRAM result for JSON output (S17.8)."""
    per_stratum = []
    for s in result.get("per_stratum", []):
        per_stratum.append(
            {
                "dimension": s.get("dimension"),
                "stratum": s.get("stratum"),
                "metric": s.get("metric_name"),
                "delta_observed": s.get("delta_observed"),
                "posterior_mu": s.get("posterior_mu"),
                "posterior_sigma": s.get("posterior_sigma"),
                "tail_probability": s.get("tail_probability"),
                "risk_weight": s.get("risk_weight"),
                "weighted_risk": s.get("weighted_risk"),
                "flagged": s.get("flagged", False),
                "tier": s.get("tier"),
            }
        )

    return {
        "enabled": True,
        "verdict": result.get("verdict", "BRAM_NOT_RUN"),
        "aggregate_risk_score": result.get("aggregate_risk_score", 0.0),
        "mean_risk_score": result.get("mean_risk_score", 0.0),
        "flagged_count": result.get("flagged_stratum_count", 0),
        "alert_threshold": result.get("alert_threshold"),
        "degradation_threshold": result.get("degradation_threshold"),
        "per_stratum": per_stratum,
    }


def serialize_report(report: dict) -> str:
    """Serialize a report dict to a JSON string.

    Args:
        report: Report dict

    Returns:
        Pretty-printed JSON string
    """
    return json.dumps(report, indent=2, default=str)


def compute_report_checksum(report_json: str) -> str:
    """Compute SHA-256 checksum of a serialized report.

    Args:
        report_json: JSON string

    Returns:
        Hex digest
    """
    return hashlib.sha256(report_json.encode()).hexdigest()
