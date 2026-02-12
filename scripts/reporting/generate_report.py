"""Snakemake script for generating validation reports."""

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


def generate_reports(
    database: str,
    output_json: str,
    output_html: str,
    sample: str,
    mode: str = "germline",
    pipeline_version: str = "",
    caller: str = "",
    comparison_tool: str = "",
    acceptance_config: dict | None = None,
    qc_trending_config: dict | None = None,
) -> None:
    """Generate JSON and HTML validation reports."""
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker

    from models.base import Base
    from models.run import Run
    from models.stratified_metric import StratifiedMetric
    from models.verification import VerificationResult
    from analysis.acceptance import AcceptanceEngine
    from reporting.json_report import generate_json_report, serialize_report
    from reporting.html_report import generate_html_report

    engine = create_engine(f"sqlite:///{database}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        # Get the most recent run
        run = session.query(Run).order_by(Run.id.desc()).first()
        if not run:
            raise ValueError("No runs found in database")

        # Get stratified metrics
        strat_metrics = (
            session.query(StratifiedMetric)
            .filter(StratifiedMetric.run_id == run.id)
            .all()
        )

        metrics_list = [
            {
                "dimension": sm.dimension,
                "stratum": sm.stratum,
                "variant_type": sm.variant_type,
                "precision": sm.precision,
                "recall": sm.recall,
                "f1": sm.f1,
                "tp_count": sm.tp_count,
                "fp_count": sm.fp_count,
                "fn_count": sm.fn_count,
                "total_variants": sm.total_variants,
                "low_confidence": sm.low_confidence,
                "precision_ci_lower": sm.precision_ci_lower,
                "precision_ci_upper": sm.precision_ci_upper,
                "recall_ci_lower": sm.recall_ci_lower,
                "recall_ci_upper": sm.recall_ci_upper,
                "f1_ci_lower": sm.f1_ci_lower,
                "f1_ci_upper": sm.f1_ci_upper,
                "mcc": sm.mcc,
                "genotype_concordance": sm.genotype_concordance,
            }
            for sm in strat_metrics
        ]

        # Run acceptance evaluation
        ac_engine = AcceptanceEngine(acceptance_config or {})
        acceptance_result = ac_engine.evaluate_all(metrics_list)

        # Query verification data if available
        variant_diff_data = None
        rca_data = None
        baseline_name = None
        vr = (
            session.query(VerificationResult)
            .filter(VerificationResult.run_id == run.id)
            .first()
        )
        if vr:
            baseline_name = None  # Will be populated if baseline queried
            variant_diff_data = []
            rca_data = []
            for t in vr.transitions:
                variant_diff_data.append({
                    "chrom": t.chrom,
                    "pos": t.pos,
                    "ref": t.ref,
                    "alt": t.alt,
                    "transition": t.transition,
                    "diff_class": t.diff_class,
                })
                for e in t.evidence:
                    rca_data.append({
                        "chrom": t.chrom,
                        "pos": t.pos,
                        "category": e.category,
                        "score": e.score,
                        "description": e.change_description,
                    })

        # Compute QC trending if configured
        qc_trending_data = None
        trending_cfg = qc_trending_config or {}
        if trending_cfg.get("enabled", False):
            from analysis.qc_trending import compute_trending_report
            qc_trending_data = compute_trending_report(
                session, run.id, config=trending_cfg
            )

        # Generate JSON report
        report_type = "continuous_verification" if vr else "initial_validation"
        report_data = generate_json_report(
            report_type=report_type,
            mode=mode,
            run_metadata={
                "sample": sample,
                "pipeline_version": pipeline_version,
                "caller": caller,
                "comparison_tool": comparison_tool,
            },
            stratified_metrics=metrics_list,
            acceptance_result=acceptance_result,
            variant_diff=variant_diff_data,
            root_cause_evidence=rca_data,
            baseline_name=baseline_name,
            qc_trending=qc_trending_data,
        )

        # Write JSON
        Path(output_json).parent.mkdir(parents=True, exist_ok=True)
        Path(output_json).write_text(serialize_report(report_data))
        logging.info("Wrote JSON report: %s", output_json)

        # Generate and write HTML
        html = generate_html_report(report_data)
        Path(output_html).write_text(html)
        logging.info("Wrote HTML report: %s", output_html)

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
    generate_reports(
        database=params.database,
        output_json=snakemake.output.json,
        output_html=snakemake.output.html,
        sample=params.sample,
        mode=params.execution_mode,
        pipeline_version=params.pipeline_version,
        caller=params.caller,
        comparison_tool=params.comparison_tool,
        acceptance_config=params.acceptance_config,
    )

except NameError:
    pass  # Not running via Snakemake
