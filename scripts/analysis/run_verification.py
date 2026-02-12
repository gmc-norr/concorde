"""Snakemake script for running continuous verification against a baseline.

Orchestrates the full verification workflow:
1. Load baseline + current run variants from DB
2. Compute variant diff
3. Stratify and annotate transitions
4. Check baseline envelopes
5. Classify transitions (drift vs biological)
6. Collect root cause evidence
7. Determine overall verdict
8. Write VerificationResult + VariantTransition + RootCauseEvidence to DB
9. Optionally run BRAM assessment
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent)
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

log = logging.getLogger(__name__)


def run_verification(
    database: str,
    baseline_name: str,
    mode: str = "germline",
    stratification_config: dict | None = None,
    bram_config: dict | None = None,
    acceptance_config: dict | None = None,
) -> dict:
    """Run continuous verification against a locked baseline.

    Args:
        database: Path to SQLite database
        baseline_name: Name of the locked baseline to verify against
        mode: Execution mode (germline or somatic)
        stratification_config: Stratification config dict
        bram_config: Optional BRAM configuration
        acceptance_config: Optional acceptance criteria config

    Returns:
        Dict with verdict, transition counts, and summary
    """
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker

    from models.base import Base
    from models.baseline import Baseline
    from models.evidence import RootCauseEvidence
    from models.run import Run
    from models.variant import Variant
    from models.verification import VariantTransition, VerificationResult

    from analysis.baseline_manager import BaselineManager
    from analysis.diff_classifier import classify_transitions, determine_verdict
    from analysis.root_cause import collect_evidence, determine_likely_cause
    from analysis.stratification import StratificationEngine
    from analysis.variant_diff import assign_strata_to_transitions, compute_variant_diff

    engine = create_engine(f"sqlite:///{database}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        # 1. Load baseline and current run
        baseline = session.query(Baseline).filter(
            Baseline.name == baseline_name,
            Baseline.locked.is_(True),
        ).first()

        if not baseline:
            raise ValueError(f"No locked baseline found with name: {baseline_name}")

        current_run = session.query(Run).order_by(Run.id.desc()).first()
        if not current_run:
            raise ValueError("No runs found in database")

        baseline_run = session.query(Run).filter(Run.id == baseline.run_id).first()
        if not baseline_run:
            raise ValueError(f"Baseline run id={baseline.run_id} not found")

        log.info(
            "Verifying run id=%d against baseline '%s' (run id=%d)",
            current_run.id, baseline_name, baseline_run.id,
        )

        # 2. Load variants
        baseline_variants = session.query(Variant).filter(
            Variant.run_id == baseline_run.id
        ).all()
        current_variants = session.query(Variant).filter(
            Variant.run_id == current_run.id
        ).all()

        # 3. Compute variant diff
        transitions = compute_variant_diff(baseline_variants, current_variants)

        # 4. Stratify current variants and annotate transitions
        strat_engine = StratificationEngine(
            config=stratification_config or {},
            mode=mode,
        )
        strata = strat_engine.stratify(current_variants)
        transitions = assign_strata_to_transitions(transitions, current_variants, strata)

        # 5. Check baseline envelopes
        baseline_mgr = BaselineManager(session)
        strat_metrics = strat_engine.compute_metrics(strata)
        envelope_violations = baseline_mgr.check_envelopes(baseline.id, strat_metrics)

        # 6. Build high-impact and regulated sets
        high_impact_keys = set()
        for v in current_variants:
            if getattr(v, "impact", None) in ("HIGH", "MODERATE"):
                high_impact_keys.add((v.chrom, v.pos, v.ref, v.alt))

        # 7. Classify transitions
        transitions = classify_transitions(
            transitions,
            envelope_violations,
            high_impact_variants=high_impact_keys,
        )

        # 8. Collect evidence for non-drift transitions
        variant_map_base = {
            (v.chrom, v.pos, v.ref, v.alt): {
                "dp": v.dp, "mq": v.mq, "qual": v.qual, "filter_status": v.filter_status
            }
            for v in baseline_variants
        }
        variant_map_curr = {
            (v.chrom, v.pos, v.ref, v.alt): {
                "dp": v.dp, "mq": v.mq, "qual": v.qual, "filter_status": v.filter_status
            }
            for v in current_variants
        }

        for t in transitions:
            key = (t["chrom"], t["pos"], t["ref"], t["alt"])
            t["evidence"] = collect_evidence(
                t,
                baseline_variant=variant_map_base.get(key),
                verification_variant=variant_map_curr.get(key),
            )
            t["likely_cause"] = determine_likely_cause(t["evidence"])

        # 9. Determine verdict
        verdict = determine_verdict(transitions, envelope_violations)

        # 10. Write results to DB
        drift_count = sum(1 for t in transitions if t.get("diff_class") == "numerical_drift")
        bio_count = sum(1 for t in transitions if t.get("diff_class") == "biological_difference")

        vr = VerificationResult(
            run_id=current_run.id,
            baseline_id=baseline.id,
            verdict=verdict,
            summary=f"{len(transitions)} transitions: {drift_count} drift, {bio_count} biological",
            total_transitions=len(transitions),
            equivalent_count=0,
            drift_count=drift_count,
            biological_count=bio_count,
            envelope_breaches=len(envelope_violations),
        )
        session.add(vr)
        session.flush()

        for t in transitions:
            vt = VariantTransition(
                verification_id=vr.id,
                chrom=t["chrom"],
                pos=t["pos"],
                ref=t["ref"],
                alt=t["alt"],
                baseline_classification=t.get("baseline_classification"),
                verification_classification=t.get("verification_classification"),
                transition=t["transition"],
                diff_class=t.get("diff_class", "equivalent"),
                strata_json=json.dumps(t.get("strata", [])),
            )
            session.add(vt)
            session.flush()

            for e in t.get("evidence", []):
                rce = RootCauseEvidence(
                    transition_id=vt.id,
                    category=e["category"],
                    baseline_value=e.get("baseline_value"),
                    verification_value=e.get("verification_value"),
                    change_description=e.get("change_description", ""),
                    change_magnitude=e.get("change_magnitude"),
                    score=e.get("score", "NONE"),
                )
                session.add(rce)

        # 11. Optional BRAM assessment
        bram_result = None
        if bram_config and bram_config.get("enabled", False):
            try:
                from analysis.bayesian_risk import BRAMEngine

                bram_engine = BRAMEngine(config=bram_config)

                # Build strata_data from stratified metrics
                bram_strata_data = []
                baseline_metrics = baseline_mgr.get_baseline_metrics(baseline.id)
                for sm in strat_metrics:
                    dim = sm.get("dimension", "")
                    strat = sm.get("stratum", "")
                    # Find matching baseline metric
                    bl_metric = next(
                        (
                            bm
                            for bm in baseline_metrics
                            if bm.get("dimension") == dim
                            and bm.get("stratum") == strat
                        ),
                        None,
                    )
                    if bl_metric is None:
                        continue
                    tp = sm.get("tp_count", 0)
                    fp = sm.get("fp_count", 0)
                    fn = sm.get("fn_count", 0)
                    for metric_name in ("recall", "precision"):
                        bl_val = bl_metric.get(metric_name)
                        ver_val = sm.get(metric_name)
                        n = sm.get("total_variants", 0)
                        if bl_val is not None and ver_val is not None and n > 0:
                            # Pass raw counts to avoid rounding from proportions
                            if metric_name == "recall":
                                k, t = tp, tp + fn
                            else:
                                k, t = tp, tp + fp
                            entry = {
                                "dimension": dim,
                                "stratum": strat,
                                "metric_name": metric_name,
                                "baseline_value": bl_val,
                                "verification_value": ver_val,
                                "sample_size": n,
                                "tier": "tier_3",
                            }
                            if t > 0:
                                entry["successes"] = k
                                entry["trials"] = t
                            bram_strata_data.append(entry)

                bram_result = bram_engine.assess_all(bram_strata_data)
                log.info("BRAM assessment complete: verdict=%s", bram_result.get("verdict"))
            except Exception:
                log.exception("BRAM assessment failed, continuing without it")

        session.commit()
        log.info("Verification complete: verdict=%s", verdict)

        return {
            "verdict": verdict,
            "total_transitions": len(transitions),
            "drift_count": drift_count,
            "biological_count": bio_count,
            "envelope_breaches": len(envelope_violations),
            "bram_result": bram_result,
        }

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

    import yaml

    strat_config = {}
    strat_config_path = snakemake.params.get("stratification_config", "")
    if strat_config_path and Path(strat_config_path).exists():
        with open(strat_config_path) as f:
            strat_config = yaml.safe_load(f) or {}

    bram_config = {}
    bram_config_path = snakemake.params.get("bram_config", "")
    if bram_config_path and Path(bram_config_path).exists():
        with open(bram_config_path) as f:
            bram_config = yaml.safe_load(f) or {}

    result = run_verification(
        database=snakemake.params.database,
        baseline_name=snakemake.params.baseline_name,
        mode=snakemake.params.get("execution_mode", "germline"),
        stratification_config=strat_config,
        bram_config=bram_config,
    )

    # Touch sentinel
    Path(snakemake.output.done).touch()
    log.info("Verification done: %s", result)

except NameError:
    pass  # Not running via Snakemake
