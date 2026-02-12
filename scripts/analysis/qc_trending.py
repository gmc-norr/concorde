"""QC metric trending: z-score analysis against historical runs.

Compares current run metrics against historical data to detect
outliers. Flags metrics with |z-score| > threshold (default 2.0).

Uses scipy.stats for z-score computation.
"""

from __future__ import annotations

import logging

from scipy.stats import zscore as _scipy_zscore  # noqa: F401 (validates import)

log = logging.getLogger(__name__)


def compute_historical_stats(
    session,
    metric_name: str,
    current_run_id: int,
    min_runs: int = 5,
) -> dict | None:
    """Compute mean and stdev from historical runs for a given metric.

    Args:
        session: SQLAlchemy session
        metric_name: Name of the QC metric
        current_run_id: Current run ID (excluded from historical stats)
        min_runs: Minimum number of historical runs required

    Returns:
        Dict with mean, stdev, count, or None if insufficient history
    """
    from models.qc_metric import QCMetric

    historical = (
        session.query(QCMetric.metric_value_float)
        .filter(
            QCMetric.metric_name == metric_name,
            QCMetric.run_id != current_run_id,
            QCMetric.metric_value_float.isnot(None),
        )
        .all()
    )

    values = [row[0] for row in historical if row[0] is not None]

    if len(values) < min_runs:
        log.info(
            "Insufficient history for %s: %d runs (min %d)",
            metric_name, len(values), min_runs,
        )
        return None

    import numpy as np

    arr = np.array(values, dtype=np.float64)
    mean = float(arr.mean())
    stdev = float(arr.std(ddof=1))  # Sample stdev (Bessel's correction)

    return {
        "mean": mean,
        "stdev": stdev,
        "count": len(values),
    }


def compute_trending_report(
    session,
    run_id: int,
    config: dict | None = None,
) -> list[dict]:
    """Compute QC trending z-scores for all metrics in a run.

    Args:
        session: SQLAlchemy session
        run_id: Current run ID
        config: Optional config with min_historical_runs, z_score_threshold

    Returns:
        List of dicts with metric_name, current_value, historical_mean,
        historical_stdev, z_score, flagged
    """
    from models.qc_metric import QCMetric

    config = config or {}
    min_runs = config.get("min_historical_runs", 5)
    z_threshold = config.get("z_score_threshold", 2.0)

    # Get all metrics for the current run
    current_metrics = (
        session.query(QCMetric)
        .filter(
            QCMetric.run_id == run_id,
            QCMetric.metric_value_float.isnot(None),
        )
        .all()
    )

    report = []
    for qm in current_metrics:
        stats = compute_historical_stats(
            session, qm.metric_name, run_id, min_runs=min_runs
        )

        if stats is None:
            report.append({
                "metric_name": qm.metric_name,
                "current_value": qm.metric_value_float,
                "historical_mean": None,
                "historical_stdev": None,
                "z_score": None,
                "flagged": False,
                "insufficient_history": True,
            })
            continue

        # Compute z-score
        if stats["stdev"] > 0:
            z_score = (qm.metric_value_float - stats["mean"]) / stats["stdev"]
        else:
            z_score = 0.0 if qm.metric_value_float == stats["mean"] else float("inf")

        flagged = abs(z_score) > z_threshold

        report.append({
            "metric_name": qm.metric_name,
            "current_value": qm.metric_value_float,
            "historical_mean": stats["mean"],
            "historical_stdev": stats["stdev"],
            "z_score": z_score,
            "flagged": flagged,
            "insufficient_history": False,
        })

        if flagged:
            log.warning(
                "QC ALERT: %s z-score=%.2f (current=%.4f, mean=%.4f, stdev=%.4f)",
                qm.metric_name, z_score, qm.metric_value_float,
                stats["mean"], stats["stdev"],
            )

    log.info(
        "QC trending: %d metrics analyzed, %d flagged",
        len(report),
        sum(1 for r in report if r["flagged"]),
    )
    return report
