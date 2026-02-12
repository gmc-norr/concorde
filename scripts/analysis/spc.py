"""Statistical Process Control (SPC) for metric drift detection.

Computes control limits (mean +/- N*sigma) from historical runs and
evaluates whether the current run's metrics are in statistical control.

Uses numpy for mean/stdev computation with Bessel's correction.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


class SPCEngine:
    """Statistical Process Control engine for tracking metric drift.

    Args:
        session: SQLAlchemy session
        config: SPC config dict with sigma_multiplier, min_data_points
    """

    def __init__(self, session, config=None):
        self.session = session
        self.config = config or {}
        self.sigma_multiplier = self.config.get("sigma_multiplier", 3.0)
        self.min_data_points = self.config.get("min_data_points", 10)

    def compute_control_limits(
        self,
        metric_name: str,
        dimension: str | None = None,
        stratum: str | None = None,
    ) -> dict | None:
        """Compute control limits for a metric from historical data.

        Args:
            metric_name: Name of the metric (e.g., "precision", "recall")
            dimension: Optional stratification dimension filter
            stratum: Optional stratum value filter

        Returns:
            Dict with mean, stdev, ucl, lcl, data_points, or None if insufficient data
        """
        from models.stratified_metric import StratifiedMetric

        query = self.session.query(StratifiedMetric)

        if dimension:
            query = query.filter(StratifiedMetric.dimension == dimension)
        if stratum:
            query = query.filter(StratifiedMetric.stratum == stratum)

        # Get the metric values
        metrics = query.all()
        values = []
        for m in metrics:
            val = getattr(m, metric_name, None)
            if val is not None:
                values.append(val)

        if len(values) < self.min_data_points:
            return None

        import numpy as np

        arr = np.array(values, dtype=np.float64)
        mean = float(arr.mean())
        stdev = float(arr.std(ddof=1))  # Sample stdev (Bessel's correction)

        ucl = mean + self.sigma_multiplier * stdev
        lcl = mean - self.sigma_multiplier * stdev

        return {
            "mean": mean,
            "stdev": stdev,
            "ucl": ucl,
            "lcl": lcl,
            "data_points": len(values),
        }

    def evaluate_run(self, run_id: int) -> list[dict]:
        """Evaluate all stratified metrics for a run against control limits.

        Args:
            run_id: Run ID to evaluate

        Returns:
            List of dicts with metric_name, dimension, stratum, value,
            mean, ucl, lcl, in_control
        """
        from models.stratified_metric import StratifiedMetric

        current_metrics = (
            self.session.query(StratifiedMetric)
            .filter(StratifiedMetric.run_id == run_id)
            .all()
        )

        results = []
        for sm in current_metrics:
            for metric_name in ("precision", "recall", "f1"):
                value = getattr(sm, metric_name, None)
                if value is None:
                    continue

                limits = self.compute_control_limits(
                    metric_name,
                    dimension=sm.dimension,
                    stratum=sm.stratum,
                )

                if limits is None:
                    results.append({
                        "metric_name": metric_name,
                        "dimension": sm.dimension,
                        "stratum": sm.stratum,
                        "value": value,
                        "mean": None,
                        "ucl": None,
                        "lcl": None,
                        "in_control": True,
                        "insufficient_data": True,
                    })
                    continue

                in_control = limits["lcl"] <= value <= limits["ucl"]

                results.append({
                    "metric_name": metric_name,
                    "dimension": sm.dimension,
                    "stratum": sm.stratum,
                    "value": value,
                    "mean": limits["mean"],
                    "ucl": limits["ucl"],
                    "lcl": limits["lcl"],
                    "in_control": in_control,
                    "insufficient_data": False,
                })

                if not in_control:
                    log.warning(
                        "SPC OUT OF CONTROL: %s [%s/%s] = %.4f "
                        "(UCL=%.4f, LCL=%.4f)",
                        metric_name, sm.dimension, sm.stratum, value,
                        limits["ucl"], limits["lcl"],
                    )

        log.info(
            "SPC evaluation: %d metrics checked, %d out of control",
            len(results),
            sum(1 for r in results if not r["in_control"]),
        )
        return results

    def generate_chart_data(
        self,
        metric_name: str,
        dimension: str | None = None,
        stratum: str | None = None,
        n_runs: int = 50,
    ) -> dict:
        """Generate data for SPC control chart visualization.

        Args:
            metric_name: Name of the metric
            dimension: Optional dimension filter
            stratum: Optional stratum filter
            n_runs: Maximum number of recent runs to include

        Returns:
            Dict with run_ids, values, mean, ucl, lcl, flagged_runs
        """
        from models.stratified_metric import StratifiedMetric

        query = self.session.query(StratifiedMetric)

        if dimension:
            query = query.filter(StratifiedMetric.dimension == dimension)
        if stratum:
            query = query.filter(StratifiedMetric.stratum == stratum)

        metrics = query.order_by(StratifiedMetric.run_id.desc()).limit(n_runs).all()
        metrics.reverse()  # Chronological order

        run_ids = []
        values = []
        for m in metrics:
            val = getattr(m, metric_name, None)
            if val is not None:
                run_ids.append(m.run_id)
                values.append(val)

        limits = self.compute_control_limits(metric_name, dimension, stratum)
        mean = limits["mean"] if limits else None
        ucl = limits["ucl"] if limits else None
        lcl = limits["lcl"] if limits else None

        flagged = []
        if limits:
            for run_id, val in zip(run_ids, values):
                if val < limits["lcl"] or val > limits["ucl"]:
                    flagged.append(run_id)

        return {
            "run_ids": run_ids,
            "values": values,
            "mean": mean,
            "ucl": ucl,
            "lcl": lcl,
            "flagged_runs": flagged,
        }
