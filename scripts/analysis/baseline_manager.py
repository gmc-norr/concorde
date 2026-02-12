"""Baseline management for initial validation (Spec S12).

Creates, locks, and queries baselines. Auto-derives performance envelopes
from stratified metrics with configurable tolerance.
"""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import UTC, datetime

from sqlalchemy.orm import Session

log = logging.getLogger(__name__)

DEFAULT_TOLERANCE = 0.005  # ±0.5% for auto-derived envelopes
ENVELOPE_METRICS = ("precision", "recall", "f1")


class BaselineManager:
    """Manages baseline creation, envelope computation, and signing.

    Args:
        session: SQLAlchemy database session
    """

    def __init__(self, session: Session):
        self.session = session

    def create_baseline(
        self,
        name: str,
        run_id: int,
        config_snapshot: str | None = None,
        stratification_snapshot: str | None = None,
        pipeline_version: str | None = None,
        mode: str = "germline",
        tolerance: float = DEFAULT_TOLERANCE,
    ) -> "Baseline":
        """Create a new baseline from a completed run.

        Args:
            name: Unique baseline name (e.g. "v1.0_germline_GRCh38")
            run_id: ID of the run that produced the baseline metrics
            config_snapshot: YAML config snapshot (serialized string)
            stratification_snapshot: YAML stratification config snapshot
            pipeline_version: Pipeline version string
            mode: Execution mode ("germline" or "somatic")
            tolerance: Default tolerance for auto-derived envelopes

        Returns:
            Created Baseline object

        Raises:
            ValueError: If name already exists or run not found
        """
        from models.baseline import Baseline

        # Check uniqueness
        existing = (
            self.session.query(Baseline).filter(Baseline.name == name).first()
        )
        if existing:
            raise ValueError(f"Baseline with name '{name}' already exists")

        # Verify run exists
        from models.run import Run

        run = self.session.query(Run).filter(Run.id == run_id).first()
        if not run:
            raise ValueError(f"Run with id {run_id} not found")

        baseline = Baseline(
            name=name,
            run_id=run_id,
            config_snapshot=config_snapshot,
            stratification_snapshot=stratification_snapshot,
            pipeline_version=pipeline_version or run.pipeline_version,
            mode=mode,
        )
        self.session.add(baseline)
        self.session.flush()

        # Auto-derive envelopes from stratified metrics
        envelope_count = self._auto_derive_envelopes(
            baseline.id, run_id, tolerance
        )

        log.info(
            "Created baseline '%s' (id=%d) with %d envelopes from run %d",
            name,
            baseline.id,
            envelope_count,
            run_id,
        )
        return baseline

    def _auto_derive_envelopes(
        self, baseline_id: int, run_id: int, tolerance: float
    ) -> int:
        """Auto-derive performance envelopes from stratified metrics.

        For each metric in each stratum, creates an envelope with:
        - expected_value = observed value
        - lower_bound = expected_value - tolerance
        - upper_bound = expected_value + tolerance (clamped to 1.0)

        Args:
            baseline_id: Baseline to attach envelopes to
            run_id: Run containing stratified metrics
            tolerance: Absolute tolerance for bounds

        Returns:
            Number of envelopes created
        """
        from models.baseline import BaselineEnvelope
        from models.stratified_metric import StratifiedMetric

        strat_metrics = (
            self.session.query(StratifiedMetric)
            .filter(StratifiedMetric.run_id == run_id)
            .all()
        )

        envelopes = []
        for sm in strat_metrics:
            for metric_name in ENVELOPE_METRICS:
                value = getattr(sm, metric_name, None)
                if value is None:
                    continue

                envelope = BaselineEnvelope(
                    baseline_id=baseline_id,
                    dimension=sm.dimension,
                    stratum=sm.stratum,
                    metric_name=metric_name,
                    expected_value=value,
                    lower_bound=max(0.0, value - tolerance),
                    upper_bound=min(1.0, value + tolerance),
                    manually_set=False,
                )
                envelopes.append(envelope)

        if envelopes:
            self.session.add_all(envelopes)

        return len(envelopes)

    def set_envelope(
        self,
        baseline_id: int,
        dimension: str,
        stratum: str,
        metric_name: str,
        expected_value: float,
        lower_bound: float | None = None,
        upper_bound: float | None = None,
    ) -> "BaselineEnvelope":
        """Manually set or override an envelope for a specific stratum/metric.

        Args:
            baseline_id: Baseline to update
            dimension: Stratification dimension
            stratum: Stratum value
            metric_name: Metric name (precision, recall, f1)
            expected_value: Expected metric value
            lower_bound: Minimum acceptable value
            upper_bound: Maximum acceptable value

        Returns:
            Created or updated BaselineEnvelope

        Raises:
            ValueError: If baseline is locked or not found
        """
        from models.baseline import Baseline, BaselineEnvelope

        baseline = (
            self.session.query(Baseline)
            .filter(Baseline.id == baseline_id)
            .first()
        )
        if not baseline:
            raise ValueError(f"Baseline {baseline_id} not found")
        if baseline.locked:
            raise ValueError(
                f"Baseline '{baseline.name}' is locked and cannot be modified"
            )

        # Find existing or create new
        existing = (
            self.session.query(BaselineEnvelope)
            .filter(
                BaselineEnvelope.baseline_id == baseline_id,
                BaselineEnvelope.dimension == dimension,
                BaselineEnvelope.stratum == stratum,
                BaselineEnvelope.metric_name == metric_name,
            )
            .first()
        )

        if existing:
            existing.expected_value = expected_value
            existing.lower_bound = lower_bound
            existing.upper_bound = upper_bound
            existing.manually_set = True
            return existing

        envelope = BaselineEnvelope(
            baseline_id=baseline_id,
            dimension=dimension,
            stratum=stratum,
            metric_name=metric_name,
            expected_value=expected_value,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            manually_set=True,
        )
        self.session.add(envelope)
        return envelope

    def lock_baseline(
        self,
        baseline_id: int,
        approver: str,
        comment: str | None = None,
    ) -> "Baseline":
        """Lock and sign a baseline (S12.3).

        Args:
            baseline_id: Baseline to lock
            approver: Approver identity
            comment: Approval justification

        Returns:
            Locked Baseline object

        Raises:
            ValueError: If baseline not found or already locked
        """
        from models.baseline import Baseline

        baseline = (
            self.session.query(Baseline)
            .filter(Baseline.id == baseline_id)
            .first()
        )
        if not baseline:
            raise ValueError(f"Baseline {baseline_id} not found")
        if baseline.locked:
            raise ValueError(f"Baseline '{baseline.name}' is already locked")

        # Compute artifact hash over the baseline content
        artifact_hash = self._compute_artifact_hash(baseline)

        baseline.locked = True
        baseline.locked_at = datetime.now(UTC)
        baseline.approver = approver
        baseline.approval_comment = comment
        baseline.artifact_hash = artifact_hash

        log.info(
            "Locked baseline '%s' (hash=%s, approver=%s)",
            baseline.name,
            artifact_hash[:12],
            approver,
        )
        return baseline

    def _compute_artifact_hash(self, baseline: "Baseline") -> str:
        """Compute SHA-256 hash of the baseline artifact content.

        Includes: name, config, envelopes, pipeline version.
        """
        content = {
            "name": baseline.name,
            "run_id": baseline.run_id,
            "mode": baseline.mode,
            "pipeline_version": baseline.pipeline_version,
            "config_snapshot": baseline.config_snapshot,
            "stratification_snapshot": baseline.stratification_snapshot,
            "envelopes": [
                {
                    "dimension": e.dimension,
                    "stratum": e.stratum,
                    "metric_name": e.metric_name,
                    "expected_value": e.expected_value,
                    "lower_bound": e.lower_bound,
                    "upper_bound": e.upper_bound,
                }
                for e in sorted(
                    baseline.envelopes,
                    key=lambda e: (e.dimension, e.stratum, e.metric_name),
                )
            ],
        }
        serialized = json.dumps(content, sort_keys=True, default=str)
        return hashlib.sha256(serialized.encode()).hexdigest()

    def get_baseline(self, name: str) -> "Baseline | None":
        """Look up a baseline by name.

        Args:
            name: Baseline name

        Returns:
            Baseline object or None
        """
        from models.baseline import Baseline

        return (
            self.session.query(Baseline).filter(Baseline.name == name).first()
        )

    def check_envelopes(
        self, baseline_id: int, current_metrics: list[dict]
    ) -> list[dict]:
        """Check current metrics against baseline envelopes.

        Args:
            baseline_id: Baseline to compare against
            current_metrics: List of metric dicts with dimension, stratum, and metric values

        Returns:
            List of violation dicts with dimension, stratum, metric_name,
            current_value, expected_value, lower_bound, upper_bound
        """
        from models.baseline import BaselineEnvelope

        envelopes = (
            self.session.query(BaselineEnvelope)
            .filter(BaselineEnvelope.baseline_id == baseline_id)
            .all()
        )

        # Build lookup
        env_map = {}
        for e in envelopes:
            key = (e.dimension, e.stratum, e.metric_name)
            env_map[key] = e

        violations = []
        for m in current_metrics:
            dim = m.get("dimension")
            stratum = m.get("stratum")
            for metric_name in ENVELOPE_METRICS:
                value = m.get(metric_name)
                if value is None:
                    continue

                key = (dim, stratum, metric_name)
                envelope = env_map.get(key)
                if envelope is None:
                    continue

                below = envelope.lower_bound is not None and value < envelope.lower_bound
                above = envelope.upper_bound is not None and value > envelope.upper_bound

                if below or above:
                    violations.append(
                        {
                            "dimension": dim,
                            "stratum": stratum,
                            "metric_name": metric_name,
                            "current_value": value,
                            "expected_value": envelope.expected_value,
                            "lower_bound": envelope.lower_bound,
                            "upper_bound": envelope.upper_bound,
                            "violation": "below_lower_bound" if below else "above_upper_bound",
                        }
                    )

        return violations
