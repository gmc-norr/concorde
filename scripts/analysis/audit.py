"""Audit trail recording for pipeline operations.

Records key events (ingestion, baseline creation/locking, verification)
with actor, timestamp, and detail for compliance and traceability.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


def record_event(
    session,
    event_type: str,
    run_id: int | None = None,
    baseline_id: int | None = None,
    actor: str | None = None,
    detail: str | None = None,
):
    """Record an audit event.

    Args:
        session: SQLAlchemy session
        event_type: Type of event (e.g., "ingestion_complete", "baseline_created",
            "baseline_locked", "verification_complete")
        run_id: Optional run ID
        baseline_id: Optional baseline ID
        actor: Optional actor name (user or system)
        detail: Optional detail string

    Returns:
        The created AuditEvent object
    """
    from models.audit_event import AuditEvent

    event = AuditEvent(
        event_type=event_type,
        run_id=run_id,
        baseline_id=baseline_id,
        actor=actor or "system",
        detail=detail,
    )
    session.add(event)
    log.info(
        "Audit event: %s (run_id=%s, baseline_id=%s, actor=%s)",
        event_type, run_id, baseline_id, actor,
    )
    return event


def get_audit_trail(
    session,
    run_id: int | None = None,
    baseline_id: int | None = None,
    event_type: str | None = None,
    limit: int = 100,
) -> list:
    """Query the audit trail.

    Args:
        session: SQLAlchemy session
        run_id: Optional filter by run ID
        baseline_id: Optional filter by baseline ID
        event_type: Optional filter by event type
        limit: Maximum number of events to return

    Returns:
        List of AuditEvent objects, most recent first
    """
    from models.audit_event import AuditEvent

    query = session.query(AuditEvent)

    if run_id is not None:
        query = query.filter(AuditEvent.run_id == run_id)
    if baseline_id is not None:
        query = query.filter(AuditEvent.baseline_id == baseline_id)
    if event_type is not None:
        query = query.filter(AuditEvent.event_type == event_type)

    return query.order_by(AuditEvent.created_at.desc()).limit(limit).all()
