"""Report storage and naming convention (Spec S16.4).

Provides naming conventions, checksum computation, and storage
utilities for validation reports.
"""

from __future__ import annotations

import hashlib
import logging
from datetime import UTC, datetime
from pathlib import Path

log = logging.getLogger(__name__)


def generate_report_name(
    sample: str,
    mode: str,
    report_type: str,
    extension: str = "json",
    date: datetime | None = None,
) -> str:
    """Generate a report filename following the naming convention.

    Format: {sample}_{mode}_{date}_{report_type}.{ext}

    Args:
        sample: Sample identifier
        mode: Execution mode (germline/somatic)
        report_type: Report type (initial_validation/continuous_verification)
        extension: File extension (json/html/pdf)
        date: Report date (defaults to now)

    Returns:
        Formatted filename string
    """
    dt = date or datetime.now(UTC)
    date_str = dt.strftime("%Y%m%d_%H%M%S")

    # Sanitize sample name
    safe_sample = sample.replace(" ", "_").replace("/", "_")

    return f"{safe_sample}_{mode}_{date_str}_{report_type}.{extension}"


def compute_file_checksum(file_path: str) -> str:
    """Compute SHA-256 checksum of a file.

    Args:
        file_path: Path to the file

    Returns:
        Hex digest string
    """
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        while True:
            chunk = f.read(8192)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def save_report(
    content: str,
    output_dir: str,
    filename: str,
) -> str:
    """Save a report to disk and return its path.

    Args:
        content: Report content (JSON or HTML string)
        output_dir: Directory to save the report in
        filename: Report filename

    Returns:
        Full path to the saved report
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    file_path = out_dir / filename
    file_path.write_text(content)

    checksum = compute_file_checksum(str(file_path))
    log.info("Saved report: %s (sha256=%s)", file_path, checksum[:12])

    return str(file_path)
