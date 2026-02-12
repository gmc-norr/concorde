"""PDF report generation (Spec S16.3).

Converts HTML reports to PDF using weasyprint (if available).
Falls back gracefully if weasyprint is not installed.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)


def generate_pdf_report(html_content: str, output_path: str) -> bool:
    """Generate a PDF from an HTML report string.

    Args:
        html_content: Rendered HTML string
        output_path: Path to write the PDF file

    Returns:
        True if PDF was generated, False if weasyprint is unavailable
    """
    try:
        from weasyprint import HTML
    except ImportError:
        log.warning(
            "weasyprint not installed; PDF generation skipped. "
            "Install with: pip install weasyprint"
        )
        return False

    try:
        HTML(string=html_content).write_pdf(output_path)
        log.info("Generated PDF report: %s", output_path)
        return True
    except Exception as e:
        log.error("Failed to generate PDF: %s", e)
        return False
