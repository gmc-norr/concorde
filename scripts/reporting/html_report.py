"""HTML report generation (Spec S16.2).

Renders a Jinja2-based HTML report with 7 sections: header, executive
summary, per-stratum metrics, variant details, root cause analysis,
traceability, and appendices.
"""

from __future__ import annotations

import logging
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

log = logging.getLogger(__name__)

TEMPLATES_DIR = Path(__file__).parent / "templates"


def generate_html_report(
    report_data: dict,
    template_name: str = "report.html.j2",
    templates_dir: str | Path | None = None,
) -> str:
    """Generate an HTML report from a JSON report dict.

    Args:
        report_data: JSON report dict from generate_json_report()
        template_name: Jinja2 template file name
        templates_dir: Override templates directory path

    Returns:
        Rendered HTML string
    """
    tpl_dir = Path(templates_dir) if templates_dir else TEMPLATES_DIR
    env = Environment(
        loader=FileSystemLoader(str(tpl_dir)),
        autoescape=True,
    )

    # Add custom filters
    env.filters["format_pct"] = _format_pct
    env.filters["format_float"] = _format_float
    env.filters["status_color"] = _status_color

    template = env.get_template(template_name)
    html = template.render(report=report_data)

    log.info("Generated HTML report (%d chars)", len(html))
    return html


def _format_pct(value: float | None) -> str:
    """Format a float as a percentage."""
    if value is None:
        return "N/A"
    return f"{value * 100:.2f}%"


def _format_float(value: float | None, decimals: int = 4) -> str:
    """Format a float to N decimal places."""
    if value is None:
        return "N/A"
    return f"{value:.{decimals}f}"


def _status_color(status: str) -> str:
    """Map status to CSS color class."""
    return {
        "PASS": "green",
        "FAIL": "red",
        "REVIEW_REQUIRED": "amber",
        "CONDITIONAL_PASS": "amber",
        "WARN": "amber",
        "EXCLUDED": "gray",
        "BRAM_PASS": "green",
        "BRAM_FLAG": "red",
        "BRAM_NOT_RUN": "gray",
    }.get(status, "gray")
