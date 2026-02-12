"""Setup import paths for Snakemake workflow scripts.

This module ensures that all workflow scripts can import from workflow/lib
regardless of where they are executed from.
"""

import sys
from pathlib import Path


def setup_workflow_paths():
    """Add workflow/lib to sys.path for imports."""
    # Get the workflow/lib directory
    lib_dir = Path(__file__).resolve().parent

    # Add to sys.path if not already present
    lib_str = str(lib_dir)
    if lib_str not in sys.path:
        sys.path.insert(0, lib_str)

    return lib_dir


# Auto-setup when imported
setup_workflow_paths()
