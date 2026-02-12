"""Shared fixtures for parser tests."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable

import pytest

# Add scripts directory to path
SCRIPTS_DIR = str(Path(__file__).resolve().parent.parent.parent / "scripts")
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)


@pytest.fixture
def temp_tsv_file(tmp_path) -> Callable:
    """Factory fixture for creating temporary TSV files.

    Returns a function that creates a TSV file with specified data.

    Example:
        >>> tsv_path = temp_tsv_file(
        ...     headers=["chrom", "pos", "type"],
        ...     rows=[["chr1", "100", "SNP"], ["chr2", "200", "INDEL"]]
        ... )
    """
    def _create_tsv(headers: list[str], rows: list[list[str]]) -> Path:
        tsv_path = tmp_path / "test.tsv"

        with open(tsv_path, 'w') as f:
            # Write header
            f.write('\t'.join(headers) + '\n')
            # Write rows
            for row in rows:
                f.write('\t'.join(str(val) for val in row) + '\n')

        return tsv_path

    return _create_tsv


@pytest.fixture
def temp_json_file(tmp_path) -> Callable:
    """Factory fixture for creating temporary JSON files.

    Returns a function that creates a JSON file with specified data.

    Example:
        >>> json_path = temp_json_file({"key": "value", "num": 42})
    """
    def _create_json(data: dict) -> Path:
        json_path = tmp_path / "test.json"

        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)

        return json_path

    return _create_json


@pytest.fixture
def sample_multiqc_data() -> dict:
    """Sample nf-core MultiQC JSON data structure.

    Returns a realistic but minimal MultiQC data structure
    for testing parse_nfcore_qc.py.
    """
    return {
        "report_general_stats_data": [
            {
                "NA12878": {
                    "mosdepth_mean_coverage": 35.2,
                    "mosdepth_median_coverage": 34.8,
                    "mosdepth_10x": 99.1,
                    "mosdepth_30x": 95.3,
                    "picard_PCT_PF_READS_ALIGNED": 99.8,
                    "picard_MEDIAN_INSERT_SIZE": 350.0,
                    "picard_MEAN_INSERT_SIZE": 365.2,
                    "picard_PERCENT_DUPLICATION": 12.5,
                    "samtools_mapped_passed": 45000000,
                    "samtools_raw_total_sequences": 50000000,
                    "fastqc_avg_sequence_length": 150,
                    "fastqc_percent_gc": 42.3
                }
            }
        ],
        "report_plot_data": {
            "picard_gcbias": {
                "NA12878": {
                    "NORMALIZED_COVERAGE": [0.95, 1.0, 1.05, 1.0, 0.98],
                    "GC": [20, 30, 40, 50, 60]
                }
            },
            "mosdepth-coverage-dist-cov": {
                "NA12878": {
                    "0": 0.5,
                    "1": 2.1,
                    "10": 99.1,
                    "30": 95.3,
                    "50": 85.2
                }
            }
        }
    }


@pytest.fixture
def sample_happy_summary_data() -> list[list[str]]:
    """Sample hap.py summary metrics data.

    Returns realistic hap.py summary.csv data as rows.
    """
    headers = [
        "Type", "Filter", "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
        "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK", "FP.gt", "FP.al",
        "METRIC.Recall", "METRIC.Precision", "METRIC.Frac_NA",
        "METRIC.F1_Score", "TRUTH.TOTAL.TiTv_ratio", "QUERY.TOTAL.TiTv_ratio"
    ]

    rows = [
        ["INDEL", "ALL", "500", "475", "25", "520", "30", "15", "10", "20",
         "0.95", "0.94", "0.029", "0.945", ".", "."],
        ["SNP", "ALL", "3000", "2950", "50", "3100", "80", "70", "30", "50",
         "0.983", "0.974", "0.023", "0.978", "2.1", "2.05"]
    ]

    return [headers] + rows


@pytest.fixture
def sample_rtg_summary_data() -> list[list[str]]:
    """Sample RTG vcfeval weighted_roc.tsv data.

    Returns realistic RTG metrics data as rows.
    """
    headers = [
        "#score", "true_positives_baseline", "false_positives",
        "false_negatives", "precision", "sensitivity", "f_measure"
    ]

    rows = [
        ["0.0", "2950", "80", "50", "0.974", "0.983", "0.978"],
        ["10.0", "2920", "45", "80", "0.985", "0.973", "0.979"],
        ["20.0", "2850", "20", "150", "0.993", "0.950", "0.971"]
    ]

    return [headers] + rows
