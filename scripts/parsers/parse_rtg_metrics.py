"""Parse rtg vcfeval summary.txt into a metrics TSV.

Called via Snakemake `script:` directive.
Reads the rtg vcfeval summary.txt file and produces a tab-separated file with
precision, recall, F1, and TP/FP/FN counts per variant type.

rtg vcfeval summary.txt format example:
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
--------   -----------------  -------------  ---------  ---------  ---------  -----------  ---------
None       7043               7043           123        456        0.9828     0.9392       0.9605
"""

import logging
import re
from typing import TYPE_CHECKING

import pandas as pd

from utils import safe_float, safe_int

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
log = logging.getLogger(__name__)


def parse_rtg_summary(summary_path):
    """Parse rtg vcfeval summary.txt file into metric records."""
    records = []

    with open(summary_path) as f:
        content = f.read()

    # rtg vcfeval produces sections for different variant types
    # Look for lines with variant type headers (e.g., "SNP", "INDEL", "None")
    # The summary typically has a table format

    # Split by lines and look for data rows
    lines = content.strip().split("\n")

    # Find header line (contains "Threshold", "True-pos", etc.)
    header_idx = None
    for i, line in enumerate(lines):
        if "Threshold" in line and "True-pos" in line:
            header_idx = i
            break

    if header_idx is None:
        log.warning("Could not find header line in summary.txt")
        return records

    # Skip the separator line (dashes)
    data_start_idx = header_idx + 2

    # Parse data lines
    for line in lines[data_start_idx:]:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Split by whitespace
        parts = line.split()
        if len(parts) < 8:
            continue

        threshold = parts[0]
        # tp_baseline = safe_int(parts[1])  # Not used, tp_call is sufficient
        tp_call = safe_int(parts[2])
        fp = safe_int(parts[3])
        fn = safe_int(parts[4])
        precision = safe_float(parts[5])
        sensitivity = safe_float(parts[6])  # recall
        f_measure = safe_float(parts[7])  # F1

        # Determine variant type from threshold
        # rtg uses "None" for overall, or may include variant type prefixes
        variant_type = "ALL"
        stratification = None

        # Check if this is a filtered or stratified result
        if threshold != "None":
            stratification = f"threshold={threshold}"

        records.append(
            {
                "variant_type": variant_type,
                "stratification": stratification,
                "precision": precision,
                "recall": sensitivity,
                "f1": f_measure,
                "tp_count": tp_call,
                "fp_count": fp,
                "fn_count": fn,
            }
        )

    # rtg vcfeval may also produce separate sections for SNP, INDEL
    # Try to extract variant-type-specific metrics
    snp_match = re.search(
        r"SNP\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        content,
    )
    if snp_match:
        records.append(
            {
                "variant_type": "SNP",
                "stratification": None,
                "precision": safe_float(snp_match.group(5)),
                "recall": safe_float(snp_match.group(6)),
                "f1": safe_float(snp_match.group(7)),
                "tp_count": safe_int(snp_match.group(2)),
                "fp_count": safe_int(snp_match.group(3)),
                "fn_count": safe_int(snp_match.group(4)),
            }
        )

    indel_match = re.search(
        r"INDEL\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        content,
    )
    if indel_match:
        records.append(
            {
                "variant_type": "INDEL",
                "stratification": None,
                "precision": safe_float(indel_match.group(5)),
                "recall": safe_float(indel_match.group(6)),
                "f1": safe_float(indel_match.group(7)),
                "tp_count": safe_int(indel_match.group(2)),
                "fp_count": safe_int(indel_match.group(3)),
                "fn_count": safe_int(indel_match.group(4)),
            }
        )

    return records


# --- Main execution via Snakemake ---
log.info("Parsing rtg vcfeval summary: %s", snakemake.input.summary)

records = parse_rtg_summary(snakemake.input.summary)

df = pd.DataFrame(records)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

log.info("Wrote %d metric records to %s", len(df), snakemake.output.tsv)
