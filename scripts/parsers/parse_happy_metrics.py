"""Parse hap.py extended.csv into a metrics TSV.

Called via Snakemake `script:` directive.
Reads the hap.py extended CSV and produces a tab-separated file with
precision, recall, F1, and TP/FP/FN counts per variant type and stratification.
"""

import logging
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


def parse_extended_csv(csv_path):
    """Parse hap.py extended.csv into metric records."""
    df = pd.read_csv(csv_path)
    records = []

    for _, row in df.iterrows():
        filter_val = row.get("Filter", "ALL")
        variant_type = row.get("Type", "UNKNOWN")
        subset = row.get("Subset", "*")
        subtype = row.get("Subtype", "*")

        # Only keep ALL and PASS filter rows
        if filter_val not in ("ALL", "PASS"):
            continue

        # Build stratification string
        strat_parts = []
        if filter_val == "PASS":
            strat_parts.append("filter=PASS")
        if pd.notna(subset) and subset != "*":
            strat_parts.append(f"subset={subset}")
        if pd.notna(subtype) and subtype != "*":
            strat_parts.append(f"subtype={subtype}")
        stratification = ";".join(strat_parts) if strat_parts else None

        records.append(
            {
                "variant_type": variant_type,
                "stratification": stratification,
                "precision": safe_float(row.get("METRIC.Precision", 0.0)),
                "recall": safe_float(row.get("METRIC.Recall", 0.0)),
                "f1": safe_float(row.get("METRIC.F1_Score", 0.0)),
                "tp_count": safe_int(row.get("TRUTH.TP")),
                "fp_count": safe_int(row.get("QUERY.FP")),
                "fn_count": safe_int(row.get("TRUTH.FN")),
            }
        )

    return records


# --- Main execution via Snakemake ---
log.info("Parsing hap.py extended CSV: %s", snakemake.input.csv)

records = parse_extended_csv(snakemake.input.csv)

df = pd.DataFrame(records)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

log.info("Wrote %d metric records to %s", len(df), snakemake.output.tsv)
