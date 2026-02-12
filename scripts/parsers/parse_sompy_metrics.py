"""Parse som.py (hap.py somatic module) summary metrics into a metrics TSV.

Called via Snakemake `script:` directive.
Reads the som.py stats CSV output and writes a tab-separated file with
one row per variant_type/stratification combination.
"""

import logging

import pandas as pd

from utils import safe_float, safe_int

log = logging.getLogger(__name__)


def parse_sompy_stats(stats_csv_path):
    """Parse som.py stats CSV and extract metrics.

    Args:
        stats_csv_path: Path to som.py *.stats.csv file

    Returns:
        List of metric record dictionaries
    """
    try:
        df = pd.read_csv(stats_csv_path)
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        log.warning("Could not read som.py stats CSV %s: %s", stats_csv_path, e)
        return []

    records = []

    for _, row in df.iterrows():
        variant_type = str(row.get("type", "ALL"))

        # Build stratification string
        strat_parts = []
        subset = str(row.get("subtype", "*"))
        filter_val = str(row.get("filter", "ALL"))
        if filter_val == "PASS":
            strat_parts.append("filter=PASS")
        if subset != "*":
            strat_parts.append(f"subset={subset}")
        stratification = ";".join(strat_parts) if strat_parts else None

        tp = safe_int(row.get("tp"))
        fp = safe_int(row.get("fp"))
        fn = safe_int(row.get("fn"))

        recall = safe_float(row.get("recall"))
        precision = safe_float(row.get("precision"))

        # Compute F1 using shared utility for consistency
        from utils import compute_classification_metrics
        prf = compute_classification_metrics(tp or 0, fp or 0, fn or 0)
        # Use som.py's precision/recall if available, else computed values
        final_precision = precision if precision is not None else (prf["precision"] or 0.0)
        final_recall = recall if recall is not None else (prf["recall"] or 0.0)
        f1 = prf["f1"] or 0.0

        records.append(
            {
                "variant_type": variant_type,
                "stratification": stratification,
                "precision": final_precision,
                "recall": final_recall,
                "f1": f1,
                "tp_count": tp,
                "fp_count": fp,
                "fn_count": fn,
            }
        )

    return records


# --- Main execution via Snakemake ---
try:
    from typing import TYPE_CHECKING

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

    log.info("Parsing som.py stats: %s", snakemake.input.csv)

    records = parse_sompy_stats(snakemake.input.csv)

    df = pd.DataFrame(records)
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)

    log.info("Wrote %d somatic metric records to %s", len(df), snakemake.output.tsv)
except NameError:
    pass  # Not running via Snakemake (e.g., imported for testing)
