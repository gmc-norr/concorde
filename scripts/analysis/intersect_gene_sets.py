"""Intersect variant positions with gene set BED files using bedtools.

Called via Snakemake `script:` directive.
Reads the parsed variants TSV, converts to BED format, runs bedtools intersect
against each gene set BED, and writes a TSV mapping variant indices to gene sets.
"""

import logging
import subprocess  # nosec B404
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

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


def variants_to_bed(variants_df):
    """Convert variants DataFrame to BED format string (0-based half-open)."""
    lines = []
    for idx, row in variants_df.iterrows():
        # BED is 0-based half-open; VCF POS is 1-based
        start = int(row["pos"]) - 1
        end = int(row["pos"])
        lines.append(f"{row['chrom']}\t{start}\t{end}\t{idx}")
    return "\n".join(lines) + "\n"


def intersect_with_bedtools(variants_bed_path, gene_set_bed_path):
    """Run bedtools intersect and return list of (variant_idx, gene_name)."""
    result = subprocess.run(  # nosec B603 B607
        [
            "bedtools",
            "intersect",
            "-a",
            variants_bed_path,
            "-b",
            gene_set_bed_path,
            "-wa",
            "-wb",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.warning("bedtools intersect failed: %s", result.stderr)
        return []

    assignments = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        variant_idx = int(cols[3])
        # Gene set BED: chrom(4), start(5), end(6), gene_name(7) if present
        gene_name = cols[7] if len(cols) > 7 else None
        assignments.append((variant_idx, gene_name))
    return assignments


# --- Main execution via Snakemake ---
gene_sets_config = snakemake.params.gene_sets

variants_df = pd.read_csv(snakemake.input.variants_tsv, sep="\t")
log.info("Loaded %d variants for gene set intersection", len(variants_df))

all_assignments = []

if gene_sets_config and len(variants_df) > 0:
    bed_content = variants_to_bed(variants_df)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
        tmp.write(bed_content)
        variants_bed_path = tmp.name

    try:
        for gs in gene_sets_config:
            gs_name = gs["name"]
            gs_bed = gs["bed"]

            if not Path(gs_bed).exists():
                log.warning(
                    "Gene set BED not found, skipping: %s (%s)", gs_name, gs_bed
                )
                continue

            log.info("Intersecting with gene set: %s (%s)", gs_name, gs_bed)
            assignments = intersect_with_bedtools(variants_bed_path, gs_bed)

            for variant_idx, gene in assignments:
                all_assignments.append(
                    {
                        "variant_idx": variant_idx,
                        "gene_set_name": gs_name,
                        "gene": gene,
                    }
                )
            log.info("  Found %d overlapping variants", len(assignments))
    finally:
        Path(variants_bed_path).unlink(missing_ok=True)

result_df = pd.DataFrame(
    all_assignments,
    columns=["variant_idx", "gene_set_name", "gene"],
)
result_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

log.info("Wrote %d gene set assignments to %s", len(result_df), snakemake.output.tsv)
