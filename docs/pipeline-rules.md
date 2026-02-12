# Pipeline Rules

The Concorde pipeline consists of 11 Snakemake rules that execute in dependency order. Each rule delegates to a Python script containing testable pure functions.

## Rule 0: validate_inputs

**Purpose**: Pre-flight validation of all inputs (31 checks).

**Script**: `scripts/validation/validate_inputs.py`

Runs automatically before any processing. The pipeline will not proceed if validation fails. Uses a collect-all-errors strategy to report every issue at once.

---

## Rule 1: compare_variants

**Purpose**: Run variant benchmarking against truth sets.

**Tools**:

| Mode | Tool | Container |
|------|------|-----------|
| Germline | hap.py | `docker://pkrusche/hap.py` |
| Germline | RTG vcfeval | `docker://realtimegenomics/rtg-tools` |
| Somatic | som.py | `docker://pkrusche/hap.py` |
| Somatic | RTG vcfeval | `docker://realtimegenomics/rtg-tools` |
| Somatic | Internal matcher | No container needed |

**Output**: Comparison VCF with TP/FP/FN annotations and summary metrics.

---

## Rule 2: parse_comparison

**Purpose**: Extract TP/FP/FN variants and performance metrics from comparison output.

**Scripts** (selected based on comparison_tool and mode):

- `scripts/parsers/parse_happy_vcf.py` + `parse_happy_metrics.py`
- `scripts/parsers/parse_rtg_vcf.py` + `parse_rtg_metrics.py`
- `scripts/parsers/parse_sompy_vcf.py` + `parse_sompy_metrics.py`

Parsers extract genotype fields (`truth_gt`, `query_gt`) and compute `indel_size` for INDEL variants.

**Output**: `variants.tsv` (classified variants) and `metrics.tsv` (aggregate metrics).

---

## Rule 3: parse_qc_metrics

**Purpose**: Parse QC metrics from nf-core/raredisease outputs.

**Script**: `scripts/parsers/parse_nfcore_qc.py`

**Output**: `qc_summary.tsv`, `qc_metrics.tsv`, `pipeline_metadata.tsv`, `software_versions.tsv`

Only runs when `nfcore_qc_dir` is configured.

---

## Rule 4: assign_gene_sets

**Purpose**: Annotate variants with gene names and gene set memberships.

**Script**: `scripts/analysis/intersect_gene_sets.py`

**Output**: `gene_set_assignments.tsv` mapping variants to gene panels.

---

## Rule 5: run_vep (optional)

**Purpose**: Annotate variants with VEP consequences, gene symbols, and impact tiers.

**Script**: `scripts/parsers/parse_vep_annotations.py`

Only runs when `vep.enabled: true`. Adds: `consequence`, `gene_symbol`, `gene_id`, `transcript_id`, `hgvsc`, `hgvsp`, `impact`, `sift_prediction`, `polyphen_prediction`.

---

## Rule 6: stratify_variants

**Purpose**: Assign variants to strata across 11 dimensions and compute per-stratum metrics.

**Script**: `scripts/analysis/stratification.py`

Computes precision, recall, F1, Wilson score CIs, bootstrap F1 CI, MCC, and genotype concordance independently for each enabled dimension. Strata with fewer than `min_variants_per_stratum` variants are flagged as low-confidence.

See [Stratification](stratification.md) for details.

---

## Rule 7: evaluate_acceptance

**Purpose**: Apply risk-based acceptance criteria with optional BRAM augmentation.

**Scripts**: `scripts/analysis/acceptance.py`, `scripts/analysis/bayesian_risk.py`

Produces a per-stratum pass/fail result and an overall verdict. When BRAM is enabled, a two-stage decision combines deterministic and probabilistic results.

See [Acceptance Criteria](acceptance.md) and [BRAM](bram.md) for details.

---

## Rule 8: generate_report

**Purpose**: Generate JSON, HTML, and PDF validation reports.

**Scripts**: `scripts/reporting/generate_report.py`, `json_report.py`, `html_report.py`, `pdf_report.py`

Reports include: header, executive summary, per-stratum metrics (with CIs and MCC), acceptance results, BRAM assessment, QC trending, variant diffs (verification), and traceability.

See [Reporting](reporting.md) for details.

---

## Rule 9: ingest_to_db

**Purpose**: Load all data into the SQLite database.

**Scripts**: `scripts/ingestion/ingest.py`, `ingest_helpers.py`

Populates 19 ORM models. During ingestion:
- Region annotation is applied (GC content, LCR, segdup, mappability) when BED files are configured
- Zygosity is derived from genotype when not explicitly provided
- Ti/Tv ratio and Het/Hom ratio are computed and stored as QCMetric records
- Audit events are recorded

---

## Rule 10: llm_analysis (optional)

**Purpose**: Generate AI-powered QC insights using local LLM.

**Script**: `scripts/analysis/llm_qc_analyzer.py`

Only runs when `llm_analysis.enabled: true`. Requires Ollama running locally.

---

## Output Directory Structure

```
results/
  {sample}_{caller}_{version}_{mode}/
    comparison/              # Comparison tool output
    variants.tsv             # Parsed variants with classifications
    metrics.tsv              # Aggregate performance metrics
    qc_summary.tsv           # Key QC metrics (optional)
    qc_metrics.tsv           # Detailed QC metrics (optional)
    gene_set_assignments.tsv # Gene set annotations
    report.json              # JSON validation report
    report.html              # HTML validation report
    report.pdf               # PDF validation report
    ingest_done.txt          # Ingestion complete marker
```

---

## Rule 11: verify (conditional)

**Purpose**: Run continuous verification against a locked baseline.

**Script**: `scripts/analysis/run_verification.py`

Only runs when `verification_baseline` is configured. Orchestrates:
1. Variant-level diff between baseline and current run
2. Drift vs. biological classification of transitions
3. Root cause evidence collection
4. Baseline envelope checking
5. BRAM assessment (when enabled)
6. Audit event recording
