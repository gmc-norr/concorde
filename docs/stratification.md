# Stratification Engine

The stratification engine assigns each variant to one or more strata across 11 independent dimensions, then computes per-stratum classification counts and derived metrics. This decomposition reveals performance variation across variant classes, genomic contexts, and clinical significance categories that global metrics obscure.

## Dimensions

| Dimension | Mode | Strata | Data Source |
|-----------|------|--------|-------------|
| `variant_class` | Both | SNP, INDEL, COMPLEX | Variant type field |
| `indel_size` | Both | 1bp, 2-5bp, 6-15bp, 16-50bp, >50bp | Allele length difference |
| `functional_impact` | Both | HIGH, MODERATE, LOW, MODIFIER | VEP `impact` field |
| `gene_panel` | Both | Auto-derived from configured gene sets + off_panel | BED intersection |
| `zygosity` | Germline | HET, HOM_ALT | Genotype field (derived from GT) |
| `vaf_bins` | Somatic | <0.05, 0.05-0.10, 0.10-0.20, 0.20-0.50, >0.50 | `tumor_af` field |
| `low_complexity` | Both | low_complexity, non_lcr | BED annotation (SimpleRepeat regions) |
| `gc_content` | Both | 0-25%, 25-30%, 30-55%, 55-65%, 65-100% | Reference FASTA (100bp window) |
| `segdup` | Both | in_segdup, non_segdup | BED annotation (segmental duplications) |
| `mappability` | Both | low (<0.5), medium (0.5-0.9), high (>0.9) | BED annotation (mappability scores) |
| `coverage_depth` | Both | <10x, 10-30x, 30-50x, >50x | Variant DP field |

A variant may belong to multiple strata across dimensions simultaneously. For example, a variant could be assigned to the `SNP` stratum (variant_class), the `HIGH` stratum (functional_impact), the `ACMG59` stratum (gene_panel), and the `30-50x` stratum (coverage_depth).

## Workflow

1. **Region annotation**: The `RegionAnnotator` annotates each variant with GC content, low-complexity status, segmental duplication overlap, and mappability score using tabix-indexed BED files and `pysam.FastaFile` for reference lookups. Variants retain `None` for any dimension whose BED file is not configured.

2. **Zygosity derivation**: For germline mode, zygosity (HET, HOM_ALT, HOM_REF) is derived from the genotype string (e.g., `0/1` -> HET, `1/1` -> HOM_ALT) during ingestion when not explicitly provided by the parser.

3. **Stratum assignment**: Each variant is assigned to strata in each enabled dimension by the `StratificationEngine.stratify()` method. Mode-specific dimensions are automatically skipped when inapplicable.

4. **Metric computation**: For each (dimension, stratum) pair, the engine computes TP, FP, and FN counts from variant classifications, then derives the full metric set.

5. **Low-confidence flagging**: Strata with fewer than `min_variants_per_stratum` (default: 20) variants are flagged.

## Per-Stratum Metrics

Each stratum produces the following metrics:

### Classification Metrics

| Metric | Formula | Implementation |
|--------|---------|----------------|
| Precision | TP / (TP + FP) | Shared `compute_classification_metrics()` utility |
| Recall (Sensitivity) | TP / (TP + FN) | Shared `compute_classification_metrics()` utility |
| F1 | 2 * Precision * Recall / (Precision + Recall) | Shared `compute_classification_metrics()` utility |

All three metrics return `None` when the denominator is zero (undefined), ensuring consistent handling of edge cases across all computation paths.

### Confidence Intervals

| Metric | Method | Library |
|--------|--------|---------|
| Precision CI | Wilson score interval | `scipy.stats.norm.ppf` for z-quantile |
| Recall CI | Wilson score interval | `scipy.stats.norm.ppf` for z-quantile |
| F1 CI | Percentile bootstrap (n=9,999) | `scipy.stats.bootstrap` |

Wilson score intervals are preferred over normal approximation because they maintain coverage accuracy for extreme proportions (near 0 or 1) and small sample sizes. The F1 score is a nonlinear function of precision and recall; its confidence interval is computed via bootstrap resampling to avoid distributional assumptions.

### Additional Metrics

| Metric | Description | Formula |
|--------|-------------|---------|
| MCC | Matthews Correlation Coefficient | (TP\*TN - FP\*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN)) |
| Genotype concordance | GT match rate among TPs | Count(truth_gt == query_gt) / Count(TPs with both GTs) |

MCC ranges from -1 to +1 and accounts for all four confusion matrix quadrants. It returns `None` when the denominator is zero.

### Run-Level QC Metrics

Computed once per run (not per stratum) and stored as QCMetric records:

| Metric | Expected Range | Interpretation |
|--------|---------------|----------------|
| Ti/Tv ratio | ~2.0-2.1 (WGS), ~2.8-3.0 (WES) | Departure indicates systematic error |
| Het/Hom ratio | ~1.5-2.0 (diploid) | Departure may indicate contamination or relatedness |

## Mode-Specific Dimensions

- **zygosity** applies only when `mode: "germline"`. Zygosity is derived from the genotype string during ingestion.
- **vaf_bins** applies only when `mode: "somatic"`. Bins are based on `tumor_af` (preferred) or `af`.
- **functional_impact** requires VEP annotation to be enabled.
- **low_complexity**, **gc_content**, **segdup**, **mappability** require corresponding BED files or reference FASTA to be configured in `config/config.yaml` under `region_beds`. When not configured, these dimensions produce no strata (graceful degradation).
- **coverage_depth** uses the variant's `dp` field directly and requires no external file.

## Low-Count Handling

Strata with fewer than `min_variants_per_stratum` (default: 20) variants are flagged as low-confidence. The `low_count_policy` in acceptance criteria configuration controls downstream treatment:

| Policy | Behavior |
|--------|----------|
| `warn` | Include in results with a low-confidence flag (default) |
| `exclude` | Exclude from acceptance evaluation |
| `fail` | Treat as an acceptance failure |

## Configuration

See `config/stratifications.yaml` for the full dimension reference. Each dimension can be independently enabled or disabled, and numeric dimensions (`indel_size`, `vaf_bins`, `gc_content`, `mappability`, `coverage_depth`) use configurable bin boundaries.

Region annotation BED files are specified in `config/config.yaml` under `region_beds`:

```yaml
region_beds:
  low_complexity: "/path/to/SimpleRepeat.bed.gz"
  segmental_duplications: "/path/to/segdup.bed.gz"
  mappability: "/path/to/mappability.bed.gz"
```

GIAB stratification BED files can be auto-discovered when configured:

```yaml
region_beds:
  giab:
    enabled: true
    bed_dir: "/path/to/giab/stratification-bed-files/"
    categories: ["LowComplexity", "SegmentalDuplications", "GCcontent"]
```
