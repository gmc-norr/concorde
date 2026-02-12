# Configuration

Concorde uses four YAML configuration files:

| File | Purpose |
|------|---------|
| `config/config.yaml` | Main pipeline configuration |
| `config/stratifications.yaml` | 11 stratification dimension definitions |
| `config/acceptance_criteria.yaml` | Per-tier thresholds and decision logic |
| `config/bram_config.yaml` | Bayesian risk assessment parameters |

## Main Configuration (`config/config.yaml`)

```yaml
# ── Execution Mode ──
mode: "germline"  # "germline" or "somatic"

# ── Sample Identity ──
sample: "NA12878"
caller: "GATK-HC"
pipeline_version: "v4.5.0"
parameters: "GATK HaplotypeCaller default settings"

# ── Comparison Tool ──
# Germline: "happy" (hap.py) or "rtg" (RTG vcfeval)
# Somatic: "sompy" (som.py), "rtg", or "internal"
comparison_tool: "happy"

# ── Input Files ──
reference: "/data/references/GRCh38/GRCh38.fa"
truth_vcf: "/data/truth/HG001_truth.vcf.gz"
confident_regions: "/data/truth/HG001_confident.bed"
query_vcfs:
  - "/data/query/NA12878_gatk_hc.vcf.gz"

# ── Output ──
database: "data/concorde.db"
output_dir: "results"
decomposition_modes:
  - "decomposed"
  - "non-decomposed"

# ── Gene Sets ──
gene_sets:
  - name: "ACMG59"
    description: "ACMG SF v3.1 actionable genes"
    version: "v3.1"
    bed: "/data/gene_sets/ACMG59_GRCh38.bed"
```

## Region Annotation BED Files

Configure BED files for genomic context stratification dimensions (low_complexity, segdup, mappability). GC content is computed from the reference FASTA and does not require a BED file.

```yaml
region_beds:
  low_complexity: "/path/to/SimpleRepeat.bed.gz"       # tabix-indexed
  segmental_duplications: "/path/to/segdup.bed.gz"       # tabix-indexed
  mappability: "/path/to/mappability.bed.gz"             # tabix-indexed, score in column 4
  giab:
    enabled: false
    bed_dir: "/path/to/giab/stratification-bed-files/"
    categories: ["LowComplexity", "SegmentalDuplications", "GCcontent"]
```

When a BED file path is empty or not configured, the corresponding stratification dimension gracefully produces no strata.

## QC Trending

Longitudinal monitoring of metrics across historical runs using z-score flagging:

```yaml
qc_trending:
  enabled: false
  min_historical_runs: 5     # Minimum runs before trending activates
  z_score_threshold: 2.0     # Flag metrics with |z| > threshold
```

## Statistical Process Control (SPC)

Control chart computation for production monitoring:

```yaml
spc:
  enabled: false
  sigma_multiplier: 3.0      # Control limits: mean +/- N*sigma
  min_data_points: 10         # Minimum runs before SPC activates
```

## Somatic Mode

When `mode: "somatic"`, add the somatic block:

```yaml
mode: "somatic"
comparison_tool: "sompy"

somatic:
  tumor_sample: "TUMOR"
  normal_sample: "NORMAL"          # omit for tumor-only
  calling_mode: "tumor_normal"     # "tumor_only" or "tumor_normal"
  min_vaf: 0.0                     # minimum VAF for truth variants
  vaf_tolerance: 0.10              # absolute VAF difference tolerance
  expected_contamination: null     # optional upper-bound contamination
```

### Somatic Comparison Tools

| Tool | Use Case |
|------|----------|
| `sompy` | Default for somatic. Uses som.py for somatic-aware comparison. |
| `rtg` | RTG vcfeval with `--squash-ploidy` for somatic. |
| `internal` | Built-in position+allele matcher with VAF concordance checking. |

## Truth Matching (optional)

Override default matching behavior per mode:

```yaml
truth_matching:
  germline:
    tool: "happy"
    genotype_match: true
    partial_credit: false
  somatic:
    tool: "sompy"
    vaf_tolerance: 0.10
    position_window: 0
    require_allele_match: true
    track_partial_matches: true
```

## VEP Annotation (optional)

Enable VEP annotation for functional impact stratification:

```yaml
vep:
  enabled: true
  vep_version: "110"
  cache_version: "110"
  genome_build: "GRCh38"
  transcript_source: "ensembl"       # "ensembl", "refseq", or "merged"
  transcript_selection: "canonical"  # "canonical", "mane_select", or "most_severe"
  cache_dir: "/data/vep/cache"
  plugins: []                        # e.g., ["CADD", "SpliceAI"]
  container: "docker://ensemblorg/ensembl-vep:release_110.1"
```

## Stratification (`config/stratifications.yaml`)

Defines how variants are binned into strata for per-stratum metrics. See [Stratification](stratification.md) for details.

```yaml
min_variants_per_stratum: 20

dimensions:
  variant_class:
    enabled: true
    strata: [SNP, INDEL, COMPLEX]

  indel_size:
    enabled: true
    bins:
      - { name: "1bp", min: 1, max: 1 }
      - { name: "2-5bp", min: 2, max: 5 }
      - { name: "6-15bp", min: 6, max: 15 }
      - { name: "16-50bp", min: 16, max: 50 }
      - { name: ">50bp", min: 51, max: null }

  functional_impact:
    enabled: true
    strata: [HIGH, MODERATE, LOW, MODIFIER]

  gene_panel:
    enabled: true

  zygosity:
    enabled: true
    mode: "germline"
    strata: [HET, HOM_ALT]

  vaf_bins:
    enabled: true
    mode: "somatic"
    bins:
      - { name: "<0.05", min: 0.0, max: 0.05 }
      - { name: "0.05-0.10", min: 0.05, max: 0.10 }
      - { name: "0.10-0.20", min: 0.10, max: 0.20 }
      - { name: "0.20-0.50", min: 0.20, max: 0.50 }
      - { name: ">0.50", min: 0.50, max: 1.0 }

  low_complexity:
    enabled: true                    # Requires region_beds.low_complexity

  gc_content:
    enabled: true                    # Requires reference FASTA
    bins:
      - { name: "0-25%", min: 0, max: 0.25 }
      - { name: "25-30%", min: 0.25, max: 0.30 }
      - { name: "30-55%", min: 0.30, max: 0.55 }
      - { name: "55-65%", min: 0.55, max: 0.65 }
      - { name: "65-100%", min: 0.65, max: 1.0 }
    window_size: 100

  segdup:
    enabled: true                    # Requires region_beds.segmental_duplications

  mappability:
    enabled: true                    # Requires region_beds.mappability
    bins:
      - { name: "low", min: 0, max: 0.5 }
      - { name: "medium", min: 0.5, max: 0.9 }
      - { name: "high", min: 0.9, max: 1.01 }

  coverage_depth:
    enabled: true                    # Uses variant DP field directly
    bins:
      - { name: "<10x", min: 0, max: 10 }
      - { name: "10-30x", min: 10, max: 30 }
      - { name: "30-50x", min: 30, max: 50 }
      - { name: ">50x", min: 50, max: null }
```

## Acceptance Criteria (`config/acceptance_criteria.yaml`)

Risk-based tiered thresholds. See [Acceptance Criteria](acceptance.md) for details.

```yaml
acceptance_criteria:
  mode: "germline"

  global:
    min_sensitivity: 0.99
    min_precision: 0.999

  per_stratum:
    SNP:
      min_sensitivity: 0.995
      min_precision: 0.999
    INDEL:
      min_sensitivity: 0.98
      min_precision: 0.995

  per_panel:
    ACMG59:
      min_sensitivity: 0.999
      min_precision: 0.9999

  tiers:
    tier_1:
      label: "Critical"
      panels: ["ACMG59"]
      impacts: ["HIGH"]
    tier_2:
      label: "Standard"
      impacts: ["MODERATE", "LOW"]
    tier_3:
      label: "Informational"
      impacts: ["MODIFIER"]

  low_count_threshold: 20
  low_count_policy: "warn"  # warn | exclude | fail
```

## BRAM Configuration (`config/bram_config.yaml`)

Bayesian Risk Assessment Module parameters. See [BRAM](bram.md) for details.

```yaml
bram:
  enabled: true
  alert_threshold: 0.80
  aggregation_method: "max"

  default_prior:
    concentration: 20          # Effective prior sample size (kappa)

  degradation_thresholds:      # Per-tier (smaller = stricter)
    tier_1: 0.005              # 0.5% for critical panels
    tier_2: 0.01               # 1% for standard
    tier_3: 0.02               # 2% for informational

  use_historical_priors: true
  min_historical_runs: 5

  per_stratum_overrides: {}    # Can override concentration and/or threshold
```

