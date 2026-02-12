# Getting Started

## Prerequisites

### Required Software

- **pixi** - Package manager (installs all other dependencies)
- **Apptainer** or **Docker** - Container runtime for hap.py/RTG/som.py

### Input Data Requirements

1. **Query VCF** - Variant calls to validate (bgzipped + tabix indexed)
2. **Truth VCF** - Gold standard truth set (bgzipped + tabix indexed)
3. **Truth BED** - Confident regions for comparison
4. **Reference Genome** - FASTA file matching VCF coordinates
5. **nf-core/raredisease Output** (optional) - QC metrics directory

## Installation

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and install
git clone <repo-url>
cd concorde-pipeline
pixi install
```

## Running the Pipeline

### Dry Run (Preview)

Preview what the pipeline will do without executing:

```bash
pixi run dry-run
```

### Full Execution

```bash
pixi run run
```

### Targeted Execution

```bash
# Run only variant comparison
snakemake --cores all --configfile config/config.yaml --until compare_variants

# Run only data ingestion
snakemake --cores all --configfile config/config.yaml --until ingest_to_db
```

### Cluster Execution

```bash
snakemake --cores all --configfile config/config.yaml \
  --cluster "sbatch -p compute -c {threads} --mem={resources.mem_mb}" \
  --jobs 100
```

## Input Validation

The pipeline includes 31 validators that run automatically before execution.

### What Gets Validated

- **Configuration structure** - Required fields, valid values, proper types, mode-specific fields
- **VCF files** - Bgzipped format, tabix index (.tbi) present
- **Reference FASTA** - Index file (.fai) present
- **BED files** - Proper format (chr, start, end), valid coordinates
- **Chromosome consistency** - Same naming convention across all files
- **Sample names** - VCF header matches config
- **Gene sets** - Proper structure with required fields
- **Somatic config** - Tumor/normal samples, calling mode, VAF tolerance (when mode=somatic)
- **VEP config** - Version, cache, genome build, transcript selection (when enabled)

### Fixing Common Validation Errors

**Missing VCF index:**
```bash
tabix -p vcf your_file.vcf.gz
```

**Missing FASTA index:**
```bash
samtools faidx your_reference.fa
```

**Chromosome naming mismatch:**
```bash
# Check chromosome names
bcftools view -H input.vcf.gz | cut -f1 | uniq

# Rename chromosomes
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz
```
