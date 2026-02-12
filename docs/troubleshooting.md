# Troubleshooting

## Validation Errors

The validation system provides specific file, exact problem, and recommended fix for each error. Run `pixi run dry-run` to check validation before full execution.

### Missing VCF Index

```
VCF file query.vcf.gz is missing tabix index.
Expected query.vcf.gz.tbi or query.vcf.gz.csi.
```

**Fix:**
```bash
tabix -p vcf your_file.vcf.gz
```

### Missing FASTA Index

```
Reference FASTA GRCh38.fa is missing index.
Expected GRCh38.fa.fai.
```

**Fix:**
```bash
samtools faidx your_reference.fa
```

### Chromosome Naming Mismatch

```
Chromosome naming inconsistency detected:
  - Query VCF uses with_chr format but truth VCF uses without_chr format
```

**Fix:**
```bash
# Check chromosome names
bcftools view -H input.vcf.gz | cut -f1 | uniq

# Create mapping file
echo "1 chr1" > chr_map.txt
echo "2 chr2" >> chr_map.txt
# ... etc

# Rename chromosomes
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz
tabix -p vcf output.vcf.gz
```

### Somatic Config Errors

```
Somatic config requires 'tumor_sample' when mode is "somatic"
```

**Fix:** Ensure `somatic.tumor_sample` and `somatic.calling_mode` are set in `config/config.yaml` when `mode: "somatic"`.

### Invalid BED Format

**Fix:** Ensure BED file has at least 3 tab-separated columns (chr, start, end) with numeric coordinates where start < end and both >= 0.

## Runtime Errors

### Apptainer Not Found

```
Command 'apptainer' not found
```

**Fix:**
```bash
sudo apt-get install apptainer
```

### Memory Error During Comparison

```
MemoryError: Unable to allocate array
```

**Fix:** Add memory constraints to the Snakefile rule:
```python
rule compare_variants:
    resources:
        mem_mb=16000  # 16GB RAM
```

### Ollama Connection Refused

```
ConnectionRefusedError: Connection refused (localhost:11434)
```

**Fix:**
```bash
# Start Ollama server
ollama serve

# Or disable LLM analysis
# In config.yaml: llm_analysis.enabled: false
```

## Debug Mode

Run with verbose logging:

```bash
snakemake --cores all --configfile config/config.yaml --verbose --printshellcmds
```

Resume from failed state:

```bash
snakemake --cores all --configfile config/config.yaml --rerun-incomplete
```
