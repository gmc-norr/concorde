"""
Concorde NGS Variant Validation Pipeline
=========================================
Snakemake workflow for truth-vs-query concordance analysis.

Supports two execution modes:
- germline: Diploid variant calling with genotype concordance
- somatic: Tumor-normal or tumor-only calling with VAF-based evaluation

Supports multiple comparison tools:
- hap.py (germline default): Genome Analysis Toolkit variant comparison
- rtg vcfeval: Real Time Genomics vcfeval tool
- som.py (somatic default): hap.py somatic module
- internal: Position + allele matcher (somatic fallback)

Supports running both decomposed and non-decomposed variant normalization,
producing separate Run records for each mode.

Usage:
    snakemake --cores all --configfile config/config.yaml --use-apptainer
    snakemake --cores all --configfile config/config.yaml -n   # dry-run

Note: The run_happy rule uses a container (via Apptainer/Singularity) because
the hap.py bioconda package only supports Python 2.7. rtg vcfeval requires
rtg-tools to be available in PATH. All other rules run natively in the pixi
pipeline environment.
"""

import os
from pathlib import Path

# ─── Configuration ───
configfile: "config/config.yaml"

OUTDIR = config["output_dir"]
REFERENCE = config["reference"]
TRUTH_VCF = config["truth_vcf"]
CONFIDENT_BED = config.get("confident_regions", "")
QUERY_VCFS = config["query_vcfs"]
GENE_SETS = config.get("gene_sets", [])
DATABASE = config["database"]
SAMPLE = config["sample"]
CALLER = config["caller"]
PIPELINE_VERSION = config["pipeline_version"]
PARAMETERS = config.get("parameters", "")
DECOMP_MODES = config.get("decomposition_modes", ["decomposed"])
COMPARISON_TOOL = config.get("comparison_tool", "happy")
MODE = config.get("mode", "germline")
SOMATIC_CONFIG = config.get("somatic", {})
TRUTH_MATCHING = config.get("truth_matching", {})
VEP_CONFIG = config.get("vep", {})
ENABLE_VEP = VEP_CONFIG.get("enabled", False)
STRATIFICATION_CONFIG = config.get("stratification", {})
NFCORE_QC_DIR = config.get("nfcore_qc_dir", "")
NFCORE_SAMPLE_MAP = config.get("nfcore_sample_mapping", {})


def query_id(vcf_path):
    """Derive a short identifier from a query VCF filename."""
    return Path(vcf_path).name.replace(".vcf.gz", "").replace(".vcf", "")


QUERY_IDS = [query_id(q) for q in QUERY_VCFS]
QUERY_MAP = dict(zip(QUERY_IDS, QUERY_VCFS))


def norm_flags(decomp_mode):
    """Return bcftools norm flags for the given decomposition mode.

    In somatic mode, multi-allelic splitting is skipped by default
    (somatic callers typically handle multi-allelic records differently).
    """
    if MODE == "somatic" and decomp_mode == "decomposed":
        return ""  # Somatic: preserve multi-allelic records by default
    elif decomp_mode == "decomposed":
        return "-m-both"  # Germline: split multi-allelic + left-align
    else:
        return ""  # left-align only (bcftools norm default)


def get_nfcore_sample_id(qid):
    """Get nf-core/raredisease sample ID for a query ID."""
    return NFCORE_SAMPLE_MAP.get(qid, qid)


# ─── Pre-flight validation ───
rule validate_inputs:
    input:
        config="config/config.yaml",
        truth_vcf=TRUTH_VCF,
        query_vcfs=QUERY_VCFS,
        reference=REFERENCE,
        confident_bed=CONFIDENT_BED if CONFIDENT_BED else [],
        gene_set_beds=[gs.get("bed", "") for gs in GENE_SETS if gs.get("bed")],
    output:
        validation_done=os.path.join(OUTDIR, ".validation_done"),
    log:
        os.path.join(OUTDIR, "logs", "validate_inputs.log"),
    script:
        "scripts/validation/validate_inputs.py"


# ─── Top-level target ───
rule all:
    input:
        validation_done=os.path.join(OUTDIR, ".validation_done"),
        ingestion_outputs=expand(
            os.path.join(OUTDIR, "{mode}", "{qid}", "ingest.done"),
            mode=DECOMP_MODES,
            qid=QUERY_IDS,
        ),


# ─── Rule 1: Normalize truth VCF ───
rule normalize_truth:
    input:
        vcf=TRUTH_VCF,
        ref=REFERENCE,
    output:
        vcf=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        tbi=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz.tbi"),
    params:
        norm_flags=lambda wc: norm_flags(wc.mode),
    log:
        os.path.join(OUTDIR, "logs", "normalize_truth_{mode}.log"),
    shell:
        """
        (bcftools norm {params.norm_flags} -f {input.ref} {input.vcf} \
            | bcftools sort \
            | bgzip -c > {output.vcf}) 2> {log}
        tabix -p vcf {output.vcf}
        """


# ─── Rule 2: Normalize each query VCF ───
rule normalize_query:
    input:
        vcf=lambda wc: QUERY_MAP[wc.qid],
        ref=REFERENCE,
    output:
        vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        tbi=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz.tbi"),
    params:
        norm_flags=lambda wc: norm_flags(wc.mode),
    log:
        os.path.join(OUTDIR, "logs", "normalize_query_{mode}_{qid}.log"),
    shell:
        """
        (bcftools norm {params.norm_flags} -f {input.ref} {input.vcf} \
            | bcftools sort \
            | bgzip -c > {output.vcf}) 2> {log}
        tabix -p vcf {output.vcf}
        """


# ─── Rule 3a: Run hap.py ───
# Uses container because hap.py bioconda package only supports Python 2.7.
# Run with: snakemake --use-apptainer
rule run_happy:
    input:
        truth=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        query=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        ref=REFERENCE,
    output:
        vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output.vcf.gz"),
        csv=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output.extended.csv"),
        summary=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output.summary.csv"),
    params:
        prefix=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output"),
        confident=f"-f {CONFIDENT_BED}" if CONFIDENT_BED else "",
    log:
        os.path.join(OUTDIR, "logs", "happy_{mode}_{qid}.log"),
    container:
        "docker://jmcdani20/hap.py:v0.3.12"
    threads: 4
    shell:
        """
        /opt/hap.py/bin/hap.py {input.truth} {input.query} \
            -r {input.ref} \
            {params.confident} \
            -o {params.prefix} \
            --threads {threads} \
            2>&1 | tee {log}
        """


# ─── Rule 3b: Run rtg vcfeval ───
# Requires rtg-tools in PATH
rule run_rtg:
    input:
        truth=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        query=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        ref=REFERENCE,
    output:
        tp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "tp.vcf.gz"),
        fp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "fp.vcf.gz"),
        fn_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "fn.vcf.gz"),
        summary=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "summary.txt"),
    params:
        outdir=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg"),
        confident=f"--bed-regions={CONFIDENT_BED}" if CONFIDENT_BED else "",
    log:
        os.path.join(OUTDIR, "logs", "rtg_{mode}_{qid}.log"),
    threads: 4
    shell:
        """
        # rtg vcfeval requires an SDF reference (pre-formatted)
        # For simplicity, assume reference is already in SDF format or use rtg format
        # If not SDF, you may need: rtg format -o ref.sdf {input.ref}

        rtg vcfeval \
            -b {input.truth} \
            -c {input.query} \
            -t {input.ref} \
            -o {params.outdir} \
            {params.confident} \
            --threads {threads} \
            2>&1 | tee {log}
        """


# ─── Rule 3c: Run som.py (somatic hap.py module) ───
rule run_sompy:
    input:
        truth=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        query=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        ref=REFERENCE,
    output:
        vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "sompy", "output.vcf.gz"),
        csv=os.path.join(OUTDIR, "{mode}", "{qid}", "sompy", "output.stats.csv"),
    params:
        prefix=os.path.join(OUTDIR, "{mode}", "{qid}", "sompy", "output"),
        confident=f"-f {CONFIDENT_BED}" if CONFIDENT_BED else "",
    log:
        os.path.join(OUTDIR, "logs", "sompy_{mode}_{qid}.log"),
    container:
        "docker://jmcdani20/hap.py:v0.3.12"
    threads: 4
    shell:
        """
        /opt/hap.py/bin/som.py {input.truth} {input.query} \
            -r {input.ref} \
            {params.confident} \
            -o {params.prefix} \
            2>&1 | tee {log}
        """


# ─── Rule 3d: Run rtg vcfeval for somatic (--squash-ploidy) ───
rule run_rtg_somatic:
    input:
        truth=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        query=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        ref=REFERENCE,
    output:
        tp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg_somatic", "tp.vcf.gz"),
        fp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg_somatic", "fp.vcf.gz"),
        fn_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg_somatic", "fn.vcf.gz"),
        summary=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg_somatic", "summary.txt"),
    params:
        outdir=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg_somatic"),
        confident=f"--bed-regions={CONFIDENT_BED}" if CONFIDENT_BED else "",
    log:
        os.path.join(OUTDIR, "logs", "rtg_somatic_{mode}_{qid}.log"),
    threads: 4
    shell:
        """
        rtg vcfeval \
            -b {input.truth} \
            -c {input.query} \
            -t {input.ref} \
            -o {params.outdir} \
            {params.confident} \
            --squash-ploidy \
            --threads {threads} \
            2>&1 | tee {log}
        """


# ─── Rule 3e: Run internal matcher (somatic fallback) ───
rule run_internal_matcher:
    input:
        truth_vcf=os.path.join(OUTDIR, "{mode}", "truth_normalized.vcf.gz"),
        query_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
    output:
        variants_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "internal", "variants.tsv"),
        metrics_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "internal", "metrics.tsv"),
    params:
        vaf_tolerance=SOMATIC_CONFIG.get("vaf_tolerance", 0.10),
        position_window=TRUTH_MATCHING.get("somatic", {}).get("position_window", 0),
        require_allele_match=TRUTH_MATCHING.get("somatic", {}).get("require_allele_match", True),
        track_partial_matches=TRUTH_MATCHING.get("somatic", {}).get("track_partial_matches", True),
    log:
        os.path.join(OUTDIR, "logs", "internal_matcher_{mode}_{qid}.log"),
    script:
        "scripts/matching/internal_matcher.py"


# ─── Rule 4: Parse comparison tool output -> variants TSV ───
# Branches on MODE and COMPARISON_TOOL
if MODE == "somatic" and COMPARISON_TOOL == "internal":
    rule parse_variants:
        input:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "internal", "variants.tsv"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_variants_{mode}_{qid}.log"),
        shell:
            """
            cp {input.tsv} {output.tsv}
            echo "Copied internal matcher variants" >> {log}
            """
elif MODE == "somatic" and COMPARISON_TOOL == "sompy":
    rule parse_variants:
        input:
            vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "sompy", "output.vcf.gz"),
            query_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
        params:
            tumor_sample=SOMATIC_CONFIG.get("tumor_sample", ""),
            normal_sample=SOMATIC_CONFIG.get("normal_sample", ""),
        log:
            os.path.join(OUTDIR, "logs", "parse_variants_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_sompy_vcf.py"
elif COMPARISON_TOOL == "rtg":
    rule parse_variants:
        input:
            tp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "tp.vcf.gz"),
            fp_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "fp.vcf.gz"),
            fn_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "fn.vcf.gz"),
            query_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_variants_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_rtg_vcf.py"
else:
    rule parse_variants:
        input:
            vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output.vcf.gz"),
            query_vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_variants_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_happy_vcf.py"


# ─── Rule 5: Parse comparison tool output -> metrics TSV ───
if MODE == "somatic" and COMPARISON_TOOL == "internal":
    rule parse_metrics:
        input:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "internal", "metrics.tsv"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "metrics.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_metrics_{mode}_{qid}.log"),
        shell:
            """
            cp {input.tsv} {output.tsv}
            echo "Copied internal matcher metrics" >> {log}
            """
elif MODE == "somatic" and COMPARISON_TOOL == "sompy":
    rule parse_metrics:
        input:
            csv=os.path.join(OUTDIR, "{mode}", "{qid}", "sompy", "output.stats.csv"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "metrics.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_metrics_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_sompy_metrics.py"
elif COMPARISON_TOOL == "rtg":
    rule parse_metrics:
        input:
            summary=os.path.join(OUTDIR, "{mode}", "{qid}", "rtg", "summary.txt"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "metrics.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_metrics_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_rtg_metrics.py"
else:
    rule parse_metrics:
        input:
            csv=os.path.join(OUTDIR, "{mode}", "{qid}", "happy", "output.extended.csv"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "metrics.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_metrics_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_happy_metrics.py"


# ─── Rule 6: Intersect variants with gene set BEDs ───
rule intersect_gene_sets:
    input:
        variants_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
    output:
        tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "gene_set_assignments.tsv"),
    params:
        gene_sets=GENE_SETS,
    log:
        os.path.join(OUTDIR, "logs", "intersect_gene_sets_{mode}_{qid}.log"),
    script:
        "scripts/analysis/intersect_gene_sets.py"


# ─── Rule 6a: VEP annotation (optional) ───
if ENABLE_VEP:
    rule run_vep:
        input:
            vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "query_normalized.vcf.gz"),
            ref=REFERENCE,
        output:
            vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "vep", "annotated.vcf.gz"),
        params:
            cache_dir=VEP_CONFIG.get("cache_dir", ""),
            genome_build=VEP_CONFIG.get("genome_build", "GRCh38"),
            transcript_source=VEP_CONFIG.get("transcript_source", "ensembl"),
            plugins=",".join(VEP_CONFIG.get("plugins", [])),
        log:
            os.path.join(OUTDIR, "logs", "vep_{mode}_{qid}.log"),
        container:
            VEP_CONFIG.get("container", "docker://ensemblorg/ensembl-vep:release_110.1")
        threads: 4
        shell:
            """
            vep --input_file {input.vcf} \
                --output_file {output.vcf} \
                --vcf --compress_output bgzip \
                --fasta {input.ref} \
                --dir_cache {params.cache_dir} \
                --assembly {params.genome_build} \
                --{params.transcript_source} \
                --cache --offline \
                --canonical --symbol --numbers --hgvs \
                --sift b --polyphen b \
                --force_overwrite \
                --fork {threads} \
                2>&1 | tee {log}
            tabix -p vcf {output.vcf}
            """

    rule parse_vep_annotations:
        input:
            vcf=os.path.join(OUTDIR, "{mode}", "{qid}", "vep", "annotated.vcf.gz"),
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "vep_annotations.tsv"),
        params:
            transcript_selection=VEP_CONFIG.get("transcript_selection", "canonical"),
        log:
            os.path.join(OUTDIR, "logs", "parse_vep_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_vep_annotations.py"
else:
    rule parse_vep_annotations:
        output:
            tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "vep_annotations.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_vep_{mode}_{qid}.log"),
        shell:
            """
            # Create empty VEP annotations file when VEP is disabled
            echo -e "chrom\tpos\tref\talt\tconsequence\timpact\tgene_symbol\tgene_id\ttranscript_id\thgvsc\thgvsp\texon\tprotein_position\tsift_prediction\tpolyphen_prediction" > {output.tsv}
            echo "VEP annotation disabled in config" >> {log}
            """


# ─── Rule 6b: Parse nf-core/raredisease QC outputs (optional) ───
if NFCORE_QC_DIR:
    rule parse_nfcore_qc:
        input:
            multiqc_json=lambda wc: os.path.join(
                NFCORE_QC_DIR, "multiqc", "multiqc_data", "multiqc_data.json"
            )
            if os.path.exists(os.path.join(NFCORE_QC_DIR, "multiqc", "multiqc_data", "multiqc_data.json"))
            else [],
            mosdepth_summary=lambda wc: os.path.join(
                NFCORE_QC_DIR,
                "qc_bam",
                f"{get_nfcore_sample_id(wc.qid)}_mosdepth.summary.txt",
            )
            if os.path.exists(
                os.path.join(
                    NFCORE_QC_DIR,
                    "qc_bam",
                    f"{get_nfcore_sample_id(wc.qid)}_mosdepth.summary.txt",
                )
            )
            else [],
            mosdepth_dist=lambda wc: os.path.join(
                NFCORE_QC_DIR,
                "qc_bam",
                f"{get_nfcore_sample_id(wc.qid)}_mosdepth.global.dist.txt",
            )
            if os.path.exists(
                os.path.join(
                    NFCORE_QC_DIR,
                    "qc_bam",
                    f"{get_nfcore_sample_id(wc.qid)}_mosdepth.global.dist.txt",
                )
            )
            else [],
            verifybamid_selfsm=lambda wc: os.path.join(
                NFCORE_QC_DIR, "qc_bam", f"{get_nfcore_sample_id(wc.qid)}.selfSM"
            )
            if os.path.exists(
                os.path.join(
                    NFCORE_QC_DIR, "qc_bam", f"{get_nfcore_sample_id(wc.qid)}.selfSM"
                )
            )
            else [],
            software_versions_yml=lambda wc: os.path.join(
                NFCORE_QC_DIR, "pipeline_info", "software_versions.yml"
            )
            if os.path.exists(
                os.path.join(NFCORE_QC_DIR, "pipeline_info", "software_versions.yml")
            )
            else [],
            pipeline_params_json=lambda wc: os.path.join(
                NFCORE_QC_DIR, "pipeline_info", "params.json"
            )
            if os.path.exists(
                os.path.join(NFCORE_QC_DIR, "pipeline_info", "params.json")
            )
            else [],
        output:
            summary_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_summary.tsv"),
            metrics_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_metrics.tsv"),
            versions_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "software_versions.tsv"),
            metadata_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "pipeline_metadata.tsv"),
        params:
            sample_id=lambda wc: get_nfcore_sample_id(wc.qid),
        log:
            os.path.join(OUTDIR, "logs", "parse_nfcore_qc_{mode}_{qid}.log"),
        script:
            "scripts/parsers/parse_nfcore_qc.py"
else:
    rule parse_nfcore_qc:
        output:
            summary_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_summary.tsv"),
            metrics_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_metrics.tsv"),
            versions_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "software_versions.tsv"),
            metadata_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "pipeline_metadata.tsv"),
        log:
            os.path.join(OUTDIR, "logs", "parse_nfcore_qc_{mode}_{qid}.log"),
        shell:
            """
            # Create empty QC files when nf-core QC is not available
            echo -e "mean_coverage\\tmedian_coverage\\tpct_bases_10x\\tpct_bases_20x\\tpct_bases_30x\\tmapping_rate\\tduplicate_rate\\tmedian_insert_size\\tcontamination_pct\\tsex_inferred" > {output.summary_tsv}
            echo -e "metric_name\\tmetric_value_float\\tmetric_value_str\\tmetric_source\\tmetric_category\\tnotes" > {output.metrics_tsv}
            echo -e "tool_name\\tversion\\tcategory\\tnotes" > {output.versions_tsv}
            echo -e "nfcore_pipeline_version\\tpipeline_params_json" > {output.metadata_tsv}
            echo "nf-core QC directory not configured, created empty QC files" >> {log}
            """


# ─── Rule 7: Ingest everything into the database ───
rule ingest:
    input:
        variants_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "variants.tsv"),
        metrics_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "metrics.tsv"),
        gene_set_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "gene_set_assignments.tsv"),
        qc_summary_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_summary.tsv"),
        qc_metrics_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "qc_metrics.tsv"),
        software_versions_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "software_versions.tsv"),
        pipeline_metadata_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "pipeline_metadata.tsv"),
        vep_annotations_tsv=os.path.join(OUTDIR, "{mode}", "{qid}", "vep_annotations.tsv"),
    output:
        done=os.path.join(OUTDIR, "{mode}", "{qid}", "ingest.done"),
    params:
        database=DATABASE,
        sample=SAMPLE,
        caller=CALLER,
        pipeline_version=PIPELINE_VERSION,
        parameters=PARAMETERS,
        decomposition_mode=lambda wc: wc.mode,
        comparison_tool=COMPARISON_TOOL,
        execution_mode=MODE,
        somatic_config=SOMATIC_CONFIG,
        gene_sets_config=GENE_SETS,
        stratification_config=STRATIFICATION_CONFIG,
    log:
        os.path.join(OUTDIR, "logs", "ingest_{mode}_{qid}.log"),
    script:
        "scripts/ingestion/ingest.py"


# ─── Rule 8: Create baseline (optional, on-demand) ───
# Run with: snakemake create_baseline --config baseline_name="v1.0_germline"
BASELINE_NAME = config.get("baseline_name", "")
BASELINE_TOLERANCE = config.get("baseline_tolerance", 0.005)

if BASELINE_NAME:
    rule create_baseline:
        input:
            ingest_done=expand(
                os.path.join(OUTDIR, "{mode}", "{qid}", "ingest.done"),
                mode=DECOMP_MODES,
                qid=QUERY_IDS,
            ),
        output:
            done=os.path.join(OUTDIR, "baseline", f"{BASELINE_NAME}.done"),
        params:
            database=DATABASE,
            baseline_name=BASELINE_NAME,
            tolerance=BASELINE_TOLERANCE,
            execution_mode=MODE,
            pipeline_version=PIPELINE_VERSION,
        log:
            os.path.join(OUTDIR, "logs", f"create_baseline_{BASELINE_NAME}.log"),
        script:
            "scripts/analysis/create_baseline.py"


# ─── Rule 9: Generate validation report (optional) ───
ENABLE_REPORT = config.get("generate_report", False)
ACCEPTANCE_CONFIG = config.get("acceptance_criteria", {})

if ENABLE_REPORT:
    rule generate_report:
        input:
            ingest_done=expand(
                os.path.join(OUTDIR, "{mode}", "{qid}", "ingest.done"),
                mode=DECOMP_MODES,
                qid=QUERY_IDS,
            ),
        output:
            json=os.path.join(OUTDIR, "reports", "report.json"),
            html=os.path.join(OUTDIR, "reports", "report.html"),
        params:
            database=DATABASE,
            sample=SAMPLE,
            execution_mode=MODE,
            pipeline_version=PIPELINE_VERSION,
            caller=CALLER,
            comparison_tool=COMPARISON_TOOL,
            acceptance_config=ACCEPTANCE_CONFIG,
        log:
            os.path.join(OUTDIR, "logs", "generate_report.log"),
        script:
            "scripts/reporting/generate_report.py"
