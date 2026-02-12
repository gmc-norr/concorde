"""Tests for SQLAlchemy database models."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pytest

# Add project root to path
PROJECT_ROOT = str(Path(__file__).resolve().parent.parent.parent)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from models import (
    Baseline,
    BaselineEnvelope,
    BayesianRiskAssessment,
    FixtureRecord,
    RootCauseEvidence,
    Run,
    StratumPosterior,
    Variant,
    VariantTransition,
    VerificationResult,
    Metric,
    QCMetric,
    SoftwareVersion,
    GeneSet,
)


class TestRunModel:
    """Tests for Run model."""

    def test_create_run(self, test_db):
        """Test creating a Run instance."""
        run = Run(
            sample="NA12878",
            caller="GATK-HC",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
        )
        test_db.add(run)
        test_db.commit()

        assert run.id is not None
        assert run.sample == "NA12878"
        assert run.date is not None
        assert isinstance(run.date, datetime)

    def test_run_with_qc_metrics(self, test_db):
        """Test Run with QC metrics."""
        run = Run(
            sample="NA12878",
            caller="GATK-HC",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
            mean_coverage=35.2,
            median_coverage=34.8,
            pct_bases_10x=99.1,
        )
        test_db.add(run)
        test_db.commit()

        assert run.mean_coverage == 35.2
        assert run.median_coverage == 34.8
        assert run.pct_bases_10x == 99.1

    def test_run_string_representation(self, test_db):
        """Test Run __repr__ method."""
        run = Run(
            sample="NA12878",
            caller="GATK-HC",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
        )
        test_db.add(run)
        test_db.commit()

        repr_str = repr(run)
        assert "NA12878" in repr_str
        assert "GATK-HC" in repr_str


class TestVariantModel:
    """Tests for Variant model."""

    def test_create_variant(self, test_db, sample_run):
        """Test creating a Variant instance."""
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=12345,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)
        test_db.commit()

        assert variant.id is not None
        assert variant.chrom == "chr1"
        assert variant.pos == 12345
        assert variant.type == "SNP"

    def test_variant_run_relationship(self, test_db, sample_run):
        """Test Variant-Run relationship."""
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=12345,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)
        test_db.commit()

        # Access relationship
        assert variant.run.id == sample_run.id
        assert variant.run.sample == "NA12878"

    def test_variant_cascade_delete(self, test_db, sample_run):
        """Test that variants are deleted when run is deleted."""
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=12345,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)
        test_db.commit()

        variant_id = variant.id

        # Delete run
        test_db.delete(sample_run)
        test_db.commit()

        # Variant should be deleted (cascade)
        deleted_variant = test_db.query(Variant).filter(Variant.id == variant_id).first()
        assert deleted_variant is None


class TestMetricModel:
    """Tests for Metric model."""

    def test_create_metric(self, test_db, sample_run):
        """Test creating a Metric instance."""
        metric = Metric(
            run_id=sample_run.id,
            variant_type="SNP",
            stratification="PASS",
            precision=0.974,
            recall=0.983,
            f1=0.978,
            tp_count=2950,
            fp_count=80,
            fn_count=50,
        )
        test_db.add(metric)
        test_db.commit()

        assert metric.id is not None
        assert metric.variant_type == "SNP"
        assert metric.recall == 0.983


class TestQCMetricModel:
    """Tests for QCMetric model."""

    def test_create_qc_metric(self, test_db, sample_run):
        """Test creating a QCMetric instance."""
        qc_metric = QCMetric(
            run_id=sample_run.id,
            metric_name="mosdepth_mean_coverage",
            metric_value_float=35.2,
            metric_source="mosdepth",
            metric_category="coverage",
        )
        test_db.add(qc_metric)
        test_db.commit()

        assert qc_metric.id is not None
        assert qc_metric.metric_name == "mosdepth_mean_coverage"
        assert qc_metric.metric_value_float == 35.2


class TestGeneSetModel:
    """Tests for GeneSet model."""

    def test_create_gene_set(self, test_db):
        """Test creating a GeneSet instance."""
        gene_set = GeneSet(
            name="ACMG59",
            description="ACMG SF v3.1 actionable genes",
            version="3.1",
        )
        test_db.add(gene_set)
        test_db.commit()

        assert gene_set.id is not None
        assert gene_set.name == "ACMG59"

    def test_gene_set_variant_relationship(self, test_db, sample_run):
        """Test GeneSet-Variant many-to-many relationship."""
        gene_set = GeneSet(
            name="ACMG59",
            description="Test gene set",
            version="1.0",
        )
        test_db.add(gene_set)
        test_db.commit()

        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=12345,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)
        test_db.commit()

        # Add variant to gene set
        variant.gene_sets.append(gene_set)
        test_db.commit()

        # Verify relationship
        assert gene_set in variant.gene_sets
        assert variant in gene_set.variants


class TestSoftwareVersionModel:
    """Tests for SoftwareVersion model."""

    def test_create_software_version(self, test_db, sample_run):
        """Test creating a SoftwareVersion instance."""
        sw_version = SoftwareVersion(
            run_id=sample_run.id,
            tool_name="bcftools",
            version="1.17",
        )
        test_db.add(sw_version)
        test_db.commit()

        assert sw_version.id is not None
        assert sw_version.tool_name == "bcftools"
        assert sw_version.version == "1.17"


class TestSomaticRunModel:
    """Tests for Run model somatic fields."""

    def test_run_with_somatic_mode(self, test_db):
        run = Run(
            sample="TUMOR_001",
            caller="Mutect2",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="sompy",
            mode="somatic",
            tumor_sample="TUMOR",
            normal_sample="NORMAL",
            calling_mode="tumor_normal",
        )
        test_db.add(run)
        test_db.commit()

        assert run.mode == "somatic"
        assert run.tumor_sample == "TUMOR"
        assert run.normal_sample == "NORMAL"
        assert run.calling_mode == "tumor_normal"

    def test_run_defaults_to_germline(self, test_db):
        run = Run(
            sample="NA12878",
            caller="GATK-HC",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
        )
        test_db.add(run)
        test_db.commit()

        assert run.mode == "germline"
        assert run.tumor_sample is None
        assert run.normal_sample is None
        assert run.calling_mode is None

    def test_tumor_only_run(self, test_db):
        run = Run(
            sample="TUMOR_001",
            caller="Mutect2",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="internal",
            mode="somatic",
            tumor_sample="TUMOR",
            calling_mode="tumor_only",
        )
        test_db.add(run)
        test_db.commit()

        assert run.calling_mode == "tumor_only"
        assert run.normal_sample is None


class TestSomaticVariantModel:
    """Tests for Variant model somatic and annotation fields."""

    def test_variant_with_somatic_fields(self, test_db, sample_run):
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr7",
            pos=55249071,
            ref="C",
            alt="T",
            type="SNP",
            classification="TP",
            tumor_dp=150,
            tumor_af=0.35,
            normal_dp=80,
            normal_af=0.001,
            somatic_quality=95.5,
            caller_filter="PASS",
        )
        test_db.add(variant)
        test_db.commit()

        assert variant.tumor_dp == 150
        assert variant.tumor_af == 0.35
        assert variant.normal_dp == 80
        assert variant.normal_af == 0.001
        assert variant.somatic_quality == 95.5
        assert variant.caller_filter == "PASS"

    def test_variant_with_stratification_fields(self, test_db, sample_run):
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=100,
            ref="AT",
            alt="A",
            type="INDEL",
            classification="TP",
            zygosity="HET",
            indel_size=1,
        )
        test_db.add(variant)
        test_db.commit()

        assert variant.zygosity == "HET"
        assert variant.indel_size == 1

    def test_variant_with_vep_annotations(self, test_db, sample_run):
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr17",
            pos=7577538,
            ref="C",
            alt="T",
            type="SNP",
            classification="TP",
            consequence="missense_variant",
            gene_symbol="TP53",
            gene_id="ENSG00000141510",
            transcript_id="ENST00000269305",
            hgvsc="c.743G>A",
            hgvsp="p.Arg248Gln",
            exon="7/11",
            protein_position=248,
            sift_prediction="deleterious(0.0)",
            polyphen_prediction="probably_damaging(1.0)",
            impact="HIGH",
        )
        test_db.add(variant)
        test_db.commit()

        assert variant.consequence == "missense_variant"
        assert variant.gene_symbol == "TP53"
        assert variant.gene_id == "ENSG00000141510"
        assert variant.transcript_id == "ENST00000269305"
        assert variant.hgvsc == "c.743G>A"
        assert variant.hgvsp == "p.Arg248Gln"
        assert variant.exon == "7/11"
        assert variant.protein_position == 248
        assert variant.sift_prediction == "deleterious(0.0)"
        assert variant.polyphen_prediction == "probably_damaging(1.0)"
        assert variant.impact == "HIGH"

    def test_somatic_fields_nullable(self, test_db, sample_run):
        variant = Variant(
            run_id=sample_run.id,
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)
        test_db.commit()

        assert variant.tumor_dp is None
        assert variant.tumor_af is None
        assert variant.normal_dp is None
        assert variant.normal_af is None
        assert variant.somatic_quality is None
        assert variant.caller_filter is None
        assert variant.zygosity is None
        assert variant.indel_size is None
        assert variant.consequence is None
        assert variant.gene_symbol is None


class TestBaselineModel:
    """Tests for Baseline model."""

    def test_create_baseline(self, test_db, sample_run):
        baseline = Baseline(
            name="v1.0_germline_GRCh38",
            run_id=sample_run.id,
            mode="germline",
            pipeline_version="v4.5.0",
        )
        test_db.add(baseline)
        test_db.commit()

        assert baseline.id is not None
        assert baseline.name == "v1.0_germline_GRCh38"
        assert baseline.locked is False
        assert baseline.created_at is not None

    def test_baseline_with_config_snapshot(self, test_db, sample_run):
        baseline = Baseline(
            name="v1.0_with_config",
            run_id=sample_run.id,
            config_snapshot="mode: germline\nreference: GRCh38",
            stratification_snapshot="dimensions:\n  variant_class: true",
        )
        test_db.add(baseline)
        test_db.commit()

        assert baseline.config_snapshot is not None
        assert baseline.stratification_snapshot is not None

    def test_baseline_locking_fields(self, test_db, sample_run):
        from datetime import UTC, datetime

        baseline = Baseline(
            name="locked_baseline",
            run_id=sample_run.id,
            locked=True,
            locked_at=datetime.now(UTC),
            approver="Dr. Smith",
            approval_comment="Approved for production",
            artifact_hash="a" * 64,
        )
        test_db.add(baseline)
        test_db.commit()

        assert baseline.locked is True
        assert baseline.approver == "Dr. Smith"
        assert baseline.artifact_hash is not None

    def test_baseline_run_relationship(self, test_db, sample_run):
        baseline = Baseline(
            name="rel_test",
            run_id=sample_run.id,
        )
        test_db.add(baseline)
        test_db.commit()

        assert baseline.run.id == sample_run.id
        assert baseline in sample_run.baselines

    def test_baseline_cascade_delete(self, test_db, sample_run):
        baseline = Baseline(
            name="cascade_test",
            run_id=sample_run.id,
        )
        test_db.add(baseline)
        test_db.commit()
        baseline_id = baseline.id

        test_db.delete(sample_run)
        test_db.commit()

        deleted = test_db.query(Baseline).filter(Baseline.id == baseline_id).first()
        assert deleted is None


class TestBaselineEnvelopeModel:
    """Tests for BaselineEnvelope model."""

    def test_create_envelope(self, test_db, sample_run):
        baseline = Baseline(name="env_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        envelope = BaselineEnvelope(
            baseline_id=baseline.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            expected_value=0.9975,
            lower_bound=0.9925,
            upper_bound=1.0,
            manually_set=False,
        )
        test_db.add(envelope)
        test_db.commit()

        assert envelope.id is not None
        assert envelope.expected_value == 0.9975
        assert envelope.manually_set is False

    def test_envelope_cascade_delete(self, test_db, sample_run):
        baseline = Baseline(name="env_cascade", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        envelope = BaselineEnvelope(
            baseline_id=baseline.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="recall",
            expected_value=0.99,
        )
        test_db.add(envelope)
        test_db.commit()
        env_id = envelope.id

        test_db.delete(baseline)
        test_db.commit()

        deleted = test_db.query(BaselineEnvelope).filter(BaselineEnvelope.id == env_id).first()
        assert deleted is None


class TestFixtureRecordModel:
    """Tests for FixtureRecord model."""

    def test_create_fixture_record(self, test_db):
        record = FixtureRecord(
            file_path="/data/truth_sets/NA12878_v4.2.1.vcf.gz",
            checksum_type="sha256",
            checksum="a" * 64,
            source="GIAB",
            genome_build="GRCh38",
            sample_identifier="NA12878",
            truth_set_version="4.2.1",
        )
        test_db.add(record)
        test_db.commit()

        assert record.id is not None
        assert record.source == "GIAB"
        assert record.date_registered is not None

    def test_fixture_record_minimal(self, test_db):
        record = FixtureRecord(
            file_path="/data/ref.fa",
            checksum="b" * 64,
        )
        test_db.add(record)
        test_db.commit()

        assert record.id is not None
        assert record.checksum_type == "sha256"
        assert record.source is None


class TestVerificationResultModel:
    """Tests for VerificationResult model."""

    def test_create_verification_result(self, test_db, sample_run):
        baseline = Baseline(name="vr_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="PASS",
            total_transitions=0,
            equivalent_count=0,
            drift_count=0,
            biological_count=0,
            envelope_breaches=0,
        )
        test_db.add(vr)
        test_db.commit()

        assert vr.id is not None
        assert vr.verdict == "PASS"
        assert vr.created_at is not None

    def test_verification_with_transitions(self, test_db, sample_run):
        baseline = Baseline(name="vr_trans", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="FAIL",
            total_transitions=2,
            biological_count=1,
            drift_count=1,
        )
        test_db.add(vr)
        test_db.flush()

        t = VariantTransition(
            verification_id=vr.id,
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            baseline_classification="TP",
            verification_classification="FN",
            transition="TP→FN",
            diff_class="biological_difference",
        )
        test_db.add(t)
        test_db.commit()

        assert len(vr.transitions) == 1
        assert vr.transitions[0].transition == "TP→FN"


class TestVariantTransitionModel:
    """Tests for VariantTransition model."""

    def test_create_transition(self, test_db, sample_run):
        baseline = Baseline(name="trans_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="REVIEW_REQUIRED",
        )
        test_db.add(vr)
        test_db.flush()

        t = VariantTransition(
            verification_id=vr.id,
            chrom="chr7",
            pos=55249071,
            ref="C",
            alt="T",
            baseline_classification="TP",
            verification_classification="FP",
            transition="TP→FP",
            diff_class="numerical_drift",
            strata_json='[["variant_class", "SNP"]]',
        )
        test_db.add(t)
        test_db.commit()

        assert t.id is not None
        assert t.diff_class == "numerical_drift"

    def test_transition_cascade_delete(self, test_db, sample_run):
        baseline = Baseline(name="trans_cascade", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="PASS",
        )
        test_db.add(vr)
        test_db.flush()

        t = VariantTransition(
            verification_id=vr.id,
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            transition="TP→FN",
        )
        test_db.add(t)
        test_db.commit()
        t_id = t.id

        test_db.delete(vr)
        test_db.commit()

        deleted = test_db.query(VariantTransition).filter(VariantTransition.id == t_id).first()
        assert deleted is None


class TestRootCauseEvidenceModel:
    """Tests for RootCauseEvidence model."""

    def test_create_evidence(self, test_db, sample_run):
        baseline = Baseline(name="ev_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="FAIL",
        )
        test_db.add(vr)
        test_db.flush()

        t = VariantTransition(
            verification_id=vr.id,
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            transition="TP→FN",
        )
        test_db.add(t)
        test_db.flush()

        ev = RootCauseEvidence(
            transition_id=t.id,
            category="coverage",
            baseline_value="30",
            verification_value="5",
            change_description="Depth dropped from 30 to 5",
            change_magnitude=0.83,
            score="STRONG",
        )
        test_db.add(ev)
        test_db.commit()

        assert ev.id is not None
        assert ev.score == "STRONG"
        assert ev.category == "coverage"

    def test_evidence_cascade_delete(self, test_db, sample_run):
        baseline = Baseline(name="ev_cascade", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        vr = VerificationResult(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            verdict="FAIL",
        )
        test_db.add(vr)
        test_db.flush()

        t = VariantTransition(
            verification_id=vr.id,
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            transition="TP→FN",
        )
        test_db.add(t)
        test_db.flush()

        ev = RootCauseEvidence(
            transition_id=t.id,
            category="filter",
            score="MODERATE",
        )
        test_db.add(ev)
        test_db.commit()
        ev_id = ev.id

        test_db.delete(t)
        test_db.commit()

        deleted = test_db.query(RootCauseEvidence).filter(RootCauseEvidence.id == ev_id).first()
        assert deleted is None


class TestModelConstraints:
    """Tests for model constraints and validations."""

    def test_run_requires_sample(self, test_db):
        """Test that Run requires sample field."""
        run = Run(
            caller="GATK-HC",
            pipeline_version="v4.5.0",
            decomposition_mode="decomposed",
            comparison_tool="happy",
            # Missing sample
        )
        test_db.add(run)

        with pytest.raises(Exception):  # SQLAlchemy IntegrityError
            test_db.commit()
        test_db.rollback()

    def test_variant_requires_foreign_key(self, test_db):
        """Test that Variant requires valid run_id."""
        variant = Variant(
            run_id=99999,  # Non-existent run_id
            chrom="chr1",
            pos=12345,
            ref="A",
            alt="T",
            type="SNP",
            classification="TP",
        )
        test_db.add(variant)

        with pytest.raises(Exception):  # Foreign key constraint
            test_db.commit()
        test_db.rollback()


class TestBayesianRiskAssessmentModel:
    """Tests for BayesianRiskAssessment model (S17.7.1)."""

    def test_create_assessment(self, test_db, sample_run):
        baseline = Baseline(name="bram_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        assessment = BayesianRiskAssessment(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            aggregate_risk_score=0.654,
            mean_risk_score=0.327,
            flagged_stratum_count=0,
            alert_threshold=0.80,
            degradation_threshold=0.01,
            aggregation_method="max",
            verdict="BRAM_PASS",
        )
        test_db.add(assessment)
        test_db.commit()

        assert assessment.id is not None
        assert assessment.verdict == "BRAM_PASS"
        assert assessment.aggregate_risk_score == 0.654
        assert assessment.created_at is not None

    def test_assessment_with_posteriors(self, test_db, sample_run):
        baseline = Baseline(name="bram_post", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        assessment = BayesianRiskAssessment(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            aggregate_risk_score=0.85,
            mean_risk_score=0.42,
            flagged_stratum_count=1,
            alert_threshold=0.80,
            degradation_threshold=0.01,
            verdict="BRAM_FLAG",
        )
        test_db.add(assessment)
        test_db.flush()

        posterior = StratumPosterior(
            assessment_id=assessment.id,
            dimension="gene_panel",
            stratum="ACMG59",
            metric_name="sensitivity",
            delta_observed=-0.02,
            prior_mu=0.0,
            prior_sigma=0.005,
            sigma_obs=0.0101,
            posterior_mu=-0.00798,
            posterior_sigma=0.00451,
            tail_probability=0.327,
            risk_weight=2.0,
            weighted_risk=0.654,
            flagged=False,
            tier="tier_1",
        )
        test_db.add(posterior)
        test_db.commit()

        assert len(assessment.stratum_posteriors) == 1
        assert assessment.stratum_posteriors[0].stratum == "ACMG59"

    def test_assessment_cascade_delete(self, test_db, sample_run):
        baseline = Baseline(name="bram_cascade", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        assessment = BayesianRiskAssessment(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            aggregate_risk_score=0.0,
            mean_risk_score=0.0,
            alert_threshold=0.80,
            degradation_threshold=0.01,
            verdict="BRAM_PASS",
        )
        test_db.add(assessment)
        test_db.flush()

        posterior = StratumPosterior(
            assessment_id=assessment.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            delta_observed=0.0,
            prior_mu=0.0,
            prior_sigma=0.005,
            sigma_obs=0.001,
            posterior_mu=0.0,
            posterior_sigma=0.001,
            tail_probability=0.02,
            risk_weight=0.5,
            weighted_risk=0.01,
            flagged=False,
        )
        test_db.add(posterior)
        test_db.commit()
        post_id = posterior.id

        test_db.delete(assessment)
        test_db.commit()

        deleted = test_db.query(StratumPosterior).filter(
            StratumPosterior.id == post_id
        ).first()
        assert deleted is None


class TestStratumPosteriorModel:
    """Tests for StratumPosterior model (S17.7.2)."""

    def test_create_posterior(self, test_db, sample_run):
        baseline = Baseline(name="sp_test", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        assessment = BayesianRiskAssessment(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            aggregate_risk_score=0.327,
            mean_risk_score=0.327,
            alert_threshold=0.80,
            degradation_threshold=0.01,
            verdict="BRAM_PASS",
        )
        test_db.add(assessment)
        test_db.flush()

        posterior = StratumPosterior(
            assessment_id=assessment.id,
            dimension="functional_impact",
            stratum="HIGH",
            metric_name="sensitivity",
            delta_observed=-0.005,
            prior_mu=0.0,
            prior_sigma=0.005,
            sigma_obs=0.008,
            posterior_mu=-0.0028,
            posterior_sigma=0.0042,
            tail_probability=0.15,
            risk_weight=2.0,
            weighted_risk=0.30,
            flagged=False,
            tier="tier_1",
        )
        test_db.add(posterior)
        test_db.commit()

        assert posterior.id is not None
        assert posterior.tier == "tier_1"
        assert posterior.assessment.verdict == "BRAM_PASS"

    def test_posterior_repr(self, test_db, sample_run):
        baseline = Baseline(name="sp_repr", run_id=sample_run.id)
        test_db.add(baseline)
        test_db.flush()

        assessment = BayesianRiskAssessment(
            run_id=sample_run.id,
            baseline_id=baseline.id,
            aggregate_risk_score=0.0,
            mean_risk_score=0.0,
            alert_threshold=0.80,
            degradation_threshold=0.01,
            verdict="BRAM_PASS",
        )
        test_db.add(assessment)
        test_db.flush()

        posterior = StratumPosterior(
            assessment_id=assessment.id,
            dimension="variant_class",
            stratum="SNP",
            metric_name="precision",
            delta_observed=0.0,
            prior_mu=0.0,
            prior_sigma=0.005,
            sigma_obs=0.001,
            posterior_mu=0.0,
            posterior_sigma=0.001,
            tail_probability=0.02,
            risk_weight=0.5,
            weighted_risk=0.01,
            flagged=False,
        )
        test_db.add(posterior)
        test_db.commit()

        repr_str = repr(posterior)
        assert "SNP" in repr_str
        assert "precision" in repr_str
