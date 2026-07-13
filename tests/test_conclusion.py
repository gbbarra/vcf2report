"""The deterministic, QC-aware conclusion (interpretive summary) of the report."""
from vcf2report.models import Annotation, Classification, QCSummary, SeqQuality, Variant
from vcf2report.report.assemble import ReportModel, summarize


def _cls(gene, tier, hpo=0.0):
    return Classification(
        variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene=gene),
        annotation=Annotation(hpo_match_score=hpo, hpo_best_match=hpo),
        criteria=[], tier=tier, rule_path="")


def _report(classifications, dp_median=44.0):
    return ReportModel(sample_id="x", hpo_terms=[], qc=QCSummary(),
                       classifications=classifications,
                       seq_quality=SeqQuality(dp_median=dp_median))


def test_diagnostic_finding_named():
    txt = " ".join(summarize(_report([_cls("SCN1A", "Pathogenic", hpo=1.0)])))
    assert "Likely explanatory finding" in txt and "SCN1A — Pathogenic" in txt


def test_no_finding_reports_vus_count():
    txt = " ".join(summarize(_report([_cls("G", "Uncertain Significance (VUS)", hpo=1.0)])))
    assert "No Pathogenic / Likely Pathogenic finding" in txt and "1 variant" in txt


def test_clinvar_pathogenic_surfaced_despite_vus():
    # A known ClinVar Pathogenic (>=2-star) the engine calls VUS MUST be flagged — never
    # presented as 'no finding' (clinical-safety; independent of the ACMG math).
    c = Classification(
        variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene="MECP2"),
        annotation=Annotation(hpo_match_score=1.0, hpo_best_match=1.0,
                              clinvar_significance="Pathogenic",
                              clinvar_review_status="criteria_provided,_multiple_submitters,_no_conflicts"),
        criteria=[], tier="Uncertain Significance (VUS)", rule_path="")
    txt = " ".join(summarize(_report([c])))
    assert "Classified Pathogenic/Likely Pathogenic in ClinVar" in txt and "MECP2" in txt


def test_single_submitter_clinvar_not_surfaced():
    # 1-star ClinVar P must NOT trip the >=2-star safety flag.
    c = Classification(
        variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene="G"),
        annotation=Annotation(hpo_match_score=0.0, hpo_best_match=0.0,
                              clinvar_significance="Pathogenic",
                              clinvar_review_status="criteria_provided,_single_submitter"),
        criteria=[], tier="Uncertain Significance (VUS)", rule_path="")
    assert "Classified Pathogenic" not in " ".join(summarize(_report([c])))


def test_secondary_sf_finding():
    # RB1 is an ACMG SF v3.2 gene; unrelated (hpo=0) P/LP -> reportable secondary.
    txt = " ".join(summarize(_report([_cls("RB1", "Pathogenic", hpo=0.0)])))
    assert "secondary finding" in txt and "RB1" in txt


def test_low_coverage_caveat_fires():
    txt = " ".join(summarize(_report([_cls("SCN1A", "Pathogenic", hpo=1.0)], dp_median=12.0)))
    assert "Coverage limitation" in txt and "12" in txt


def test_adequate_coverage_note():
    txt = " ".join(summarize(_report([_cls("SCN1A", "Pathogenic", hpo=1.0)], dp_median=44.0)))
    assert "adequate" in txt and "Coverage limitation" not in txt


def test_always_recommends_next_steps_and_single_proband():
    lines = summarize(_report([_cls("SCN1A", "Pathogenic", hpo=1.0)]))
    joined = " ".join(lines)
    assert "Recommended next steps" in joined and "Single-proband" in joined


def test_renders_in_report():
    from vcf2report import config, pipeline
    from vcf2report.report.render import render_markdown
    md = render_markdown(pipeline.run_pipeline(str(config.SAMPLE_VCF)))
    assert "## Conclusion (draft interpretation)" in md
