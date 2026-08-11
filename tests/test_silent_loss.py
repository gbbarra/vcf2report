"""Silent-loss defects: the pipeline runs, the suite is green, and a variant is gone.

Every test here pins a path where vcf2report used to drop or mis-attribute a variant
WITHOUT saying so — the failure mode that a passing test suite is worst at catching,
because nothing crashes and the report looks complete. Each one asserts the *observable
consequence* (what the report says, which bucket a variant lands in), not the internals.
"""
from __future__ import annotations

import gzip
import textwrap

import pytest

from vcf2report.annotate import from_vcf
from vcf2report.models import Annotation, Classification, QCSummary, SeqQuality, Variant
from vcf2report.pipeline import run_pipeline
from vcf2report.report.assemble import (ReportModel, clinvar_pathogenic_flags,
                                        is_hom_absent_artifact, is_hom_gnomad_uncovered,
                                        phenotype_compared, split_findings, summarize)
from vcf2report.vcf import annparse
from vcf2report.vcf.filter import is_impactful
from vcf2report.vcf.parse import parse_vcf

HDR = ("##fileformat=VCFv4.2\n##reference=GRCh38\n"
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")


def _write(tmp_path, name, body, header=HDR):
    p = tmp_path / name
    p.write_text(header + textwrap.dedent(body).lstrip())
    return p


def _report(classifications, qc=None):
    return ReportModel(sample_id="x", hpo_terms=[], qc=qc or QCSummary(),
                       classifications=classifications,
                       seq_quality=SeqQuality(dp_median=44.0))


def _cls(gene, tier, **ann):
    return Classification(
        variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene=gene,
                        zygosity=ann.pop("zygosity", "het")),
        annotation=Annotation(**ann), criteria=[], tier=tier, rule_path="")


# --------------------------------------------------------------------------------------
# 1-3. The QC gate is the one stage that can remove everything and still look clean.
# --------------------------------------------------------------------------------------

def test_qc_dropped_clinvar_pathogenic_is_named_not_deleted(tmp_path, monkeypatch):
    """One GQ point below threshold used to delete a known pathogenic call from the report
    while the conclusion announced 'no finding'. The report's do-not-dismiss net cannot
    reach it — QC runs before annotation — so the pipeline must flag it itself."""
    src = "data/example/SYN-001.synthetic.vcf.gz"
    lines = gzip.open(src, "rt").read().splitlines()
    out = []
    for ln in lines:
        if not ln.startswith("#") and "SCN1A" in ln:
            c = ln.split("\t")
            f, s = c[8].split(":"), c[9].split(":")
            s[f.index("GQ")] = "19"                      # threshold is 20
            c[9] = ":".join(s)
            ln = "\t".join(c)
        out.append(ln)
    p = tmp_path / "gq19.vcf"
    p.write_text("\n".join(out) + "\n")

    r = run_pipeline(p, hpo_terms=["HP:0001250", "HP:0002133", "HP:0001263"])
    assert "SCN1A" not in [c.variant.gene for c in r.classifications]   # still dropped...
    assert any("SCN1A" in s for s in r.qc.qc_rescued)                   # ...but NAMED
    txt = " ".join(summarize(r))
    assert "Removed by QC before classification" in txt and "SCN1A" in txt
    # ...and it must reach the rendered laudo, not only the model. BOTH renderers: the
    # Jinja template is used when jinja2 is installed and templates/report.md.j2 exists,
    # the built-in otherwise, and a section added to only one of them is invisible to half
    # the installs.
    from vcf2report.report.render import _render_markdown_builtin, render_markdown
    for md in (render_markdown(r), _render_markdown_builtin(r)):
        assert "Removed by QC, but known to ClinVar" in md and "GQ=19<20" in md


def test_sites_only_vcf_warns_instead_of_reporting_no_finding(tmp_path):
    """A VCF with no FORMAT/sample column loses 100% of variants at the carrier gate."""
    p = tmp_path / "sites.vcf"
    p.write_text("##fileformat=VCFv4.2\n##reference=GRCh38\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr1\t400\t.\tA\tG\t50\tPASS\tGENE=BRCA1;CSQ=stop_gained\n"
                 "chr1\t401\t.\tC\tT\t50\tPASS\tGENE=BRCA2;CSQ=frameshift_variant\n")
    w = " ".join(run_pipeline(p).qc.warnings)
    assert "ALL 2 variants were removed" in w
    assert "sites-only VCF" in w and "NO evidential weight" in w


def test_unannotated_vcf_says_the_empty_shortlist_is_not_biological(tmp_path):
    p = _write(tmp_path, "bare.vcf", """
        chr1\t400\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:40:60\t
        chr1\t401\t.\tC\tT\t50\tPASS\t.\tGT:DP:GQ\t0/1:40:60\t
        """)
    w = " ".join(run_pipeline(p).qc.warnings)
    assert "consequence annotation" in w and "NOT" in w


def test_qc_rescue_is_quiet_on_a_real_exome():
    """The rescue must not become noise: a real callset has ~100 borderline coding drops
    and none of them are ClinVar-pathogenic. Run without a phenotype, so the two
    phenotype-gated rescue paths cannot fire either — the net stays silent."""
    r = run_pipeline("data/example/SYN-004.NIPBL.annotated.vcf.gz")
    assert r.qc.qc_rescued == []


_SEIZURE_HPO = ["HP:0001250", "HP:0002133", "HP:0001263"]   # SCN1A scores high on these


def test_qc_dropped_reviewed_vus_in_phenotype_gene_is_rescued(tmp_path):
    """The MYH11 case: a well-reviewed (2★) ClinVar VUS dropped at low GQ, in a gene that
    matches the indication, must be surfaced — the P/LP-only bar missed exactly this."""
    p = _write(tmp_path, "vus.vcf", """
        chr2\t166000000\t.\tC\tT\t50\tPASS\tGENE=SCN1A;CSQ=missense_variant;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;gnomad_AF=0.000003\tGT:AD:DP:GQ\t0/1:20,20:40:14\t
        """)
    r = run_pipeline(p, hpo_terms=_SEIZURE_HPO)             # SCN1A matches seizures
    assert "SCN1A" not in [c.variant.gene for c in r.classifications]      # still dropped...
    assert any("SCN1A" in s and "well-reviewed non-benign" in s for s in r.qc.qc_rescued)

    # ...but the SAME reviewed VUS with NO phenotype supplied must stay quiet (not noise).
    r0 = run_pipeline(p, hpo_terms=[])
    assert not any("SCN1A" in s for s in r0.qc.qc_rescued)


def test_qc_dropped_rare_variant_in_phenotype_gene_is_rescued_without_clinvar(tmp_path):
    """The general case: a rare, QC-dropped coding variant in a phenotype-matched gene, with
    NO ClinVar record at all, is still worth a confirm-before-dismiss flag."""
    p = _write(tmp_path, "rare.vcf", """
        chr2\t166000000\t.\tC\tT\t50\tPASS\tGENE=SCN1A;CSQ=missense_variant;gnomad_AF=0.000003\tGT:AD:DP:GQ\t0/1:20,20:40:14\t
        """)
    r = run_pipeline(p, hpo_terms=_SEIZURE_HPO)
    assert any("SCN1A" in s and "matches the phenotype" in s for s in r.qc.qc_rescued)
    # no phenotype -> the rare-in-gene path cannot fire (would be noise on a real callset)
    assert run_pipeline(p, hpo_terms=[]).qc.qc_rescued == []


# --------------------------------------------------------------------------------------
# 4-5. Absence of data must never be rendered as evidence.
# --------------------------------------------------------------------------------------

def test_af_zero_without_an_is_not_a_vouched_absence():
    vouched = _cls("G1", "Pathogenic", zygosity="hom", gnomad_af=0.0, gnomad_ac=0,
                   gnomad_an=152000, hpo_match_score=0.9)
    uncovered = _cls("G2", "Pathogenic", zygosity="hom", gnomad_af=0.0, gnomad_ac=0,
                     gnomad_an=0, hpo_match_score=0.9)
    assert is_hom_absent_artifact(vouched) and not is_hom_gnomad_uncovered(vouched)
    assert is_hom_gnomad_uncovered(uncovered) and not is_hom_absent_artifact(uncovered)


def test_no_live_artifact_tag_rests_on_missing_data():
    """Regression: every hom-absent tag on the committed exomes had AN=0, i.e. AF was 0/0."""
    for f in ("data/example/SYN-070.RBSN.annotated.vcf.gz",
              "data/example/SYN-004.NIPBL.annotated.vcf.gz"):
        for c in run_pipeline(f).classifications:
            if is_hom_absent_artifact(c):
                assert c.annotation.gnomad_an, f"{c.variant.gene}: tagged absent with AN=0"


def test_phenotype_matched_homozygote_in_uncovered_region_is_not_demoted():
    """SYN-070's own planted candidate: hom, phenotype 1.0, ClinVar Pathogenic, AM 0.97,
    demoted to 'Not routinely reported' purely because gnomAD does not cover the site."""
    r = run_pipeline("data/example/SYN-070.RBSN.annotated.vcf.gz",
                     hpo_terms=["HP:0001510", "HP:0001903", "HP:0001873"])
    primary, _sec, other = split_findings(r.classifications)
    rbsn = [c for c in r.classifications if c.variant.gene == "RBSN"]
    assert rbsn, "the planted RBSN candidate is missing from the report"
    assert rbsn[0] in primary and rbsn[0] not in other


# --------------------------------------------------------------------------------------
# 6-7. A comparison that never ran must not be reported as a comparison that failed.
# --------------------------------------------------------------------------------------

def test_missing_hpo_score_is_not_a_non_match():
    plp = _cls("SCN1A", "Pathogenic", hpo_match_score=None)
    assert not phenotype_compared([plp])
    primary, _sec, other = split_findings([plp])
    assert plp in primary and plp not in other


def test_run_without_phenotype_does_not_claim_a_phenotype_mismatch():
    r = run_pipeline("data/sample/sample_exome.vcf", hpo_terms=[])
    txt = " ".join(summarize(r))
    assert "No phenotype was available" in txt
    assert "not matching the stated phenotype" not in txt
    assert "in phenotype-matched genes" not in txt
    primary, _sec, _other = split_findings(r.classifications)
    assert "SCN1A" in [c.variant.gene for c in primary]


def test_phenotype_claims_survive_when_a_phenotype_was_supplied():
    """The guard must not mute the real phenotype language on a normal run."""
    r = run_pipeline("data/example/SYN-001.synthetic.vcf.gz",
                     hpo_terms=["HP:0001250", "HP:0002133", "HP:0001263"])
    txt = " ".join(summarize(r))
    assert phenotype_compared(r.classifications)
    assert "No phenotype was available" not in txt
    assert "phenotype overlap NOT assessed" not in txt


# --------------------------------------------------------------------------------------
# 8. An actionable secondary finding must not be suppressed by a QC caution.
# --------------------------------------------------------------------------------------

def test_acmg_sf_gene_is_not_suppressed_by_the_artifact_gate():
    sf = _cls("RB1", "Pathogenic", zygosity="hom", gnomad_af=0.0, gnomad_ac=0,
              gnomad_an=152000, hpo_match_score=0.0)
    _p, secondary, other = split_findings([sf])
    assert sf in secondary and sf not in other
    txt = " ".join(summarize(_report([sf])))
    assert "not on the ACMG SF actionable list" not in txt
    assert "secondary finding" in txt


# --------------------------------------------------------------------------------------
# 9-10. Per-allele vs per-transcript, and the fallback that never ran.
# --------------------------------------------------------------------------------------

def test_per_allele_insilico_scores_are_not_shared_between_alleles(tmp_path):
    """am_pathogenicity="0.03,0.98" at a 2-ALT site used to give BOTH alleles 0.98 —
    the benign allele inheriting the pathogenic allele's evidence."""
    p = _write(tmp_path, "multi.vcf", """
        chr1\t100\t.\tA\tG,T\t50\tPASS\tam_pathogenicity=0.03,0.98;GENE=X;CSQ=missense_variant\tGT:AD\t1/2:0,20,20\t
        """)
    scores = {v.alt: from_vcf.extract(v).get("am_pathogenicity") for v in parse_vcf(p)[0]}
    assert scores == {"G": 0.03, "T": 0.98}


def test_per_transcript_scores_still_aggregate_to_the_max(tmp_path):
    """The disambiguator is comma arity vs ALT count — dbNSFP's per-transcript list at a
    single-ALT site must still take the max, not element 0."""
    p = _write(tmp_path, "single.vcf", """
        chr1\t100\t.\tA\tG\t50\tPASS\tREVEL=0.1,0.9,0.4;GENE=X;CSQ=missense_variant\tGT:AD\t0/1:20,20\t
        """)
    assert from_vcf.extract(parse_vcf(p)[0][0])["revel"] == 0.9


def test_short_per_allele_array_falls_back_instead_of_claiming_no_frequency(tmp_path):
    """A per-allele array that does not cover this ALT used to be recorded as
    _source='VCF INFO' with AF=None, silently skipping the local gnomAD snapshot."""
    from vcf2report.annotate import annotate_variant
    p = _write(tmp_path, "short.vcf", """
        chr1\t200\t.\tA\tG,T\t50\tPASS\tgnomad_AF=0.42,.;GENE=X;CSQ=missense_variant\tGT:AD\t1/2:0,20,20\t
        """)
    alt2 = [v for v in parse_vcf(p)[0] if v.alt == "T"][0]
    assert from_vcf.extract(alt2).get("gnomad_af") is None
    ann = annotate_variant(alt2, [])
    assert ann.source.get("gnomad") != "VCF INFO"       # the fallback actually ran
    # ALT #1, which the array DOES cover, still reads straight from INFO.
    alt1 = [v for v in parse_vcf(p)[0] if v.alt == "G"][0]
    assert annotate_variant(alt1, []).source.get("gnomad") == "VCF INFO"


# --------------------------------------------------------------------------------------
# 11-12. Annotator dialects: consequences that were parsed or filtered wrongly.
# --------------------------------------------------------------------------------------

def test_vep_written_ann_is_not_parsed_with_snpeff_offsets(tmp_path):
    """VEP --vcf_info_field ANN writes VEP's field order under the ANN key. SnpEff offset
    10 is HGVS.p but VEP offset 10 is HGVSc, so the protein column got a 'c.' string."""
    hdr = ('##fileformat=VCFv4.2\n##reference=GRCh38\n'
           '##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from '
           'Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|'
           'BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position">\n'
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    p = _write(tmp_path, "vepann.vcf", """
        chr17\t43093464\t.\tA\tG\t50\tPASS\tANN=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000357654|protein_coding|10/23|9/22|ENST00000357654.9:c.1067A>G|ENSP00000350283.3:p.Gln356Arg|1067|356|356\tGT:AD:DP:GQ\t0/1:20,20:40:60\t
        """, header=hdr)
    v = parse_vcf(p)[0][0]
    assert v.gene == "BRCA1" and v.consequence == "missense_variant"
    assert v.hgvs_c == "ENST00000357654.9:c.1067A>G"
    assert v.hgvs_p == "ENSP00000350283.3:p.Gln356Arg"
    assert v.exon == "10/23"


def test_snpeff_ann_is_still_parsed_positionally():
    """The VEP detection must not disturb the SnpEff fast path (no ANN header declared)."""
    r = annparse.parse_snpeff(
        "G|missense_variant|MODERATE|RBSN|RBSN|transcript|NM_022340.4|protein_coding|5/14|"
        "c.289G>C|p.Gly97Arg", "G", "C", 1)
    assert r["gene"] == "RBSN" and r["hgvs_c"] == "c.289G>C" and r["hgvs_p"] == "p.Gly97Arg"


@pytest.mark.parametrize("csq,first", [
    ("stop_gained&splice_region_variant", "stop_gained"),
    ("missense_variant&splice_region_variant", "missense_variant"),
    ("splice_acceptor_variant&intron_variant", "splice_acceptor_variant"),
])
def test_compound_consequence_on_the_plain_key_path_is_split(csq, first):
    """The ANN/CSQ paths take the first (most severe) term; the plain-key path passed the
    whole '&'-joined string through, so is_impactful rejected a stop_gained."""
    r = annparse.extract({"GENE": "X", "CSQ": csq}, "G", None, "A", 0, 1)
    assert r["consequence"] == first and is_impactful(r["consequence"])


@pytest.mark.parametrize("term", ["transcript_ablation", "exon_loss_variant",
                                  "protein_altering_variant"])
def test_whole_transcript_loss_terms_reach_the_shortlist(term):
    """annparse translates SnpEff EXON_DELETED to 'transcript_ablation' — the filter then
    discarded the very term the parser had just produced."""
    assert is_impactful(term)


def test_eff_exon_deleted_survives_end_to_end():
    r = annparse.parse_snpeff_eff("EXON_DELETED(HIGH||||NIPBL|protein_coding|CODING|NM_133433.4|5|1)")
    assert r["consequence"] == "transcript_ablation" and is_impactful(r["consequence"])


# --------------------------------------------------------------------------------------
# 13. The same ClinVar test guards two different safety mechanisms.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("sig", ["Likely_pathogenic", "Likely pathogenic",
                                 "Pathogenic/Likely_pathogenic", "Pathogenic"])
def test_raw_clinvar_underscores_do_not_disable_the_safety_net(sig):
    """clinvar_stars normalises underscores; the significance test beside it did not — so a
    raw CLNSIG value was dropped by the rarity filter AND invisible to the net at once."""
    from vcf2report.vcf.filter import _is_clinvar_plp
    c = _cls("G", "Uncertain Significance (VUS)", clinvar_significance=sig,
             clinvar_review_status="reviewed_by_expert_panel")
    assert clinvar_pathogenic_flags([c]) == [c]
    assert _is_clinvar_plp(c.annotation)


# --------------------------------------------------------------------------------------
# 14. Allele balance at a multiallelic site.
# --------------------------------------------------------------------------------------

def test_allele_balance_uses_the_called_alleles_not_the_whole_site(tmp_path):
    """A 1/2 call with AD=0,20,45 measured each allele against the site total (0.235 and
    0.53), so the QC gate dropped a perfectly balanced compound-heterozygote."""
    p = _write(tmp_path, "ab.vcf", """
        chr1\t100\t.\tA\tG,T\t50\tPASS\tGENE=X;CSQ=missense_variant\tGT:AD:DP:GQ\t1/2:0,20,45:65:60\t
        """)
    ab = {v.alt: v.allele_balance for v in parse_vcf(p)[0]}
    assert 0.25 <= ab["G"] <= 0.75 and 0.25 <= ab["T"] <= 0.75

    from vcf2report.vcf.qc import apply_qc
    kept, _dropped = apply_qc(parse_vcf(p)[0])
    assert len(kept) == 2


def test_biallelic_allele_balance_is_unchanged(tmp_path):
    p = _write(tmp_path, "ab2.vcf", """
        chr1\t100\t.\tA\tG\t50\tPASS\tGENE=X;CSQ=missense_variant\tGT:AD:DP:GQ\t0/1:20,20:40:60\t
        """)
    assert parse_vcf(p)[0][0].allele_balance == 0.5


# --------------------------------------------------------------------------------------
# 15-18. Judgement calls, once decided, become assertions like everything else.
# --------------------------------------------------------------------------------------

def test_near_splice_exclusion_is_disclosed_not_silent():
    """splice_region variants stay out of the shortlist by design — so the funnel must
    say how many were set aside, in both renderers."""
    from vcf2report.report.render import _render_markdown_builtin, render_markdown
    r = run_pipeline("data/example/SYN-004.NIPBL.annotated.vcf.gz")
    assert r.qc.near_splice_excluded > 0
    for md in (render_markdown(r), _render_markdown_builtin(r)):
        assert "Set aside by design" in md
        assert str(r.qc.near_splice_excluded) in md
        assert "stated sensitivity limit" in md


def test_known_pathogenic_splice_region_still_reaches_the_shortlist():
    """The exclusion is only defensible because ClinVar P/LP bypasses the impact step."""
    from vcf2report.vcf.filter import filter_variants
    v = Variant(chrom="1", pos=1, ref="A", alt="G", gene="X",
                consequence="splice_region_variant", zygosity="het")
    known, _f = filter_variants([(v, Annotation(clinvar_significance="Pathogenic",
                                                gnomad_af=0.0))])
    novel, _f = filter_variants([(v, Annotation(gnomad_af=0.0))])
    assert len(known) == 1 and len(novel) == 0


def test_per_allele_frequency_is_not_broadcast_at_a_multiallelic_site(tmp_path):
    """A lone Number=A value at a 2-ALT site resolved only one allele; giving ALT #2 ALT
    #1's frequency can flip BA1/BS1. Falls back to the snapshot instead."""
    p = _write(tmp_path, "scalar.vcf", """
        chr1\t100\t.\tA\tG,T\t50\tPASS\tgnomad_AF=0.42;GENE=X;CSQ=missense_variant\tGT:AD\t1/2:0,20,20\t
        """)
    got = {v.alt: from_vcf.extract(v).get("gnomad_af") for v in parse_vcf(p)[0]}
    assert got == {"G": 0.42, "T": None}


def test_a_site_level_faf95_is_not_attributed_to_an_individual_allele(tmp_path):
    """This test previously asserted the opposite, and the premise was wrong.

    fafmax_faf95_max is per-SITE: the MAXIMUM filtering AF across the site's alleles. But
    `_benign_af` returns faf95 ahead of the allele's own AF, so handing the site maximum to
    every ALT called a rare allele Benign on a neighbour's frequency — an allele at AF 1e-5
    met BA1 from a site faf95 of 0.29, in a trail that also reported PM2 ("absent from
    population databases"). A filtering AF is derived from one allele's AC/AN, so a site
    maximum cannot be attributed to a specific allele: leave it unset and let BA1/BS1 use
    this allele's own AF."""
    p = _write(tmp_path, "faf.vcf", """
        chr1\t100\t.\tA\tG,T\t50\tPASS\tgnomad_AF=0.4,0.01;gnomad_faf95=0.3;GENE=X;CSQ=missense_variant\tGT:AD\t1/2:0,20,20\t
        """)
    got = {v.alt: from_vcf.extract(v).get("gnomad_faf95") for v in parse_vcf(p)[0]}
    assert got == {"G": None, "T": None}


def test_a_single_alt_site_still_uses_its_faf95(tmp_path):
    """With one ALT the site value and the allele value are the same thing — the stricter
    rule must not throw away the filtering AF that BA1/BS1 are calibrated on."""
    p = _write(tmp_path, "faf1.vcf", """
        chr1\t100\t.\tA\tG\t50\tPASS\tgnomad_AF=0.4;gnomad_faf95=0.3;GENE=X;CSQ=missense_variant\tGT:AD\t0/1:20,20\t
        """)
    assert from_vcf.extract(parse_vcf(p)[0][0])["gnomad_faf95"] == 0.3


def test_a_per_allele_faf95_array_is_indexed_by_allele(tmp_path):
    """An annotator that DOES emit one faf95 per ALT should be believed, per allele."""
    p = _write(tmp_path, "faf2.vcf", """
        chr1\t100\t.\tA\tG,T\t50\tPASS\tgnomad_AF=0.4,0.01;gnomad_faf95=0.3,0.008;GENE=X;CSQ=missense_variant\tGT:AD\t1/2:0,20,20\t
        """)
    got = {v.alt: from_vcf.extract(v).get("gnomad_faf95") for v in parse_vcf(p)[0]}
    assert got == {"G": 0.3, "T": 0.008}


def test_half_call_is_kept_as_a_flagged_het_not_dropped(tmp_path):
    """'./1' certainly carries the ALT; only the zygosity is unknown. Dropping it as a
    non-carrier deleted a real variant."""
    from vcf2report.vcf.qc import apply_qc
    p = _write(tmp_path, "half.vcf", """
        chr1\t300\t.\tA\tG\t50\tPASS\tGENE=X;CSQ=stop_gained\tGT:AD:DP:GQ\t./1:0,25:25:60\t
        """)
    v = parse_vcf(p)[0][0]
    assert v.zygosity == "het" and v.partial_call is True
    kept, dropped = apply_qc([v])
    assert kept == [v] and not dropped        # the het AB window must not evict it

    from vcf2report.report.render import _render_markdown_builtin, _zyg
    assert _zyg(v) == "het (partial call)"


@pytest.mark.parametrize("gt,zyg,partial", [
    ("0/1", "het", False), ("1/1", "hom", False),
    ("./1", "het", True),         # carries ALT #1; zygosity unknown -> conservative het
    ("1|.", "het", True),
    ("./.", None, False),         # nothing called at all: not a partial call
    ("./0", None, True),          # partial GT, but no evidence for ANY alt -> non-carrier
    ("./2", None, True),          # carries ALT #2, not ALT #1 — non-carrier for this row
])
def test_zygosity_and_partial_call_across_genotypes(tmp_path, gt, zyg, partial):
    p = _write(tmp_path, f"gt_{abs(hash(gt))}.vcf", f"""
        chr1\t300\t.\tA\tG,C\t50\tPASS\tGENE=X;CSQ=stop_gained\tGT:AD:DP:GQ\t{gt}:0,25,25:50:60\t
        """)
    v = [x for x in parse_vcf(p)[0] if x.alt == "G"][0]
    assert v.zygosity == zyg and v.partial_call is partial


def test_every_demoted_hom_absent_variant_is_explained():
    """The caveat used to be gated on _PLP — the same condition that EXEMPTS a variant from
    the demotion — so a phenotype-matched hom VUS was demoted with no explanation."""
    from vcf2report.report.assemble import split_findings
    vus = _cls("GENEV", "Uncertain Significance (VUS)", zygosity="hom", gnomad_af=0.0,
               gnomad_ac=0, gnomad_an=152000, hpo_match_score=0.9, hpo_best_match=0.9)
    _p, _s, other = split_findings([vus])
    assert vus in other                                   # still demoted (healthy-exome guard)
    assert "Verify the genotype" in " ".join(summarize(_report([vus])))   # ...but explained


def test_benign_hom_absent_variants_do_not_generate_caveat_noise():
    benign = _cls("GENEB", "Benign", zygosity="hom", gnomad_af=0.0, gnomad_ac=0,
                  gnomad_an=152000, hpo_match_score=0.0)
    assert "Verify the genotype" not in " ".join(summarize(_report([benign])))


def test_a_store_vouched_absence_is_not_mistaken_for_missing_coverage():
    """Regression on the fix itself. `is_hom_absent_artifact` first required gnomad_an > 0, but
    the Parquet store's vouched-absence sentinel has no gnomAD row and therefore no AN. That
    silently disabled the calling-artifact guard on the production fast path, while the report
    simultaneously told the reader gnomAD "does not survey" a locus the store had vouched for.

    Three states must stay distinguishable:
      store covered the locus, no record   -> vouched absence  (artifact guard applies)
      annotator wrote AF=0 with AN=0       -> no coverage      (weaker caveat, no demotion)
      annotator wrote AF=0 with a real AN  -> vouched absence  (artifact guard applies)
    """
    from vcf2report.annotate import annotate_variant, gnomad_parquet
    from vcf2report.report.assemble import split_findings

    v = Variant(chrom="5", pos=157000000, ref="C", alt="CA", gene="CYFIP2",
                consequence="frameshift_variant", zygosity="hom")

    # The sentinel is taken from the store module itself, so this test tracks the store's
    # contract rather than a copy of it that can drift.
    gnomad_parquet._primed[v.key] = {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0,
                                     "pop": None, "vouched_absent": True}
    try:
        ann = annotate_variant(v, [])
        assert ann.gnomad_absence_vouched is True
        c = Classification(variant=v, annotation=ann, criteria=[], tier="Likely Pathogenic",
                           rule_path="")
        assert is_hom_absent_artifact(c) and not is_hom_gnomad_uncovered(c)
        _p, _s, other = split_findings([c])
        assert c in other                      # the artifact demotion actually happens
    finally:
        gnomad_parquet._primed.clear()

    uncovered = Classification(
        variant=v, annotation=Annotation(gnomad_af=0.0, gnomad_ac=0, gnomad_an=0),
        criteria=[], tier="Likely Pathogenic", rule_path="")
    assert not is_hom_absent_artifact(uncovered) and is_hom_gnomad_uncovered(uncovered)

    surveyed = Classification(
        variant=v, annotation=Annotation(gnomad_af=0.0, gnomad_ac=0, gnomad_an=125000),
        criteria=[], tier="Likely Pathogenic", rule_path="")
    assert is_hom_absent_artifact(surveyed) and not is_hom_gnomad_uncovered(surveyed)


def test_a_filter_dropped_known_pathogenic_is_named(tmp_path):
    """A caller FILTER is a categorical rejection, so the variant is correctly not re-admitted
    as a candidate — but it must not vanish. A 3-star expert-panel ClinVar Pathogenic dropped
    on FILTER used to leave no trace anywhere in the report."""
    p = _write(tmp_path, "filterdrop.vcf", """
        chr11\t2444000\t.\tC\tT\t50\tLowQual\tGENE=KCNQ2;CSQ=stop_gained;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT:AD:DP:GQ\t0/1:20,20:40:60\t
        """)
    r = run_pipeline(p)
    assert "KCNQ2" not in [c.variant.gene for c in r.classifications]   # not re-admitted...
    assert any("KCNQ2" in s and "FILTER=LowQual" in s for s in r.qc.qc_rescued)   # ...but named
    from vcf2report.report.render import _render_markdown_builtin, render_markdown
    for md in (render_markdown(r), _render_markdown_builtin(r)):
        assert "KCNQ2" in md


# --------------------------------------------------------------------------------------
# 19-22. The annotate-layer follow-ups: a store that cannot be asked must not answer.
# --------------------------------------------------------------------------------------

def test_a_multiallelic_site_does_not_share_its_clinvar_assertion(tmp_path):
    """CLNSIG carries literal commas so it cannot be allele-indexed. Applying the site value
    to every ALT was not a harmless approximation: a 21%-frequency synonymous allele inherited
    a 2-star Pathogenic assertion, bypassed the rarity gate through filter.py's ClinVar P/LP
    rescue, and reached the report's do-not-dismiss net."""
    p = _write(tmp_path, "cvmulti.vcf", """
        chr2\t166003360\t.\tC\tT,A\t50\tPASS\tGENE=SCN1A;CSQ=missense_variant;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;gnomad_AF=0.00001,0.21\tGT:AD\t1/2:0,20,20\t
        """)
    for v in parse_vcf(p)[0]:
        assert from_vcf.extract(v).get("clinvar_significance") is None


def test_a_single_alt_site_still_reads_clinvar_from_info(tmp_path):
    """The stricter rule must not disable the pre-annotated fast path for ordinary records."""
    p = _write(tmp_path, "cvone.vcf", """
        chr2\t166003360\t.\tC\tT\t50\tPASS\tGENE=SCN1A;CSQ=missense_variant;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;gnomad_AF=0.00001\tGT:AD\t0/1:20,20\t
        """)
    got = from_vcf.extract(parse_vcf(p)[0][0])
    assert got["clinvar_significance"] == "Pathogenic"


@pytest.mark.parametrize("key,canonical", [
    ("1-100-A-T", True),                 # SNV
    ("1-100-CAG-C", True),               # already trimmed deletion
    ("1-100-C-CAG", True),               # already trimmed insertion
    ("1-100-CAGG-CG", False),            # shared leading AND trailing base
    ("1-100-AGA-A", True),               # trimmed
    ("1-100-CCA-CTA", False),            # shared leading base, both length>1
])
def test_untrimmed_alleles_are_recognised(key, canonical):
    """A store is keyed on normalised alleles, so an un-normalised key is a different string
    for the same variant. Vouching absence for it converts a representation difference into a
    positive false observation — a 31%-frequency indel read as surveyed-zero."""
    from vcf2report.annotate.gnomad_parquet import _is_canonical
    assert _is_canonical(key) is canonical


def test_pm2_does_not_claim_local_cohort_was_consulted_when_it_was_not():
    """README promised a variant absent from the local cohort table is treated as unknown. The
    decision converted the honest None to 0.0, and the criterion name advertised both
    databases — across the real inputs, 100% of PM2 grants rested on an unsurveyed leg.
    PM2 still fires (gnomAD absence is real evidence) but no longer overstates its basis."""
    from vcf2report.acmg.criteria import all_criteria

    pm2 = all_criteria()["PM2"]
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="SCN1A")

    unchecked = pm2(v, Annotation(gnomad_af=0.0, gnomad_an=152000, local_cohort_af=None))
    assert unchecked.met is True                       # still fires on gnomAD alone
    assert "local cohort not consulted" in unchecked.name
    assert unchecked.confidence == "moderate"
    assert unchecked.evidence["local_cohort_checked"] is False

    checked = pm2(v, Annotation(gnomad_af=0.0, gnomad_an=152000, local_cohort_af=0.0))
    assert checked.met is True and "gnomAD + local cohort" in checked.name
    assert checked.confidence == "high" and checked.evidence["local_cohort_checked"] is True

    common_in_brazil = pm2(v, Annotation(gnomad_af=0.0, gnomad_an=152000, local_cohort_af=0.032))
    assert common_in_brazil.met is False               # the differentiator still works
