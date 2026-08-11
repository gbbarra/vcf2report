"""The benchmark's --withhold-clinvar re-classification (run_benchmark.withhold_own_clinvar_tier).

The hpo-spiked-exomes benchmark plants KNOWN ClinVar pathogenics, so PP5 (+ the ClinVar safety
flag) recover them no matter what the residue/constraint criteria do — masking whether PS1/PM5/PP2
can carry a variant alone. Nulling the plant's own ClinVar assertion turns it into the novel-variant
scenario those criteria exist for. This locks the re-classification logic: with ClinVar removed, a
variant at an established pathogenic residue is recovered to P/LP by PS1 (which the own-record guard
had suppressed while the ClinVar assertion was present).
"""

import importlib.util
from pathlib import Path

from vcf2report.acmg.engine import classify
from vcf2report.models import Annotation, Variant

_RB = Path(__file__).resolve().parent.parent / "scripts" / "run_benchmark.py"
_spec = importlib.util.spec_from_file_location("run_benchmark", _RB)
rb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rb)


def _kcnq2_arg213trp():
    # A missense at an established pathogenic residue in a missense-constrained gene, rare in
    # gnomAD, and itself catalogued by ClinVar as Pathogenic — the typical benchmark plant.
    v = Variant(
        chrom="20",
        pos=63446204,
        ref="G",
        alt="A",
        gene="KCNQ2",
        consequence="missense_variant",
        hgvs_p="p.Arg213Trp",
    )
    ps1 = {
        "alt_aa": "Trp",
        "ref_aa": "Arg",
        "stars": 2,
        "genomic_key": "20-63444712-G-A",
        "accession": "VCV000021795",
    }
    a = Annotation(
        gnomad_af=0.0,
        gnomad_faf95=0.0,
        clinvar_significance="Pathogenic",
        clinvar_review_status="criteria provided, single submitter",
        clinvar_accession="VCV1",
        clinvar_ps1=ps1,
        clinvar_residue_available=True,
        gene_mis_z=5.2,
        gene_missense_constrained=True,
    )  # PP2
    return classify(v, a)


def test_normal_run_is_pp5_driven_ps1_suppressed():
    # With the ClinVar assertion present: PP5 fires, PS1 is withheld (double-count guard).
    c = _kcnq2_arg213trp()
    codes = {x.code: x for x in c.criteria}
    assert codes["PP5"].met is True
    assert codes["PS1"].met is False


def test_withhold_flips_pp5_to_ps1():
    # The mechanism, tier-independent: nulling the own ClinVar assertion turns PP5 off and lets
    # PS1 fire from the residue index (no longer suppressed as a double-count).
    from dataclasses import replace

    c = _kcnq2_arg213trp()
    ann = replace(
        c.annotation,
        clinvar_significance=None,
        clinvar_review_status=None,
        clinvar_accession=None,
        clinvar_condition=None,
    )
    codes = {x.code: x for x in classify(c.variant, ann).criteria}
    assert codes["PP5"].met is False
    assert codes["PS1"].met is True


def test_withhold_recovers_plp_via_ps1():
    # Null the own ClinVar assertion -> PS1(strong) + PM2(supporting) + PP2(supporting) reaches
    # Likely Pathogenic on the residue/constraint/rarity evidence alone. This is what the mode measures.
    c = _kcnq2_arg213trp()
    withheld = rb.withhold_own_clinvar_tier(c)
    assert withheld.startswith(("Pathogenic", "Likely Pathogenic")), withheld


def test_withhold_does_not_mutate_original():
    # dataclasses.replace must leave the original classification's annotation untouched.
    c = _kcnq2_arg213trp()
    rb.withhold_own_clinvar_tier(c)
    assert c.annotation.clinvar_significance == "Pathogenic"


# --------------------------------- the scorer must describe ONE variant, not two spliced together


class _FakeReport:
    def __init__(self, classifications):
        self.classifications = classifications


def _c(gene, tier, zyg="het", hpo=1.0, sig=None):
    from vcf2report.models import Classification

    return Classification(
        variant=Variant(
            chrom="1",
            pos=1,
            ref="A",
            alt="T",
            gene=gene,
            zygosity=zyg,
            consequence="missense_variant",
        ),
        annotation=Annotation(
            hpo_match_score=hpo, hpo_best_match=hpo, clinvar_significance=sig
        ),
        criteria=[],
        tier=tier,
        rule_path="",
    )


def test_the_reported_tier_belongs_to_the_variant_in_the_reported_bucket():
    """SYN-001 was scored `primary·Pathogenic`. What it actually has is a Pathogenic variant in
    the CARRIER bucket and a separate VUS in `primary` — the scorer took the bucket from one
    variant and the tier from another, describing a pair that exists in no single row.

    30 of the 200 cohort cases classify more than one variant in the planted gene, so this is
    not a corner case; it is 15% of the per-case TSV that `--compare` diffs against.
    """
    vus = _c("TM2D3", "Uncertain Significance (VUS)")
    path = _c("TM2D3", "Pathogenic")
    bucket, tier, hit = rb._bucket_of("TM2D3", _FakeReport([path, vus]))
    assert hit is not None
    assert (bucket, tier) in {
        ("primary", vus.tier),
        ("carrier", path.tier),
        ("other", path.tier),
        ("probable_vus", vus.tier),
        ("secondary", path.tier),
    }, (bucket, tier)
    assert tier == hit.tier, (
        "the reported tier is not the tier of the reported classification"
    )


def test_within_one_bucket_the_most_severe_variant_is_scored():
    """Two classifications in the planted gene, both in the same bucket (SYN-026 BORCS8 has
    Pathogenic and Likely Pathogenic together). The report leads with the severe one."""
    lp, p = _c("BORCS8", "Likely Pathogenic"), _c("BORCS8", "Pathogenic")
    _, tier, hit = rb._bucket_of("BORCS8", _FakeReport([lp, p]))
    assert tier == "Pathogenic" and hit is p


def test_a_gene_with_no_classification_is_absent_with_no_tier():
    bucket, tier, hit = rb._bucket_of("NOPE", _FakeReport([_c("OTHER", "Pathogenic")]))
    assert (bucket, tier, hit) == ("absent", None, None)


def test_a_single_classification_is_unchanged():
    """The common case (170 of 200) must behave exactly as before."""
    only = _c("SCN1A", "Pathogenic")
    bucket, tier, hit = rb._bucket_of("SCN1A", _FakeReport([only]))
    assert (
        hit is only
        and tier == "Pathogenic"
        and bucket in {"primary", "secondary", "other"}
    )


# ------------------------------------- the planted VARIANT is identified by coordinate, not gene


def _write_manifest(tmp_path, rows):
    m = tmp_path / "manifest"
    m.mkdir(parents=True, exist_ok=True)
    with open(m / "planted_variants.tsv", "w", newline="") as fh:
        fh.write("syn_id\tchrom\tpos\tref\talt\tgene\tzygosity\tallele\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return tmp_path


def test_the_primary_planted_allele_is_the_one_scored(tmp_path):
    """40 of the 200 cases plant TWO alleles (compound het — the manifest's `allele` column is
    primary/second). "The classification for this gene" then has no single answer; both are
    real plants. The coordinate does."""
    bench = _write_manifest(
        tmp_path,
        [
            ("SYN-001", "chr15", "101645162", "C", "T", "TM2D3", "het", "second"),
            ("SYN-001", "chr15", "101652355", "C", "A", "TM2D3", "het", "primary"),
        ],
    )
    assert rb._load_planted_loci(bench) == {"SYN-001": ("chr15", 101652355, "C", "A")}


def test_a_missing_manifest_degrades_instead_of_crashing(tmp_path):
    assert rb._load_planted_loci(tmp_path) == {}


def test_the_locus_match_ignores_the_chr_prefix_and_allele_case():
    """The VCF, the manifest and the stores do not agree on `chr` prefixing, and VCF alleles are
    case-insensitive. A miss here would silently fall back to the gene-level pick."""
    c = _c("TM2D3", "Pathogenic")
    c.variant.chrom, c.variant.pos, c.variant.ref, c.variant.alt = (
        "15",
        101652355,
        "c",
        "a",
    )
    rep = _FakeReport([c])
    assert rb._at_locus(rep, ("chr15", 101652355, "C", "A")) is c
    assert rb._at_locus(rep, ("chr15", 101652356, "C", "A")) is None
    assert rb._at_locus(rep, None) is None


def test_the_planted_locus_wins_over_the_gene_level_pick():
    """SYN-001's shape: a Pathogenic variant and a separate VUS in the same gene. Whichever
    bucket logic says, `tier` must describe the allele the benchmark actually planted."""
    vus = _c("TM2D3", "Uncertain Significance (VUS)")
    vus.variant.pos = 101652355
    path = _c("TM2D3", "Pathogenic")
    path.variant.pos = 999
    rep = _FakeReport([path, vus])
    got = rb._at_locus(rep, ("1", 101652355, "A", "T"))
    assert got is vus, "the coordinate match picked the wrong variant"


def test_gene_matching_is_case_insensitive():
    """HGNC symbols are mixed-case ("C10orf71") and the benchmark's answer key upper-cases them.
    An exact string comparison scored SYN-182 `absent` — a frameshift the engine had in fact
    recovered into the PRIMARY bucket. A full point off the published headline, lost to letter
    case rather than to anything the engine did or failed to do.
    """
    hit = _c("C10orf71", "Uncertain Significance (VUS)")
    bucket, tier, got = rb._bucket_of("C10ORF71", _FakeReport([hit]))
    assert got is hit and bucket != "absent", "the gene was missed on letter case alone"
    assert tier == "Uncertain Significance (VUS)"


def test_a_genuinely_different_gene_is_still_absent():
    """The relaxation must be case only. SYN-042 plants in RNU4-2 and the annotator calls the
    locus RNU4-1 — an adjacent paralogue, a real annotation limitation, and it must keep
    counting as a miss."""
    bucket, tier, got = rb._bucket_of(
        "RNU4-2", _FakeReport([_c("RNU4-1", "Pathogenic")])
    )
    assert (bucket, tier, got) == ("absent", None, None)
