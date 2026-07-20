"""Deterministic PVS1 strength modulation (ClinGen SVI / Abou Tayoun 2018).

The NMD/exon decision tree only fires when the VCF is annotated with an exon rank
(VEP EXON / SnpEff rank). Un-annotated variants keep PVS1 at Very Strong, so the
synthetic demos and the frozen concordance panel are unaffected.
"""
import pytest

from vcf2report import config
from vcf2report.acmg import criteria
from vcf2report.acmg.criteria import _is_last_exon, _pvs1_strength
from vcf2report.acmg.engine import classify
from vcf2report.acmg.rules import LIKELY_PATHOGENIC, VUS
from vcf2report.annotate import inheritance
from vcf2report.models import Annotation, Variant


# ---------------------------------------------------------------------------
# Exon-rank parsing
# ---------------------------------------------------------------------------
def test_is_last_exon():
    assert _is_last_exon("12/12") is True
    assert _is_last_exon("1/1") is True        # single-exon transcript escapes NMD
    assert _is_last_exon("11/12") is False
    assert _is_last_exon("5") is False          # no denominator
    assert _is_last_exon("0/0") is False
    assert _is_last_exon(None) is False
    assert _is_last_exon("") is False
    assert _is_last_exon("x/y") is False


# ---------------------------------------------------------------------------
# Decision tree
# ---------------------------------------------------------------------------
def _v(consequence, exon=None):
    return Variant(chrom="1", pos=1, ref="A", alt="T", gene="G",
                   consequence=consequence, exon=exon)


def test_pvs1_strength_tree():
    assert _pvs1_strength(_v("start_lost", "1/12")) == "moderate"          # annotated start-loss
    assert _pvs1_strength(_v("stop_gained", "10/10")) == "strong"          # NMD-escaping
    assert _pvs1_strength(_v("frameshift_variant", "8/8")) == "strong"
    assert _pvs1_strength(_v("stop_gained", "5/10")) == "very_strong"      # NMD-triggering
    assert _pvs1_strength(_v("splice_donor_variant", "1/1")) == "very_strong"  # not in downgrade set


def test_pvs1_tree_gated_on_exon():
    # THE invariant: no exon rank -> the tree never engages, PVS1 stays Very Strong.
    # (start_lost must be gated too, not just the NMD-escape branch.)
    for cons in ("start_lost", "stop_gained", "frameshift_variant", "splice_donor_variant"):
        assert _pvs1_strength(_v(cons, None)) == "very_strong"
        assert _pvs1_strength(_v(cons, "")) == "very_strong"


# ---------------------------------------------------------------------------
# Criterion wiring (applied_strength + met)
# ---------------------------------------------------------------------------
def _ann(**kw):
    base = dict(gene_lof_intolerant=True, abraom_af=0.0, gnomad_af=0.0,
                gnomad_faf95=0.0, source={})
    base.update(kw)
    return Annotation(**base)


def test_pvs1_criterion_applied_strength():
    last = criteria.pvs1(_v("stop_gained", "10/10"), _ann())
    assert last.met and last.applied_strength == "strong"

    mid = criteria.pvs1(_v("stop_gained", "5/10"), _ann())
    assert mid.met and mid.applied_strength == "very_strong"

    start = criteria.pvs1(_v("start_lost", "1/12"), _ann())
    assert start.met and start.applied_strength == "moderate"

    # Not LoF-intolerant AND no recessive phenotype ("G" is not a real gene) -> never fires.
    off = criteria.pvs1(_v("stop_gained", "5/10"), _ann(gene_lof_intolerant=False))
    assert not off.met and off.applied_strength is None


# ---------------------------------------------------------------------------
# The LoF-mechanism gate: constraint OR an established recessive phenotype
# ---------------------------------------------------------------------------
@pytest.fixture
def moi(monkeypatch):
    """Stub the HPO-derived inheritance table (criteria imports the function by value,
    so patch the cache the function reads, not the name)."""
    def _set(mapping):
        monkeypatch.setattr(inheritance, "_gene_moi",
                            {g.upper(): frozenset(v) for g, v in mapping.items()})
    return _set


def test_pvs1_fires_for_recessive_gene_without_constraint(moi):
    """THE recessive fix: pLI/LOEUF measure selection against HETEROZYGOUS LoF, so a
    recessive disease gene scores as tolerant (carriers are healthy) and the old
    constraint-only gate rejected the very variants whose mechanism IS loss of function."""
    moi({"RECGENE": ["AR"]})
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="RECGENE",
                consequence="stop_gained", exon="5/10")
    r = criteria.pvs1(v, _ann(gene_lof_intolerant=False))
    assert r.met and r.applied_strength == "very_strong"
    assert "autosomal-recessive" in r.evidence["lof_mechanism_basis"]


def test_pvs1_does_not_fire_for_unconstrained_dominant_gene(moi):
    """A dominant phenotype can be gain-of-function or dominant-negative, where a null is
    NOT the mechanism — so dominant genes still have to earn PVS1 via constraint."""
    moi({"DOMGENE": ["AD"]})
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="DOMGENE",
                consequence="stop_gained", exon="5/10")
    r = criteria.pvs1(v, _ann(gene_lof_intolerant=False))
    assert not r.met and r.evidence["lof_mechanism_basis"] is None


def test_pvs1_constraint_still_wins_when_inheritance_unknown(moi):
    moi({})
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="ANY",
                consequence="stop_gained", exon="5/10")
    assert criteria.pvs1(v, _ann(gene_lof_intolerant=True)).met
    assert not criteria.pvs1(v, _ann(gene_lof_intolerant=False)).met


def test_pvs1_basis_names_which_route_fired(moi):
    """The laudo must show WHY PVS1 fired — constraint and recessive-phenotype are
    different strengths of evidence and a reviewer has to be able to tell them apart."""
    moi({"BOTH": ["AR"]})
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="BOTH",
                consequence="stop_gained", exon="5/10")
    assert "constraint" in criteria.pvs1(v, _ann(gene_lof_intolerant=True)).evidence["lof_mechanism_basis"]
    assert "recessive" in criteria.pvs1(v, _ann(gene_lof_intolerant=False)).evidence["lof_mechanism_basis"]


def test_multimode_gene_takes_a_different_pole_per_criterion(moi):
    """THE bug this guards: "conservative" has a DIRECTION, and the two disagree.

    PM2 is PATHOGENIC evidence -> its conservative pole is the STRICTEST (lowest) ceiling.
    BS1 is BENIGN evidence -> its conservative pole is the MOST PERMISSIVE (highest) cutoff.
    Collapsing both onto one "stricter" value inverts BS1: an ATM variant at an ordinary
    ataxia-telangiectasia CARRIER frequency clears the dominant cutoff and is called Benign,
    overriding real pathogenic evidence.
    """
    moi({"BOTHGENE": ["AD", "AR"]})
    assert config.pm2_af_ceiling("BOTHGENE") == (config.PM2_AF_DOMINANT, "AD")
    assert config.bs1_af_cutoff("BOTHGENE") == (config.BS1_AF_RECESSIVE, "AR")
    # a carrier-frequency allele must NOT clear BS1 on a gene with recessive disease
    cutoff, _ = config.bs1_af_cutoff("BOTHGENE")
    assert 0.005 < cutoff
    # LoF is still an established mechanism there, because the recessive disease is real
    assert inheritance.lof_is_disease_mechanism("BOTHGENE") is True
    # ...and the label stays honest about the ambiguity rather than hiding it
    assert config.gene_inheritance("BOTHGENE") == "AD+AR"


def test_single_mode_gene_thresholds_unchanged(moi):
    moi({"REC": ["AR"], "DOM": ["AD"]})
    assert config.pm2_af_ceiling("REC") == (config.PM2_AF_RECESSIVE, "AR")
    assert config.bs1_af_cutoff("REC") == (config.BS1_AF_RECESSIVE, "AR")
    assert config.pm2_af_ceiling("DOM") == (config.PM2_AF_DOMINANT, "AD")
    assert config.bs1_af_cutoff("DOM") == (config.BS1_AF_DOMINANT, "AD")
    assert config.pm2_af_ceiling("UNKNOWN") == (config.PM2_AF_DEFAULT, None)
    assert config.bs1_af_cutoff("UNKNOWN") == (config.BS1_AF_DEFAULT, None)


def test_curated_inheritance_beats_hpo(moi):
    moi({"CFTR": ["AD"]})            # HPO disagrees...
    assert config.gene_inheritance("CFTR") == "AR"   # ...curated map wins
    assert config.gene_inheritance("NOTAGENE") is None


@pytest.mark.skipif(not config.HPO_GENES_LOCAL.exists(), reason="HPO table not installed")
def test_real_multimode_genes_never_call_a_carrier_frequency_benign():
    """Named regression, on the REAL HPO table — the variants a directional slip destroys.

    HBB chr11:5226925 (sickle-cell, ClinVar Pathogenic 3-star) sits at faf95=0.00433: HbS is
    common where malaria was endemic. HBB is AD+AR in HPO, so collapsing it to the dominant
    BS1 cutoff (0.001) calls sickle-cell disease Benign. An audit found 38 ClinVar P/LP
    >=2-star variants in this exposed band (ABCA4, ANO5, FLG, SPG7, GNE, MPO, GBA1, MEFV,
    PKLR, ALPL, MSH2, CDH23...). Calling a known pathogenic variant benign is the worst
    error this tool can make.
    """
    for gene in ("HBB", "ABCA4", "ANO5", "FLG", "SPG7", "GNE", "MPO", "GBA1", "MEFV",
                 "PKLR", "ALPL", "MSH2", "CDH23"):
        cutoff, moi_used = config.bs1_af_cutoff(gene)
        assert moi_used == "AR", f"{gene}: BS1 must take the recessive pole, got {moi_used}"
        assert cutoff == config.BS1_AF_RECESSIVE
        assert 0.00433 < cutoff, f"{gene}: a carrier-frequency allele would be called benign"


@pytest.mark.skipif(not config.HPO_GENES_LOCAL.exists(), reason="HPO table not installed")
def test_every_multimode_gene_keeps_both_poles():
    """The property, over EVERY gene HPO gives two modes (~590), not just the named ones.

    The named-gene test above only covers the 13 the audit happened to surface; the bug was a
    whole CLASS. The invariant: a gene with any recessive disease never takes a BS1 cutoff
    below the recessive one (else carrier frequencies read as benign), and a gene with any
    dominant disease never takes a PM2 ceiling above the dominant one (else rarity is granted
    too cheaply). The two poles are independent — that independence IS the fix.
    """
    from vcf2report.annotate.inheritance import _load
    multi = [g for g, m in _load().items() if len(m) > 1]
    assert len(multi) > 100, f"expected HPO to give many multi-mode genes, got {len(multi)}"
    for gene in multi:
        m = config.gene_inheritance_modes(gene)
        if "AR" in m:
            assert config.bs1_af_cutoff(gene)[0] == config.BS1_AF_RECESSIVE, gene
        if "AD" in m or "XL" in m:
            assert config.pm2_af_ceiling(gene)[0] == config.PM2_AF_DOMINANT, gene


# ---------------------------------------------------------------------------
# End-to-end tier flip (visible under the ClinGen points model)
# ---------------------------------------------------------------------------
def test_pvs1_downgrade_flips_tier_clingen(monkeypatch):
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    # Mid-exon PTC: PVS1_VeryStrong(8) + PM2_Supporting(1) = 9 -> Likely Pathogenic.
    mid = classify(_v("stop_gained", "5/10"), _ann())
    assert mid.tier == LIKELY_PATHOGENIC
    # Last-exon PTC (NMD escape): PVS1_Strong(4) + PM2_Supporting(1) = 5 -> VUS.
    last = classify(_v("stop_gained", "10/10"), _ann())
    assert last.tier == VUS


def test_pvs1_unannotated_unchanged(monkeypatch):
    # No exon rank -> Very Strong under both models (panel/demo behaviour preserved).
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    c = classify(_v("stop_gained", None), _ann())
    pvs1 = next(x for x in c.criteria if x.code == "PVS1")
    assert pvs1.applied_strength == "very_strong" and c.tier == LIKELY_PATHOGENIC


def test_pvs1_clingen_haploinsufficiency_route(moi):
    """ClinGen Haploinsufficiency=3 is a curated 'LoF causes disease' statement — it opens PVS1 for
    late-onset / incompletely-penetrant dominants that population constraint misses (TP53: LOEUF
    0.469, 'not intolerant', yet a textbook haploinsufficient tumour suppressor)."""
    from vcf2report.annotate import dosage
    moi({})  # no HPO inheritance
    monkey = dosage._hi
    dosage._hi = {"TP53"}
    try:
        v = Variant(chrom="17", pos=1, ref="C", alt="A", gene="TP53",
                    consequence="stop_gained", exon="4/11")
        r = criteria.pvs1(v, _ann(gene_lof_intolerant=False))
        assert r.met and "Haploinsufficiency=3" in r.evidence["lof_mechanism_basis"]
        # a gene that is neither HI, constrained, nor recessive still does not fire
        v2 = Variant(chrom="1", pos=1, ref="A", alt="T", gene="RANDO",
                     consequence="stop_gained", exon="2/5")
        assert not criteria.pvs1(v2, _ann(gene_lof_intolerant=False)).met
    finally:
        dosage._hi = monkey


def test_clingen_hi_store_data_integrity():
    """Lock the committed ClinGen HI store against a silent data regression — the wiring test
    monkeypatches the set, so nothing else reads the real file. A bad fetch/merge that empties it,
    drops TP53, or leaks a score-30 (recessive) / score-40 (dosage-unlikely) gene would otherwise
    stay green. Skips cleanly if the file is not installed."""
    from vcf2report.annotate import dosage
    if not config.CLINGEN_HI_LOCAL.exists():
        pytest.skip("ClinGen HI store not installed")
    dosage._hi = None                       # force a real load from the committed file
    try:
        genes = dosage._load()
        assert len(genes) >= 400, f"HI store shrank to {len(genes)} — suspect a truncated fetch"
        for canary in ("TP53", "BRCA1", "BRCA2", "NF1", "PTEN"):   # classic HI=3
            assert dosage.haploinsufficient(canary), f"{canary} missing — HI=3 filter/fetch broke"
        for leaked in ("A4GALT", "AARS2"):    # ClinGen score 30 (autosomal recessive) — must NOT leak
            assert not dosage.haploinsufficient(leaked), f"{leaked} (score-30) leaked past the HI=3 filter"
    finally:
        dosage._hi = None
