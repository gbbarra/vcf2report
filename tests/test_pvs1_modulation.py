"""Deterministic PVS1 strength modulation (ClinGen SVI / Abou Tayoun 2018).

The NMD/exon decision tree only fires when the VCF is annotated with an exon rank
(VEP EXON / SnpEff rank). Un-annotated variants keep PVS1 at Very Strong, so the
synthetic demos and the frozen concordance panel are unaffected.
"""
import pytest

from vcf2report.acmg import criteria
from vcf2report.acmg.criteria import _is_last_exon, _pvs1_strength
from vcf2report.acmg.engine import classify
from vcf2report.acmg.rules import LIKELY_PATHOGENIC, VUS
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
    assert _pvs1_strength(_v("start_lost")) == "moderate"
    assert _pvs1_strength(_v("start_lost", "5/10")) == "moderate"          # exon irrelevant
    assert _pvs1_strength(_v("stop_gained", "10/10")) == "strong"          # NMD-escaping
    assert _pvs1_strength(_v("frameshift_variant", "8/8")) == "strong"
    assert _pvs1_strength(_v("stop_gained", "5/10")) == "very_strong"      # NMD-triggering
    assert _pvs1_strength(_v("stop_gained", None)) == "very_strong"        # unannotated
    assert _pvs1_strength(_v("splice_donor_variant", "1/1")) == "very_strong"  # not in downgrade set


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

    # Not LoF-intolerant -> PVS1 never fires, strength stays None.
    off = criteria.pvs1(_v("stop_gained", "5/10"), _ann(gene_lof_intolerant=False))
    assert not off.met and off.applied_strength is None


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
