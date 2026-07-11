"""Regression tests for the stability/adversarial-audit fixes."""
from vcf2report.acmg.engine import classify
from vcf2report.acmg.rules import VUS
from vcf2report.annotate import abraom, gnomad_remote
from vcf2report.models import Annotation, Variant


def _rare(**kw):
    base = dict(gene_lof_intolerant=True, gnomad_af=0.0, gnomad_faf95=0.0,
                abraom_af=0.0, source={})
    base.update(kw)
    return Annotation(**base)


# H1 — lowercase alleles must not defeat snapshot-key matching.
def test_key_uppercases_alleles():
    assert Variant(chrom="chr1", pos=100, ref="a", alt="t").key == "1-100-A-T"
    assert (Variant(chrom="1", pos=100, ref="A", alt="T").key
            == Variant(chrom="1", pos=100, ref="a", alt="t").key)


# H3 — stop_lost is PM4-only, never PVS1 (would double-count and over-call).
def test_stop_lost_not_lof_no_pvs1():
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="G", consequence="stop_lost")
    assert v.is_lof is False
    c = classify(v, _rare())
    assert "PVS1" not in c.met_codes and "PM4" in c.met_codes
    assert c.tier == VUS               # PM2 + PM4 (2 moderate) is insufficient


def test_stop_gained_still_pvs1():
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="G", consequence="stop_gained")
    assert v.is_lof is True and "PVS1" in classify(v, _rare()).met_codes


# H2 — a single queryable callset must NOT assert absence (af 0.0).
def test_remote_single_callset_no_absence(monkeypatch):
    class _Empty:
        def fetch(self, *a, **k):
            return []

    monkeypatch.setattr(gnomad_remote, "_get_pysam", lambda: object())
    monkeypatch.setattr(gnomad_remote, "_open",
                        lambda kind, chrom: _Empty() if kind == "exomes" else None)
    # only exomes opened, no allele -> opened==1 -> unknown (None), not fabricated 0.0.
    assert gnomad_remote.query(Variant(chrom="1", pos=1, ref="A", alt="T")) is None


def test_remote_both_callsets_absent_is_zero(monkeypatch):
    class _Empty:
        def fetch(self, *a, **k):
            return []

    monkeypatch.setattr(gnomad_remote, "_get_pysam", lambda: object())
    monkeypatch.setattr(gnomad_remote, "_open", lambda kind, chrom: _Empty())
    r = gnomad_remote.query(Variant(chrom="1", pos=1, ref="A", alt="T"))
    assert r is not None and r["af"] == 0.0   # both opened, absent -> genuine 0.0


# L1 — ABraOM miss is UNKNOWN (None), not a fabricated checked 0.0.
def test_abraom_miss_is_unknown():
    r = abraom.lookup(Variant(chrom="99", pos=1, ref="A", alt="T"))
    assert r["af"] is None
