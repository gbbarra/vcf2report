"""gnomAD remote-tabix client: record matching, grpmax, and fallback.

The htslib/GCS layer is mocked with fake VCF records — no network. These lock the
grpmax-vs-global logic, ref/alt matching, exome-vs-genome max, and the
covered-but-absent (AF 0) vs not-runnable (None) distinction.
"""
import pytest

from vcf2report.annotate import gnomad, gnomad_remote
from vcf2report.models import Variant


class _Rec:
    def __init__(self, pos, ref, alts, info):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self._info = info

    class _Info:
        def __init__(self, d):
            self._d = d

        def get(self, k):
            return self._d.get(k)

    @property
    def info(self):
        return _Rec._Info(self._info)


class _VF:
    def __init__(self, recs):
        self._recs = recs

    def fetch(self, chrom, start, end):
        return [r for r in self._recs if start < r.pos <= end]


def _reset():
    gnomad_remote._handles.clear()
    gnomad_remote._failed.clear()


def _patch(monkeypatch, exomes=None, genomes=None):
    _reset()
    monkeypatch.setattr(gnomad_remote, "_get_pysam", lambda: object())  # non-None

    def fake_open(kind, chrom):
        recs = exomes if kind == "exomes" else genomes
        return _VF(recs) if recs is not None else None

    monkeypatch.setattr(gnomad_remote, "_open", fake_open)


def test_grpmax_used(monkeypatch):
    rec = _Rec(100, "G", ("T",), {
        "grpmax": ("nfe",), "AF_grpmax": (1.1e-05,), "AC_grpmax": (1,),
        "AN_grpmax": (90000,), "nhomalt_grpmax": (0,), "AF": (5e-06,)})
    _patch(monkeypatch, exomes=[rec], genomes=[])
    r = gnomad_remote.query(Variant(chrom="chr21", pos=100, ref="G", alt="T"))
    assert r["pop"] == "nfe" and abs(r["af"] - 1.1e-05) < 1e-12 and r["an"] == 90000


def test_excluded_grpmax_falls_back_to_global(monkeypatch):
    rec = _Rec(100, "G", ("T",), {
        "grpmax": ("fin",), "AF_grpmax": (0.02,),   # fin excluded from popmax
        "AF": (0.001,), "AC": (30,), "AN": (30000,), "nhomalt": (0,)})
    _patch(monkeypatch, exomes=[rec], genomes=[])
    r = gnomad_remote.query(Variant(chrom="21", pos=100, ref="G", alt="T"))
    assert r["pop"] is None and abs(r["af"] - 0.001) < 1e-12


def test_ref_alt_must_match(monkeypatch):
    rec = _Rec(100, "G", ("A",), {"AF": (0.5,)})   # different ALT
    _patch(monkeypatch, exomes=[rec], genomes=[])
    r = gnomad_remote.query(Variant(chrom="21", pos=100, ref="G", alt="T"))
    assert r["af"] == 0.0            # covered position, allele absent


def test_takes_higher_of_exome_genome(monkeypatch):
    ex = _Rec(100, "G", ("T",), {"AF": (0.001,), "AC": (10,), "AN": (10000,), "nhomalt": (0,)})
    ge = _Rec(100, "G", ("T",), {"AF": (0.004,), "AC": (8,), "AN": (2000,), "nhomalt": (0,)})
    _patch(monkeypatch, exomes=[ex], genomes=[ge])
    r = gnomad_remote.query(Variant(chrom="21", pos=100, ref="G", alt="T"))
    assert abs(r["af"] - 0.004) < 1e-12    # genome AF is higher


def test_not_runnable_returns_none(monkeypatch):
    _reset()
    monkeypatch.setattr(gnomad_remote, "_get_pysam", lambda: object())
    monkeypatch.setattr(gnomad_remote, "_open", lambda kind, chrom: None)  # nothing opens
    assert gnomad_remote.query(Variant(chrom="21", pos=100, ref="G", alt="T")) is None


def test_no_pysam_returns_none(monkeypatch):
    _reset()
    monkeypatch.setattr(gnomad_remote, "_get_pysam", lambda: None)
    assert gnomad_remote.query(Variant(chrom="21", pos=1, ref="A", alt="G")) is None


def test_lookup_prefers_remote_tabix(monkeypatch):
    """gnomad.lookup should use the remote-tabix result and label its source."""
    monkeypatch.setenv("OFFLINE", "")
    monkeypatch.setenv("VCF2REPORT_ALLOW_NETWORK", "1")
    rec = _Rec(166003360, "C", ("T",), {
        "grpmax": ("amr",), "AF_grpmax": (3e-04,), "AC_grpmax": (5,),
        "AN_grpmax": (16000,), "nhomalt_grpmax": (0,)})
    _patch(monkeypatch, exomes=[rec], genomes=[])
    r = gnomad.lookup(Variant(chrom="2", pos=166003360, ref="C", alt="T"))
    assert r["pop"] == "amr" and "remote tabix" in r["_source"]


def test_lookup_offline_skips_remote(monkeypatch):
    monkeypatch.setenv("OFFLINE", "1")

    def _boom(*a, **k):
        raise AssertionError("remote tabix must not run in OFFLINE mode")

    monkeypatch.setattr(gnomad_remote, "query", _boom)
    # TTN in local snapshot -> served offline without touching remote
    r = gnomad.lookup(Variant(chrom="2", pos=178562809, ref="G", alt="A"))
    assert r["af"] == 0.081
