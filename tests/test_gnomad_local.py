"""Local gnomAD reduced-tabix client + its integration into gnomad.lookup."""
import json

import pytest

from vcf2report import config
from vcf2report.annotate import gnomad, gnomad_local
from vcf2report.models import Variant

pysam = pytest.importorskip("pysam")

_ROWS = [
    "#chrom\tpos\tref\talt\taf\tac\tan\thom\tfaf95\tpop",
    "1\t100\tA\tT\t0.30\t300\t1000\t20\t0.28\tnfe",     # common
    "1\t100\tA\tG\t0.001\t1\t1000\t0\t0.0005\tafr",     # same site, other allele
    "1\t200\tC\tG\t0.0\t0\t1000\t0\t0.0\t",             # covered, absent allele (empty pop)
]


def _make_table(tmp_path, mode, contigs=("1",)):
    tsv = tmp_path / "g.tsv"
    tsv.write_text("\n".join(_ROWS) + "\n")
    out = tmp_path / "g.tsv.gz"
    pysam.tabix_compress(str(tsv), str(out), force=True)
    pysam.tabix_index(str(out), seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    meta = {"mode": mode}
    if mode == "full":
        meta["contigs"] = list(contigs)
    (tmp_path / "g.tsv.gz.meta").write_text(json.dumps(meta))
    return out


@pytest.fixture
def local_table(tmp_path, monkeypatch):
    def _setup(mode="partial", contigs=("1",)):
        out = _make_table(tmp_path, mode, contigs)
        monkeypatch.setattr(config, "GNOMAD_LOCAL_TABIX", out)
        gnomad_local._reset_for_tests()
        return out
    yield _setup
    gnomad_local._reset_for_tests()


def _v(pos, ref, alt, chrom="1"):
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


def test_exact_match_returns_grpmax_fields(local_table):
    local_table("partial")
    r = gnomad_local.query(_v(100, "A", "T"))
    assert r["af"] == 0.30 and r["faf95"] == 0.28
    assert r["pop"] == "nfe" and r["hom"] == 20 and r["an"] == 1000


def test_partial_other_allele_does_not_fabricate_absence(local_table):
    # Site 100 carries A>T and A>G, but a partial table may not enumerate every allele
    # (a --from-vcf table holds only the sample's), so a non-exact hit must be None,
    # never a fabricated absence.
    local_table("partial")
    assert gnomad_local.query(_v(100, "A", "C")) is None


def test_full_other_allele_is_absent(local_table):
    # A full table has every gnomAD allele, so a non-exact hit at a covered site is
    # a genuine absence.
    local_table("full")
    assert gnomad_local.query(_v(100, "A", "C"))["af"] == 0.0


def test_partial_miss_returns_none(local_table):
    # A partial table cannot assert absence off its covered sites -> None -> fall back.
    local_table("partial")
    assert gnomad_local.query(_v(999, "A", "T")) is None


def test_full_miss_asserts_absence(local_table):
    # A full table treats a miss as genuinely absent from gnomAD.
    local_table("full")
    assert gnomad_local.query(_v(999, "A", "T"))["af"] == 0.0


def test_chrom_prefix_tolerated(local_table):
    local_table("partial")
    assert gnomad_local.query(_v(100, "A", "T", chrom="chr1"))["af"] == 0.30


def test_case_insensitive_alleles(local_table):
    # Lowercase input alleles (soft-masked callers) must still match — a byte-exact
    # compare would fabricate an absence in full mode.
    local_table("full")
    assert gnomad_local.query(_v(100, "a", "t"))["af"] == 0.30


def test_full_uncovered_contig_returns_none(local_table):
    # A full table only vouches for absence on contigs it covered. chr2 / chrM are not
    # in the sidecar's `contigs`, so a miss there is None (fallback), not a fake 0.0.
    local_table("full", contigs=("1",))
    assert gnomad_local.query(_v(500, "A", "T", chrom="2")) is None
    assert gnomad_local.query(_v(1, "A", "G", chrom="MT")) is None
    # ...but a covered-contig miss is still a genuine absence.
    assert gnomad_local.query(_v(999, "A", "T", chrom="1"))["af"] == 0.0


def test_no_table_returns_none(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "GNOMAD_LOCAL_TABIX", tmp_path / "absent.tsv.gz")
    gnomad_local._reset_for_tests()
    assert gnomad_local.query(_v(100, "A", "T")) is None
    gnomad_local._reset_for_tests()


def test_lookup_prefers_local(local_table, monkeypatch):
    from vcf2report.annotate import cache
    local_table("partial")
    monkeypatch.setattr(cache, "get", lambda *a, **k: None)
    monkeypatch.setattr(cache, "put", lambda *a, **k: None)
    r = gnomad.lookup(_v(100, "A", "T"))
    assert r["af"] == 0.30 and "local tabix" in r["_source"]


def test_lookup_partial_miss_falls_through_offline(local_table, monkeypatch):
    # Local partial miss -> None -> offline with no cache/bundled -> 'unavailable' (af None).
    from vcf2report.annotate import cache
    local_table("partial")
    monkeypatch.setattr(cache, "get", lambda *a, **k: None)
    monkeypatch.setattr(config, "offline", lambda: True)
    r = gnomad.lookup(_v(999, "A", "T"))
    assert r["af"] is None
