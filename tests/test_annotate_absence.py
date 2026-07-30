"""Six ways the annotation layer turned "we do not know" into "we looked, and it is not there".

All six were reported by adversarial hunts over the annotate package and left unfixed. They are
one defect wearing six costumes: a lookup path that cannot distinguish *no answer* from *the
answer is no*, and then reports the second. That direction over-calls — PM2 ("absent from
controls") is the criterion most of them feed, and it fires on absence.

Each test here is written from the dangerous direction: what would have to be true for a
variant gnomAD/ClinVar/AlphaMissense actually knows about to be reported as unobserved.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from vcf2report import config
from vcf2report.annotate import alphamissense, cache, clinvar, gnomad
from vcf2report.models import Variant
from vcf2report.vcf.parse import detect_build

REPO = config.REPO_ROOT


def _v(chrom="chr1", pos=100, ref="A", alt="G", **kw):
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, **kw)


# ----------------------------------------------------- 1 · ClinVar consulted the cache first

def test_clinvar_prefers_the_local_store_over_a_stale_cache(monkeypatch):
    """The store is the versioned, refreshed source; data/cache/ is flat unversioned JSON.

    gnomad.lookup already documents this order. ClinVar had it inverted — and ClinVar is the
    one store with a weekly staleness window, so a cache entry could serve a months-old
    assertion under a _source that says nothing about its age.
    """
    v = _v()
    monkeypatch.setattr(cache, "get", lambda src, key: {
        "significance": "Benign", "review_status": "no assertion criteria provided"})
    from vcf2report.annotate import clinvar_parquet
    monkeypatch.setattr(clinvar_parquet, "get", lambda key: {
        "significance": "Pathogenic", "review_status": "reviewed by expert panel"})

    got = clinvar.lookup(v)
    assert got["significance"] == "Pathogenic", "the cache shadowed a fresh store answer"
    assert got["_source"] == "ClinVar (local)"


def test_clinvar_cache_is_still_used_when_no_store_answers(monkeypatch):
    """Reordering must not throw the cache away — it is the offline fallback that earns its
    keep when neither store has the variant."""
    v = _v()
    monkeypatch.setattr(cache, "get", lambda src, key: {
        "significance": "Likely pathogenic", "review_status": "criteria provided"})
    from vcf2report.annotate import clinvar_parquet
    monkeypatch.setattr(clinvar_parquet, "get", lambda key: None)
    monkeypatch.setattr(clinvar, "_tabix_lookup", lambda variant: None)
    monkeypatch.setattr(config, "offline", lambda: True)
    monkeypatch.setattr(clinvar, "_load_local", lambda: {})

    got = clinvar.lookup(v)
    assert got["significance"] == "Likely pathogenic"
    assert got["_source"] == "ClinVar (cache)"


def test_a_legacy_poisoned_cache_entry_cannot_outrank_a_store(monkeypatch):
    """Caches written by the old warm_cache still exist on operators' disks, carrying
    {significance: None}. They must not read as "ClinVar says nothing"."""
    v = _v()
    monkeypatch.setattr(cache, "get", lambda src, key: {
        "significance": None, "review_status": None, "accession": None})
    from vcf2report.annotate import clinvar_parquet
    monkeypatch.setattr(clinvar_parquet, "get", lambda key: None)
    monkeypatch.setattr(clinvar, "_tabix_lookup", lambda variant: None)
    monkeypatch.setattr(config, "offline", lambda: True)
    monkeypatch.setattr(clinvar, "_load_local", lambda: {})

    got = clinvar.lookup(v)
    # It falls through to the honest "no record", NOT a cache hit dressed as an answer.
    assert got["_source"] == "ClinVar (no record)"


# --------------------------------------------------- 2 · warm_cache persisted the sentinels

def test_warm_cache_does_not_persist_a_negative():
    """A cached miss is an absence of data stored where data is read from.

    Both lookups return a dict whether or not they found anything, so the old script wrote
    {significance: None} to disk for every unmatched variant — permanently, through every
    later run, including ones where a correctly-built store DOES have the record.
    """
    src = (REPO / "scripts/warm_cache.py").read_text()
    # The guard must be on the value, before the put.
    assert 'if cv.get("significance"):' in src
    assert 'if g.get("af") is not None:' in src
    # And the two cache.put calls must be the guarded ones, not bare.
    for line in src.splitlines():
        if "cache.put(" in line:
            assert line.startswith(" " * 12), (
                f"cache.put at top-level indent — outside the guard: {line!r}")


def test_warm_cache_still_runs_and_reports_what_it_skipped(tmp_path):
    """It must stay usable: warm the bundled sample, and say out loud how many were left
    uncached rather than silently writing fewer entries."""
    env = {"VCF2REPORT_CACHE": str(tmp_path / "cache"), "OFFLINE": "1",
           "PATH": "/usr/bin:/bin", "HOME": str(tmp_path)}
    p = subprocess.run([sys.executable, "scripts/warm_cache.py"], cwd=REPO,
                       capture_output=True, text=True, timeout=600, env=env)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "Warmed cache" in p.stdout


# ------------------------------------------ 3 · an AlphaMissense store miss silenced tabix

def test_alphamissense_store_miss_falls_through_to_tabix(monkeypatch):
    """The store answers by LEFT JOIN, so it returns an entry for EVERY key — None for the
    rows it lacks. Writing those Nones into _primed and returning meant a store that merely
    lacks a contig permanently silenced the tabix file that carries the score.

    Unlike gnomad_parquet, this store has no vouched-absence machinery: it cannot tell "no
    score exists for this substitution" from "this store does not cover it".
    """
    a, b = _v(pos=1), _v(pos=2)
    alphamissense.clear_primed()
    from vcf2report.annotate import alphamissense_parquet
    monkeypatch.setattr(alphamissense_parquet, "available", lambda: True)
    # Store knows `a`, has no row for `b`.
    monkeypatch.setattr(alphamissense_parquet, "prime",
                        lambda variants: {a.key: {"am_pathogenicity": 0.9}, b.key: None})
    # Tabix DOES have `b`.
    monkeypatch.setattr(alphamissense, "_open", lambda: object())
    monkeypatch.setattr(alphamissense, "_fetch", lambda t, c, p: ["row"])
    monkeypatch.setattr(alphamissense, "_best",
                        lambda rows, v: {"am_pathogenicity": 0.77})

    alphamissense.prime([a, b])
    assert alphamissense._primed[a.key]["am_pathogenicity"] == 0.9
    assert alphamissense._primed[b.key] is not None, "store miss silenced the tabix fallback"
    assert alphamissense._primed[b.key]["am_pathogenicity"] == 0.77
    alphamissense.clear_primed()


def test_alphamissense_skips_tabix_when_the_store_answered_everything(monkeypatch):
    """The fall-through must not cost a tabix sweep on the happy path."""
    a = _v(pos=1)
    alphamissense.clear_primed()
    from vcf2report.annotate import alphamissense_parquet
    monkeypatch.setattr(alphamissense_parquet, "available", lambda: True)
    monkeypatch.setattr(alphamissense_parquet, "prime",
                        lambda variants: {a.key: {"am_pathogenicity": 0.5}})

    def _boom():
        raise AssertionError("tabix opened even though the store answered every key")

    monkeypatch.setattr(alphamissense, "_open", _boom)
    alphamissense.prime([a])
    assert alphamissense._primed[a.key]["am_pathogenicity"] == 0.5
    alphamissense.clear_primed()


# ------------------------------------------------- 4 · case-sensitive remote allele compare

@pytest.mark.parametrize("ref,alt", [("a", "g"), ("A", "g"), ("a", "G")])
def test_remote_matches_lowercase_alleles(ref, alt):
    """Soft-masked references make lowercase calls ordinary, the VCF spec treats bases as
    case-insensitive, and vcf.parse deliberately does not rewrite the operator's alleles.

    A case-sensitive compare turned "a" vs "A" into a gnomAD MISS — so the variant looked
    unobserved and PM2 fired on an allele gnomAD knows.
    """
    src = (REPO / "src/vcf2report/annotate/gnomad_remote.py").read_text()
    assert "variant.ref.upper()" in src and "variant.alt.upper()" in src
    assert 'rec.ref or ""' in src, "record REF must be upper-cased defensively too"
    # And the equality must not be against the raw attribute any more.
    assert "rec.ref != variant.ref" not in src
    assert "variant.alt not in alts" not in src or "v_alt not in alts" in src


# ----------------------------------- 5 · _from_payload invented a confident zero-frequency

def _payload(pops, hom=0):
    return {"exome": {"homozygote_count": hom,
                      "populations": [{"id": k, "ac": ac, "an": an}
                                      for k, (ac, an) in pops.items()]}}


def test_founder_only_allele_is_unknown_not_absent():
    """The Finnish case. A founder allele can be percent-level in FIN and invisible to
    popmax — fin/asj/ami/mid are excluded on purpose (ClinGen SVI excludes them too, a
    founder allele inflates popmax). But excluding it from the MAX is not the same as
    reporting the allele absent, and af=0.0 satisfied PM2 on a variant gnomAD has seen.
    """
    got = gnomad._from_payload(_payload({"fin": (250, 12000)}))
    assert got["af"] is None, "a founder-only allele was reported as frequency 0.0"
    assert "founder" in got["_freq_unknown"] or "excluded" in got["_freq_unknown"]


def test_underpowered_site_is_unknown_not_absent():
    """Below _MIN_AN the frequency estimate is noise, not a small number."""
    got = gnomad._from_payload(_payload({"nfe": (0, 10), "afr": (0, 4)}))
    assert got["af"] is None
    assert str(gnomad._MIN_AN) in got["_freq_unknown"]


def test_a_well_powered_zero_is_still_reported_as_zero():
    """The fix must not make every absence unknown — a deeply-surveyed zero is real
    evidence and PM2 is entitled to it."""
    got = gnomad._from_payload(_payload({"nfe": (0, 60000), "afr": (0, 20000)}))
    assert got["af"] == 0.0
    # It names the population that surveyed it most deeply, so the trail shows the basis.
    assert got["pop"] == "nfe"
    assert got["an"] == 60000


def test_a_real_frequency_still_wins_popmax():
    got = gnomad._from_payload(_payload({"nfe": (5, 50000), "afr": (300, 20000)}))
    assert got["pop"] == "afr"
    assert got["af"] == pytest.approx(300 / 20000)


def test_homozygote_count_survives_an_unknown_frequency():
    """hom IS surveyed even when popmax is not, and BS2 may still use it."""
    got = gnomad._from_payload(_payload({"fin": (250, 12000)}, hom=7))
    assert got["af"] is None
    assert got["hom"] == 7


# --------------------------------------------------- 6 · build guard missed real GRCh37 names

@pytest.mark.parametrize("token", [
    "hs37d5",                                   # 1000 Genomes phase-2 decoy reference
    "Homo_sapiens_assembly19.fasta",            # the Broad/GATK bundle name
    "HS37D5",                                   # case must not matter
])
def test_common_grch37_reference_names_are_detected(token):
    """A header naming only these fell through to None — "build not declared" instead of
    "this is GRCh37". A GRCh37 callset that slips the guard is classified against GRCh38
    coordinates, and every store lookup silently mismatches."""
    assert detect_build([f"##reference=file:///refs/{token}"]) == "GRCh37"


def test_grch38_still_wins_and_unknown_stays_unknown():
    assert detect_build(["##reference=GRCh38_full_analysis_set.fa"]) == "GRCh38"
    assert detect_build(["##fileformat=VCFv4.2"]) is None
    # The point of returning None is that a loose guess is worse than admitting ignorance.
    assert detect_build(["##source=someCaller_v38.1"]) is None
