"""AlphaMissense batch prime: an in-memory fast path that avoids the on-disk cache's
per-key whole-file rewrite (O(n^2) across a candidate list), mirroring gnomad_parquet.prime.

The scores are identical to the per-variant lookup (same _fetch/_best), so no
classification changes — this only removes the disk churn that dominated the phase.
"""
from vcf2report import config, pipeline
from vcf2report.annotate import alphamissense, alphamissense_parquet, cache
from vcf2report.models import Variant


def _v(pos, ref="A", alt="T"):
    return Variant(chrom="1", pos=pos, ref=ref, alt=alt, gene="G",
                   consequence="missense_variant")


def _fake_store(monkeypatch, scores):
    """Resolve prime/lookup from an in-test store (pos -> score or None) via _best,
    with NO real tabix file."""
    # Isolate the tabix path: prime() now prefers the Parquet store when one is present.
    monkeypatch.setattr(alphamissense_parquet, "available", lambda: False)
    monkeypatch.setattr(alphamissense, "_open", lambda: object())      # non-None handle
    monkeypatch.setattr(alphamissense, "_fetch", lambda t, c, p: ["row"])

    def fake_best(rows, v):
        s = scores.get(v.pos)
        return None if s is None else {"am_pathogenicity": s, "am_class": "x"}
    monkeypatch.setattr(alphamissense, "_best", fake_best)


def _boom(*a, **k):
    raise AssertionError("cache must NOT be touched for a primed variant")


def test_prime_populates_and_lookup_uses_it(monkeypatch):
    alphamissense._reset_for_tests()
    _fake_store(monkeypatch, {100: 0.9, 200: None})
    monkeypatch.setattr(cache, "get", _boom)
    monkeypatch.setattr(cache, "put", _boom)
    assert alphamissense.prime([_v(100), _v(200)]) == 2
    hit = alphamissense.lookup(_v(100))
    assert hit["am_pathogenicity"] == 0.9 and "primed" in hit["_source"]
    miss = alphamissense.lookup(_v(200))       # primed None -> definitive no-score
    assert miss["am_pathogenicity"] is None and "primed" in miss["_source"]
    alphamissense._reset_for_tests()


def test_prime_is_idempotent(monkeypatch):
    alphamissense._reset_for_tests()
    _fake_store(monkeypatch, {100: 0.9})
    assert alphamissense.prime([_v(100)]) == 1
    assert alphamissense.prime([_v(100)]) == 0   # already primed -> not recounted
    alphamissense._reset_for_tests()


def test_unprimed_variant_falls_through(monkeypatch):
    alphamissense._reset_for_tests()
    _fake_store(monkeypatch, {100: 0.9})
    alphamissense.prime([_v(100)])
    seen = []
    monkeypatch.setattr(cache, "get", lambda s, k: (seen.append(k), None)[1])
    alphamissense.lookup(_v(300))              # not primed -> normal cache/tabix path
    assert seen == ["1-300-A-T"]               # cache WAS consulted for the unprimed one
    alphamissense._reset_for_tests()


def test_prime_noop_without_file(monkeypatch):
    alphamissense._reset_for_tests()
    monkeypatch.setattr(alphamissense_parquet, "available", lambda: False)  # no parquet either
    monkeypatch.setattr(alphamissense, "_open", lambda: None)   # no local file/pysam
    assert alphamissense.prime([_v(100)]) == 0
    alphamissense._reset_for_tests()


def test_pipeline_primes_am_no_disk_churn(monkeypatch):
    # In the pipeline, candidates are primed -> add_alphamissense must NOT hit cache.put
    # for AlphaMissense (the O(n^2) per-key disk write we removed).
    alphamissense._reset_for_tests()
    _fake_store(monkeypatch, {})               # every candidate -> no score, via prime
    puts = []
    monkeypatch.setattr(cache, "put", lambda source, k, val: puts.append(source))
    report = pipeline.run_pipeline(config.SAMPLE_VCF)
    assert report.classifications
    assert "alphamissense" not in puts
    alphamissense._reset_for_tests()
