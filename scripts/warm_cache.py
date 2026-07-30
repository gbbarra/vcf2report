#!/usr/bin/env python3
"""Pre-fetch annotations for a VCF into the on-disk cache.

Run this once with network access before a demo so the whole run is then
network-independent (set OFFLINE=1 for the demo itself). For each variant in the
VCF it resolves gnomAD + ClinVar (live if online, else the bundled snapshot) and
writes the result into data/cache/, so subsequent lookups are instant and offline.

    python scripts/warm_cache.py [VCF]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import config  # noqa: E402
from vcf2report.annotate import cache, clinvar, gnomad  # noqa: E402
from vcf2report.vcf.parse import parse_vcf  # noqa: E402


def main() -> int:
    vcf = sys.argv[1] if len(sys.argv) > 1 else str(config.SAMPLE_VCF)
    variants, _build, _ = parse_vcf(vcf)
    n = 0
    skipped = 0
    for v in variants:
        # ONLY persist positive answers. Caching a "not found" sentinel writes an absence of
        # data into a store that reads back as data: gnomad.lookup / clinvar.lookup return
        # a dict either way, so a cached {significance: None} is indistinguishable from a
        # real "ClinVar says nothing about this variant" — and it survives every later run,
        # including ones where the store DOES have the record. Warming with a store
        # unmounted, or before a weekly ClinVar refresh, poisoned those keys permanently.
        # A miss simply stays uncached: the next lookup re-resolves it, which is the correct
        # cost of not knowing.
        g = gnomad.lookup(v)
        if g.get("af") is not None:
            cache.put("gnomad", v.key,
                      {k: g[k] for k in ("af", "ac", "an", "hom", "pop") if k in g})
        else:
            skipped += 1
        cv = clinvar.lookup(v)
        if cv.get("significance"):
            cache.put("clinvar", v.key, {k: cv.get(k) for k in
                      ("significance", "review_status", "accession", "condition", "date")})
        n += 1
    print(f"Warmed cache for {n} variants into {config.CACHE_DIR}")
    if skipped:
        print(f"  {skipped} variant(s) had no gnomAD frequency and were left UNCACHED "
              f"(an unknown must not be persisted as an answer).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
