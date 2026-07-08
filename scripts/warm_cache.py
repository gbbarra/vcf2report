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
    for v in variants:
        g = gnomad.lookup(v)
        cache.put("gnomad", v.key, {k: g[k] for k in ("af", "ac", "an", "hom", "pop") if k in g})
        cv = clinvar.lookup(v)
        cache.put("clinvar", v.key, {k: cv.get(k) for k in
                  ("significance", "review_status", "accession", "condition", "date")})
        n += 1
    print(f"Warmed cache for {n} variants into {config.CACHE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
