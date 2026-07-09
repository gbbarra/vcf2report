#!/usr/bin/env python3
"""Annotate the spiked variants of a synthetic VCF with LIVE gnomAD frequency.

The background exome is unannotated in-sandbox (no SnpEff/vcfanno), so it is
filtered at the impact stage and only the spiked candidates need frequency data.
Rather than let the pipeline hit gnomAD live for all ~25k variants, this queries
gnomAD v4.1 (grpmax, remote tabix) for just the ``SPIKED=1`` records and writes
``gnomad_AF``/``AC``/``AN`` into their INFO — exactly the field vcfanno would add,
but only where it matters. Honest: the frequency is the real live gnomAD value.

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/bake_spiked_gnomad.py in.vcf
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def main() -> int:
    from vcf2report import config
    if config.offline():
        raise SystemExit("network egress required: set VCF2REPORT_ALLOW_NETWORK=1")
    from vcf2report.annotate import gnomad_remote
    from vcf2report.models import Variant

    gnomad_hdr = [
        '##INFO=<ID=gnomad_AF,Number=A,Type=Float,Description="gnomAD v4.1 grpmax allele frequency (live remote tabix; baked for spiked variants)">',
        '##INFO=<ID=gnomad_AC,Number=A,Type=Integer,Description="gnomAD v4.1 allele count">',
        '##INFO=<ID=gnomad_AN,Number=1,Type=Integer,Description="gnomAD v4.1 allele number">',
    ]

    path = Path(sys.argv[1])
    out = []
    n = 0
    hdr_written = False
    for line in path.read_text().splitlines():
        # Declare the gnomad_* INFO fields once, just before the #CHROM line, so
        # the baked values are self-describing and parsers don't warn.
        if line.startswith("#CHROM") and not hdr_written:
            if not any("ID=gnomad_AF" in l for l in out):
                out.extend(gnomad_hdr)
            hdr_written = True
        if line.startswith("#") or "SPIKED=1" not in line:
            out.append(line)
            continue
        f = line.split("\t")
        chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
        r = gnomad_remote.query(Variant(chrom=chrom, pos=pos, ref=ref, alt=alt))
        af = (r or {}).get("af", 0.0)
        ac = (r or {}).get("ac", 0) or 0
        an = (r or {}).get("an", 0) or 0
        f[7] = f"gnomad_AF={af};gnomad_AC={ac};gnomad_AN={an};" + f[7]
        out.append("\t".join(f))
        n += 1
        print(f"  {chrom}:{pos} {ref}>{alt[:12]} -> live gnomAD AF={af}", file=sys.stderr)
    path.write_text("\n".join(out) + "\n")
    print(f"baked live gnomAD into {n} spiked variants of {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
