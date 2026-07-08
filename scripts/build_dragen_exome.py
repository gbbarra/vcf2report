#!/usr/bin/env python3
"""In-sandbox real exome from a 1000G DRAGEN 4.4.7 single-sample VCF.

A pure-Python alternative to make_synthetic_exomes.sh's download+bcftools steps:
reads a per-sample DRAGEN small-variant VCF over remote tabix (htslib HTTPS range
queries against the public s3://1000genomes-dragen-v4-4-7 bucket), keeps only the
PASS variants that fall inside the IDT xGen Exome v2 target BED, and writes a
single-sample exome VCF with the participant's real genotype/DP/GQ/AD. No AWS CLI,
no bcftools, no reference FASTA, no 440 MB download of the whole file.

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_dragen_exome.py \
        --sample NA12878 --bed scratch/idt_exome_v2.bed --out data/real/NA12878_exome.vcf
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DRAGEN = ("https://1000genomes-dragen-v4-4-7.s3.amazonaws.com/data/individuals/"
          "hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2/{s}/{s}.hard-filtered.vcf.gz")
AUTOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def load_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    by: dict[str, list[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        by.setdefault(f[0], []).append((int(f[1]), int(f[2])))  # BED is 0-based half-open
    for c in by:
        by[c].sort()
    return by


def _member(starts, ivs, pos1):
    import bisect
    i = bisect.bisect_right(starts, pos1 - 1) - 1   # BED start is 0-based; VCF pos 1-based
    return i >= 0 and (pos1 - 1) < ivs[i][1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="NA12878")
    ap.add_argument("--bed", default=str(REPO / "scratch" / "idt_exome_v2.bed"))
    ap.add_argument("--out", default=str(REPO / "data" / "real" / "NA12878_exome.vcf"))
    ap.add_argument("--chroms", default=",".join(AUTOSOMES))
    args = ap.parse_args()

    from vcf2report import config
    if config.offline():
        raise SystemExit("network egress required: set VCF2REPORT_ALLOW_NETWORK=1")
    import pysam

    bed = load_bed(Path(args.bed))
    url = DRAGEN.format(s=args.sample)
    vf = pysam.VariantFile(url)
    if args.sample not in vf.header.samples:
        raise SystemExit(f"{args.sample} not in {url}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    chroms = [c for c in args.chroms.split(",") if c in bed]
    t0 = time.time()
    rows: list[str] = []
    for c in chroms:
        ivs = bed[c]
        starts = [s for s, _ in ivs]
        lo, hi = ivs[0][0], ivs[-1][1]
        n = kept = 0
        for rec in vf.fetch(c, lo, hi):          # one streamed pass over the exome span
            n += 1
            if list(rec.filter.keys()) not in ([], ["PASS"]):
                continue
            if not _member(starts, ivs, rec.pos):
                continue
            smp = rec.samples[args.sample]
            gt = smp.get("GT")
            if not gt or all(a in (0, None) for a in gt):
                continue
            for ai in sorted({a for a in gt if a not in (0, None)}):
                alt = rec.alts[ai - 1]
                if not alt or alt.startswith("<") or alt == "*":
                    continue
                cnt = sum(1 for x in gt if x == ai)
                gts = "1/1" if cnt == 2 else "0/1"
                dp = smp.get("DP") or "."
                gq = smp.get("GQ") or "."
                ad = smp.get("AD")
                ads = f"{ad[0]},{ad[ai]}" if ad and len(ad) > ai else "."
                rows.append("\t".join([c, str(rec.pos), ".", rec.ref, alt, "999", "PASS",
                                       ".", "GT:DP:GQ:AD", f"{gts}:{dp}:{gq}:{ads}"]))
                kept += 1
        print(f"  {c}: {kept} exome variants ({n} scanned) [{time.time()-t0:.0f}s]", flush=True)

    contigs = "\n".join(f"##contig=<ID={c}>" for c in chroms)
    header = (f"##fileformat=VCFv4.2\n##reference=GRCh38\n"
              f"##source=1000G DRAGEN v4.4.7 {args.sample}, subset to IDT xGen Exome v2 targets\n"
              f"##note=Real per-sample genotypes (DRAGEN 4.4.7); PASS variants inside the exome BED.\n"
              f"{contigs}\n"
              "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
              "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n"
              "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
              "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
              f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
    out.write_text(header + "\n".join(rows) + "\n")
    print(f"\nWrote {len(rows)} real exome variants for {args.sample} to {out} "
          f"in {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
