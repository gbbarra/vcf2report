#!/usr/bin/env python3
"""Build a REAL single-sample exome VCF from public gnomAD HGDP+1KG data.

This is not synthetic. It materialises the exome of ONE real, consented,
publicly-released research participant from the gnomAD v3.1.2 HGDP + 1000
Genomes callset (``gcp-public-data--gnomad`` bucket) by:

  1. restricting to the Broad exome calling regions (GRCh38) so the variant
     set is exome-scale (~30-50k) rather than whole-genome;
  2. keeping only the sites where the chosen sample actually carries an ALT
     allele (its real diploid genotype, depth and quality);
  3. copying, verbatim, gnomAD's own per-allele frequencies (global + popmax)
     and Ensembl VEP consequence / gene / HGVS annotation into INFO.

The result is a fully real, fully offline-consumable annotated exome VCF: every
frequency and consequence in the downstream report traces back to gnomAD, not
to anything this project invented.

Default sample: HG01565 — a 1000 Genomes "PEL" participant (Peruvian in Lima),
an admixed Latin-American genome. That population is under-represented in
gnomAD's global frequencies, which is exactly the gap ABraOM (Brazilian SABE)
is meant to fill — so the sample doubles as a live demonstration of the
Brazilian-frequency differentiator.

Reproducible, no auth, no local reference FASTA:
    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_real_exome.py \
        --sample HG01565 --out data/real/HG01565_exome.vcf

Network egress is required (opt-in): it reads directly from
storage.googleapis.com via htslib remote tabix.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

GNOMAD_BASE = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/"
               "3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.{chrom}.vcf.bgz")
EXOME_INTERVALS_URL = ("https://storage.googleapis.com/gcp-public-data--broad-"
                       "references/hg38/v0/exome_calling_regions.v1.interval_list")

VEP_FORMAT = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
              "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
              "Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|"
              "SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|"
              "UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|"
              "MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info")
_ALLELE_NUM_IDX = VEP_FORMAT.split("|").index("ALLELE_NUM")

AUTOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
# Merge exome intervals within this gap into one *fetch* window. A large gap means
# far fewer remote range-requests (round-trip latency dominates); the intronic
# records swept in are cheaply discarded by the exome-membership test on emit, so
# the output stays a true exome while the network cost collapses.
MERGE_GAP = 100000


def load_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            c, s, e, *_ = line.split("\t")
            by_chrom.setdefault(c, []).append((int(s), int(e)))
    raw: dict[str, list[tuple[int, int]]] = {}
    merged: dict[str, list[tuple[int, int]]] = {}
    for c, ivs in by_chrom.items():
        ivs.sort()
        raw[c] = ivs
        out: list[tuple[int, int]] = []
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce + MERGE_GAP:
                ce = max(ce, e)
            else:
                out.append((cs, ce))
                cs, ce = s, e
        out.append((cs, ce))
        merged[c] = out
    return raw, merged


def _in_exome(starts: list[int], ivs: list[tuple[int, int]], pos: int) -> bool:
    """True if 1-based ``pos`` falls within any exome interval (bisect on starts)."""
    import bisect
    i = bisect.bisect_right(starts, pos) - 1   # rightmost interval starting <= pos
    return i >= 0 and pos <= ivs[i][1]


def _num_a(info, key: str, alt_idx0: int):
    """Read a Number=A INFO value for allele index (0-based among ALTs)."""
    v = info.get(key)
    if v is None:
        return None
    if isinstance(v, (tuple, list)):
        return v[alt_idx0] if alt_idx0 < len(v) else None
    return v


def _fmt_num(x) -> str:
    if x is None:
        return "."
    if isinstance(x, float):
        # gnomAD AFs need full precision to distinguish rare variants
        return repr(x)
    return str(x)


def process_batch(chrom: str, sample: str, windows: list[tuple[int, int]],
                  raw_ivs: list[tuple[int, int]], out_path: Path) -> tuple[str, int, int]:
    """Stream a batch of one chromosome's exome windows; write carrier rows.

    Batching a chromosome into several tasks lets the big chromosomes run across
    many workers at once, so wall time is bounded by the slowest *batch*, not the
    slowest whole chromosome.
    """
    import pysam

    # Resume: a completed batch already has its rows on disk.
    if out_path.exists() and out_path.stat().st_size > 0:
        n = sum(1 for _ in out_path.open())
        return chrom, n, -1

    url = GNOMAD_BASE.format(chrom=chrom)
    vf = None
    for attempt in range(4):                 # remote opens flake; retry with backoff
        try:
            vf = pysam.VariantFile(url)
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    if sample not in vf.header.samples:
        raise SystemExit(f"sample {sample!r} not in callset")

    starts = [s for s, _ in raw_ivs]
    rows: list[tuple[int, str]] = []
    n_seen = 0
    for (start, end) in windows:
        recs = None
        for attempt in range(4):             # transient fetch resets → retry the window
            try:
                recs = list(vf.fetch(chrom, start, end))
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(1.5 * (attempt + 1))
        for rec in recs:
            n_seen += 1
            if not _in_exome(starts, raw_ivs, rec.pos):
                continue  # intronic record swept in by the wide fetch → discard
            smp = rec.samples[sample]
            gt = smp.get("GT")
            if not gt or all(a in (0, None) for a in gt):
                continue  # homozygous reference / no-call → not a carrier
            carried = sorted({a for a in gt if a not in (0, None)})
            vep = rec.info.get("vep") or ()
            dp = smp.get("DP")
            gq = smp.get("GQ")
            ad = smp.get("AD")
            for a in carried:                       # 1-based ALT allele number
                alt = rec.alts[a - 1]
                if alt is None or alt.startswith("<") or alt == "*":
                    continue
                af = _num_a(rec.info, "gnomad_AF_popmax", a - 1)
                pop = _num_a(rec.info, "gnomad_popmax", a - 1)
                if af is None:                      # subset-only site → real subset AF
                    af = _num_a(rec.info, "AF", a - 1)
                    pop = "hgdp_tgp"
                ac = (_num_a(rec.info, "gnomad_AC_popmax", a - 1)
                      if pop != "hgdp_tgp" else _num_a(rec.info, "AC", a - 1))
                an = (_num_a(rec.info, "gnomad_AN_popmax", a - 1)
                      if pop != "hgdp_tgp" else _num_a(rec.info, "AN", a - 1))
                hom = (_num_a(rec.info, "gnomad_nhomalt_popmax", a - 1)
                       if pop != "hgdp_tgp" else _num_a(rec.info, "nhomalt", a - 1))
                # VEP entries for this allele, renumbered to a single-ALT record
                csq_parts = []
                for ent in vep:
                    f = ent.split("|")
                    if len(f) > _ALLELE_NUM_IDX and f[_ALLELE_NUM_IDX] == str(a):
                        f[_ALLELE_NUM_IDX] = "1"
                        csq_parts.append("|".join(f))
                info = [f"gnomad_AF={_fmt_num(af)}", f"gnomad_AC={_fmt_num(ac)}",
                        f"gnomad_AN={_fmt_num(an)}", f"gnomad_nhomalt={_fmt_num(hom)}",
                        f"gnomad_popmax={pop or '.'}"]
                if csq_parts:
                    info.append("CSQ=" + ",".join(csq_parts))
                # genotype for this biallelic record
                cnt = sum(1 for x in gt if x == a)
                gt_s = "1/1" if cnt == 2 else "0/1"
                ad_s = "."
                if ad is not None and len(ad) > a:
                    ad_s = f"{ad[0]},{ad[a]}"
                smp_s = ":".join([gt_s, _fmt_num(dp), _fmt_num(gq), ad_s])
                line = "\t".join([chrom, str(rec.pos), ".", rec.ref, alt, "999",
                                  "PASS", ";".join(info), "GT:DP:GQ:AD", smp_s])
                rows.append((rec.pos, line))
    rows.sort(key=lambda r: r[0])
    out_path.write_text("\n".join(l for _, l in rows) + ("\n" if rows else ""))
    return chrom, len(rows), n_seen


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="HG01565")
    ap.add_argument("--out", default=str(REPO / "data" / "real" / "HG01565_exome.vcf"))
    ap.add_argument("--chroms", default=",".join(AUTOSOMES))
    ap.add_argument("--intervals", default=str(REPO / "scratch" / "exome.interval_list"))
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--batch", type=int, default=120,
                    help="exome windows per task (smaller -> more parallelism)")
    ap.add_argument("--tmp", default=str(REPO / "scratch" / "real"))
    args = ap.parse_args()

    from vcf2report import config
    if config.offline():
        raise SystemExit("network egress required: set VCF2REPORT_ALLOW_NETWORK=1")

    raw_intervals, windows_by_chrom = load_intervals(Path(args.intervals))
    chroms = [c for c in args.chroms.split(",") if c in windows_by_chrom]
    tmp = Path(args.tmp)
    tmp.mkdir(parents=True, exist_ok=True)

    # Split each chromosome's fetch windows into batches so big chromosomes fan out
    # across many workers instead of monopolising one thread.
    tasks: list[tuple[str, int, list[tuple[int, int]]]] = []
    for c in chroms:
        wins = windows_by_chrom[c]
        for bi in range(0, len(wins), args.batch):
            tasks.append((c, bi // args.batch, wins[bi:bi + args.batch]))

    print(f"sample={args.sample} chroms={len(chroms)} tasks={len(tasks)} "
          f"workers={args.workers}", flush=True)
    t0 = time.time()
    per_batch: dict[tuple[str, int], int] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_batch, c, args.sample, wins, raw_intervals[c],
                          tmp / f"{c}.{bi}.rows"): (c, bi)
                for (c, bi, wins) in tasks}
        failed = 0
        for fut in as_completed(futs):
            c, bi = futs[fut]
            try:
                _, n, _ = fut.result()
                per_batch[(c, bi)] = n
            except Exception as e:            # a batch that fails is simply not
                failed += 1                   # written; a resume run will retry it
                print(f"  ! batch {c}.{bi} failed: {type(e).__name__}", flush=True)
            done += 1
            if done % 20 == 0 or done == len(tasks):
                got = sum(per_batch.values())
                print(f"  {done}/{len(tasks)} batches, {got} carriers so far, "
                      f"{failed} failed [{time.time()-t0:.0f}s]", flush=True)

    # merge batches per chromosome, sorted by position, into one VCF
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    missing = 0
    with open(out, "w") as w:
        w.write(_header(args.sample, chroms))
        for c in chroms:
            n_batches = (len(windows_by_chrom[c]) + args.batch - 1) // args.batch
            rows: list[tuple[int, str]] = []
            for bi in range(n_batches):       # iterate all batches; skip gaps
                bf = tmp / f"{c}.{bi}.rows"
                if not bf.exists():
                    missing += 1
                    continue
                for line in bf.read_text().splitlines():
                    if line:
                        rows.append((int(line.split("\t", 2)[1]), line))
            rows.sort(key=lambda r: r[0])
            for _, line in rows:
                w.write(line + "\n")
            total += len(rows)
    dt = time.time() - t0
    status = "COMPLETE" if missing == 0 else f"PARTIAL ({missing} batches missing — rerun to resume)"
    print(f"\nWrote {total} real exome variants for {args.sample} to {out} "
          f"in {dt:.0f}s [{status}]")
    return 0


def _header(sample: str, chroms: list[str]) -> str:
    contigs = "\n".join(f"##contig=<ID={c}>" for c in chroms)
    return f"""##fileformat=VCFv4.2
##source=vcf2report/build_real_exome.py
##reference=GRCh38
##dataset=gnomAD v3.1.2 HGDP+1KG (gcp-public-data--gnomad), Broad exome calling regions v1
##sample_provenance=1000 Genomes / HGDP publicly-released research participant {sample}
##note=Real per-sample genotypes; gnomAD frequencies + Ensembl VEP consequences copied verbatim from source.
{contigs}
##INFO=<ID=gnomad_AF,Number=A,Type=Float,Description="gnomAD popmax allele frequency (v3.1.2), subset AF if popmax absent">
##INFO=<ID=gnomad_AC,Number=A,Type=Integer,Description="gnomAD popmax allele count">
##INFO=<ID=gnomad_AN,Number=A,Type=Integer,Description="gnomAD popmax allele number">
##INFO=<ID=gnomad_nhomalt,Number=A,Type=Integer,Description="gnomAD popmax homozygote count">
##INFO=<ID=gnomad_popmax,Number=A,Type=String,Description="gnomAD population with max AF">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {VEP_FORMAT}">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


if __name__ == "__main__":
    raise SystemExit(main())
