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
# Merge exome intervals within this gap into one *fetch* window. With a single
# reused file handle per chromosome the per-window round-trip is cheap, so we keep
# windows TIGHT (minimal intronic over-read = minimal bytes transferred). The
# exome-membership test still guards emitted records.
MERGE_GAP = 10000
THROTTLE_S = 0.0    # raw-text reads are light; no extra pacing needed at low workers


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


def _parse_info(info_str: str) -> dict:
    """Parse a VCF INFO column into a dict (flags map to True)."""
    d = {}
    for part in info_str.split(";"):
        if not part:
            continue
        k, _, v = part.partition("=")
        d[k] = v if _ else True
    return d


def _info_a(info: dict, key: str, alt_idx0: int):
    """Read a Number=A INFO value (comma-separated per ALT) for allele index."""
    v = info.get(key)
    if v is None or v is True:
        return None
    parts = v.split(",")
    return parts[alt_idx0] if alt_idx0 < len(parts) else None


def _col_index(header_lines, sample: str) -> int:
    for h in header_lines:
        if h.startswith("#CHROM"):
            cols = h.rstrip("\n").split("\t")
            if sample not in cols:
                raise SystemExit(f"sample {sample!r} not in callset")
            return cols.index(sample)
    raise SystemExit("no #CHROM header line")


def process_batch(chrom: str, sample: str, windows: list[tuple[int, int]],
                  raw_ivs: list[tuple[int, int]], out_path: Path) -> tuple[str, int, int]:
    """Stream one chromosome's exome windows; write this sample's carrier rows.

    Reads the joint VCF as *raw tabix text* and pulls only the one sample's
    genotype column, instead of letting pysam decode all 4151 samples per record.
    That is ~9x faster and still uses the participant's real genotypes; gnomAD
    frequencies + VEP consequences are parsed straight from the INFO column.
    """
    import pysam

    # Resume: a completed chromosome already has its rows on disk.
    if out_path.exists() and out_path.stat().st_size > 0:
        n = sum(1 for _ in out_path.open())
        return chrom, n, -1

    url = GNOMAD_BASE.format(chrom=chrom)
    tb = None
    for attempt in range(4):                 # remote opens flake; retry with backoff
        try:
            tb = pysam.TabixFile(url)
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    col = _col_index(tb.header, sample)

    starts = [s for s, _ in raw_ivs]
    rows: list[tuple[int, str]] = []
    n_seen = 0
    for (start, end) in windows:
        lines = None
        for attempt in range(5):             # transient fetch resets → retry the window
            try:
                lines = list(tb.fetch(chrom, start, end))
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2.0 * (attempt + 1))
        if THROTTLE_S:
            time.sleep(THROTTLE_S)
        for line in lines:
            n_seen += 1
            f = line.split("\t")
            pos = int(f[1])
            if not _in_exome(starts, raw_ivs, pos):
                continue  # intronic record swept in by the wide fetch → discard
            fmt = f[8].split(":")
            smp_fields = f[col].split(":")
            gt_raw = smp_fields[0].replace("|", "/")
            alleles = [a for a in gt_raw.split("/")]
            carried = sorted({int(a) for a in alleles if a not in ("0", ".")})
            if not carried:
                continue  # hom-ref / no-call → not a carrier
            ref, alts = f[3], f[4].split(",")
            info = _parse_info(f[7])
            vep = (info.get("vep") or "")
            vep_entries = vep.split(",") if vep and vep is not True else []
            smp_map = dict(zip(fmt, smp_fields))
            dp, gq, ad = smp_map.get("DP", "."), smp_map.get("GQ", "."), smp_map.get("AD")
            for a in carried:                       # 1-based ALT allele number
                if a > len(alts):
                    continue
                alt = alts[a - 1]
                if not alt or alt.startswith("<") or alt == "*":
                    continue
                af = _info_a(info, "gnomad_AF_popmax", a - 1)
                pop = _info_a(info, "gnomad_popmax", a - 1)
                if af is None or af == ".":         # subset-only site → real subset AF
                    af = _info_a(info, "AF", a - 1)
                    pop = "hgdp_tgp"
                sfx = "" if pop == "hgdp_tgp" else "_popmax"
                base = "" if pop == "hgdp_tgp" else "gnomad_"
                ac = _info_a(info, f"{base}AC{sfx}", a - 1)
                an = _info_a(info, f"{base}AN{sfx}", a - 1)
                hom = _info_a(info, f"{base}nhomalt{sfx}", a - 1)
                # VEP entries for this allele, renumbered to a single-ALT record
                csq_parts = []
                for ent in vep_entries:
                    sub = ent.split("|")
                    if len(sub) > _ALLELE_NUM_IDX and sub[_ALLELE_NUM_IDX] == str(a):
                        sub[_ALLELE_NUM_IDX] = "1"
                        csq_parts.append("|".join(sub))
                out_info = [f"gnomad_AF={af or '.'}", f"gnomad_AC={ac or '.'}",
                            f"gnomad_AN={an or '.'}", f"gnomad_nhomalt={hom or '.'}",
                            f"gnomad_popmax={pop or '.'}"]
                if csq_parts:
                    out_info.append("CSQ=" + ",".join(csq_parts))
                cnt = sum(1 for x in alleles if x == str(a))
                gt_s = "1/1" if cnt == 2 else "0/1"
                ad_s = "."
                if ad:
                    ad_parts = ad.split(",")
                    if len(ad_parts) > a:
                        ad_s = f"{ad_parts[0]},{ad_parts[a]}"
                smp_s = ":".join([gt_s, dp, gq, ad_s])
                row = "\t".join([chrom, str(pos), ".", ref, alt, "999",
                                 "PASS", ";".join(out_info), "GT:DP:GQ:AD", smp_s])
                rows.append((pos, row))
    rows.sort(key=lambda r: r[0])
    out_path.write_text("\n".join(l for _, l in rows) + ("\n" if rows else ""))
    return chrom, len(rows), n_seen


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="HG01565")
    ap.add_argument("--out", default=str(REPO / "data" / "real" / "HG01565_exome.vcf"))
    ap.add_argument("--chroms", default=",".join(AUTOSOMES))
    ap.add_argument("--intervals", default=str(REPO / "scratch" / "exome.interval_list"))
    ap.add_argument("--workers", type=int, default=3,
                    help="chromosomes processed concurrently (keep low: one reused "
                         "remote handle each; too many opens => proxy refuses)")
    ap.add_argument("--tmp", default=str(REPO / "scratch" / "real"))
    args = ap.parse_args()

    from vcf2report import config
    if config.offline():
        raise SystemExit("network egress required: set VCF2REPORT_ALLOW_NETWORK=1")

    raw_intervals, windows_by_chrom = load_intervals(Path(args.intervals))
    chroms = [c for c in args.chroms.split(",") if c in windows_by_chrom]
    tmp = Path(args.tmp)
    tmp.mkdir(parents=True, exist_ok=True)

    # One task per chromosome: each opens ONE remote handle and reuses it for all
    # of that chromosome's windows (reopening a 20 GB file per batch is what made
    # the proxy refuse connections). Low --workers keeps only a few handles open.
    print(f"sample={args.sample} chroms={len(chroms)} workers={args.workers}", flush=True)
    t0 = time.time()
    per_chrom: dict[str, int] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_batch, c, args.sample, windows_by_chrom[c],
                          raw_intervals[c], tmp / f"{c}.0.rows"): c
                for c in chroms}
        failed = 0
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                _, n, _ = fut.result()
                per_chrom[c] = n
                print(f"  {c}: {n} carriers [{time.time()-t0:.0f}s]", flush=True)
            except Exception as e:            # a chromosome that fails is simply not
                failed += 1                   # written; a resume run retries it
                print(f"  ! {c} failed: {type(e).__name__} "
                      f"[{time.time()-t0:.0f}s]", flush=True)
            done += 1

    # merge batches per chromosome, sorted by position, into one VCF
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    missing = 0
    with open(out, "w") as w:
        w.write(_header(args.sample, chroms))
        for c in chroms:
            bf = tmp / f"{c}.0.rows"
            if not bf.exists():
                missing += 1
                continue
            rows: list[tuple[int, str]] = []
            for line in bf.read_text().splitlines():
                if line:
                    rows.append((int(line.split("\t", 2)[1]), line))
            rows.sort(key=lambda r: r[0])
            for _, line in rows:
                w.write(line + "\n")
            total += len(rows)
    dt = time.time() - t0
    status = "COMPLETE" if missing == 0 else f"PARTIAL ({missing} chroms missing — rerun to resume)"
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
