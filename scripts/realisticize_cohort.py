#!/usr/bin/env python3
"""Rewrite an already-spiked synthetic VCF so the planted variant is indistinguishable from a real
call, and emit a traceability sidecar.

The default synthetic plant is trivially identifiable (SPIKED/GENE/CSQ/CLN* INFO + a minimal
GT:DP:GQ:AD FORMAT). This removes every marker from the planted record(s) and re-inserts a tell-free
record that BORROWS a real background call's INFO/FORMAT (full DRAGEN fields) of the same zygosity,
relocated to the plant's coordinate/alleles — statistically and visually a genuine call. The engine
provably ignores the markers (byte-identical classification), so this changes no result; it removes
the cosmetic tell only.

Truth is NOT lost: every planted variant is written to a per-VCF sidecar `<stem>.planted.tsv` and,
optionally, appended to a cohort-wide manifest. The truth was extracted from the SPIKED record's own
INFO before it was stripped.

    python3 scripts/realisticize_cohort.py IN.vcf.gz OUT.vcf [--sidecar OUT.planted.tsv] \
        [--syn-id SYN-101] [--manifest planted_variants.tsv]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spike_pathogenic import _CHROM_ORDER, load_exome  # noqa: E402
from spike_variant import spiked_line_realistic  # noqa: E402

_SIDE_COLS = ["syn_id", "chrom", "pos", "ref", "alt", "gene", "zygosity", "allele",
              "consequence", "clnsig", "clnrevstat", "clnvid"]


def _info(info: str) -> dict[str, str]:
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        elif kv:
            d[kv] = "1"
    return d


def _gt(fmt: str, sample: str) -> str:
    keys = fmt.split(":")
    if "GT" not in keys:
        return "0/1"
    return sample.split(":")[keys.index("GT")].replace("|", "/")


def realisticize(in_vcf: str, out_vcf: str, syn_id: str):
    meta, col_line, records, style = load_exome(in_vcf)
    planted, background = [], []
    for f in records:
        info = _info(f[7])
        if "SPIKED" in info or "SPIKED2" in info:
            planted.append(f)
        else:
            background.append(f)
    if not planted:
        raise SystemExit(f"ERROR: no SPIKED/SPIKED2 record in {in_vcf}")

    ncols = len(col_line.split("\t"))
    new_planted, rows = [], []
    for f in planted:
        info = _info(f[7])
        gt = _gt(f[8], f[9])
        zyg = "hom" if gt == "1/1" else "het"
        rec = {"chrom": f[0].replace("chr", ""), "pos": f[1], "ref": f[3], "alt": f[4],
               "gene": info.get("GENE", "")}
        r = spiked_line_realistic(rec, background, style, zyg, ncols)
        if r is None:
            raise SystemExit(f"ERROR: no {zyg} template to borrow in {in_vcf} for {f[0]}:{f[1]}")
        new_planted.append(r)
        rows.append({"syn_id": syn_id, "chrom": f[0], "pos": f[1], "ref": f[3], "alt": f[4],
                     "gene": info.get("GENE", ""), "zygosity": zyg,
                     "allele": "second" if "SPIKED2" in info else "primary",
                     "consequence": info.get("CSQ", ""), "clnsig": info.get("CLNSIG", ""),
                     "clnrevstat": info.get("CLNREVSTAT", ""), "clnvid": info.get("CLNVID", "")})

    all_rows = background + new_planted
    all_rows.sort(key=lambda x: (_CHROM_ORDER.get(x[0].replace("chr", ""), 99), int(x[1])))
    # strip the now-unused spike INFO definitions + keep everything else; header comment kept honest
    keep_meta = [m for m in meta if not any(t in m for t in (
        "ID=SPIKED", "ID=SPIKED2", "ID=GENE,", "ID=CSQ,",
        "ID=CLNSIG,", "ID=CLNREVSTAT,", "ID=CLNDN,", "ID=CLNVID,"))]
    with open(out_vcf, "w") as out:
        out.write("\n".join(keep_meta) + "\n")
        out.write(col_line + "\n")
        for x in all_rows:
            out.write("\t".join(x) + "\n")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("in_vcf")
    ap.add_argument("out_vcf")
    ap.add_argument("--syn-id", required=True)
    ap.add_argument("--sidecar")
    ap.add_argument("--manifest")
    a = ap.parse_args()
    rows = realisticize(a.in_vcf, a.out_vcf, a.syn_id)
    if a.sidecar:
        with open(a.sidecar, "w") as s:
            s.write("\t".join(_SIDE_COLS) + "\n")
            for r in rows:
                s.write("\t".join(str(r[c]) for c in _SIDE_COLS) + "\n")
    if a.manifest:
        new = not Path(a.manifest).exists()
        with open(a.manifest, "a") as m:
            if new:
                m.write("\t".join(_SIDE_COLS) + "\n")
            for r in rows:
                m.write("\t".join(str(r[c]) for c in _SIDE_COLS) + "\n")
    print(f"  {a.syn_id}: {len(rows)} planted allele(s) realisticised -> {a.out_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
