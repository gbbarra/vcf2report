#!/usr/bin/env python3
"""Cohort v2 — make the lone-het autosomal-recessive cases BIALLELIC, so AR diagnostic
recovery is testable (a lone het in an AR-only gene is a carrier, not a diagnosis).

For each AR-only carrier gene in the SYN cohort:
  * **compound het** — plant a SECOND, distinct ClinVar Pathogenic/Likely-Pathogenic variant
    from the gene's MANE region (opposite allele), as a second heterozygous call. This is the
    most common real AR presentation and exercises the engine's "2+ P/LP hits in the gene ->
    not a carrier -> diagnostic" path.
  * **homozygous** — for the 7 genes with no second ClinVar variant, set the existing spike's
    genotype to 1/1 (consanguinity/founder pattern; the engine's hom exclusion routes it out
    of carrier).

Reads the v2 plan (scratch/v2_plan.json from select_second.py) + the existing raw synthetic
VCFs, writes modified raw VCFs to <out>/SYN-NNN.v2.vcf.gz for re-annotation. The 64 non-AR
cases are unchanged (copied), so the v2 cohort is a drop-in for the routing measurement.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
COHORT = REPO / "data" / "synthetic_cohort"


def _spike2_line(chrom, pos, ref, alt, gene, sample_cols):
    """A second heterozygous spiked call. SnpEff adds ANN on re-annotation; the engine's
    ClinVar lookup fires by coordinate — so we carry only gene + the SPIKED2 flag, not a
    hand-written consequence (which would risk disagreeing with the real annotation)."""
    info = f"GENE={gene};SPIKED2=1"
    fmt = "GT:DP:GQ:AD"
    gt = "0/1:40:99:20,20"
    row = [chrom, str(pos), ".", ref, alt, "800", "PASS", info, fmt, gt]
    return row + ["." for _ in sample_cols[10:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default=str(REPO / "data" / "synthetic_cohort" / "v2_plan.json"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    plan = json.load(open(a.plan))          # keyed by syn_id (faithful phenopacket plan)
    rows = list(csv.DictReader(open(COHORT / "cohort.tsv"), delimiter="\t"))
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    from spike_pathogenic import _CHROM_ORDER  # noqa

    made = {"compound_het": 0, "hom": 0, "copied": 0}
    v2_manifest = []
    for r in rows:
        gene, syn = r["gene"], r["syn_id"]
        src = COHORT / f"{syn}.synthetic.vcf.gz"
        dst = out / f"{syn}.v2.vcf.gz"
        p = plan.get(syn)
        if p and p["mode"] not in ("hom", "compound_het"):
            p = None                        # single_het / nomatch → unchanged (faithful: source had one allele)
        with gzip.open(src, "rt") as fh:
            lines = fh.read().splitlines()
        meta = [l for l in lines if l.startswith("##")]
        header = [l for l in lines if l.startswith("#CHROM")][0]
        cols = header.split("\t")
        body = [l.split("\t") for l in lines if l and not l.startswith("#")]

        if not p:                                   # non-AR case: copy unchanged
            made["copied"] += 1
            mode = "unchanged"
        elif p["mode"] == "hom":
            for f in body:                          # flip the spike genotype to 1/1
                if f[1] == r["pos"] and "SPIKED=1" in f[7]:
                    g = f[9].split(":"); g[0] = "1/1"
                    if len(g) >= 4:
                        g[3] = "2,40"
                    f[9] = ":".join(g)
            made["hom"] += 1
            mode = "hom"
        else:                                       # compound het: add the 2nd ClinVar variant
            body.append(_spike2_line(p["chrom"], p["pos2"], p["ref2"], p["alt2"], gene, cols))
            body.sort(key=lambda f: (_CHROM_ORDER.get(f[0].replace("chr", ""), 99), int(f[1])))
            made["compound_het"] += 1
            mode = "compound_het"

        extra = ['##INFO=<ID=SPIKED2,Number=0,Type=Flag,Description="Second biallelic spiked variant (v2)">']
        meta_out = meta + [m for m in extra if m not in meta]
        tmp = out / f"{syn}.v2.vcf"
        with open(tmp, "w") as w:
            w.write("\n".join(meta_out) + "\n" + header + "\n")
            for f in body:
                w.write("\t".join(f) + "\n")
        subprocess.run(["bgzip", "-f", str(tmp)], check=True)
        v2_manifest.append({"syn_id": syn, "gene": gene, "mode": mode,
                            **({"pos2": p["pos2"], "ref2": p["ref2"], "alt2": p["alt2"],
                                "zyg2": p.get("zyg2", "het")} if p and p["mode"] == "compound_het" else {})})

    json.dump(v2_manifest, open(out / "v2_manifest.json", "w"), indent=1)
    print(f"compound_het: {made['compound_het']} | hom: {made['hom']} | unchanged: {made['copied']}")
    print(f"-> {out}")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
