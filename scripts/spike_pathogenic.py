#!/usr/bin/env python3
"""Spike REAL ClinVar pathogenic variants into a (BED-filtered) real exome VCF.

Turns a real, healthy 1000G exome into a controlled *synthetic* case: real
background + known pathogenic variants planted at their true GRCh38 coordinates
(pulled straight from the ClinVar VCF, so consequence/HGVS annotate correctly and
the calls reproduce reliably). De-identifies the sample column.

Usage:
  python spike_pathogenic.py \
      --exome  sample.exome.vcf.gz \
      --clinvar clinvar_GRCh38.vcf.gz \
      --targets spike_targets.tsv \
      --sample-id SYN-001 \
      --out sample.synthetic.vcf

`spike_targets.tsv` (tab-separated, '#'=comment):
  gene   category    zygosity
  SCN1A  primary     het
  RB1    secondary   het

Notes on downstream ACMG (vcf2report):
  * CLNSIG + CLNREVSTAT are carried so the engine's PP5 (reputable ClinVar,
    >=1-star) fires; CLNDN carries the disease name into the report.
  * For a spiked LoF to reach Pathogenic/Likely-Pathogenic via PVS1, the gene must
    be LoF-intolerant in data/constraint/gene_constraint.tsv (ClinVar alone is only
    PP5/supporting in this engine, by design).
"""
from __future__ import annotations
import argparse, gzip, sys
from pathlib import Path

# ClinVar MC (SO term names) -> canonical consequence vcf2report understands.
MC_MAP = {
    "nonsense": "stop_gained", "frameshift_variant": "frameshift_variant",
    "splice_donor_variant": "splice_donor_variant",
    "splice_acceptor_variant": "splice_acceptor_variant",
    "initiator_codon_variant": "start_lost", "stop_lost": "stop_lost",
    "missense_variant": "missense_variant",
    "synonymous_variant": "synonymous_variant",
}
LOF = {"stop_gained", "frameshift_variant", "splice_donor_variant",
       "splice_acceptor_variant", "start_lost", "stop_lost"}
_CHROM_ORDER = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "M": 25, "MT": 25}


def _open(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def parse_info(s: str) -> dict:
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); d[k] = v
        elif kv:
            d[kv] = True
    return d


def clnsig_pathogenic(clnsig: str) -> bool:
    c = clnsig.lower()
    return "pathogenic" in c and "conflicting" not in c  # P + LP, not conflicting


def mc_consequence(info: dict):
    terms = [p.split("|", 1)[1] for p in info.get("MC", "").split(",") if "|" in p]
    mapped = [MC_MAP.get(t, t) for t in terms]
    for t in mapped:
        if t in LOF:
            return t, True
    return (mapped[0], False) if mapped else (None, False)


def collect_from_clinvar(clinvar_path, genes):
    """Stream ClinVar once; keep pathogenic records for the target genes."""
    want = set(genes)
    hits = {g: [] for g in want}
    with _open(clinvar_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            if alt in (".", "") or "," in alt or ref in (".", ""):
                continue
            info = parse_info(f[7])
            gene = info.get("GENEINFO", "").split(":")[0] or None
            if gene not in want or not clnsig_pathogenic(info.get("CLNSIG", "")):
                continue
            csq, is_lof = mc_consequence(info)
            rs = info.get("RS", "")
            hits[gene].append({
                "chrom": chrom.replace("chr", ""), "pos": int(pos), "ref": ref, "alt": alt,
                "gene": gene, "clnsig": info.get("CLNSIG", ""), "csq": csq or "missense_variant",
                "is_lof": is_lof, "rs": ("rs" + rs) if rs else ".", "vid": vid,
                "clndn": info.get("CLNDN", ""),
                # review status drives PP5 (reputable-source, >=1-star) in the engine
                "clnrevstat": info.get("CLNREVSTAT", ""),
            })
    picked = {}
    for g, recs in hits.items():
        if recs:
            recs.sort(key=lambda r: (0 if r["is_lof"] else 1, r["pos"]))  # LoF first, deterministic
            picked[g] = recs[0]
    return picked


def load_exome(exome_path):
    meta, col_line, records, style = [], None, [], None
    with _open(exome_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                meta.append(line)
            elif line.startswith("#CHROM"):
                col_line = line
            elif line:
                f = line.split("\t")
                records.append(f)
                if style is None:
                    style = "chr" if f[0].startswith("chr") else "plain"
    if col_line is None:
        sys.exit("ERROR: no #CHROM header line found in exome VCF")
    return meta, col_line, records, (style or "plain")


def read_targets(path):
    out = []
    for line in _open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        gene = parts[0]
        if gene.lower() == "gene":  # skip optional header row
            continue
        category = parts[1] if len(parts) > 1 else "primary"
        zyg = parts[2] if len(parts) > 2 else "het"
        out.append((gene, category, zyg))
    return out


def spiked_line(rec, style, zyg, sample_id, col_count):
    chrom = ("chr" + rec["chrom"]) if style == "chr" else rec["chrom"]
    gt = "1/1" if zyg == "hom" else "0/1"
    dp, ad = (44, "0,44") if zyg == "hom" else (44, "22,22")
    info = (f"GENE={rec['gene']};CSQ={rec['csq']};CLNSIG={rec['clnsig']};"
            f"CLNREVSTAT={rec['clnrevstat'] or '.'};CLNDN={rec['clndn'] or '.'};"
            f"CLNVID={rec['vid']};SPIKED=1")
    fixed = [chrom, str(rec["pos"]), rec["rs"], rec["ref"], rec["alt"], "800", "PASS",
             info, "GT:DP:GQ:AD", f"{gt}:{dp}:99:{ad}"]
    # pad to the exome's column count (single-sample expected; extra sample cols = '.')
    while len(fixed) < col_count:
        fixed.append(".")
    return fixed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exome", required=True)
    ap.add_argument("--clinvar", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--sample-id", default="SYN-001")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    targets = read_targets(a.targets)
    genes = [g for g, _, _ in targets]
    picked = collect_from_clinvar(a.clinvar, genes)

    missing = [g for g in genes if g not in picked]
    if missing:
        print(f"WARNING: no pathogenic ClinVar record found for: {', '.join(missing)}",
              file=sys.stderr)

    meta, col_line, records, style = load_exome(a.exome)
    cols = col_line.split("\t")
    col_count = len(cols)

    spikes = []
    for gene, category, zyg in targets:
        if gene in picked:
            r = picked[gene]
            spikes.append(spiked_line(r, style, zyg, a.sample_id, col_count))
            print(f"  spiked {gene:8s} {style}:{r['pos']} {r['ref']}>{r['alt']} "
                  f"[{r['csq']}, {r['clnsig']}, {category}]  {r['clndn'][:40]}",
                  file=sys.stderr)

    all_rows = records + spikes
    all_rows.sort(key=lambda f: (_CHROM_ORDER.get(f[0].replace("chr", ""), 99), int(f[1])))

    # De-identify: single opaque sample id in the genotype column(s).
    new_cols = cols[:9] + [a.sample_id] + ["." for _ in cols[10:]]

    extra_meta = [
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol (spiked)">',
        '##INFO=<ID=CSQ,Number=1,Type=String,Description="Molecular consequence (spiked)">',
        '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar significance (spiked)">',
        '##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description="ClinVar review status (spiked)">',
        '##INFO=<ID=CLNDN,Number=1,Type=String,Description="ClinVar disease name (spiked)">',
        '##INFO=<ID=CLNVID,Number=1,Type=String,Description="ClinVar Variation ID (spiked)">',
        '##INFO=<ID=SPIKED,Number=0,Type=Flag,Description="Synthetically spiked pathogenic variant">',
        "##comment=SYNTHETIC: real 1000G background with ClinVar pathogenic variants "
        "spiked in for demonstration. De-identified. Not real patient data. Not for clinical use.",
    ]
    have = set(meta)
    meta_out = meta + [m for m in extra_meta if m not in have]

    with open(a.out, "w") as out:
        out.write("\n".join(meta_out) + "\n")
        out.write("\t".join(new_cols) + "\n")
        for f in all_rows:
            out.write("\t".join(f) + "\n")

    print(f"Wrote {a.out}: {len(records)} background + {len(spikes)} spiked "
          f"= {len(all_rows)} variants (sample '{a.sample_id}')", file=sys.stderr)


if __name__ == "__main__":
    main()
