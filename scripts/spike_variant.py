#!/usr/bin/env python3
"""Spike ONE exact pathogenic variant (by coordinate) into a real exome VCF.

Unlike spike_pathogenic.py (which picks *a* pathogenic ClinVar record for a gene), this plants the
EXACT chrom:pos:ref:alt supplied — the variant from the source phenopacket case — and looks that
exact coordinate up in the ClinVar VCF to carry its CLNSIG/CLNREVSTAT/CLNDN (so PP5 fires + the
disease name reaches the report). If the coordinate is not in this ClinVar release, it still plants
the variant with the supplied consequence/disease and a synthetic Pathogenic CLNSIG. De-identifies.

  python spike_variant.py --exome sample.exome.vcf.gz --clinvar clinvar_GRCh38.vcf.gz \
      --chrom chr15 --pos 101652355 --ref C --alt A --gene TM2D3 --consequence stop_gained \
      --disease "Neurocardiorenal malformation syndrome" --sample-id SYN-001 --out out.vcf
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spike_pathogenic import (_CHROM_ORDER, _open, load_exome,  # noqa: E402
                              mc_consequence, parse_info, spiked_line)


def find_exact(clinvar_path, chrom, pos, ref, alt):
    """The ClinVar record at exactly chrom:pos:ref:alt, or None."""
    bare = str(chrom).replace("chr", "").replace("CHR", "")
    pos, ref, alt = str(pos), ref.upper(), alt.upper()
    with _open(clinvar_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[1] != pos:
                continue
            if f[0].replace("chr", "") != bare or f[3].upper() != ref or f[4].upper() != alt:
                continue
            info = parse_info(f[7])
            csq, _ = mc_consequence(info)
            return {"vid": f[2] or ".", "clnsig": info.get("CLNSIG", ""),
                    "clnrevstat": info.get("CLNREVSTAT", ""), "clndn": info.get("CLNDN", ""),
                    "csq": csq}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exome", required=True)
    ap.add_argument("--clinvar", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--pos", required=True, type=int)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--alt", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--consequence", default="")
    ap.add_argument("--disease", default="")
    ap.add_argument("--zygosity", default="het")
    ap.add_argument("--sample-id", default="SYN-001")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cv = find_exact(a.clinvar, a.chrom, a.pos, a.ref, a.alt) or {}

    def _clean(x):
        # A newline or tab in an INFO value breaks the VCF record into pieces / shifts its columns.
        # A ClinVar CLNDN can carry a stray one, so strip it. (Spaces are left as-is — CLNREVSTAT
        # already uses them and the engine parses them.)
        return str(x).replace("\n", "").replace("\r", "").replace("\t", " ") if x else x

    rec = {
        "chrom": str(a.chrom).replace("chr", ""), "pos": a.pos, "ref": a.ref, "alt": a.alt,
        "gene": a.gene, "rs": ".", "vid": cv.get("vid") or ".",
        "csq": _clean(cv.get("csq") or a.consequence or "missense_variant"),
        # exact ClinVar record when the coord matches this release; else a synthetic Pathogenic label
        "clnsig": _clean(cv.get("clnsig") or "Pathogenic"),
        "clnrevstat": cv.get("clnrevstat") or ("criteria provided, single submitter" if not cv else ""),
        "clndn": _clean(cv.get("clndn") or a.disease or ""),
    }
    src = "exact ClinVar match" if cv else "coord planted; synthetic CLNSIG (not in this ClinVar)"
    print(f"  spike {a.gene} {a.chrom}:{a.pos} {a.ref}>{a.alt} [{rec['csq']}, {rec['clnsig']}] — {src}",
          file=sys.stderr)

    meta, col_line, records, style = load_exome(a.exome)
    cols = col_line.split("\t")
    spike = spiked_line(rec, style, a.zygosity, a.sample_id, len(cols))

    all_rows = records + [spike]
    all_rows.sort(key=lambda f: (_CHROM_ORDER.get(f[0].replace("chr", ""), 99), int(f[1])))
    new_cols = cols[:9] + [a.sample_id] + ["." for _ in cols[10:]]

    extra_meta = [
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol (spiked)">',
        '##INFO=<ID=CSQ,Number=1,Type=String,Description="Molecular consequence (spiked)">',
        '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar significance (spiked)">',
        '##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description="ClinVar review status (spiked)">',
        '##INFO=<ID=CLNDN,Number=1,Type=String,Description="ClinVar disease name (spiked)">',
        '##INFO=<ID=CLNVID,Number=1,Type=String,Description="ClinVar Variation ID (spiked)">',
        '##INFO=<ID=SPIKED,Number=0,Type=Flag,Description="Synthetically spiked pathogenic variant">',
        "##comment=SYNTHETIC: real 1000G background with ONE exact ClinVar/phenopacket pathogenic "
        "variant spiked in. De-identified. Not real patient data. Not for clinical use.",
    ]
    have = set(meta)
    meta_out = meta + [m for m in extra_meta if m not in have]

    with open(a.out, "w") as out:
        out.write("\n".join(meta_out) + "\n")
        out.write("\t".join(new_cols) + "\n")
        for f in all_rows:
            out.write("\t".join(f) + "\n")
    print(f"Wrote {a.out}: {len(records)} background + 1 spiked = {len(all_rows)} variants "
          f"(sample '{a.sample_id}')", file=sys.stderr)


if __name__ == "__main__":
    main()
