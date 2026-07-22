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


_ANNOT_KEYS = {"ANN", "EFF", "LOF", "NMD"}


def _strip_annot(info: str) -> str:
    """Drop position-specific functional annotation from a borrowed INFO — SnpEff re-derives it
    for the spike's real coordinate. Keeps the caller's quality stats (AC/AF/DP/MQ/QD/FS/SOR/…)."""
    kept = [kv for kv in info.split(";") if kv.split("=", 1)[0] not in _ANNOT_KEYS]
    return ";".join(kept) or "."


def _pick_template(records, zyg):
    """A REAL background call of the target zygosity whose INFO/FORMAT the spike will borrow, so the
    planted record is statistically indistinguishable from a genuine call (a different template per
    spike → values vary). Prefers a PASS SNV with the full DRAGEN FORMAT nearest 44x depth.

    The borrowed call must ITSELF pass the engine's QC comfortably (GQ>=30, DP>=25, well-balanced het) —
    otherwise the planted variant inherits a QC-failing GQ/AB and is silently dropped before
    classification (the original minimal spike hard-coded GQ=99, so this never mattered before)."""
    best, best_d = None, 1e18
    for f in records:
        if len(f) < 10 or len(f[3]) != 1 or len(f[4]) != 1 or f[6] != "PASS":
            continue
        fmt = f[8].split(":")
        if "GT" not in fmt or "F1R2" not in fmt:        # require a full DRAGEN FORMAT template
            continue
        smp = f[9].split(":")
        gt = smp[fmt.index("GT")].replace("|", "/")
        z = "hom" if gt == "1/1" else "het" if gt in ("0/1", "1/0") else None
        if z != zyg:
            continue
        try:
            gq = int(smp[fmt.index("GQ")]) if "GQ" in fmt else 0
            dp = int(smp[fmt.index("DP")]) if "DP" in fmt else 0
        except (ValueError, IndexError):
            continue
        if gq < 30 or dp < 25:                           # comfortably above QC (GQ>=20, DP>=10)
            continue
        if zyg == "het" and "AD" in fmt:                 # a well-balanced het (QC AB is 0.25-0.75)
            try:
                ad = [int(x) for x in smp[fmt.index("AD")].split(",")]
                if not (0.35 <= ad[1] / max(1, sum(ad)) <= 0.65):
                    continue
            except (ValueError, IndexError):
                continue
        d = abs(dp - 44)
        if d < best_d:
            best, best_d = f, d
    return best


def spiked_line_realistic(rec, records, style, zyg, col_count):
    """A tell-free spiked record: a real background call of the same zygosity, relocated to the
    spike's coordinate/alleles, annotation stripped (SnpEff re-adds it). No GENE/CSQ/CLN*/SPIKED —
    truth is tracked externally by coordinate (cohort.tsv / truth.tsv). Requires later SnpEff
    annotation for the engine to derive consequence. Returns None if no suitable template exists."""
    tmpl = _pick_template(records, zyg)
    if tmpl is None:
        return None
    chrom = ("chr" + rec["chrom"]) if style == "chr" else rec["chrom"]
    fixed = [chrom, str(rec["pos"]), ".", rec["ref"], rec["alt"],
             tmpl[5], "PASS", _strip_annot(tmpl[7]), tmpl[8], tmpl[9]]
    while len(fixed) < col_count:
        fixed.append(".")
    return fixed


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
    ap.add_argument("--realistic", action="store_true",
                    help="Plant a tell-free record indistinguishable from a real call (borrows a "
                         "background call's INFO/FORMAT; no GENE/CSQ/CLN*/SPIKED). Needs SnpEff "
                         "annotation before classification; truth is the external cohort/truth TSV.")
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
    used_realistic = False
    spike = None
    if a.realistic:
        spike = spiked_line_realistic(rec, records, style, a.zygosity, len(cols))
        if spike is None:
            print("  WARN: no full-FORMAT background call of matching zygosity to borrow — "
                  "falling back to the minimal (tell-bearing) record", file=sys.stderr)
        else:
            used_realistic = True
            print("  realistic: borrowed a real background call's INFO/FORMAT; no spike markers "
                  "(truth = external cohort/truth TSV; annotate with SnpEff before classifying)",
                  file=sys.stderr)
    if spike is None:
        spike = spiked_line(rec, style, a.zygosity, a.sample_id, len(cols))

    all_rows = records + [spike]
    all_rows.sort(key=lambda f: (_CHROM_ORDER.get(f[0].replace("chr", ""), 99), int(f[1])))
    new_cols = cols[:9] + [a.sample_id] + ["." for _ in cols[10:]]

    comment = ("##comment=SYNTHETIC: real 1000G background with ONE exact ClinVar/phenopacket "
               "pathogenic variant spiked in. De-identified. Not real patient data. Not for clinical use.")
    # Realistic records carry no per-variant markers, so they need no extra INFO definitions — only
    # the file-level provenance comment. The minimal record declares its GENE/CSQ/CLN*/SPIKED tags.
    if used_realistic:
        extra_meta = [comment]
    else:
        extra_meta = [
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol (spiked)">',
            '##INFO=<ID=CSQ,Number=1,Type=String,Description="Molecular consequence (spiked)">',
            '##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar significance (spiked)">',
            '##INFO=<ID=CLNREVSTAT,Number=1,Type=String,Description="ClinVar review status (spiked)">',
            '##INFO=<ID=CLNDN,Number=1,Type=String,Description="ClinVar disease name (spiked)">',
            '##INFO=<ID=CLNVID,Number=1,Type=String,Description="ClinVar Variation ID (spiked)">',
            '##INFO=<ID=SPIKED,Number=0,Type=Flag,Description="Synthetically spiked pathogenic variant">',
            comment,
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
