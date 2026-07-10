#!/usr/bin/env python3
"""Freeze a ground-truth concordance panel from REAL ClinVar + live gnomAD.

Selects a balanced, deterministic set of expert-reviewed ClinVar variants
(Pathogenic LoF, Pathogenic missense, Benign) and attaches each variant's REAL
live gnomAD frequency (grpmax AF + faf95). The frozen panel
(``data/truth/clinvar_panel.json``) is then classified OFFLINE by
``tests/test_ground_truth_concordance.py`` — with ClinVar withheld from the
engine (PP5 off) so the comparison is honest, not circular.

Input is the compact TSV distilled from the ClinVar VCF (see the header of the
concordance test / SYNTHETIC_CASES docs for the streaming filter command):

    allele_id  chrom  pos  ref  alt  gene  clnsig  clnrevstat  mc

Run (network required for the live gnomAD lookups):

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_truth_panel.py COMPACT.tsv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# ClinVar molecular-consequence (MC) SO terms -> the tokens vcf2report's engine
# understands (Variant.is_lof / PVS1). ClinVar says "nonsense"; we say "stop_gained".
MC_MAP = {
    "nonsense": "stop_gained",
    "frameshift_variant": "frameshift_variant",
    "splice_donor_variant": "splice_donor_variant",
    "splice_acceptor_variant": "splice_acceptor_variant",
    "missense_variant": "missense_variant",
    "synonymous_variant": "synonymous_variant",
    "stop_lost": "stop_lost",
    "initiator_codon_variant": "start_lost",
    "initiatior_codon_variant": "start_lost",  # ClinVar's spelling
}
# Severity order for picking one consequence out of a multi-term MC field.
SEVERITY = ["stop_gained", "splice_donor_variant", "splice_acceptor_variant",
            "frameshift_variant", "start_lost", "stop_lost", "missense_variant",
            "synonymous_variant"]

PATHOGENIC = {"Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"}
BENIGN = {"Benign", "Benign/Likely_benign", "Likely_benign"}
LOF = {"stop_gained", "splice_donor_variant", "splice_acceptor_variant",
       "frameshift_variant", "start_lost", "stop_lost"}

# Panel sizes per bucket, and the max variants kept per gene (diversity).
N_PATH_LOF, N_PATH_MIS, N_BENIGN = 30, 20, 50
PER_GENE_CAP = 2


def _consequence(mc: str) -> str | None:
    terms = []
    for part in mc.split(","):
        tok = part.split("|")[-1]
        if tok in MC_MAP:
            terms.append(MC_MAP[tok])
    for s in SEVERITY:
        if s in terms:
            return s
    return None


def _rows(tsv: Path):
    for line in tsv.read_text().splitlines():
        f = line.split("\t")
        if len(f) != 9:
            continue
        aid, chrom, pos, ref, alt, gene, clnsig, rev, mc = f
        con = _consequence(mc)
        if con is None or not gene:
            continue
        yield {"allele_id": int(aid), "chrom": chrom, "pos": int(pos),
               "ref": ref, "alt": alt, "gene": gene, "consequence": con,
               "clinvar_sig": clnsig, "clinvar_review": rev}


def _buckets(rows):
    """Partition rows into (pathogenic LoF, pathogenic missense, benign)."""
    path_lof, path_mis, benign = [], [], []
    for r in rows:
        if r["clinvar_sig"] in PATHOGENIC:
            (path_lof if r["consequence"] in LOF else path_mis).append(r)
        elif r["clinvar_sig"] in BENIGN:
            benign.append(r)
    return path_lof, path_mis, benign


class _LocalGnomad:
    """Query AF/faf95 from a LOCAL gnomAD sites VCF (bgzipped + tabixed).

    The offline, robust path — the same local file vcfanno uses (no flaky remote
    API). Reuses gnomad_remote's field extraction so grpmax/faf95 semantics match
    the live path exactly. Callable like ``gnomad_remote.query``.
    """

    def __init__(self, path):
        import pysam
        from vcf2report.annotate.gnomad_remote import _best_from_record
        self._vf = pysam.VariantFile(path)
        self._extract = _best_from_record
        # gnomAD v4.1 contigs are "chr"-prefixed; ClinVar rows are not. Detect.
        self._chr_prefixed = any(c.startswith("chr") for c in self._vf.header.contigs)

    def _contig(self, chrom):
        c = str(chrom)
        if self._chr_prefixed and not c.startswith("chr"):
            return f"chr{c}"
        if not self._chr_prefixed and c.startswith("chr"):
            return c[3:]
        return c

    def __call__(self, variant):
        chrom = self._contig(variant.chrom)
        best = None
        try:
            for rec in self._vf.fetch(chrom, variant.pos - 1, variant.pos):
                if rec.pos != variant.pos or rec.ref != variant.ref:
                    continue
                if variant.alt not in (rec.alts or ()):
                    continue
                cand = self._extract(rec)
                if cand and (best is None or (cand["af"] or 0) > (best["af"] or 0)):
                    best = cand
        except (ValueError, OSError):
            return None
        if best is not None:
            return best
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}


def _take(bucket, n):
    """First n rows by allele_id, at most PER_GENE_CAP per gene (deterministic)."""
    bucket.sort(key=lambda r: r["allele_id"])
    out, seen = [], {}
    for r in bucket:
        if seen.get(r["gene"], 0) >= PER_GENE_CAP:
            continue
        seen[r["gene"]] = seen.get(r["gene"], 0) + 1
        out.append(r)
        if len(out) >= n:
            break
    return out


def main() -> int:
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("compact_tsv", help="compact ClinVar TSV (see module docstring)")
    p.add_argument("--gnomad-vcf", default=None,
                   help="LOCAL gnomAD sites VCF (bgzipped + tabixed). Preferred: no "
                        "network, no flaky remote API. Omit to use the live remote path.")
    p.add_argument("--out", default="data/truth/clinvar_panel.json")
    args = p.parse_args()

    from vcf2report.models import Variant

    if args.gnomad_vcf:
        gq = _LocalGnomad(args.gnomad_vcf)
        src_label = f"gnomAD local sites VCF ({Path(args.gnomad_vcf).name})"
        print(f"reading gnomAD from LOCAL {args.gnomad_vcf}", file=sys.stderr)
    else:
        from vcf2report import config
        if config.offline():
            raise SystemExit("no --gnomad-vcf and network is off: pass a local gnomAD "
                             "sites VCF (--gnomad-vcf) or set VCF2REPORT_ALLOW_NETWORK=1")
        from vcf2report.annotate import gnomad_remote
        gq = gnomad_remote.query
        src_label = "gnomAD v4.1 live grpmax (remote tabix)"
        print("reading gnomAD from LIVE remote tabix", file=sys.stderr)

    tsv = Path(args.compact_tsv)
    path_lof, path_mis, benign = _buckets(list(_rows(tsv)))
    selected = (_take(path_lof, N_PATH_LOF) + _take(path_mis, N_PATH_MIS)
                + _take(benign, N_BENIGN))
    print(f"selected {len(selected)} candidates; querying gnomAD ...", file=sys.stderr)

    panel = []
    for i, r in enumerate(selected, 1):
        q = gq(Variant(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"]))
        if q is None:
            print(f"  [{i}] {r['gene']} {r['chrom']}:{r['pos']} gnomAD lookup FAILED — skip",
                  file=sys.stderr)
            continue
        r["gnomad_af"] = q.get("af")
        r["gnomad_faf95"] = q.get("faf95")
        panel.append(r)
        print(f"  [{i}/{len(selected)}] {r['gene']:8s} {r['clinvar_sig']:24s} "
              f"{r['consequence']:22s} gnomAD AF={r['gnomad_af']}", file=sys.stderr)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "_source": f"ClinVar 2023-07-17 (GRCh38, expert-reviewed SNVs) + {src_label}",
        "_note": "Ground-truth concordance panel. De-identified public data. Not for clinical use.",
        "_engine_eval": "classified WITH ClinVar withheld (PP5 off) to avoid circularity",
        "count": len(panel),
    }
    out.write_text(json.dumps({"meta": meta, "variants": panel}, indent=2) + "\n")
    print(f"wrote {out} with {len(panel)} variants", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
