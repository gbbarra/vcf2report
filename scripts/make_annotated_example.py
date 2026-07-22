#!/usr/bin/env python3
"""Turn a benchmark SnpEff-annotated exome into a small, SELF-CONTAINED annotated example VCF.

The hpo-spiked-exomes ``realistic_annotated/`` VCFs already carry the SnpEff ``ANN`` (gene /
consequence / HGVS), but are otherwise tell-free — no gnomAD/ClinVar in the INFO. This bakes the
gnomAD, ClinVar and AlphaMissense values from the LOCAL stores into the INFO (the exact fields the
engine reads back via its ``from_vcf`` path), and slims the ~100k-variant exome to a manageable size,
so the result runs through ``run_headless.py`` and produces a full laudo with **no stores needed** —
ideal as a committed repo example.

Mirrors the pipeline's batch annotation: ``gnomad_parquet.prime`` + ``clinvar_parquet.prime`` once,
then per-variant reads (no per-variant DuckDB query). Run where the stores are present
(``check_stores.py`` green). The planted variant (``--keep``) is always retained.

  python3 scripts/make_annotated_example.py \
      --vcf  ~/hpo-spiked-exomes/realistic_annotated/SYN-004.annotated.vcf.gz \
      --keep chr5:37022325:C:T \
      --out  data/example/SYN-004.NIPBL.annotated.vcf.gz \
      --max-background 8000
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report.models import Variant  # noqa: E402

# INFO fields we add (names must match config.INFO_ALIASES so the engine reads them back).
_HEADER_LINES = [
    '##INFO=<ID=gnomad_AF,Number=A,Type=Float,Description="gnomAD v4.1 grpmax allele frequency (baked from local store)">',
    '##INFO=<ID=gnomad_AC,Number=A,Type=Integer,Description="gnomAD v4.1 allele count">',
    '##INFO=<ID=gnomad_AN,Number=1,Type=Integer,Description="gnomAD v4.1 allele number">',
    '##INFO=<ID=gnomad_nhomalt,Number=A,Type=Integer,Description="gnomAD v4.1 homozygote count">',
    '##INFO=<ID=gnomad_faf95,Number=A,Type=Float,Description="gnomAD v4.1 grpmax filtering AF (faf95)">',
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance (baked from local store)">',
    '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">',
    '##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar disease name">',
    '##INFO=<ID=CLNVI,Number=.,Type=String,Description="ClinVar accession">',
    '##INFO=<ID=am_pathogenicity,Number=A,Type=Float,Description="AlphaMissense pathogenicity (baked from local store)">',
]


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def _norm_chrom(c: str) -> str:
    return c if c.lower().startswith("chr") else "chr" + c


def _ann_gene(info: str) -> str | None:
    for f in info.split(";"):
        if f.startswith("ANN="):
            first = f[4:].split(",")[0].split("|")
            return first[3] if len(first) > 3 else None
    return None


def _esc(v) -> str:
    return str(v).replace(" ", "_").replace(";", "_").replace("=", "_")


def _baked_info(ann) -> str:
    """The INFO snippet to prepend, from a populated Annotation (only present fields)."""
    parts = []
    if ann.gnomad_af is not None:
        parts.append(f"gnomad_AF={ann.gnomad_af:g}")
    if ann.gnomad_ac is not None:
        parts.append(f"gnomad_AC={ann.gnomad_ac}")
    if ann.gnomad_an is not None:
        parts.append(f"gnomad_AN={ann.gnomad_an}")
    if ann.gnomad_homozygotes is not None:
        parts.append(f"gnomad_nhomalt={ann.gnomad_homozygotes}")
    if ann.gnomad_faf95 is not None:
        parts.append(f"gnomad_faf95={ann.gnomad_faf95:g}")
    if ann.clinvar_significance:
        parts.append(f"CLNSIG={_esc(ann.clinvar_significance)}")
        if ann.clinvar_review_status:
            parts.append(f"CLNREVSTAT={_esc(ann.clinvar_review_status)}")
        if ann.clinvar_condition:
            parts.append(f"CLNDN={_esc(ann.clinvar_condition)}")
        if ann.clinvar_accession:
            parts.append(f"CLNVI={_esc(ann.clinvar_accession)}")
    if ann.am_pathogenicity is not None:
        parts.append(f"am_pathogenicity={ann.am_pathogenicity:g}")
    return ";".join(parts)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Bake store annotations into a SnpEff VCF -> example.")
    ap.add_argument("--vcf", required=True, help="SnpEff-annotated (tell-free) VCF")
    ap.add_argument("--out", required=True, help="output annotated example VCF (.vcf.gz)")
    ap.add_argument("--keep", required=True, help="planted variant chrom:pos:ref:alt (always kept)")
    ap.add_argument("--max-background", type=int, default=8000,
                    help="downsample the background to ~this many variants (default 8000)")
    args = ap.parse_args(argv)

    from vcf2report.annotate import annotate_variant, add_alphamissense, clinvar_parquet, gnomad_parquet

    kc, kp, kr, ka = args.keep.split(":")
    keep_key = (_norm_chrom(kc), kp, kr, ka)

    header, data = [], []
    with _open(args.vcf) as fh:
        for line in fh:
            (header if line.startswith("#") else data).append(line.rstrip("\n"))

    # slim: always keep the plant; stride-sample the rest to ~max-background across the genome
    stride = max(1, len(data) // max(1, args.max_background))
    kept_lines, plant_found = [], False
    for i, l in enumerate(data):
        f = l.split("\t")
        key = (_norm_chrom(f[0]), f[1], f[3], f[4].split(",")[0])
        if key == keep_key:
            kept_lines.append(l); plant_found = True
        elif i % stride == 0:
            kept_lines.append(l)
    if not plant_found:
        print(f"WARNING: planted variant {args.keep} not found in {args.vcf}", file=sys.stderr)

    # build Variants for the batch store lookup
    variants = []
    for l in kept_lines:
        f = l.split("\t")
        variants.append(Variant(chrom=_norm_chrom(f[0]), pos=int(f[1]), ref=f[3],
                                 alt=f[4].split(",")[0], gene=_ann_gene(f[7])))

    # mirror the pipeline: prime gnomAD + ClinVar once (batch), then per-variant read
    primed = gnomad_parquet.prime(variants)
    clinvar_parquet.prime(variants)
    print(f"primed gnomAD for {primed}/{len(variants)} variants", file=sys.stderr)
    if primed == 0:
        print("WARNING: 0 gnomAD hits — are the stores present? (check_stores.py)", file=sys.stderr)
    anns = [annotate_variant(v, [], with_alphamissense=False) for v in variants]
    for v, a in zip(variants, anns):
        add_alphamissense(v, a)  # AlphaMissense for all kept (per-variant tabix/parquet)

    # rewrite each kept line's INFO with the baked fields prepended
    out_data = []
    baked = 0
    for l, ann in zip(kept_lines, anns):
        f = l.split("\t")
        snippet = _baked_info(ann)
        if snippet:
            f[7] = snippet + ";" + f[7] if f[7] not in (".", "") else snippet
            baked += 1
        out_data.append("\t".join(f))

    # header: inject the INFO declarations just before #CHROM
    out_header = []
    for h in header:
        if h.startswith("#CHROM"):
            out_header += [x for x in _HEADER_LINES if not any(x.split(",")[0] in hh for hh in header)]
        out_header.append(h)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt") as out:
        out.write("\n".join(out_header) + "\n")
        out.write("\n".join(out_data) + "\n")
    print(f"wrote {len(out_data)} variants ({baked} store-annotated, plant {'kept' if plant_found else 'MISSING'}) -> {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
