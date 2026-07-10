#!/usr/bin/env python3
"""Convert a GA4GH Phenopacket into vcf2report inputs (a VCF + an HPO file).

    python scripts/phenopacket_to_inputs.py case.json --out-prefix data/out/case

Produces <prefix>.vcf and <prefix>.hpo.txt. The VCF carries coordinates + gene +
HGVS but NOT molecular consequence — annotate it (docs/ANNOTATION.md) before
classification unless your source already provides consequence in INFO.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report.phenopacket import load_phenopacket, write_inputs  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("phenopacket", help="Path to a GA4GH Phenopacket v2 JSON.")
    ap.add_argument("--out-prefix", default=None,
                    help="Output prefix (default: alongside the input).")
    args = ap.parse_args(argv)

    data = load_phenopacket(args.phenopacket)
    prefix = Path(args.out_prefix) if args.out_prefix else Path(args.phenopacket).with_suffix("")
    vcf, hpo = f"{prefix}.vcf", f"{prefix}.hpo.txt"
    write_inputs(data, vcf, hpo)
    print(f"subject: {data['subject_id']}")
    print(f"HPO terms ({len(data['hpo_terms'])}): {', '.join(data['hpo_terms'])}")
    print(f"variants: {len(data['variants'])}")
    if data.get("skipped_variants"):
        print(f"WARNING: {data['skipped_variants']} variant(s) skipped (HGVS-only, "
              f"no VCF coordinates) — resolve HGVS to coordinates to include them.")
    print(f"wrote {vcf} and {hpo}")
    print("NOTE: this VCF has no molecular consequence — annotate it "
          "(docs/ANNOTATION.md) before classification.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
