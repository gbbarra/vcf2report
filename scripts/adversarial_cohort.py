#!/usr/bin/env python3
"""Adversarial / anti-circularity analysis of the SYN cohort.

For each case, runs the engine twice: with the case's REAL HPO and with a DECOY HPO (another
case's phenotype, permuted +50). Per run it records how the planted gene surfaces:
  - engine_plp     : the gene reached Pathogenic/Likely-Pathogenic by the ACMG engine
  - pheno_primary  : the gene routed to the phenotype-matched PRIMARY findings (P/LP)
  - clinvar_flag   : surfaced only via the >=2-star ClinVar safety flag (phenotype-independent)
  - n_plp          : total P/LP calls in the ~100k-variant background (over-call proxy)
The honest, non-circular metric is ENGINE-ONLY (engine_plp OR pheno_primary — excludes the
ClinVar read-back). The DECOY pheno-primary rate measures phenotype specificity (should be low).

Run it on the ANNOTATED cohort (--annotated) for the honest numbers. On the raw cohort the
~100k-variant background carries no consequence, so the engine CANNOT call it P/LP whatever it
does — a low over-call there is an artifact of the missing annotation, not evidence of
specificity. Only the annotated run measures specificity for real.
"""
import argparse
import csv
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "src")
from vcf2report import pipeline  # noqa: E402
from vcf2report.report.assemble import clinvar_pathogenic_flags, split_findings  # noqa: E402

DIR = Path("data/synthetic_cohort")
_PLP = {"Pathogenic", "Likely Pathogenic"}
ANNOTATED = False


def _analyze(vcf, hpo, gene):
    rep = pipeline.run_pipeline(str(vcf), hpo_terms=hpo)
    cls = rep.classifications
    engine_plp = any(c.variant.gene == gene and c.tier in _PLP for c in cls)
    primary, _sec, _other = split_findings(cls)
    pheno_primary = any(c.variant.gene == gene and c.tier in _PLP for c in primary)
    clinvar_flag = any(c.variant.gene == gene for c in clinvar_pathogenic_flags(cls))
    n_plp = sum(1 for c in cls if c.tier in _PLP)
    return dict(engine_plp=engine_plp, pheno_primary=pheno_primary, clinvar_flag=clinvar_flag,
                n_plp=n_plp, surfaced=engine_plp or pheno_primary or clinvar_flag,
                engine_only=engine_plp or pheno_primary)


def _vcf_for(syn):
    return (DIR / "annotated" / f"{syn}.annotated.vcf.gz") if ANNOTATED \
        else (DIR / f"{syn}.synthetic.vcf.gz")


def _one(args):
    syn, sample, gene, real_hpo, decoy_hpo, annotated = args
    global ANNOTATED
    ANNOTATED = annotated
    vcf = _vcf_for(syn)
    try:
        r = _analyze(vcf, real_hpo, gene)
        d = _analyze(vcf, decoy_hpo, gene)
    except Exception as exc:
        return {"syn": syn, "gene": gene, "error": str(exc)[:120]}
    return {"syn": syn, "sample": sample, "gene": gene,
            "surfaced": r["surfaced"], "engine_only": r["engine_only"],
            "engine_plp": r["engine_plp"], "pheno_primary": r["pheno_primary"],
            "clinvar_flag": r["clinvar_flag"], "n_plp": r["n_plp"],
            "decoy_pheno_primary": d["pheno_primary"], "decoy_surfaced": d["surfaced"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotated", action="store_true",
                    help="analyze data/synthetic_cohort/annotated/*.annotated.vcf.gz (SnpEff MANE) "
                         "instead of the raw calls — required for a real specificity number")
    ap.add_argument("--limit", type=int, default=0, help="only the first N cases (smoke test)")
    ap.add_argument("--jobs", type=int, default=2, help="parallel workers (default 2)")
    args = ap.parse_args()
    global ANNOTATED
    ANNOTATED = args.annotated

    rows = list(csv.DictReader(open(DIR / "cohort.tsv"), delimiter="\t"))
    hpos = [r["hpo"].split(",") for r in rows]
    tasks = []
    for i, r in enumerate(rows):
        decoy = hpos[(i + 50) % len(hpos)]           # a real but MISMATCHED phenotype
        tasks.append((r["syn_id"], r["sample"], r["gene"], hpos[i], decoy, ANNOTATED))
    if args.limit:
        tasks = tasks[:args.limit]

    missing = [t[0] for t in tasks if not _vcf_for(t[0]).exists()]
    if missing:
        sys.exit(f"ERROR: {len(missing)} VCF(s) not found (e.g. {_vcf_for(missing[0])}). "
                 + ("Run: bash scripts/annotate_syn_cohort.sh" if ANNOTATED else ""))
    label = ("ANNOTATED (SnpEff MANE)" if ANNOTATED
             else "RAW calls (no consequence — specificity NOT measurable)")
    print(f"cohort: {label} · {len(tasks)} cases · {args.jobs} workers")

    with Pool(args.jobs) as pool:                     # conservative for a laptop's RAM
        results = []
        for k, res in enumerate(pool.imap_unordered(_one, tasks), 1):
            results.append(res)
            print(f"  [{k}/{len(tasks)}] {res['syn']} {res['gene']} "
                  f"{'ok' if 'error' not in res else 'ERR'}", flush=True)

    out = DIR / ("adversarial_results_annotated.tsv" if ANNOTATED else "adversarial_results.tsv")
    cols = ["syn", "sample", "gene", "surfaced", "engine_only", "engine_plp", "pheno_primary",
            "clinvar_flag", "n_plp", "decoy_pheno_primary", "decoy_surfaced"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow(r)

    good = [r for r in results if "error" not in r]
    n = len(good)

    def pct(key):
        return f"{sum(bool(r[key]) for r in good)}/{n} ({100*sum(bool(r[key]) for r in good)//max(1,n)}%)"

    nplp = sorted(r["n_plp"] for r in good)
    print("\n================ ADVERSARIAL SUMMARY ================")
    print(f"cases analyzed: {n}/100")
    print(f"surfaced (any path):         {pct('surfaced')}")
    print(f"ENGINE-ONLY (no ClinVar flag): {pct('engine_only')}   <- anti-circular metric")
    print(f"  via engine P/LP tier:      {pct('engine_plp')}")
    print(f"  via phenotype-primary:     {pct('pheno_primary')}")
    print(f"surfaced ONLY by ClinVar flag: "
          f"{sum(r['clinvar_flag'] and not r['engine_only'] for r in good)}/{n}")
    print(f"DECOY pheno-primary (wrong HPO still routes gene): {pct('decoy_pheno_primary')}   "
          f"<- should be LOW (phenotype specificity)")
    print(f"over-call — P/LP calls per case (~100k background): "
          f"median {nplp[n//2]} / max {nplp[-1]} / min {nplp[0]}")
    if not ANNOTATED:
        print("  ^ ARTIFACT: the raw background has no consequence, so it CANNOT be called P/LP.\n"
              "    This number measures the missing annotation, not the engine. Re-run --annotated.")
    print(f"results -> {out}")


if __name__ == "__main__":
    main()
