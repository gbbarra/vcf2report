#!/usr/bin/env python3
"""Score vcf2report against the hpo-spiked-exomes benchmark (200 tell-free spiked exomes).

For each case, run the pipeline on the SnpEff-annotated VCF + the case's HPO terms, and check whether
the planted gene (the answer key) is recovered as a **primary** (diagnostic) finding. Reports the
headline ``X/200`` primary recovery plus an honest breakdown of the misses by where the planted gene
landed instead — secondary / carrier / probable-pathogenic VUS / other-tier / absent — so a
regression shows up as a case moving buckets, not just the count changing.

The dataset: https://github.com/gbbarra/hpo-spiked-exomes (fetch its release with its ``fetch.sh``).
The answer key (``manifest/cohort.tsv``) and per-case HPO (``sidecars/SYN-NNN.hpo.txt``) live in that repo.

⚠️ Requires the full data stores (gnomAD + ClinVar + AlphaMissense parquet). The tell-free VCFs carry
NO baked frequency/ClinVar INFO, so the engine must look them up — run this where the stores live
(``scripts/check_stores.py`` should be green). Without them PM2/PP5/PP3 silently disable and the score
is not comparable to the published 178/200.

  python3 scripts/run_benchmark.py \
      --annotated /path/hpo-spiked-exomes/realistic_annotated \
      --bench     /path/hpo-spiked-exomes \
      --out       benchmark_results.tsv [--jobs 4] [--limit N]

Before/after an ACMG change (e.g. PP2/BP1/PS1/PM5): score the OLD tree, then the NEW one with
``--compare`` pointed at the old TSV. It prints net PRIMARY recovery before→after (overall and for
the missense subset those criteria move) and lists every case whose bucket/tier changed (★ = missense):

  git checkout main         && python3 scripts/run_benchmark.py ... --out before.tsv
  git checkout acmg-branch  && python3 scripts/run_benchmark.py ... --out after.tsv --compare before.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _load_answer_key(bench: Path) -> dict[str, str]:
    """syn_id -> planted gene, from the benchmark's cohort.tsv answer key."""
    coh = bench / "manifest" / "cohort.tsv"
    with open(coh, newline="") as fh:
        return {r["syn_id"]: r["gene"] for r in csv.DictReader(fh, delimiter="\t")}


def _bucket_of(gene: str, report) -> tuple[str, str | None]:
    """Where did the planted gene land? Returns (bucket, tier-of-that-gene-or-None)."""
    from vcf2report.report.assemble import carrier_findings, split_findings
    from vcf2report.report.vus_triage import probable_pathogenic_vus
    primary, secondary, other = split_findings(report.classifications)
    carriers = carrier_findings(report.classifications)
    vus = [e["classification"] for e in probable_pathogenic_vus(report.classifications)]
    hit = next((c for c in report.classifications if c.variant.gene == gene), None)
    tier = hit.tier if hit else None
    for name, members in (("primary", primary), ("secondary", secondary),
                          ("carrier", carriers), ("probable_vus", vus)):
        if any(c.variant.gene == gene for c in members):
            return name, tier
    return ("other" if hit else "absent"), tier


def _score_one(args: tuple[str, str, str, str]) -> dict:
    """Worker: run the pipeline for one case and locate the planted gene. Top-level for pickling."""
    sid, vcf, hpo_path, gene = args
    from vcf2report.cli import read_hpo_file
    from vcf2report.pipeline import run_pipeline
    try:
        hpo_terms = read_hpo_file(hpo_path) if Path(hpo_path).exists() else []
        report = run_pipeline(vcf, hpo_terms=hpo_terms, sample_id=sid)
        bucket, tier = _bucket_of(gene, report)
        hit = next((c for c in report.classifications if c.variant.gene == gene), None)
        return {"syn_id": sid, "gene": gene, "outcome": bucket, "tier": tier or "",
                "consequence": (hit.variant.consequence if hit else "") or "",
                "candidates": report.qc.candidates, "error": ""}
    except Exception as e:  # never let one case abort the sweep
        return {"syn_id": sid, "gene": gene, "outcome": "ERROR", "tier": "",
                "consequence": "", "candidates": 0, "error": f"{type(e).__name__}: {e}"}


def _is_missense(consequence: str | None) -> bool:
    return (consequence or "").startswith("missense")


def _compare(results: list[dict], prev_path: str) -> None:
    """Diff this run against a saved per-case TSV, spotlighting the planted variants whose
    outcome/tier moved — grouped so missense cases (the ones PP2/BP1/PS1/PM5 touch) are
    obvious. Prints net PRIMARY recovery before→after, overall and for missense only."""
    with open(prev_path, newline="") as fh:
        prev = {r["syn_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    common = [r for r in results if r["syn_id"] in prev]
    if not common:
        print(f"\n(compare: no overlapping syn_ids with {prev_path})")
        return

    def _primary(rows, key):
        return sum(1 for r in rows if key(r) == "primary")

    changed = []
    for r in common:
        p = prev[r["syn_id"]]
        if (p.get("outcome"), p.get("tier", "")) != (r["outcome"], r["tier"]):
            changed.append((r, p))

    n_before = _primary([prev[r["syn_id"]] for r in common], lambda r: r.get("outcome"))
    n_after = _primary(common, lambda r: r["outcome"])
    mis = [r for r in common if _is_missense(r.get("consequence"))]
    mis_before = _primary([prev[r["syn_id"]] for r in mis], lambda r: r.get("outcome"))
    mis_after = _primary(mis, lambda r: r["outcome"])

    print(f"\n=== compare vs {prev_path} ({len(common)} common cases) ===")
    print(f"PRIMARY recovery:  before {n_before}  →  after {n_after}  ({n_after - n_before:+d})")
    print(f"  of which missense: before {mis_before}  →  after {mis_after}  "
          f"({mis_after - mis_before:+d})  [{len(mis)} missense cases]")
    if not changed:
        print("no case changed outcome or tier.")
        return
    print(f"\n{len(changed)} case(s) changed (★ = missense):")
    for r, p in sorted(changed, key=lambda x: x[0]["syn_id"]):
        star = "★" if _is_missense(r.get("consequence")) else " "
        po = f"{p.get('outcome')}·{p.get('tier','') or '—'}"
        no = f"{r['outcome']}·{r['tier'] or '—'}"
        print(f"  {star} {r['syn_id']}  {r['gene']:10} {po:24} → {no}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Score vcf2report on the hpo-spiked-exomes benchmark.")
    ap.add_argument("--annotated", required=True, help="dir of SYN-NNN.annotated.vcf.gz")
    ap.add_argument("--bench", required=True, help="hpo-spiked-exomes repo root (manifest/ + sidecars/)")
    ap.add_argument("--out", default="benchmark_results.tsv", help="per-case TSV output")
    ap.add_argument("--jobs", type=int, default=1, help="parallel worker processes")
    ap.add_argument("--limit", type=int, default=0, help="score only the first N cases (debug)")
    ap.add_argument("--compare", default="", metavar="PREV.tsv",
                    help="diff this run against a previous --out TSV (before/after; flags missense moves)")
    args = ap.parse_args(argv)

    bench = Path(args.bench)
    ann = Path(args.annotated)
    key = _load_answer_key(bench)

    tasks: list[tuple[str, str, str, str]] = []
    missing = 0
    for sid, gene in sorted(key.items()):
        vcf = ann / f"{sid}.annotated.vcf.gz"
        if not vcf.exists():
            missing += 1
            continue
        hpo = bench / "sidecars" / f"{sid}.hpo.txt"
        tasks.append((sid, str(vcf), str(hpo), gene))
    if args.limit:
        tasks = tasks[: args.limit]
    if not tasks:
        print(f"No annotated VCFs found under {ann} matching the answer key.", file=sys.stderr)
        return 1
    if missing:
        print(f"note: {missing} answer-key cases have no annotated VCF under {ann} — skipped.",
              file=sys.stderr)

    results: list[dict] = []
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(_score_one, t) for t in tasks]
            for i, f in enumerate(as_completed(futs), 1):
                results.append(f.result())
                print(f"\r  scored {i}/{len(tasks)}", end="", file=sys.stderr, flush=True)
    else:
        for i, t in enumerate(tasks, 1):
            results.append(_score_one(t))
            print(f"\r  scored {i}/{len(tasks)}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)

    results.sort(key=lambda r: r["syn_id"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["syn_id", "gene", "outcome", "tier", "consequence",
                                       "candidates", "error"])
        w.writeheader()
        w.writerows(results)

    n = len(results)
    from collections import Counter
    by = Counter(r["outcome"] for r in results)
    primary = by.get("primary", 0)
    mis = [r for r in results if _is_missense(r.get("consequence"))]
    mis_primary = sum(1 for r in mis if r["outcome"] == "primary")
    print(f"\n=== hpo-spiked-exomes benchmark — vcf2report ===")
    print(f"PRIMARY (diagnostic) recovery: {primary}/{n}  ({100*primary/n:.1f}%)")
    if mis:
        print(f"  missense planted variants: {mis_primary}/{len(mis)} primary "
              f"({100*mis_primary/len(mis):.1f}%)  — the subset PP2/BP1/PS1/PM5 move")
    print("breakdown by where the planted gene landed:")
    for k in ("primary", "secondary", "carrier", "probable_vus", "other", "absent", "ERROR"):
        if by.get(k):
            print(f"  {k:13}: {by[k]}")
    misses = [r for r in results if r["outcome"] not in ("primary",)]
    if misses:
        print("\nnon-primary cases (gene → bucket · tier):")
        for r in sorted(misses, key=lambda r: (r["outcome"], r["syn_id"])):
            extra = f" · {r['error']}" if r["error"] else ""
            print(f"  {r['syn_id']}  {r['gene']:10} → {r['outcome']}"
                  f"{(' · ' + r['tier']) if r['tier'] else ''}{extra}")
    if args.compare:
        _compare(results, args.compare)
    print(f"\nper-case TSV: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
