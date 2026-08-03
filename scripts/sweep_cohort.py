#!/usr/bin/env python3
"""Dump EVERY classification the engine makes across a cohort — not just the planted gene.

``run_benchmark.py`` answers "was the answer found?". It says nothing about what else came out
of the same run, and that is the half a laboratory actually feels: the report a clinician opens
contains the planted variant *and* everything else the engine chose to call, and each of those
costs review time or, worse, an unwarranted result.

This sweep writes one row per classification (200 cases x ~1.2k candidates) so two questions can
be answered from one pass, offline, without re-running the pipeline per question:

**1. Over-call / review burden.** How many Pathogenic / Likely Pathogenic calls does a case
produce *besides* the planted one, and in which bucket do they land? These are not all false
positives — the backgrounds are real 1000G exomes and genuinely carry pathogenic alleles
(secondary findings, carrier status), which is why the output keeps the bucket and the gene
rather than reporting a bare "false positive rate" that would be a lie in either direction.

**2. ClinVar disagreement at scale.** The engine classifies independently. Where ClinVar has an
assertion, the two can be compared over every candidate rather than over the 200 planted loci —
and the interesting cell is not agreement but *gross* discordance: engine P/LP where ClinVar says
benign, or engine benign where ClinVar says pathogenic with review stars behind it.

    python3 scripts/sweep_cohort.py --annotated <bench>/realistic_annotated \\
        --bench <bench> --out sweep.tsv [--jobs 4] [--limit N]

Needs the full stores, like ``run_benchmark.py`` — a sweep on absent stores measures the gate,
not the engine.
"""
from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

FIELDS = ["syn_id", "chrom", "pos", "ref", "alt", "gene", "tier", "consequence", "bucket",
          "is_planted", "clinvar_significance", "clinvar_stars", "zygosity", "error"]


def _load_answer_key(bench: Path) -> dict[str, str]:
    with open(bench / "manifest" / "cohort.tsv", newline="") as fh:
        return {r["syn_id"]: r["gene"] for r in csv.DictReader(fh, delimiter="\t")}


def _buckets(report) -> dict[int, str]:
    """id(classification) -> bucket name, using the same split the report itself uses.

    Keyed on identity rather than on gene: a case can classify two variants in one gene and they
    do not have to land in the same bucket, and collapsing them would silently drop one.
    """
    from vcf2report.report.assemble import carrier_findings, split_findings
    from vcf2report.report.vus_triage import probable_pathogenic_vus
    primary, secondary, other = split_findings(report.classifications)
    out: dict[int, str] = {}
    for name, members in (("other", other), ("secondary", secondary),
                          ("carrier", carrier_findings(report.classifications)),
                          ("probable_vus", [e["classification"]
                                            for e in probable_pathogenic_vus(report.classifications)]),
                          ("primary", primary)):
        for c in members:
            out[id(c)] = name          # later entries win; primary is the most specific
    return out


def _sweep_one(args: tuple) -> list[dict]:
    sid, vcf, hpo_path, planted = args
    from vcf2report.cli import read_hpo_file
    from vcf2report.pipeline import run_pipeline
    # Stars are DERIVED from the review status, not stored on Annotation — use the same
    # helper the report does, so "2 stars" here means what it means in the laudo.
    from vcf2report.report.assemble import clinvar_stars
    try:
        hpo = read_hpo_file(hpo_path) if Path(hpo_path).exists() else []
        report = run_pipeline(vcf, hpo_terms=hpo, sample_id=sid)
        buckets = _buckets(report)
        rows = []
        for c in report.classifications:
            a, v = c.annotation, c.variant
            rows.append({
                "syn_id": sid,
                "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt,
                "gene": v.gene or "",
                "tier": c.tier or "",
                "consequence": (v.consequence or "").split("&")[0],
                "bucket": buckets.get(id(c), "unranked"),
                "is_planted": "1" if v.gene == planted else "0",
                "clinvar_significance": a.clinvar_significance or "",
                "clinvar_stars": clinvar_stars(a.clinvar_review_status)
                                 if a.clinvar_review_status else "",
                "zygosity": v.zygosity or "",
                "error": "",
            })
        return rows
    except Exception as e:                       # one bad case must not abort the sweep
        return [{**{k: "" for k in FIELDS}, "syn_id": sid, "gene": planted,
                 "error": f"{type(e).__name__}: {e}"}]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Dump every classification across the cohort.")
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--bench", required=True)
    ap.add_argument("--out", default="sweep.tsv")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)

    bench, ann = Path(a.bench), Path(a.annotated)
    key = _load_answer_key(bench)
    tasks = [(sid, str(ann / f"{sid}.annotated.vcf.gz"),
              str(bench / "sidecars" / f"{sid}.hpo.txt"), gene)
             for sid, gene in sorted(key.items())
             if (ann / f"{sid}.annotated.vcf.gz").exists()]
    if a.limit:
        tasks = tasks[: a.limit]
    if not tasks:
        print(f"No annotated VCFs under {ann}.", file=sys.stderr)
        return 1

    rows: list[dict] = []
    if a.jobs > 1:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            futs = [ex.submit(_sweep_one, t) for t in tasks]
            for i, f in enumerate(as_completed(futs), 1):
                rows.extend(f.result())
                print(f"\r  swept {i}/{len(tasks)}", end="", file=sys.stderr, flush=True)
    else:
        for i, t in enumerate(tasks, 1):
            rows.extend(_sweep_one(t))
            print(f"\r  swept {i}/{len(tasks)}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)

    rows.sort(key=lambda r: (r["syn_id"], r["gene"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    errs = [r for r in rows if r["error"]]
    print(f"\n{len(rows):,} classifications across {len(tasks)} cases -> {a.out}")
    if errs:
        print(f"⚠️  {len(errs)} case(s) errored: {[r['syn_id'] for r in errs]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
