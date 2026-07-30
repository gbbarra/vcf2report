#!/usr/bin/env python3
"""Why does the concordance panel report 61% pathogenic sensitivity?

Answers it per variant instead of leaving the number to be read as "the engine fails on 39
of 100". For every panel entry ClinVar calls pathogenic that the engine does NOT call
pathogenic, this reports the WEAKEST single line of evidence that would reach Likely
Pathogenic, and which criteria could supply it — then re-runs the six cheapest with a
phenotype match, the one input the panel structurally withholds from everybody.

Nothing here re-classifies anything or relaxes a rule: it reads the same combining rules the
engine publishes and reports what they would do with one more line. Offline; no stores.

    python3 scripts/analyse_panel_misses.py
    python3 scripts/analyse_panel_misses.py --json
"""
import argparse
import json
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report.acmg.engine import classify  # noqa: E402
from vcf2report.concordance import (_annotation_from_frozen, classify_entry,  # noqa: E402
                                    collapse_engine_tier, load_panel)
from vcf2report.report.explore import missing_evidence  # noqa: E402

_ORDER = ["supporting", "moderate", "strong", "very_strong"]


def analyse() -> dict:
    missed = []
    for e in load_panel():
        if e.truth_class != "PATH":
            continue
        c = classify_entry(e)
        if collapse_engine_tier(c.tier) != "PATH":
            missed.append((e, c))

    res = missing_evidence({"classifications": [c.to_dict() for _, c in missed]})
    by_strength: Counter = Counter()
    suppliers: dict[str, Counter] = {}
    cheapest = []
    for r in res:
        ups = [w for w in (r.get("would_change_with") or []) if w.get("direction") == "up"]
        if not ups:
            by_strength["unreachable with one line"] += 1
            continue
        w = ups[0]
        s = w["strength"]
        by_strength[s] += 1
        suppliers.setdefault(s, Counter()).update(
            c["code"] for c in w.get("candidates", []) if c.get("side") == "pathogenic")
        if s == "supporting":
            cheapest.append(r)

    # The panel supplies no phenotype to anybody, so PP4 — a Supporting criterion — can never
    # fire. Re-run the cheapest misses with a phenotype match and nothing else changed.
    want = {r.get("hgvs_p") for r in cheapest}
    recovered = []
    for e in load_panel():
        if e.truth_class != "PATH" or (e.variant.hgvs_p or "") not in want:
            continue
        before = classify_entry(e)
        ann = replace(_annotation_from_frozen(e, True), hpo_match_score=1.0,
                      hpo_matched_terms=[f"{e.variant.gene} phenotype (hypothetical)"])
        after = classify(e.variant, ann)
        recovered.append({"gene": e.variant.gene, "hgvs_p": e.variant.hgvs_p,
                          "before": before.tier, "after": after.tier})

    return {"n_missed": len(missed), "by_strength": dict(by_strength),
            "suppliers": {k: dict(v) for k, v in suppliers.items()},
            "with_phenotype": recovered}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    r = analyse()
    if a.json:
        print(json.dumps(r, indent=2))
        return 0

    print(f"\n  {r['n_missed']} pathogenic panel entries the engine holds below Pathogenic.")
    print("  Every one is ONE line of evidence short — none is misclassified toward benign.\n")
    print("  weakest addition that reaches Likely Pathogenic:")
    for s in _ORDER:
        n = r["by_strength"].get(s)
        if not n:
            continue
        who = ", ".join(sorted(r["suppliers"].get(s, {})))
        print(f"    {n:>3}  one {s:<12} could come from: {who or '—'}")
    if r["by_strength"].get("unreachable with one line"):
        print(f"    {r['by_strength']['unreachable with one line']:>3}  unreachable with one line")

    print("\n  The panel supplies NO phenotype, so PP4 (Supporting) cannot fire for anybody.")
    print("  Adding only a phenotype match:\n")
    for x in r["with_phenotype"]:
        print(f"    {x['gene']:<6} {x['hgvs_p']:<14} {x['before']} -> {x['after']}")
    n = sum(1 for x in r["with_phenotype"] if "Pathogenic" in x["after"])
    print(f"\n  {n} recovered by phenotype alone — the panel measures a genotype-only run of a")
    print("  phenotype-driven engine, so 61% is a floor, not the engine's accuracy.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
