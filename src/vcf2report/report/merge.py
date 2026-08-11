"""Merge per-chromosome ``results.json`` files into one whole-genome result.

The chunked WGS path (docs/WGS.md) annotates and classifies one chromosome at a
time to stay inside a small RAM budget, then stitches the per-chromosome
``results.json`` files back together here. ACMG classification is per-variant and
carries no cross-chromosome context, so concatenating the findings is exact; the
only care needed is the sequencing-quality panel, whose averages/medians must be
**weighted by each chunk's variant count**, never averaged naively.

Pure data-in/data-out (no I/O) so it is unit-testable without a real WGS run.
"""
from __future__ import annotations

from typing import Any

from ..vcf.seqqc import _assay_guess

# seq_quality fields that are additive counts across chunks.
_SQ_SUM = ("n_variants", "n_with_dp", "n_with_gq", "n_snv", "n_het", "n_hom",
           "n_indel", "n_sites", "n_multiallelic_sites", "n_with_rsid", "n_het_ab")
# seq_quality averages/medians/percentages -> (field, weight-count). A median cannot
# be recombined exactly from per-chunk medians; weighting by the count is the honest
# approximation (and far better than a naive mean of means).
_SQ_WAVG = {
    "dp_mean": "n_with_dp", "dp_median": "n_with_dp",
    "dp_pct_ge10": "n_with_dp", "dp_pct_ge20": "n_with_dp",
    "gq_median": "n_with_gq", "gq_pct_ge20": "n_with_gq",
    "titv": "n_snv",
    "pct_het_ab_balanced": "n_het_ab",
    "pct_pass": "n_variants", "pct_novel": "n_variants",
}


def _weighted(rows: list[dict], field: str, weight_key: str):
    num = den = 0.0
    for r in rows:
        v, w = r.get(field), r.get(weight_key) or 0
        if v is not None and w:
            num += float(v) * w
            den += w
    return round(num / den, 4) if den else None


def _ratio(a: float, b: float, nd: int = 3):
    return round(a / b, nd) if b else None


def merge_seq_quality(rows: list[dict]) -> dict:
    """Combine per-chunk seq_quality dicts, weighting rates by variant/site counts."""
    rows = [r for r in rows if r]
    if not rows:
        return {}
    out: dict[str, Any] = {}
    for k in _SQ_SUM:
        out[k] = sum(int(r.get(k) or 0) for r in rows)
    for field, wkey in _SQ_WAVG.items():
        out[field] = _weighted(rows, field, wkey)
    # Ratios/percentages that derive exactly from the summed counts — recompute them
    # (exact) instead of weighting (approximate).
    out["het_hom_ratio"] = _ratio(out["n_het"], out["n_hom"], 2)
    out["indel_snv_ratio"] = _ratio(out["n_indel"], out["n_snv"], 3)
    out["pct_multiallelic"] = (round(100 * out["n_multiallelic_sites"] / out["n_sites"], 4)
                               if out["n_sites"] else None)
    out["assay_guess"] = _assay_guess(out["n_variants"])
    notes: list[str] = []
    for r in rows:
        for n in r.get("notes", []) or []:
            if n not in notes:
                notes.append(n)
    out["notes"] = notes
    return out


def _merge_qc(rows: list[dict]) -> dict:
    """QC funnel: sum integer counters, concatenate any list fields (e.g. qc_rescued)."""
    out: dict[str, Any] = {}
    for r in rows:
        for k, v in (r or {}).items():
            if isinstance(v, list):
                out.setdefault(k, [])
                out[k].extend(v)
            elif isinstance(v, (int, float)) and not isinstance(v, bool):
                out[k] = out.get(k, 0) + v
            else:
                out.setdefault(k, v)
    return out


def merge_results(results: list[dict]) -> dict:
    """Stitch per-chromosome ``results.json`` dicts into one whole-genome result.

    Scalars are taken from the first chunk; the QC funnel counters are summed;
    every finding list (classifications, each bucket, clinvar_do_not_dismiss,
    conclusion) is concatenated; sequencing quality is recombined with count
    weighting; timings are summed. Order follows the input list — sort the inputs
    by chromosome first if a stable ordering matters.
    """
    results = [r for r in results if r]
    if not results:
        return {}
    if len(results) == 1:
        return dict(results[0])
    first = results[0]
    out: dict[str, Any] = dict(first)  # sample_id, build, hpo_terms, methods, provenance, ...

    out["qc"] = _merge_qc([r.get("qc", {}) for r in results])
    out["seq_quality"] = merge_seq_quality([r.get("seq_quality", {}) for r in results])

    for key in ("classifications", "conclusion", "clinvar_do_not_dismiss"):
        merged: list = []
        for r in results:
            merged.extend(r.get(key, []) or [])
        out[key] = merged

    buckets: dict[str, list] = {}
    for r in results:
        for name, items in (r.get("buckets", {}) or {}).items():
            buckets.setdefault(name, []).extend(items or [])
    out["buckets"] = buckets

    timings: dict[str, float] = {}
    for r in results:
        for k, v in (r.get("timings", {}) or {}).items():
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                timings[k] = round(timings.get(k, 0.0) + v, 4)
    out["timings"] = timings

    # newest generated timestamp across chunks (lexicographic ISO-8601 sorts chronologically)
    stamps = [r.get("generated") for r in results if r.get("generated")]
    if stamps:
        out["generated"] = max(stamps)
    out["chunks_merged"] = len(results)
    return out
