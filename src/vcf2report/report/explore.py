"""Persist a run's full structured result as JSON, so an operator can EXPLORE it conversationally
(with Claude) after the laudo is rendered — "show the VUS in gene X", "why did Y get PM2", "which
findings rest on ClinVar", "open case Z" — all answered from this file, with no re-run.

The laudo stays the deliverable; this keeps the underlying data live and queryable. `build_explore`
returns the whole report (every variant with its full ACMG criterion trail + evidence + reasoning),
the routed buckets, the conclusion, and the ClinVar do-not-dismiss list.

The read side (``load_explore`` + the query helpers below) answers the conversational questions
straight off that file — no engine, no network, no re-run — so a follow-up chat about a case is a
cheap lookup, not a re-analysis. Every helper takes the loaded ``data`` dict and returns plain
JSON-able structures.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .assemble import (carrier_findings, clinvar_pathogenic_flags, clinvar_stars,
                       split_findings, summarize)
from .vus_triage import probable_pathogenic_vus

# The routed buckets, in report order. Also the valid names for ``variants_in_bucket``.
BUCKETS = ("primary", "secondary", "carrier", "probable_pathogenic_vus", "other")


# ---------------------------------------------------------------------------
# Write side — persist the run
# ---------------------------------------------------------------------------
def _flag_ref(c) -> dict[str, Any]:
    """A compact, actionable record for the ClinVar do-not-dismiss safety list.

    Self-contained (does not just point into ``classifications``) because this list must never be
    lost or need a second lookup: it is the "a well-reviewed ClinVar-pathogenic variant is sitting
    below the engine's tier — do not report it as no-finding" flag.
    """
    v, a = c.variant, c.annotation
    return {
        "gene": v.gene,
        "variant": v.key,
        "hgvs": v.hgvs_p or v.hgvs_c,
        "engine_tier": c.tier,
        "clinvar_significance": a.clinvar_significance,
        "clinvar_review_status": a.clinvar_review_status,
        "clinvar_stars": clinvar_stars(a.clinvar_review_status),
        "clinvar_accession": a.clinvar_accession,
    }


def build_explore(report) -> dict[str, Any]:
    primary, secondary, other = split_findings(report.classifications)
    carriers = carrier_findings(report.classifications)
    vus = probable_pathogenic_vus(report.classifications)

    data = report.to_dict()  # sample, build, hpo, qc funnel, seq_quality, every classification (+criteria)
    data["conclusion"] = summarize(report)
    data["buckets"] = {
        "primary": [c.variant.gene for c in primary],
        "secondary": [c.variant.gene for c in secondary],
        "carrier": [c.variant.gene for c in carriers],
        "probable_pathogenic_vus": [e["classification"].variant.gene for e in vus],
        "other": [c.variant.gene for c in other],
    }
    # Serialise as compact dicts — NOT raw Classification objects. Left as objects they would fall to
    # json.dump's ``default=str`` and land as a dataclass repr string: unqueryable, and silently so
    # (the demo case has an empty list, hiding it). This keeps the safety list machine-readable.
    data["clinvar_do_not_dismiss"] = [_flag_ref(c) for c in clinvar_pathogenic_flags(report.classifications)]
    return data


def write_explore(report, path: str | Path) -> str:
    path = str(path)
    with open(path, "w") as fh:
        json.dump(build_explore(report), fh, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# Read side — explore a persisted run without re-running the pipeline
# ---------------------------------------------------------------------------
def load_explore(path: str | Path) -> dict[str, Any]:
    """Load a persisted ``*_results.json`` back into a dict. This is "open case Z"."""
    with open(str(path)) as fh:
        return json.load(fh)


def overview(data: dict[str, Any]) -> dict[str, Any]:
    """Case-level digest — the answer to "summarise this case / what did we find".

    Bucket *counts* (not the gene lists) plus the sample header and the deterministic conclusion,
    so a follow-up chat can orient without walking every classification.
    """
    buckets = data.get("buckets", {})
    return {
        "sample_id": data.get("sample_id"),
        "build": data.get("build"),
        "generated": data.get("generated"),
        "hpo_terms": data.get("hpo_terms", []),
        "n_candidates": len(data.get("classifications", [])),
        "bucket_counts": {b: len(buckets.get(b, [])) for b in BUCKETS},
        "clinvar_do_not_dismiss": data.get("clinvar_do_not_dismiss", []),
        "conclusion": data.get("conclusion", []),
    }


def findings_for_gene(data: dict[str, Any], gene: str) -> list[dict[str, Any]]:
    """Every classified variant in ``gene`` (case-insensitive) — "show me gene X"."""
    g = (gene or "").upper()
    return [c for c in data.get("classifications", [])
            if (c.get("variant", {}).get("gene") or "").upper() == g]


def variants_in_bucket(data: dict[str, Any], bucket: str) -> list[dict[str, Any]]:
    """The classifications routed into ``bucket`` — e.g. ``variants_in_bucket(d, "probable_pathogenic_vus")``
    is "show the VUS worth review". Raises ValueError on an unknown bucket name (naming the valid ones)
    so a typo surfaces instead of silently returning nothing."""
    buckets = data.get("buckets", {})
    if bucket not in buckets:
        raise ValueError(f"unknown bucket {bucket!r}; valid: {', '.join(BUCKETS)}")
    genes = set(buckets.get(bucket, []))
    return [c for c in data.get("classifications", [])
            if c.get("variant", {}).get("gene") in genes]


def criterion_basis(data: dict[str, Any], gene: str, code: str) -> list[dict[str, Any]]:
    """Why did ``gene`` get ``code`` (e.g. PM2)? Returns the criterion's full trail per matching variant —
    whether it applied/was met, at what strength, the concrete evidence, the citation, the one-line
    reasoning, and who adjudicated it (engine vs. left for the model). This is the audit answer, verbatim
    from the persisted run."""
    code_u = (code or "").upper()
    out: list[dict[str, Any]] = []
    for c in findings_for_gene(data, gene):
        v = c.get("variant", {})
        for cr in c.get("criteria", []):
            if (cr.get("code") or "").upper() != code_u:
                continue
            out.append({
                "gene": v.get("gene"),
                "variant": v.get("key"),
                "hgvs": v.get("hgvs_p") or v.get("hgvs_c"),
                "tier": c.get("tier"),
                "code": cr.get("code"),
                "name": cr.get("name"),
                "applies": cr.get("applies"),
                "met": cr.get("met"),
                "strength": cr.get("applied_strength") or cr.get("default_strength"),
                "evidence": cr.get("evidence"),
                "reasoning": cr.get("reasoning"),
                "citation": cr.get("citation"),
                "adjudicated_by": cr.get("adjudicated_by"),
            })
    return out


def _cites_clinvar(cr: dict[str, Any]) -> bool:
    """Does a criterion rest on a ClinVar assertion? PP5/BP6 are the ClinVar-assertion criteria by
    definition; otherwise look for a ClinVar accession in the citation — the trail cites the ``VCV…``
    accession (or the word ClinVar), NOT necessarily the literal string "ClinVar", so matching text
    alone would miss PP5."""
    if (cr.get("code") or "").upper() in {"PP5", "BP6"}:
        return True
    for s in (cr.get("citation") or []):
        s = (s or "").lower()
        if "clinvar" in s or s.startswith("vcv"):
            return True
    return False


def findings_citing_clinvar(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Which findings rest on ClinVar — every classification with a MET criterion that draws on a
    ClinVar assertion (PP5/BP6 or a cited ``VCV`` accession). Returns the gene, the engine tier, the
    ClinVar significance, and just the ClinVar-citing criteria, so the reader can weigh the
    circularity/anti-circularity themselves."""
    out: list[dict[str, Any]] = []
    for c in data.get("classifications", []):
        cites = [cr for cr in c.get("criteria", []) if cr.get("met") and _cites_clinvar(cr)]
        if not cites:
            continue
        v = c.get("variant", {})
        out.append({
            "gene": v.get("gene"),
            "variant": v.get("key"),
            "tier": c.get("tier"),
            "clinvar_significance": c.get("annotation", {}).get("clinvar_significance"),
            "criteria": [{"code": cr.get("code"), "met": cr.get("met"),
                          "reasoning": cr.get("reasoning"), "citation": cr.get("citation")}
                         for cr in cites],
        })
    return out


def explain(data: dict[str, Any], gene: str) -> dict[str, Any]:
    """A gene-level digest — "tell me about gene X". Per classified variant: tier, the combining-rule
    path, the met criteria, HGVS/coordinate, the headline clinical annotation, and which routed
    bucket(s) it fell into. A single call that orients a follow-up conversation on one gene."""
    buckets = data.get("buckets", {})
    g = (gene or "").upper()
    in_buckets = [b for b in BUCKETS if g in [x.upper() for x in buckets.get(b, [])]]
    variants: list[dict[str, Any]] = []
    for c in findings_for_gene(data, gene):
        v, a = c.get("variant", {}), c.get("annotation", {})
        variants.append({
            "variant": v.get("key"),
            "hgvs_c": v.get("hgvs_c"),
            "hgvs_p": v.get("hgvs_p"),
            "consequence": v.get("consequence"),
            "zygosity": v.get("zygosity"),
            "tier": c.get("tier"),
            "rule_path": c.get("rule_path"),
            "met_codes": c.get("met_codes", []),
            "clinvar_significance": a.get("clinvar_significance"),
            "gnomad_af": a.get("gnomad_af"),
            "am_pathogenicity": a.get("am_pathogenicity"),
            "hpo_match_score": a.get("hpo_match_score"),
        })
    return {"gene": gene, "buckets": in_buckets, "variants": variants}


# ---------------------------------------------------------------------------
# Tiny terminal helper — isolated from the pipeline CLI (no risk to run_headless)
# ---------------------------------------------------------------------------
def _cli(argv: list[str] | None = None) -> int:  # pragma: no cover - thin arg plumbing
    import argparse
    p = argparse.ArgumentParser(
        prog="vcf2report-explore",
        description="Query a persisted *_results.json without re-running the pipeline.")
    p.add_argument("results_json", help="Path to a <case>_results.json written by the pipeline.")
    p.add_argument("--gene", help="Show the classified variant(s) in this gene.")
    p.add_argument("--criterion", help="With --gene: why the gene got this ACMG code (e.g. PM2).")
    p.add_argument("--bucket", choices=BUCKETS, help="Show the variants routed into this bucket.")
    p.add_argument("--clinvar", action="store_true", help="Show findings that rest on ClinVar.")
    args = p.parse_args(argv)

    data = load_explore(args.results_json)
    if args.criterion:
        if not args.gene:
            p.error("--criterion requires --gene")
        result: Any = criterion_basis(data, args.gene, args.criterion)
    elif args.gene:
        result = explain(data, args.gene)
    elif args.bucket:
        result = variants_in_bucket(data, args.bucket)
    elif args.clinvar:
        result = findings_citing_clinvar(data)
    else:
        result = overview(data)
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_cli())
