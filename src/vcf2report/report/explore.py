"""Persist a run's full structured result as JSON, so an operator can EXPLORE it conversationally
(with Claude) after the laudo is rendered — "show the VUS in gene X", "why did Y get PM2", "which
findings rest on ClinVar", "open case Z" — all answered from this file, with no re-run.

The laudo stays the deliverable; this keeps the underlying data live and queryable. `build_explore`
returns the whole report (every variant with its full ACMG criterion trail + evidence + reasoning),
the routed buckets, the conclusion, and the ClinVar do-not-dismiss list.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .assemble import carrier_findings, clinvar_pathogenic_flags, split_findings, summarize
from .vus_triage import probable_pathogenic_vus


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
    data["clinvar_do_not_dismiss"] = clinvar_pathogenic_flags(report.classifications)
    return data


def write_explore(report, path: str | Path) -> str:
    path = str(path)
    with open(path, "w") as fh:
        json.dump(build_explore(report), fh, indent=2, default=str)
    return path
