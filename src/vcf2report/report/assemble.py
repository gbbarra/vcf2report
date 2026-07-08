"""Assemble the end-to-end result into a single report model."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from .. import __version__
from ..config import (AF_BA1, AF_RECESSIVE_MAX, GENOME_BUILD, QC_MIN_DP,
                      QC_MIN_GQ)
from ..models import Classification, QCSummary


@dataclass
class ReportModel:
    sample_id: str
    hpo_terms: list[str]
    qc: QCSummary
    classifications: list[Classification]  # ranked candidates, classified
    build: str = GENOME_BUILD
    tool_version: str = __version__
    generated: str = ""
    methods: dict[str, Any] = field(default_factory=dict)
    timings: dict[str, float] = field(default_factory=dict)  # per-stage seconds

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "hpo_terms": self.hpo_terms,
            "build": self.build,
            "tool_version": self.tool_version,
            "generated": self.generated,
            "qc": self.qc.to_dict(),
            "methods": self.methods,
            "timings": self.timings,
            "classifications": [c.to_dict() for c in self.classifications],
        }


def build_report(sample_id: str, hpo_terms: list[str], qc: QCSummary,
                 classifications: list[Classification]) -> ReportModel:
    methods = {
        "genome_build": GENOME_BUILD,
        "qc_thresholds": {"min_DP": QC_MIN_DP, "min_GQ": QC_MIN_GQ},
        "rarity_cutoff_popmax_af": AF_RECESSIVE_MAX,
        "ba1_cutoff": AF_BA1,
        "databases": ["ClinVar", "gnomAD r4", "ABraOM (SABE)", "HPO", "gnomAD constraint"],
        "guideline": "ACMG/AMP (Richards et al., Genet Med 2015)",
    }
    # reportable = anything not benign, ordered by clinical relevance
    order = {"Pathogenic": 0, "Likely Pathogenic": 1,
             "Uncertain Significance (VUS)": 2, "Likely Benign": 3, "Benign": 4}
    ranked = sorted(classifications, key=lambda c: order.get(c.tier, 9))
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return ReportModel(sample_id=sample_id, hpo_terms=hpo_terms, qc=qc,
                       classifications=ranked, methods=methods, generated=generated)
