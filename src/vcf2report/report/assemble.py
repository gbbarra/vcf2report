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
        "standards": [
            "ACMG/AMP variant classification (Richards et al., Genet Med 2015)",
            "ClinGen SVI criteria refinements",
            "ACMG secondary-findings list (SF v3.2, Miller et al. 2023)",
            "HGVS nomenclature",
            "GA4GH Phenopackets (phenotype exchange)",
        ],
    }
    # reportable = anything not benign, ordered by clinical relevance
    order = {"Pathogenic": 0, "Likely Pathogenic": 1,
             "Uncertain Significance (VUS)": 2, "Likely Benign": 3, "Benign": 4}
    ranked = sorted(classifications, key=lambda c: order.get(c.tier, 9))
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return ReportModel(sample_id=sample_id, hpo_terms=hpo_terms, qc=qc,
                       classifications=ranked, methods=methods, generated=generated)


_PLP = {"Pathogenic", "Likely Pathogenic"}
_BENIGN = {"Benign", "Likely Benign"}


def split_findings(classifications):
    """Partition reported candidates into (primary, secondary, other).

    * **primary** — phenotype-related (HPO overlap) AND not benign: diagnostic.
    * **secondary** — an unrelated P/LP variant in an **ACMG SF v3.2 gene**: a
      reportable, actionable secondary finding (subject to patient opt-in).
    * **other** — everything else, incl. unrelated P/LP in a non-SF gene (an
      incidental finding that is not on the actionable SF list), phenotype-matched
      benign, and unrelated VUS/benign.
    """
    from ..config import ACMG_SF_GENES
    primary, secondary, other = [], [], []
    for c in classifications:
        related = (c.annotation.hpo_match_score or 0) > 0
        is_sf = c.variant.gene in ACMG_SF_GENES
        if related and c.tier not in _BENIGN:
            primary.append(c)
        elif not related and c.tier in _PLP and is_sf:
            secondary.append(c)
        else:
            other.append(c)
    return primary, secondary, other
