"""Assemble the end-to-end result into a single report model."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from .. import __version__
from ..config import (AF_BA1, AF_RECESSIVE_MAX, GENOME_BUILD, QC_MIN_DP,
                      QC_MIN_GQ)
from ..models import Classification, QCSummary, SeqQuality


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
    seq_quality: Optional[SeqQuality] = None  # estimated from the VCF's variant sites

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "hpo_terms": self.hpo_terms,
            "build": self.build,
            "tool_version": self.tool_version,
            "generated": self.generated,
            "qc": self.qc.to_dict(),
            "seq_quality": self.seq_quality.to_dict() if self.seq_quality else None,
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
    # Report the DETECTED build (qc.build), not the assumed default, so the header
    # can't disagree with the build-mismatch warning.
    return ReportModel(sample_id=sample_id, hpo_terms=hpo_terms, qc=qc, build=qc.build,
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
    from ..config import ACMG_SF_GENES, HPO_RELATED_MIN
    primary, secondary, other = [], [], []
    for c in classifications:
        related = (c.annotation.hpo_match_score or 0) >= HPO_RELATED_MIN
        is_sf = c.variant.gene in ACMG_SF_GENES
        if related and c.tier not in _BENIGN:
            primary.append(c)
        elif not related and c.tier in _PLP and is_sf:
            secondary.append(c)
        else:
            other.append(c)
    return primary, secondary, other


def summarize(report: "ReportModel") -> list[str]:
    """A deterministic, QC-aware interpretive conclusion for the report.

    Bottom-line-up-front: the likely explanatory finding (or its absence), any
    reportable secondary finding, an honest coverage caveat, the single-proband
    limitation, and recommended next steps. Derived only from the classifications
    and the sequencing-quality estimate — no model judgment.
    """
    primary, secondary, _other = split_findings(report.classifications)
    lines: list[str] = []

    diag = [c for c in primary if c.tier in _PLP]
    if diag:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in diag)
        lines.append(f"Likely explanatory finding for the clinical indication: **{g}** "
                     "(in a gene overlapping the patient's phenotype) — confirm and review.")
    else:
        vus = [c for c in primary if c.tier == "Uncertain Significance (VUS)"]
        msg = ("**No Pathogenic / Likely Pathogenic finding** was identified in "
               "phenotype-matched genes")
        if vus:
            msg += f"; {len(vus)} variant(s) of uncertain significance need further evaluation"
        lines.append(msg + ".")

    sec = [c for c in secondary if c.tier in _PLP]
    if sec:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in sec)
        lines.append(f"Reportable **secondary finding** (ACMG SF v3.2 — actionable, subject to "
                     f"the patient's opt-in policy): {g}.")

    sq = report.seq_quality
    if sq and sq.dp_median is not None and sq.dp_median < 20:
        lines.append(f"⚠️ **Coverage limitation:** median depth at variant sites is {sq.dp_median}x, "
                     "below a 20–30x clinical target. Findings in low-coverage regions are less "
                     "reliable and a negative result does not exclude a diagnosis — consider "
                     "higher-depth resequencing before ruling out a genetic cause.")
    elif sq and sq.dp_median is not None:
        lines.append(f"Sequencing depth at called sites is adequate (median {sq.dp_median}x); note "
                     "that a variant-only VCF conveys no breadth of coverage, so poorly-covered "
                     "regions cannot be assessed from this input.")

    lines.append("Single-proband analysis: de novo / segregation / phasing criteria (PS2, PM3, "
                 "PM6, PP1, BS4) are N/A — parental or trio testing could upgrade candidates or "
                 "resolve VUS.")
    lines.append("**Recommended next steps:** expert review and sign-out; orthogonal confirmation "
                 "(e.g. Sanger) of any reported P/LP variant; segregation / functional evidence to "
                 "resolve variants of uncertain significance.")
    return lines
