"""Core data models.

Plain dataclasses (no pydantic) so the engine stays dependency-free and runs
headless anywhere. Every model is JSON-serialisable via ``to_dict`` so the MCP
tools can hand compact structures back to Claude.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------
@dataclass
class Variant:
    """A normalized, single-allele variant on GRCh38."""

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: Optional[str] = None  # e.g. missense_variant, stop_gained
    zygosity: Optional[str] = None     # het | hom | hemi
    depth: Optional[int] = None        # DP
    gq: Optional[int] = None           # genotype quality
    allele_balance: Optional[float] = None
    filter_status: Optional[str] = None  # VCF FILTER column
    info: dict[str, str] = field(default_factory=dict)  # raw INFO (annotator fields)
    alt_index: int = 0  # 0-based index of this ALT in the original record (for Number=A INFO)

    @property
    def key(self) -> str:
        """Canonical CHROM-POS-REF-ALT key used across annotators and caches."""
        chrom = self.chrom[3:] if self.chrom.lower().startswith("chr") else self.chrom
        return f"{chrom}-{self.pos}-{self.ref}-{self.alt}"

    @property
    def is_lof(self) -> bool:
        lof = {
            "stop_gained", "frameshift_variant", "splice_donor_variant",
            "splice_acceptor_variant", "start_lost", "stop_lost",
        }
        return (self.consequence or "") in lof

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("info", None)  # keep MCP/JSON output compact
        d["key"] = self.key
        return d


@dataclass
class Annotation:
    """Everything gathered about a variant from external / local sources."""

    clinvar_significance: Optional[str] = None  # Pathogenic, Benign, VUS, Conflicting...
    clinvar_review_status: Optional[str] = None
    clinvar_accession: Optional[str] = None
    clinvar_condition: Optional[str] = None
    clinvar_date: Optional[str] = None

    gnomad_af: Optional[float] = None       # popmax AF
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None
    gnomad_homozygotes: Optional[int] = None
    gnomad_popmax_pop: Optional[str] = None
    gnomad_faf95: Optional[float] = None    # filtering AF (95% CI lower bound, grpmax) — BS1/BA1

    abraom_af: Optional[float] = None       # Brazilian (SABE) allele frequency

    # gene-level constraint (for PVS1/PP2/BP1 judgment)
    gene_lof_intolerant: Optional[bool] = None  # e.g. pLI>=0.9 / low LOEUF

    # in-silico
    revel: Optional[float] = None
    cadd_phred: Optional[float] = None
    am_pathogenicity: Optional[float] = None  # AlphaMissense score (0..1)
    am_class: Optional[str] = None            # likely_benign | ambiguous | likely_pathogenic

    # phenotype
    hpo_match_score: Optional[float] = None      # 0..1 overlap patient<->gene
    hpo_matched_terms: list[str] = field(default_factory=list)

    source: dict[str, str] = field(default_factory=dict)  # field -> "db@date"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# ACMG
# ---------------------------------------------------------------------------
@dataclass
class CriterionResult:
    """One ACMG/AMP criterion evaluated with a citable trail.

    ``met`` is the on/off state; ``evidence`` holds the concrete values that
    drove it; ``citation`` names the source(s); ``reasoning`` is a one-liner.
    ``adjudicated_by`` is "engine" for deterministic criteria or "model" when a
    judgment criterion is left for Claude to decide — this separation is what
    makes the classification auditable rather than a black box.
    """

    code: str                 # PVS1, PM2, BA1, ...
    name: str
    default_strength: str     # very_strong | strong | moderate | supporting | stand_alone
    applies: bool             # False => N/A (e.g. needs a trio we don't have)
    met: bool = False
    applied_strength: Optional[str] = None
    evidence: dict[str, Any] = field(default_factory=dict)
    citation: list[str] = field(default_factory=list)
    reasoning: str = ""
    confidence: str = "high"  # high | moderate | low
    adjudicated_by: str = "engine"  # engine | model

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Classification:
    """The 5-tier call for one variant plus the full criterion trail."""

    variant: Variant
    annotation: Annotation
    criteria: list[CriterionResult]
    tier: str                 # Pathogenic | Likely Pathogenic | VUS | Likely Benign | Benign
    rule_path: str            # e.g. "PVS1 + PM2 + PP3 => Likely Pathogenic (LP-1)"

    @property
    def met_codes(self) -> list[str]:
        return [c.code for c in self.criteria if c.applies and c.met]

    def to_dict(self) -> dict[str, Any]:
        return {
            "variant": self.variant.to_dict(),
            "annotation": self.annotation.to_dict(),
            "tier": self.tier,
            "rule_path": self.rule_path,
            "met_codes": self.met_codes,
            "criteria": [c.to_dict() for c in self.criteria],
        }


# ---------------------------------------------------------------------------
# QC + report container
# ---------------------------------------------------------------------------
@dataclass
class QCSummary:
    total_variants: int = 0
    pass_filter: int = 0
    after_qc: int = 0
    after_rarity: int = 0
    after_impact: int = 0
    candidates: int = 0
    build: str = "GRCh38"
    warnings: list[str] = field(default_factory=list)
    # Spurious candidates removed thanks to ABraOM (Brazilian) frequencies.
    abraom_filtered: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
