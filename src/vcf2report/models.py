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
    exon: Optional[str] = None         # SnpEff rank / VEP EXON as "N/M" (for PVS1 NMD tree)
    transcript: Optional[str] = None   # SnpEff feature_id / VEP Feature — reported with HGVS
    zygosity: Optional[str] = None     # het | hom | hemi
    # Half-call ("./1"): the sample certainly carries this ALT but the second allele was
    # not called, so `zygosity` above is the conservative "het" and NOT an observation.
    # Rendered as "het (partial call)" so a het/hom question is confirmed, not assumed.
    partial_call: bool = False
    depth: Optional[int] = None        # DP
    gq: Optional[int] = None           # genotype quality
    allele_balance: Optional[float] = None
    filter_status: Optional[str] = None  # VCF FILTER column
    variant_id: Optional[str] = None   # VCF ID column (dbSNP rsID if present, else None)
    n_alts: int = 1                    # ALT alleles at this site (>1 => multiallelic)
    info: dict[str, str] = field(default_factory=dict)  # raw INFO (annotator fields)
    alt_index: int = 0  # 0-based index of this ALT in the original record (for Number=A INFO)

    @property
    def key(self) -> str:
        """Canonical CHROM-POS-REF-ALT key used across annotators and caches.
        REF/ALT are upper-cased: VCF alleles are case-insensitive, and lowercase
        (soft-masked / lifted / hand-edited) input must still match the uppercase
        ClinVar/local cohort/gnomAD snapshots keyed off this string."""
        chrom = self.chrom[3:] if self.chrom.lower().startswith("chr") else self.chrom
        return f"{chrom}-{self.pos}-{self.ref.upper()}-{self.alt.upper()}"

    @property
    def is_lof(self) -> bool:
        # PVS1's null set (Richards 2015 / ClinGen SVI). stop_lost is deliberately
        # excluded — a stop-loss is a C-terminal extension, not a null; it earns PM4
        # (protein length change), never PVS1 (which would double-count with PM4).
        lof = {
            "stop_gained", "frameshift_variant", "splice_donor_variant",
            "splice_acceptor_variant", "start_lost",
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

    # Residue-level ClinVar cross-match (for PS1/PM5). Each is a dict describing the
    # matched pathogenic missense, or None. clinvar_residue_available flags whether the
    # residue index was loaded at all (so PS1/PM5 can say "unavailable" vs "no match").
    clinvar_ps1: Optional[dict] = None
    clinvar_pm5: Optional[dict] = None
    clinvar_residue_available: Optional[bool] = None
    # Local density of pathogenic missense AROUND the residue (PM1 hotspot signal); the query's
    # own residue is excluded, so this never overlaps the PS1/PM5 same-residue evidence.
    clinvar_hotspot: Optional[dict] = None

    gnomad_af: Optional[float] = None       # popmax AF
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None
    gnomad_homozygotes: Optional[int] = None
    gnomad_popmax_pop: Optional[str] = None
    gnomad_faf95: Optional[float] = None    # filtering AF (95% CI lower bound, grpmax) — BS1/BA1
    # True only when a source SURVEYED this locus and observed zero alleles. AN cannot express
    # this: the Parquet store's vouched-absence sentinel has no record and therefore no AN, so
    # `gnomad_an == 0` covers both "surveyed, found nothing" and "never surveyed". Consumers that
    # must not read missing data as evidence (is_hom_absent_artifact) need the distinction.
    gnomad_absence_vouched: bool = False

    # gnomAD saw the allele but the SITE failed a gnomAD filter (InbreedingCoeff, AS_VQSR, AC0).
    # The filter is about call quality at the site, not a retraction of the observation, so these
    # travel separately from the authoritative fields above: `gnomad_af` stays None (PM2 must not
    # read a filtered record as a surveyed frequency), while the benign criteria may use these to
    # establish that an allele is COMMON. One-way by construction — a frequency can only add
    # benign evidence, never manufacture pathogenicity.
    gnomad_af_filtered: Optional[float] = None
    gnomad_homozygotes_filtered: Optional[int] = None
    gnomad_filter: Optional[str] = None     # the failing filter string, for the evidence trail

    # Allele frequency from a cohort the OPERATOR supplies (see annotate/local_cohort.py).
    # None means not consulted / not in the table — never a checked zero.
    local_cohort_af: Optional[float] = None

    # gene-level constraint (for PVS1/PP2/BP1 judgment)
    gene_lof_intolerant: Optional[bool] = None  # e.g. pLI>=0.9 / low LOEUF
    gene_mis_z: Optional[float] = None            # gnomAD missense z-score (higher = more constrained)
    gene_oe_mis_upper: Optional[float] = None     # gnomAD obs/exp missense upper CI (>=1 = tolerated)
    gene_missense_constrained: Optional[bool] = None  # mis_z >= 3.09 -> PP2
    gene_missense_tolerant: Optional[bool] = None     # oe_mis_upper >= 1.0 -> BP1 (with LoF intolerance)

    # in-silico
    revel: Optional[float] = None
    cadd_phred: Optional[float] = None
    am_pathogenicity: Optional[float] = None  # AlphaMissense score (0..1)
    am_class: Optional[str] = None            # likely_benign | ambiguous | likely_pathogenic

    # phenotype
    hpo_match_score: Optional[float] = None      # 0..1 best-match-avg patient<->gene (PP4 + primary routing)
    hpo_best_match: Optional[float] = None        # 0..1 single strongest match (surfaced for display)
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
    # Spurious candidates removed thanks to an operator-supplied local cohort.
    local_cohort_filtered: list[str] = field(default_factory=list)
    # Variants dropped by the QC gate that ClinVar nevertheless classifies P/LP with
    # >=2-star review. QC runs BEFORE annotation, so the report's do-not-dismiss net
    # can never see these — they are recovered here so a marginal-quality call on a
    # known pathogenic allele is named and confirmed, not deleted in silence.
    qc_rescued: list[str] = field(default_factory=list)
    # Rare splice-ADJACENT variants the impact step set aside by design (see
    # vcf.filter.NEAR_SPLICE). Surfaced so the shortlist's coverage is stated, not assumed.
    near_splice_excluded: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SeqQuality:
    """Sequencing-quality estimate derived from the VCF's variant sites.

    A variants-only VCF carries data ONLY at called variant sites, so these are a
    proxy for how well the sample sequenced *at those sites* — NOT genome-wide
    breadth of coverage (that needs a gVCF or BAM). Depth (DP) at called sites is
    the closest coverage proxy; Ti/Tv and het:hom are orthogonal quality signals.
    """

    n_variants: int = 0
    assay_guess: str = "unknown"          # whole-genome / exome / panel / demo, by count
    n_with_dp: int = 0
    dp_mean: Optional[float] = None       # mean read depth at variant sites
    dp_median: Optional[float] = None
    dp_pct_ge10: Optional[float] = None   # % of sites with DP >= 10
    dp_pct_ge20: Optional[float] = None
    n_with_gq: int = 0
    gq_median: Optional[float] = None
    gq_pct_ge20: Optional[float] = None
    n_snv: int = 0
    titv: Optional[float] = None          # transition/transversion ratio (SNVs)
    n_het: int = 0
    n_hom: int = 0
    het_hom_ratio: Optional[float] = None
    n_indel: int = 0
    indel_snv_ratio: Optional[float] = None
    n_sites: int = 0
    n_multiallelic_sites: int = 0
    pct_multiallelic: Optional[float] = None
    n_with_rsid: int = 0
    pct_novel: Optional[float] = None          # % without a dbSNP rsID (only if annotated)
    n_het_ab: int = 0
    pct_het_ab_balanced: Optional[float] = None
    pct_pass: Optional[float] = None
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
