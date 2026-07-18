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
    plp_hits = _plp_hits_by_gene(classifications)
    for c in classifications:
        related = (c.annotation.hpo_match_score or 0) >= HPO_RELATED_MIN
        # QC caution: a homozygous genotype for a variant the store vouches is absent from gnomAD
        # (AC=0) is genotype-implausible — a homozygote needs the allele to exist — and a classic
        # calling-artifact signature in difficult regions. On a HEALTHY exome that is noise. But
        # hom + gnomAD-absent + P/LP + phenotype-matched is ALSO the textbook signature of a
        # recessive DIAGNOSIS (a rare pathogenic allele, homozygous in an affected proband). The
        # discriminator is the phenotype: demote it UNLESS it is a phenotype-matched P/LP finding,
        # where the imperative to surface the candidate wins and the report carries the
        # "confirm the genotype" caveat (surfaced via is_hom_absent_artifact) instead of hiding it.
        if is_hom_absent_artifact(c) and not (related and c.tier in _PLP):
            other.append(c)
            continue
        # Carrier caution: a lone heterozygous null in a recessive-only gene is a PATHOGENIC
        # VARIANT but a NON-DIAGNOSTIC genotype — the carrier is healthy. Phenotype routing
        # cannot catch this, because recessive disease genes have exactly the phenotypes a
        # proband presents with, so the carrier clears HPO_RELATED_MIN and lands in `primary`
        # next to (or instead of) the real answer. Everyone carries a few of these. Keep the
        # ACMG tier — it is correct — but report it as carrier status, not as a diagnosis.
        if _is_carrier(c, plp_hits):
            other.append(c)
            continue
        # `related` (best-match-average, computed above): a random, unrelated phenotype clears the
        # max on one broad term far too often, so the max is not specific. The average requires the
        # phenotype as a whole to fit the gene — measured 2-3x more discriminative vs a decoy.
        is_sf = c.variant.gene in ACMG_SF_GENES
        if related and c.tier not in _BENIGN:
            primary.append(c)
        elif not related and c.tier in _PLP and is_sf:
            secondary.append(c)
        else:
            other.append(c)
    return primary, secondary, other


def _plp_hits_by_gene(classifications) -> dict:
    """gene -> number of P/LP calls, counted ONCE per report.

    The carrier test needs "does this gene have a second hit?", and asking that per
    classification is O(n^2) — on a real annotated exome n is ~1200 candidates and the test
    runs from three places (split_findings, the report's carrier section, the conclusion).
    """
    hits: dict = {}
    for c in classifications:
        if c.tier in _PLP:
            hits[c.variant.gene] = hits.get(c.variant.gene, 0) + 1
    return hits


def _is_carrier(c, hits: dict) -> bool:
    from .. import config
    if c.tier not in _PLP or c.variant.zygosity != "het":
        return False
    m = config.gene_inheritance_modes(c.variant.gene)
    if "AR" not in m or "AD" in m or "XL" in m:
        return False
    return hits.get(c.variant.gene, 0) < 2


def is_unconfirmed_ar_carrier(c, classifications) -> bool:
    """A lone heterozygous P/LP in a gene whose only known disease mechanism is recessive.

    ACMG classifies the VARIANT and the Pathogenic tier is right — but a single het in a
    recessive gene does not explain a proband's phenotype, it makes them a healthy carrier.
    An average person carries 2-3 such alleles, so presenting them as diagnostic findings
    both floods the report and, worse, lets a carrier outrank the true diagnosis (measured:
    a het LIPA/SKIC2 carrier displaced the real answer into "other").

    It also keeps the ACMG SF v3.2 contract: the recessive SF genes (ATP7B, MUTYH, BTD,
    GAA, HFE, CASQ2, TRDN, RPE65) are reportable as actionable secondary findings ONLY when
    biallelic — a carrier must not be reported. Routing them out of `secondary` here honours
    that generically, via the gene's mechanism rather than a hard-coded list.

    Deliberately narrow — it must never hide a real diagnosis:
      * genes with ANY dominant/X-linked disease are excluded: a het there can be diagnostic;
      * ``hom`` is excluded: that is biallelic, i.e. exactly the diagnostic genotype;
      * a SECOND P/LP hit in the same gene is excluded: possible compound heterozygote (we
        cannot phase it, but it is a genuine candidate the clinician must see).
    """
    return _is_carrier(c, _plp_hits_by_gene(classifications))


def carrier_findings(classifications):
    """The recessive carrier alleles routed out of the diagnostic sections.

    Not noise to be discarded: carrier status carries real reproductive relevance and the
    report should show it — just not competing with the diagnosis.
    """
    hits = _plp_hits_by_gene(classifications)
    return [c for c in classifications if _is_carrier(c, hits)]


def is_hom_absent_artifact(c) -> bool:
    """Homozygous genotype for a variant the store vouches is absent from gnomAD (AC=0). A
    homozygote requires the allele to exist in the population, so AC=0 + hom is implausible for a
    real allele and a common calling-artifact signature in difficult regions (segdup / low-complexity
    / homopolymer). A QC caution only — the ACMG tier is untouched; it is just not presented as a
    confident diagnostic finding. Heterozygous variants (incl. genuine novel dominant LoF) are
    unaffected."""
    return (c.variant.zygosity == "hom") and (c.annotation.gnomad_af == 0.0)


def clinvar_stars(review_status) -> int:
    """ClinVar review status -> star count (0-4).

    Normalize underscores to spaces first: the VCF-INFO path (from_vcf) and the live
    E-utilities path both deliver a space-delimited status, so matching only underscore
    tokens would silently score every real assertion 0 and disable the safety flag.
    """
    r = (review_status or "").lower().replace("_", " ").strip()
    if "practice guideline" in r:
        return 4
    if "reviewed by expert panel" in r:
        return 3
    if "multiple submitters" in r and "no conflict" in r:
        return 2
    if r.startswith("criteria provided") or "single submitter" in r or "conflicting" in r:
        return 1
    return 0


def clinvar_pathogenic_flags(classifications):
    """Candidates ClinVar classifies P/LP with >=2-star review whose independent engine ACMG
    tier is NOT P/LP. These MUST be surfaced: never present a well-reviewed known-pathogenic
    variant as 'no finding'. This does not touch the ACMG criteria math (avoids PP5
    circularity) — it is a report-level safety flag."""
    out = []
    for c in classifications:
        sig = (c.annotation.clinvar_significance or "").lower()
        is_plp = sig.startswith("pathogenic") or sig.startswith("likely pathogenic")
        if is_plp and clinvar_stars(c.annotation.clinvar_review_status) >= 2 and c.tier not in _PLP:
            out.append(c)
    return out


def summarize(report: "ReportModel") -> list[str]:
    """A deterministic, QC-aware interpretive conclusion for the report.

    Bottom-line-up-front: the likely explanatory finding (or its absence), any
    reportable secondary finding, an honest coverage caveat, the single-proband
    limitation, and recommended next steps. Derived only from the classifications
    and the sequencing-quality estimate — no model judgment.
    """
    primary, secondary, other = split_findings(report.classifications)
    lines: list[str] = []

    diag = [c for c in primary if c.tier in _PLP]
    if diag:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in diag)
        lines.append(f"Likely explanatory finding for the clinical indication: **{g}** "
                     "(in a gene overlapping the patient's phenotype) — confirm and review.")
    else:
        vus = [c for c in primary if c.tier == "Uncertain Significance (VUS)"]
        msg = ("**No Pathogenic / Likely Pathogenic finding** by the engine's independent "
               "ACMG classification in phenotype-matched genes")
        if vus:
            msg += f"; {len(vus)} variant(s) of uncertain significance need further evaluation"
        lines.append(msg + ".")

    # Safety flag: a known, well-reviewed ClinVar-Pathogenic variant must never be hidden
    # behind a lower engine tier (surfaced independently of the ACMG math).
    flagged = clinvar_pathogenic_flags(report.classifications)
    if flagged:
        g = "; ".join(f"{c.variant.gene} ({clinvar_stars(c.annotation.clinvar_review_status)}★; "
                      f"engine: {c.tier})" for c in flagged)
        lines.append(f"⚠️ **Classified Pathogenic/Likely Pathogenic in ClinVar** (≥2-star review) — "
                     f"the engine's independent tier is lower, but DO NOT dismiss: **{g}**. Review the "
                     "ClinVar assertion and its underlying evidence.")

    artifacts = [c for c in report.classifications if is_hom_absent_artifact(c) and c.tier in _PLP]
    if artifacts:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in artifacts)
        lines.append(f"⚠️ **Verify the genotype before interpreting** — {len(artifacts)} homozygous "
                     f"variant(s) that are absent from gnomAD (AC=0), which is implausible for a real "
                     f"allele and a common calling-artifact signature in difficult regions: {g}. Confirm "
                     "the call (orthogonal / Sanger) before interpreting these.")

    sec = [c for c in secondary if c.tier in _PLP]
    if sec:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in sec)
        lines.append(f"Reportable **secondary finding** (ACMG SF v3.2 — actionable, subject to "
                     f"the patient's opt-in policy): {g}.")

    # Recessive carriers get their own sentence: "clinical relevance is uncertain" is simply
    # WRONG for them — the relevance is known and it is reproductive, not diagnostic. Lumping a
    # carrier in with genuine incidental P/LP invites the reader to weigh it as a candidate.
    carriers = carrier_findings(report.classifications)
    if carriers:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in carriers)
        lines.append(f"**Carrier finding(s)** — heterozygous {('allele' if len(carriers) == 1 else 'alleles')} "
                     f"in gene(s) whose disease mechanism is recessive: {g}. A single copy does NOT "
                     "cause the condition and does NOT explain the indication; this is carrier "
                     "status, relevant to reproductive counselling, not a diagnosis.")

    # An engine-P/LP variant that is neither phenotype-matched nor on the SF list still
    # belongs in the conclusion — it is in the ranked table but must not be silent here.
    inc = [c for c in other if c.tier in _PLP and c not in carriers]
    if inc:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in inc)
        lines.append(f"Additional **Pathogenic / Likely Pathogenic** variant(s) not matching the "
                     f"stated phenotype and not on the ACMG SF actionable list: {g}. Clinical "
                     "relevance to the indication is uncertain — review in context.")

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
