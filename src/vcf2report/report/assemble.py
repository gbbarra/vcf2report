"""Assemble the end-to-end result into a single report model."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from .. import __version__, demo
from ..config import (AF_BA1, AF_RECESSIVE_MAX, GENOME_BUILD, QC_MIN_DP,
                      QC_MIN_GQ)
from ..models import Classification, QCSummary, SeqQuality
from .vus_triage import probable_pathogenic_vus


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
    # Non-empty only when this run was NOT an ordinary full-store analysis — currently just
    # demo mode (see vcf2report.demo). Carried on the model rather than derived at render time
    # so every renderer, the JSON, and the Methods block stamp it from one source.
    provenance: dict[str, Any] = field(default_factory=dict)

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
            "provenance": self.provenance,
            "timings": self.timings,
            "classifications": [c.to_dict() for c in self.classifications],
        }


def build_report(sample_id: str, hpo_terms: list[str], qc: QCSummary,
                 classifications: list[Classification],
                 provenance: Optional[dict] = None) -> ReportModel:
    methods = {
        "genome_build": GENOME_BUILD,
        "qc_thresholds": {"min_DP": QC_MIN_DP, "min_GQ": QC_MIN_GQ},
        "rarity_cutoff_popmax_af": AF_RECESSIVE_MAX,
        "ba1_cutoff": AF_BA1,
        "databases": ["ClinVar", "gnomAD r4", "local cohort", "HPO", "gnomAD constraint"],
        "standards": [
            "ACMG/AMP variant classification (Richards et al., Genet Med 2015)",
            "ClinGen SVI criteria refinements",
            "ACMG secondary-findings list (SF v3.2, Miller et al. 2023)",
            "HGVS nomenclature",
            "GA4GH Phenopackets (phenotype exchange)",
        ],
    }
    provenance = dict(provenance or {})
    if provenance.get("mode") == "demo":
        # Stamped in Methods as well as the conclusion: a reader who skips to "which databases
        # was this run against?" must not find the ordinary answer.
        absent = provenance.get("stores_absent") or []
        methods["data_mode"] = (
            "DEMONSTRATION — committed synthetic example VCF; Parquet stores "
            + (f"absent ({', '.join(absent)}) — nothing looked up live" if absent else "present")
        )
        # Two separate keys, because "baked into the fixture from the full databases" and
        # "a slice covering only the demo genes" are different claims about the same laudo.
        if provenance.get("criteria_from_fixture_info"):
            methods["criteria_from_fixture_baked_info"] = provenance["criteria_from_fixture_info"]
        if provenance.get("criteria_from_frozen_slice"):
            methods["criteria_from_frozen_demo_slice"] = provenance["criteria_from_frozen_slice"]
    # reportable = anything not benign, ordered by clinical relevance
    order = {"Pathogenic": 0, "Likely Pathogenic": 1,
             "Uncertain Significance (VUS)": 2, "Likely Benign": 3, "Benign": 4}
    ranked = sorted(classifications, key=lambda c: order.get(c.tier, 9))
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    # Report the DETECTED build (qc.build), not the assumed default, so the header
    # can't disagree with the build-mismatch warning.
    return ReportModel(sample_id=sample_id, hpo_terms=hpo_terms, qc=qc, build=qc.build,
                       classifications=ranked, methods=methods, generated=generated,
                       provenance=provenance)


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

    A phenotype score of ``None`` means NO COMPARISON WAS POSSIBLE (no HPO terms were
    supplied, or the gene is absent from the local HPO table) — never "compared, and it
    did not match". A P/LP variant is therefore NOT demoted on an absent comparison; see
    :func:`phenotype_compared`.
    """
    from ..config import ACMG_SF_GENES, HPO_RELATED_MIN
    primary, secondary, other = [], [], []
    plp_hits = _plp_hits_by_gene(classifications)
    for c in classifications:
        score = c.annotation.hpo_match_score
        compared = score is not None
        related = compared and score >= HPO_RELATED_MIN
        # QC caution: a homozygous genotype for a variant the store vouches is absent from gnomAD
        # (AC=0) is genotype-implausible — a homozygote needs the allele to exist — and a classic
        # calling-artifact signature in difficult regions. On a HEALTHY exome that is noise. But
        # hom + gnomAD-absent + P/LP + phenotype-matched is ALSO the textbook signature of a
        # recessive DIAGNOSIS (a rare pathogenic allele, homozygous in an affected proband). The
        # discriminator is the phenotype: demote it UNLESS it is a phenotype-matched P/LP finding,
        # where the imperative to surface the candidate wins and the report carries the
        # "confirm the genotype" caveat (surfaced via is_hom_absent_artifact) instead of hiding it.
        # An ACMG-SF P/LP is likewise never demoted here: this gate used to run BEFORE the SF
        # branch below, so an actionable secondary finding was suppressed AND the conclusion
        # then asserted the gene was "not on the ACMG SF actionable list". It carries the same
        # confirm-the-genotype caveat instead.
        _sf_plp = c.variant.gene in ACMG_SF_GENES and c.tier in _PLP
        if is_hom_absent_artifact(c) and not (related and c.tier in _PLP) and not _sf_plp:
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
        elif not compared and c.tier in _PLP and not is_sf:
            # Phenotype comparison was not possible. Routing this to `other` made the report
            # state the variant does "not match the stated phenotype" — a claim about a
            # comparison that never happened, and on a run with no --hpo at all it emptied
            # `primary` while a Pathogenic variant sat under "Not routinely reported".
            primary.append(c)
        elif not related and c.tier in _PLP and is_sf:
            secondary.append(c)
        else:
            other.append(c)
    return primary, secondary, other


def phenotype_compared(classifications) -> bool:
    """True when at least one candidate had a computable phenotype score.

    ``hpo.match`` returns None both when no HPO terms were supplied and when the gene is
    absent from the local HPO table; either way the report must not describe a variant as
    "not matching the stated phenotype"."""
    return any(c.annotation.hpo_match_score is not None for c in classifications)


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
    """Homozygous genotype for a variant gnomAD **actually observed and vouches is absent**
    (AN>0 alleles surveyed, AC=0 seen). A homozygote requires the allele to exist in the
    population, so a vouched AC=0 + hom is implausible for a real allele and a common
    calling-artifact signature in difficult regions (segdup / low-complexity / homopolymer).
    A QC caution only — the ACMG tier is untouched; it is just not presented as a confident
    diagnostic finding. Heterozygous variants (incl. genuine novel dominant LoF) are unaffected.

    The absence must be VOUCHED. Annotators write ``gnomAD_AF=0`` alongside ``gnomAD_AN=0``
    wherever gnomAD has no coverage at all — AF is 0/0, undefined, not a survey that found
    nothing. Reading that as a vouched absence turns missing data into evidence, and demotes
    exactly the recessive candidates that live in gnomAD-uncovered regions. That case is
    :func:`is_hom_gnomad_uncovered` and carries a different, weaker caveat.

    ``gnomad_absence_vouched`` carries the distinction because ``AN`` cannot: the Parquet
    store's vouched-absence sentinel has no gnomAD row and therefore no AN, so testing
    ``AN > 0`` silently disabled this guard on the Parquet path — the production fast path —
    while the report simultaneously claimed gnomAD "does not survey" a locus the store had
    explicitly vouched for.
    """
    a = c.annotation
    return (c.variant.zygosity == "hom") and (a.gnomad_af == 0.0) \
        and (a.gnomad_absence_vouched or bool(a.gnomad_an))


def is_hom_gnomad_uncovered(c) -> bool:
    """Homozygous genotype at a site gnomAD does not cover (AN=0 / absent) — so the
    frequency is UNKNOWN, not zero.

    Worth a caveat (uncovered sites are typically the same difficult regions that generate
    artifacts) but NOT a demotion: a genuine novel recessive allele in a poorly-surveyed
    region looks exactly like this, and no evidence exists either way.

    A store-vouched absence is NOT uncovered, even though it also carries AN=0 — see
    :func:`is_hom_absent_artifact`."""
    a = c.annotation
    return (c.variant.zygosity == "hom") and not a.gnomad_an and not a.gnomad_absence_vouched


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
        # Normalize underscores for the SAME reason clinvar_stars does: raw ClinVar CLNSIG
        # is underscore-delimited ("Likely_pathogenic"). Every ingest path in this repo
        # normalizes, but an externally-built Parquet or a hand-made slice need not — and
        # the identical test in vcf/filter.py guards the rarity/impact RESCUE, so a raw
        # value would be dropped by the filter and invisible to this net at the same time.
        sig = (c.annotation.clinvar_significance or "").lower().replace("_", " ")
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

    # BEFORE any finding: if this run only proceeded because it is a demonstration on a
    # committed fixture, say so first. A reader who stops after one bullet must still have
    # learned that this is not a patient result.
    stamp = demo.stamp_line(report.provenance)
    if stamp:
        lines.append(stamp)

    # Every phenotype claim below is gated on a comparison having actually happened. With no
    # HPO terms (the CLI's --hpo is optional and MCP run_report defaults to None) the scores
    # are all None, and "in phenotype-matched genes" / "not matching the stated phenotype"
    # would describe a comparison that never ran.
    compared = phenotype_compared(report.classifications)
    if not compared:
        lines.append(
            "**No phenotype was available for this analysis** (no HPO terms supplied, or none "
            "of the candidate genes are in the local HPO table), so findings below are NOT "
            "phenotype-filtered and no statement about matching the indication is made. "
            "Supplying HPO terms materially sharpens the prioritisation."
        )

    diag = [c for c in primary if c.tier in _PLP]
    if diag:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in diag)
        where = (" (in a gene overlapping the patient's phenotype)" if compared else
                 " (phenotype overlap NOT assessed)")
        lines.append(f"Likely explanatory finding for the clinical indication: **{g}**"
                     f"{where} — confirm and review.")
    elif not report.classifications and report.qc.total_variants == 0:
        # "No pathogenic finding" is a RESULT. With nothing analysed there is no result — only
        # an absence of input, which is the one thing this report may never dress up as
        # evidence. The QC section already carries a warning, but the conclusion is the
        # bottom-line-up-front line a clinician reads first, and it was announcing a negative.
        #
        # Gated on there being NO classifications, not on the counter alone: a caller may build
        # a report model directly with the QC counters unset, and a report that classified
        # something plainly analysed something whatever the counters say. An existing
        # conclusion test caught exactly that.
        lines.append(
            "**No analysis was performed** — the callset contained no variants to classify "
            f"({report.qc.total_variants} in the VCF, {report.qc.after_qc} after QC). This is "
            "NOT a negative result: nothing was examined. Check that the input is a variant "
            "callset and re-run.")
    else:
        vus = [c for c in primary if c.tier == "Uncertain Significance (VUS)"]
        msg = ("**No Pathogenic / Likely Pathogenic finding** by the engine's independent "
               "ACMG classification" + (" in phenotype-matched genes" if compared else ""))
        if vus:
            msg += f"; {len(vus)} variant(s) of uncertain significance need further evaluation"
        lines.append(msg + ".")

        # A P/LP OUTSIDE the phenotype-matched set must not be silently absent from the
        # headline. Measured: with HPO terms that matched no candidate gene, the conclusion
        # opened with "No Pathogenic / Likely Pathogenic finding" while the report below listed
        # TWO Pathogenic variants. The qualifier "in phenotype-matched genes" is accurate and
        # easy to miss, and an incomplete referral phenotype — the norm, not the exception —
        # is exactly what puts a real finding outside the match.
        elsewhere = [c for c in (secondary + other) if c.tier in _PLP
                     and not is_unconfirmed_ar_carrier(c, report.classifications)]
        if elsewhere and compared:
            g = "; ".join(f"{c.variant.gene or c.variant.key} — {c.tier}" for c in elsewhere)
            lines.append(
                f"**But {len(elsewhere)} Pathogenic / Likely Pathogenic variant(s) were called "
                f"OUTSIDE the phenotype match**: {g}. The phenotype supplied does not overlap "
                "these genes — which may mean the referral phenotype is incomplete rather than "
                "that the variant is irrelevant. Review before concluding the case is negative.")

    # Safety flag: a known ClinVar-Pathogenic variant the QC gate removed BEFORE annotation.
    # The do-not-dismiss net below cannot reach these — they are not classifications at all —
    # so without this line the report deletes a known pathogenic allele and then states there
    # is no finding. Placed immediately after the headline it may be contradicting.
    if report.qc.qc_rescued:
        lines.append(
            f"⚠️ **Removed by QC before classification, but known to ClinVar** — "
            f"{len(report.qc.qc_rescued)} variant(s): "
            + "; ".join(report.qc.qc_rescued)
        )

    # Safety flag: a known, well-reviewed ClinVar-Pathogenic variant must never be hidden
    # behind a lower engine tier (surfaced independently of the ACMG math).
    flagged = clinvar_pathogenic_flags(report.classifications)
    if flagged:
        g = "; ".join(f"{c.variant.gene} ({clinvar_stars(c.annotation.clinvar_review_status)}★; "
                      f"engine: {c.tier})" for c in flagged)
        lines.append(f"⚠️ **Classified Pathogenic/Likely Pathogenic in ClinVar** (≥2-star review) — "
                     f"the engine's independent tier is lower, but DO NOT dismiss: **{g}**. Review the "
                     "ClinVar assertion and its underlying evidence.")

    # gnomAD has no coverage here, so the frequency is UNKNOWN. Not a demotion (see
    # is_hom_gnomad_uncovered) but the reader must not read a blank AF as "absent".
    uncovered = [c for c in report.classifications
                 if is_hom_gnomad_uncovered(c) and c.tier not in _BENIGN]
    if uncovered:
        g = "; ".join(f"{c.variant.gene} — {c.tier}" for c in uncovered)
        lines.append(
            f"**Population frequency unknown, not zero** — {len(uncovered)} homozygous "
            f"variant(s) at sites gnomAD does not survey (AN=0): {g}. Rarity criteria (PM2) rest "
            "on no observation here, and such sites are often the difficult regions that also "
            "generate genotyping artifacts — confirm the genotype and seek an alternative "
            "frequency source before weighing rarity."
        )

    # Every non-benign hom-absent variant, not only the P/LP ones. The caveat used to be
    # gated on _PLP — the same condition that exempts a variant from the demotion — so a
    # phenotype-matched homozygote the engine held at VUS was routed to "Not routinely
    # reported" AND given no explanation. Nothing leaves the diagnostic sections silently.
    artifacts = [c for c in report.classifications
                 if is_hom_absent_artifact(c) and c.tier not in _BENIGN]
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
        # Only claim a non-match where a phenotype score actually existed for that gene.
        why = ("not matching the stated phenotype and not on the ACMG SF actionable list"
               if compared else "not on the ACMG SF actionable list (phenotype not assessed)")
        lines.append(f"Additional **Pathogenic / Likely Pathogenic** variant(s) {why}: {g}. "
                     "Clinical relevance to the indication is uncertain — review in context.")

    # Probable-pathogenic VUS: the engine held these at Uncertain Significance (correctly — the
    # evidence is Supporting-only), but they overlap the phenotype AND carry molecular signal, so
    # they are the VUS most worth a human+model second look. Tier unchanged; this is triage.
    vus = probable_pathogenic_vus(report.classifications)
    if vus:
        g = "; ".join(f"{e['classification'].variant.gene} ("
                      f"{', '.join(s['signal'] for s in e['signals'])})" for e in vus[:3])
        lines.append(f"**Probable-pathogenic VUS for expert review** — phenotype-relevant variant(s) "
                     f"the engine held at Uncertain Significance but which carry suggestive evidence: "
                     f"{g}. Still VUS by ACMG — prioritised for exploration (literature, domain/"
                     "functional context, the ClinVar assertion), which Claude can help work through.")

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
