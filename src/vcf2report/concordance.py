"""Concordance panel: validate the engine against ClinVar ground truth.

The panel answers one question honestly: *when a real, expert-classified ClinVar
variant is fed through our ACMG engine, does the engine reach the same broad
call?* It is the quantitative counterpart to the per-variant audit trail — a
single number (plus a confusion matrix) for "how often do we agree with the
experts, and do we ever get it dangerously wrong?".

**Non-circularity.** The engine consumes ClinVar directly via PP5 (a supporting
line that fires when ClinVar already says Pathogenic). Comparing the engine to
ClinVar *while feeding it ClinVar* would be circular and inflate concordance for
free. So the panel evaluates every variant with ClinVar **withheld** from the
annotation (``withhold_clinvar=True``): the engine must recover the call from the
population-frequency, LoF-mechanic and in-silico axes alone. ClinVar is used only
as the answer key.

**What v1 can and cannot reach.** With ClinVar withheld and no in-silico scores
bundled, a pathogenic *missense* variant earns only PM2 (rare) and lands in VUS —
the engine is designed to defer missense pathogenicity to in-silico/hotspot
(``adjudicated_by="model"``) evidence it does not have here. The deterministic
axes the panel does exercise are **LoF pathogenicity** (PVS1 + PM2 -> LP/P) and
**common-benign** (BA1 / BS1). Missense-heavy pathogenic recovery is reported as a
known limitation, not hidden.

The panel data (``ground_truth.tsv`` + ``gnomad_frozen.json``) is produced once,
with network, by ``scripts/build_concordance_panel.py``; thereafter this module —
and the test that guards it — runs fully offline.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from . import config
from .acmg.engine import classify
from .acmg.rules import (
    BENIGN,
    LIKELY_BENIGN,
    LIKELY_PATHOGENIC,
    PATHOGENIC,
    VUS,
)
from .annotate import abraom, extra
from .models import Annotation, Classification, Variant

# ---------------------------------------------------------------------------
# Panel data locations (produced by scripts/build_concordance_panel.py).
# ---------------------------------------------------------------------------
CONCORDANCE_DIR = config.DATA_DIR / "concordance"
GROUND_TRUTH = CONCORDANCE_DIR / "ground_truth.tsv"
FROZEN_GNOMAD = CONCORDANCE_DIR / "gnomad_frozen.json"
FROZEN_ALPHAMISSENSE = CONCORDANCE_DIR / "alphamissense_frozen.json"

_GT_COLUMNS = [
    "key", "gene", "consequence", "hgvs_p",
    "clinvar_significance", "review_status", "accession", "condition",
]

# The three collapsed classes the panel scores on.
PATH, UNCERTAIN, BEN = "PATH", "VUS", "BEN"
CLASSES = (PATH, UNCERTAIN, BEN)


# ---------------------------------------------------------------------------
# Collapsing 5-tier calls / ClinVar text to the 3 comparison classes
# ---------------------------------------------------------------------------
def collapse_engine_tier(tier: str) -> str:
    """Map a 5-tier engine call to PATH / VUS / BEN."""
    if tier in (PATHOGENIC, LIKELY_PATHOGENIC):
        return PATH
    if tier in (BENIGN, LIKELY_BENIGN):
        return BEN
    return UNCERTAIN


def collapse_clinvar(significance: Optional[str]) -> Optional[str]:
    """Map a ClinVar significance string to PATH / VUS / BEN, or None.

    Conflicting / not-provided / unmappable strings collapse to VUS. A missing or
    empty significance returns None so the caller can drop the variant from the
    answer key rather than scoring it as an uncertain truth.
    """
    sig = (significance or "").strip().lower()
    if not sig:
        return None
    if sig.startswith("pathogenic") or sig.startswith("likely pathogenic"):
        return PATH
    if sig.startswith("benign") or sig.startswith("likely benign"):
        return BEN
    return UNCERTAIN


# ---------------------------------------------------------------------------
# Panel model
# ---------------------------------------------------------------------------
@dataclass
class PanelEntry:
    """One ground-truth variant plus its frozen gnomAD frequencies."""

    variant: Variant
    truth_significance: Optional[str]
    truth_class: Optional[str]
    review_status: Optional[str] = None
    accession: Optional[str] = None
    frozen_gnomad: dict = field(default_factory=dict)
    frozen_alphamissense: dict = field(default_factory=dict)


@dataclass
class PanelRow:
    """The engine's call for one panel variant, alongside the truth."""

    key: str
    gene: Optional[str]
    consequence: Optional[str]
    is_lof: bool
    truth_class: str
    truth_significance: Optional[str]
    engine_tier: str
    engine_class: str
    rule_path: str
    met_codes: list[str]
    concordant: bool
    gross: bool


@dataclass
class ConcordanceResult:
    """Confusion matrix + headline metrics over the whole panel."""

    rows: list[PanelRow]
    matrix: dict[str, dict[str, int]]
    metrics: dict[str, float]

    @property
    def gross_discordances(self) -> list[PanelRow]:
        return [r for r in self.rows if r.gross]

    def to_dict(self) -> dict:
        return {
            "metrics": self.metrics,
            "matrix": self.matrix,
            "rows": [vars(r) for r in self.rows],
        }

    def to_markdown(self) -> str:
        return _render_markdown(self)


# ---------------------------------------------------------------------------
# Building the engine input from frozen data (ClinVar withheld by default)
# ---------------------------------------------------------------------------
def _variant_from_row(row: dict) -> Variant:
    chrom, pos, ref, alt = row["key"].split("-")
    return Variant(
        chrom=chrom, pos=int(pos), ref=ref, alt=alt,
        gene=row.get("gene") or None,
        hgvs_p=row.get("hgvs_p") or None,
        consequence=row.get("consequence") or None,
    )


def _annotation_from_frozen(
    entry: PanelEntry, withhold_clinvar: bool
) -> Annotation:
    """Assemble an :class:`Annotation` for a panel variant, fully offline.

    gnomAD comes from the frozen snapshot; ABraOM and gene constraint from the
    bundled local datasets. ClinVar is withheld (all None) when ``withhold_clinvar``
    so PP5 cannot fire and the concordance signal stays non-circular. No in-silico
    scores are attached in v1 (missense pathogenicity is left to model adjudication).
    """
    g = entry.frozen_gnomad or {}
    am = entry.frozen_alphamissense or {}
    v = entry.variant
    ab = abraom.lookup(v)
    con = extra.gene_constraint(v.gene)

    if withhold_clinvar:
        cv_sig = cv_review = cv_acc = None
    else:
        cv_sig, cv_review, cv_acc = (
            entry.truth_significance, entry.review_status, entry.accession)

    return Annotation(
        clinvar_significance=cv_sig,
        clinvar_review_status=cv_review,
        clinvar_accession=cv_acc,
        gnomad_af=g.get("af"),
        gnomad_ac=g.get("ac"),
        gnomad_an=g.get("an"),
        gnomad_homozygotes=g.get("hom"),
        gnomad_popmax_pop=g.get("pop"),
        gnomad_faf95=g.get("faf95"),
        abraom_af=ab.get("af"),
        gene_lof_intolerant=con.get("lof_intolerant"),
        revel=None,
        cadd_phred=None,
        am_pathogenicity=am.get("am_pathogenicity"),
        am_class=am.get("am_class"),
        hpo_match_score=None,
        source={
            "gnomad": f"gnomAD v{g.get('release', '4.1')} (frozen panel)",
            "abraom": ab.get("_source", ""),
            "gene_lof_intolerant": con.get("_source", ""),
            "alphamissense": "AlphaMissense (frozen panel)" if am.get("am_pathogenicity") is not None
            else "AlphaMissense (no score)",
            "clinvar": "withheld (concordance panel)" if withhold_clinvar
            else "ClinVar (panel truth)",
        },
    )


def classify_entry(entry: PanelEntry, withhold_clinvar: bool = True) -> Classification:
    """Run the full ACMG engine on one panel variant using frozen inputs."""
    annotation = _annotation_from_frozen(entry, withhold_clinvar)
    return classify(entry.variant, annotation)


# ---------------------------------------------------------------------------
# Loading + evaluating the panel
# ---------------------------------------------------------------------------
def load_panel(
    ground_truth: str | Path = GROUND_TRUTH,
    frozen_gnomad: str | Path = FROZEN_GNOMAD,
    frozen_alphamissense: str | Path = FROZEN_ALPHAMISSENSE,
) -> list[PanelEntry]:
    """Load ground-truth variants + frozen gnomAD (+ AlphaMissense) into entries.

    Only variants whose ClinVar significance collapses to a definite PATH or BEN
    truth are kept in the answer key (an uncertain/conflicting truth is not a
    usable label). Rows missing a frozen record keep an empty snapshot, which the
    engine reads as "unavailable" (never a fabricated AF 0 / score 0). The
    AlphaMissense frozen file is optional — absent means v1 (frequency-only)
    behaviour, present means the calibrated in-silico axis is active.
    """
    gt_path, frozen_path = Path(ground_truth), Path(frozen_gnomad)
    frozen: dict = {}
    if frozen_path.exists():
        frozen = json.loads(frozen_path.read_text())
    am_path = Path(frozen_alphamissense)
    am_frozen: dict = json.loads(am_path.read_text()) if am_path.exists() else {}

    entries: list[PanelEntry] = []
    for line in gt_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        row = dict(zip(_GT_COLUMNS, parts))
        if not row.get("key"):
            continue
        truth_class = collapse_clinvar(row.get("clinvar_significance"))
        if truth_class not in (PATH, BEN):
            continue
        entries.append(PanelEntry(
            variant=_variant_from_row(row),
            truth_significance=row.get("clinvar_significance"),
            truth_class=truth_class,
            review_status=row.get("review_status") or None,
            accession=row.get("accession") or None,
            frozen_gnomad=frozen.get(row["key"], {}),
            frozen_alphamissense=am_frozen.get(row["key"], {}),
        ))
    return entries


def _empty_matrix() -> dict[str, dict[str, int]]:
    return {t: {e: 0 for e in CLASSES} for t in CLASSES}


def evaluate_panel(
    entries: Iterable[PanelEntry], withhold_clinvar: bool = True
) -> ConcordanceResult:
    """Classify every entry and reduce to a confusion matrix + metrics.

    A **gross discordance** is the clinically dangerous case: the engine calls a
    ClinVar-benign variant Pathogenic, or a ClinVar-pathogenic variant Benign.
    These must be zero — a VUS on either side is a conservative non-call, a
    PATH<->BEN flip is a real error.
    """
    rows: list[PanelRow] = []
    matrix = _empty_matrix()

    for entry in entries:
        result = classify_entry(entry, withhold_clinvar=withhold_clinvar)
        engine_class = collapse_engine_tier(result.tier)
        truth = entry.truth_class
        matrix[truth][engine_class] += 1
        gross = (truth == PATH and engine_class == BEN) or (
            truth == BEN and engine_class == PATH)
        rows.append(PanelRow(
            key=entry.variant.key,
            gene=entry.variant.gene,
            consequence=entry.variant.consequence,
            is_lof=entry.variant.is_lof,
            truth_class=truth,
            truth_significance=entry.truth_significance,
            engine_tier=result.tier,
            engine_class=engine_class,
            rule_path=result.rule_path,
            met_codes=result.met_codes,
            concordant=(engine_class == truth),
            gross=gross,
        ))

    return ConcordanceResult(rows=rows, matrix=matrix,
                             metrics=_metrics(rows, matrix))


def _rate(numerator: int, denominator: int) -> float:
    return round(numerator / denominator, 4) if denominator else 0.0


def _metrics(rows: list[PanelRow], matrix: dict[str, dict[str, int]]) -> dict:
    total = len(rows)
    truth_path = [r for r in rows if r.truth_class == PATH]
    truth_ben = [r for r in rows if r.truth_class == BEN]
    lof_path = [r for r in truth_path if r.is_lof]

    concordant = sum(1 for r in rows if r.concordant)
    gross = sum(1 for r in rows if r.gross)

    # The engine defers to VUS unless it has strong deterministic evidence, so a
    # recall-style "concordance" understates it. The honest pair is: how OFTEN does
    # it commit (decisiveness), and how RIGHT is it WHEN it commits (precision)?
    decisive = [r for r in rows if r.engine_class != UNCERTAIN]
    engine_path = [r for r in rows if r.engine_class == PATH]
    engine_ben = [r for r in rows if r.engine_class == BEN]

    return {
        "n": total,
        "n_pathogenic": len(truth_path),
        "n_benign": len(truth_ben),
        "n_lof_pathogenic": len(lof_path),
        # Overall agreement on the collapsed 3-class call (recall-flavoured).
        "concordance": _rate(concordant, total),
        # How often the engine commits to a non-VUS call at all.
        "decisiveness": _rate(len(decisive), total),
        # Of the calls the engine DID commit to, how many match ClinVar.
        "concordance_when_decisive": _rate(
            sum(1 for r in decisive if r.concordant), len(decisive)),
        # Of the engine's PATH calls, how many are truly ClinVar-pathogenic.
        "pathogenic_precision": _rate(
            sum(1 for r in engine_path if r.truth_class == PATH), len(engine_path)),
        # Of the engine's BEN calls, how many are truly ClinVar-benign.
        "benign_precision": _rate(
            sum(1 for r in engine_ben if r.truth_class == BEN), len(engine_ben)),
        # Of ClinVar-pathogenic variants, how many the engine also calls PATH.
        "pathogenic_sensitivity": _rate(
            sum(1 for r in truth_path if r.engine_class == PATH), len(truth_path)),
        # Restricted to the deterministic axis v1 actually exercises (LoF).
        "lof_pathogenic_sensitivity": _rate(
            sum(1 for r in lof_path if r.engine_class == PATH), len(lof_path)),
        # Of ClinVar-benign variants, how many the engine also calls BEN.
        "benign_agreement": _rate(
            sum(1 for r in truth_ben if r.engine_class == BEN), len(truth_ben)),
        # Of ClinVar-benign variants, how many the engine does NOT wrongly call PATH.
        "benign_specificity": _rate(
            sum(1 for r in truth_ben if r.engine_class != PATH), len(truth_ben)),
        # The clinical safety number — must be zero.
        "gross_discordances": gross,
        "gross_discordance_rate": _rate(gross, total),
    }


# ---------------------------------------------------------------------------
# Markdown rendering (the "panel")
# ---------------------------------------------------------------------------
def _render_markdown(res: ConcordanceResult) -> str:
    m = res.metrics
    lines: list[str] = []
    lines.append("# Concordance panel — ClinVar ground truth vs vcf2report engine")
    lines.append("")
    lines.append(f"- Variants scored: **{m['n']}** "
                 f"({m['n_pathogenic']} pathogenic, {m['n_benign']} benign; "
                 f"{m['n_lof_pathogenic']} of the pathogenic are LoF)")
    lines.append(f"- ClinVar withheld from the engine (non-circular)")
    lines.append("")
    lines.append("## Confusion matrix (truth rows x engine columns)")
    lines.append("")
    lines.append("| truth \\ engine | PATH | VUS | BEN | total |")
    lines.append("|---|---:|---:|---:|---:|")
    for t in CLASSES:
        r = res.matrix[t]
        total = r[PATH] + r[UNCERTAIN] + r[BEN]
        lines.append(f"| **{t}** | {r[PATH]} | {r[UNCERTAIN]} | {r[BEN]} | {total} |")
    lines.append("")
    lines.append("## Safety (the numbers that must hold)")
    lines.append("")
    lines.append(f"- **Gross discordances (PATH<->BEN flips):** "
                 f"**{m['gross_discordances']}** — must be 0")
    lines.append(f"- **Pathogenic precision (engine PATH → truly pathogenic):** "
                 f"{m['pathogenic_precision']:.1%}")
    lines.append(f"- **Benign precision (engine BEN → truly benign):** "
                 f"{m['benign_precision']:.1%}")
    lines.append("")
    lines.append("## Behaviour (a conservative engine defers to VUS)")
    lines.append("")
    lines.append(f"- **Decisiveness (non-VUS calls):** {m['decisiveness']:.1%}")
    lines.append(f"- **Concordance WHEN decisive:** {m['concordance_when_decisive']:.1%}")
    lines.append(f"- **Pathogenic sensitivity:** {m['pathogenic_sensitivity']:.1%} "
                 f"(LoF-only: {m['lof_pathogenic_sensitivity']:.1%})")
    lines.append(f"- **Benign agreement:** {m['benign_agreement']:.1%}")
    lines.append(f"- **Overall concordance (3-class, recall-flavoured):** {m['concordance']:.1%}")
    lines.append("")
    if res.gross_discordances:
        lines.append("### ⚠ Gross discordances")
        lines.append("")
        lines.append("| variant | gene | truth | engine | rule path |")
        lines.append("|---|---|---|---|---|")
        for r in res.gross_discordances:
            lines.append(f"| {r.key} | {r.gene or ''} | {r.truth_class} | "
                         f"{r.engine_tier} | {r.rule_path} |")
        lines.append("")
    return "\n".join(lines)
