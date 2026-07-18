"""Triage the VUS pile — surface the *probable-pathogenic* ones for expert review.

The engine is deliberately conservative: a variant with strong AlphaMissense, a perfect
phenotype match and a ClinVar Pathogenic assertion can still land at **Uncertain
Significance** because each of those is only *Supporting* evidence and the ACMG point model
does not reach Likely-Pathogenic on Supporting alone (and PP5 is capped for anti-circularity).
That call is correct — but a VUS carrying several suggestive signals is NOT the same as a VUS
carrying none, and a clinician should look at the first kind first.

This module does NOT change any ACMG tier. It ranks the VUS candidates by how much
suggestive-pathogenic evidence they carry, itemises *why* (so the operator sees the reasoning,
not a black-box score), and hands the top ones back to be surfaced for expert exploration —
the deterministic engine having done its part, this is where human + model judgement takes over.
"""
from __future__ import annotations

from .. import config

_VUS = "Uncertain Significance (VUS)"


def _signals(c) -> list[dict]:
    """The suggestive-pathogenic signals a VUS carries, each itemised for the operator."""
    a = c.annotation
    out: list[dict] = []

    am = a.am_pathogenicity
    if am is not None:
        if am >= config.AM_PP3_SUPPORTING:
            out.append({"signal": "AlphaMissense likely-pathogenic", "value": round(am, 3),
                        "weight": 3, "note": f"AM {am:.3f} ≥ {config.AM_PP3_SUPPORTING} (PP3 threshold)"})
        elif am >= config.AM_BP4_SUPPORTING:
            out.append({"signal": "AlphaMissense ambiguous", "value": round(am, 3),
                        "weight": 1, "note": f"AM {am:.3f} in the ambiguous band — leans deleterious"})

    hpo = a.hpo_match_score
    if hpo is not None and hpo >= config.HPO_RELATED_MIN:
        out.append({"signal": "phenotype match", "value": round(hpo, 3),
                    "weight": 2, "note": f"gene↔patient HPO match {hpo:.2f} ≥ {config.HPO_RELATED_MIN}"})

    sig = (a.clinvar_significance or "").lower()
    if "pathogenic" in sig and "conflict" not in sig:
        out.append({"signal": "ClinVar Pathogenic assertion", "value": a.clinvar_significance,
                    "weight": 2, "note": "not counted at tier (capped for anti-circularity) — "
                                         "review the assertion + its evidence directly"})
    elif "conflict" in sig:
        out.append({"signal": "ClinVar conflicting", "value": a.clinvar_significance,
                    "weight": 1, "note": "conflicting interpretations — expert adjudication needed"})

    if a.gene_lof_intolerant:
        out.append({"signal": "LoF-intolerant gene", "value": True,
                    "weight": 1, "note": "constraint suggests the gene is dosage-sensitive"})

    return out


def probable_pathogenic_vus(classifications, min_molecular: int = 2):
    """Phenotype-relevant VUS carrying suggestive-pathogenic molecular evidence — ranked for review.

    Returns {classification, signals, score}, most-suggestive first. ``score`` is the summed signal
    weight — an ordering aid, deliberately NOT a probability and NOT an ACMG tier.

    Two gates keep this SPECIFIC to the indication (a whole exome carries dozens of incidental VUS —
    surfacing them all would be noise, not triage):
      1. **phenotype-relevant** — the gene must overlap the patient's HPO (>= HPO_RELATED_MIN). An
         incidental VUS in an unrelated gene, however deleterious it looks, is not what the operator
         is triaging FOR this patient's indication.
      2. **molecular support** — beyond the phenotype, at least ``min_molecular`` weight of
         molecular signal (AlphaMissense / ClinVar / constraint), so a bare phenotype overlap on a
         benign-looking variant does not qualify.
    """
    ranked = []
    for c in classifications:
        if c.tier != _VUS:
            continue
        sig = _signals(c)
        names = {s["signal"] for s in sig}
        if "phenotype match" not in names:
            continue
        molecular = sum(s["weight"] for s in sig if s["signal"] != "phenotype match")
        if molecular < min_molecular:
            continue
        ranked.append({"classification": c, "signals": sig, "score": sum(s["weight"] for s in sig)})
    ranked.sort(key=lambda r: r["score"], reverse=True)
    return ranked


def exploration_prompt(entry) -> str:
    """A one-line, operator-facing 'why look here + what to check' for a prioritised VUS."""
    c = entry["classification"]
    v = c.variant
    names = ", ".join(s["signal"] for s in entry["signals"])
    hgvs = v.hgvs_p or v.hgvs_c or v.key
    return (f"{v.gene} {hgvs} — VUS with {names}. Worth expert exploration "
            f"(literature on this residue/gene, domain/functional context, splicing prediction, "
            f"and the ClinVar assertion's underlying evidence).")
