"""ACMG/AMP combining rules — Richards et al. 2015, Table 5.

Pure boolean logic over the counts of *met* criteria by strength, so the whole
thing is inspectable and unit-testable. Returns the 5-tier call plus a short
label of which rule fired (used verbatim in the report for auditability).
"""
from __future__ import annotations

from typing import Iterable

from ..models import CriterionResult

PATHOGENIC = "Pathogenic"
LIKELY_PATHOGENIC = "Likely Pathogenic"
VUS = "Uncertain Significance (VUS)"
LIKELY_BENIGN = "Likely Benign"
BENIGN = "Benign"

# strength -> which side it counts for
_PATHOGENIC_STRENGTHS = {"very_strong": "PVS", "strong": "PS", "moderate": "PM", "supporting": "PP"}
_BENIGN_STRENGTHS = {"stand_alone": "BA", "strong": "BS", "supporting": "BP"}

# Codes that count as benign evidence (everything else met counts pathogenic).
_BENIGN_CODES = {"BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"}


def _counts(criteria: Iterable[CriterionResult]) -> dict[str, int]:
    c = {"PVS": 0, "PS": 0, "PM": 0, "PP": 0, "BA": 0, "BS": 0, "BP": 0}
    for cr in criteria:
        if not (cr.applies and cr.met):
            continue
        strength = cr.applied_strength or cr.default_strength
        if cr.code in _BENIGN_CODES:
            bucket = _BENIGN_STRENGTHS.get(strength)
        else:
            bucket = _PATHOGENIC_STRENGTHS.get(strength)
        if bucket:
            c[bucket] += 1
    return c


def _pathogenic_rule(c: dict[str, int]) -> str | None:
    pvs, ps, pm, pp = c["PVS"], c["PS"], c["PM"], c["PP"]
    # Pathogenic
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm == 1 and pp >= 1) or pp >= 2):
        return "PATH-1 (PVS1 + strong/moderate/supporting)"
    if ps >= 2:
        return "PATH-2 (>=2 Strong)"
    if ps == 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm == 1 and pp >= 4)):
        return "PATH-3 (1 Strong + Moderate/Supporting)"
    return None


def _likely_pathogenic_rule(c: dict[str, int]) -> str | None:
    pvs, ps, pm, pp = c["PVS"], c["PS"], c["PM"], c["PP"]
    # PVS1 (very strong) + at least one other criterion -> Likely Pathogenic. On the ClinGen/Tavtigian
    # scale PVS1(8) + any Supporting(1) = 9 = LP, but the 2015 Table-5 only had "PVS1 + 1 Moderate";
    # once PM2 is Supporting (the default) a novel truncation (PVS1 + PM2) would wrongly fall to VUS.
    # This floors it to LP. It requires a SECOND criterion on purpose: PVS1 ALONE (e.g. a null variant
    # that is present in gnomAD, so PM2 does not fire) stays VUS — which is what keeps incidental het
    # LoF from flooding a healthy exome.
    if pvs >= 1 and (pm >= 1 or pp >= 1):
        return "LP-1 (PVS1 + Moderate)" if pm >= 1 else "LP-1 (PVS1 + Supporting)"
    if ps == 1 and 1 <= pm <= 2:
        return "LP-2 (1 Strong + 1-2 Moderate)"
    if ps == 1 and pp >= 2:
        return "LP-3 (1 Strong + >=2 Supporting)"
    if pm >= 3:
        return "LP-4 (>=3 Moderate)"
    if pm >= 2 and pp >= 2:
        return "LP-5 (2 Moderate + >=2 Supporting)"
    if pm == 1 and pp >= 4:
        return "LP-6 (1 Moderate + >=4 Supporting)"
    return None


def _benign_rule(c: dict[str, int]) -> str | None:
    if c["BA"] >= 1:
        return "BEN-1 (BA1 stand-alone)"
    if c["BS"] >= 2:
        return "BEN-2 (>=2 Strong benign)"
    return None


def _likely_benign_rule(c: dict[str, int]) -> str | None:
    if c["BS"] == 1 and c["BP"] >= 1:
        return "LB-1 (1 Strong + 1 Supporting benign)"
    if c["BP"] >= 2:
        return "LB-2 (>=2 Supporting benign)"
    return None


# ClinGen/Tavtigian (2020) naturally-scaled points: each strength is a power-of-two
# so the categorical Richards rules fall out of a single additive score. Benign
# evidence is negative. Thresholds (Tavtigian 2020, adopted by ClinGen SVI):
# Pathogenic >=10, Likely Pathogenic 6..9, VUS 0..5, Likely Benign -1..-6, Benign <=-7.
_PATH_POINTS = {"very_strong": 8, "strong": 4, "moderate": 2, "supporting": 1}
_BENIGN_POINTS = {"stand_alone": -8, "strong": -4, "moderate": -2, "supporting": -1}


def _combine_points(criteria: list[CriterionResult]) -> tuple[str, str]:
    """ClinGen/Tavtigian points model (used when VCF2REPORT_ACMG_MODEL=clingen)."""
    total = 0
    for cr in criteria:
        if not (cr.applies and cr.met):
            continue
        strength = cr.applied_strength or cr.default_strength
        table = _BENIGN_POINTS if cr.code in _BENIGN_CODES else _PATH_POINTS
        total += table.get(strength, 0)

    if total >= 10:
        tier = PATHOGENIC
    elif total >= 6:
        tier = LIKELY_PATHOGENIC
    elif total >= 0:
        tier = VUS
    elif total >= -6:
        tier = LIKELY_BENIGN
    else:
        tier = BENIGN

    met = [cr.code for cr in criteria if cr.applies and cr.met]
    trail = " + ".join(met) if met else "no criteria met"
    return tier, f"{trail} => {total:+d} points (ClinGen/Tavtigian) => {tier}"


def combine(criteria: list[CriterionResult]) -> tuple[str, str]:
    """Return (tier, rule_path).

    Uses Richards 2015 Table 5 by default; the ClinGen/Tavtigian points model when
    ``VCF2REPORT_ACMG_MODEL=clingen``. Under Richards, if pathogenic and benign
    evidence both fire the result is reported as VUS (conflicting).
    """
    from .. import config
    if config.acmg_model() == "clingen":
        return _combine_points(criteria)

    c = _counts(criteria)
    path = _pathogenic_rule(c) or _likely_pathogenic_rule(c)
    benign = _benign_rule(c) or _likely_benign_rule(c)

    met = [cr.code for cr in criteria if cr.applies and cr.met]
    trail = " + ".join(met) if met else "no criteria met"

    if path and benign:
        return VUS, f"{trail} => conflicting pathogenic & benign evidence => VUS"

    if path:
        tier = PATHOGENIC if _pathogenic_rule(c) else LIKELY_PATHOGENIC
        return tier, f"{trail} => {tier} [{path}]"

    if benign:
        tier = BENIGN if _benign_rule(c) else LIKELY_BENIGN
        return tier, f"{trail} => {tier} [{benign}]"

    return VUS, f"{trail} => criteria insufficient for a benign or pathogenic call => VUS"
