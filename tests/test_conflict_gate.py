"""The conflicting-evidence gate in `rules.combine`.

Richards 2015 defaults a variant to Uncertain Significance when "the evidence for benign and
pathogenic is conflicting" — EVIDENCE, not rules. The gate used to test whether a *rule* had
fired on each side, and Table 5 has no rule for a lone Strong-benign line (Likely Benign needs
1 Strong AND 1 Supporting) nor for a lone PVS1. Both therefore vanished: a nonsense variant
present in 2% of the population reached **Pathogenic** with "+ BS1" printed in the audit trail
as though it had been weighed.

The replacement is deliberately narrower than "met on both sides" — see the calibration tests
at the bottom, which are the ones that fail if the gate is widened.
"""
from __future__ import annotations

import pytest

from vcf2report.acmg import rules
from vcf2report.acmg.engine import classify
from vcf2report.models import Annotation, CriterionResult, Variant

_STRENGTH = {
    "PVS1": "very_strong", "PS1": "strong", "PS3": "strong", "PM2": "supporting",
    "PM4": "moderate", "PM5": "moderate", "PP2": "supporting", "PP3": "supporting",
    "PP4": "supporting", "PP5": "supporting",
    "BA1": "stand_alone", "BS1": "strong", "BS2": "strong",
    "BP4": "supporting", "BP6": "supporting", "BP7": "supporting",
}


def _combine(*codes):
    crs = [CriterionResult(code=c, name=c, applies=True, met=True,
                           default_strength=_STRENGTH[c], applied_strength=_STRENGTH[c])
           for c in codes]
    return rules.combine(crs)


# ---------------------------------------------------------------------------------------
# The gate must fire: Strong-or-above evidence on the losing side is never discarded.
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("codes", [
    ("PVS1", "PP3", "PP4", "BS1"),      # a 2% allele: BS1 alone reaches no Table 5 rule
    ("PVS1", "PP3", "PP4", "BS2"),      # homozygous in healthy adults
    ("PS1", "PS3", "BS1"),              # 2 Strong pathogenic vs 1 Strong benign
])
def test_strong_benign_evidence_cannot_be_discarded_into_a_pathogenic_call(codes):
    tier, path = _combine(*codes)
    assert tier == rules.VUS
    assert "conflicting" in path


@pytest.mark.parametrize("codes", [
    ("PVS1", "BS2", "BP6"),             # LB-1 fires; PVS1 reaches no rule alone
    ("PVS1", "BA1"),                    # BEN-1 fires; PVS1 reaches no rule alone
    ("PVS1", "BP4", "BP6"),             # LB-2 fires
])
def test_strong_pathogenic_evidence_cannot_be_discarded_into_a_benign_call(codes):
    tier, path = _combine(*codes)
    assert tier == rules.VUS
    assert "conflicting" in path


def test_a_common_null_variant_is_not_issued_as_a_pathogenic_finding():
    """End-to-end through the engine, and through the report's routing."""
    from vcf2report.report.assemble import split_findings
    v = Variant(chrom="14", pos=23400000, ref="C", alt="T", gene="MYH7",
                consequence="stop_gained", zygosity="het")
    a = Annotation(gnomad_af=0.020, gnomad_ac=2900, gnomad_an=145000, gnomad_faf95=0.019,
                   gene_lof_intolerant=True, cadd_phred=35.0,
                   hpo_match_score=0.90, hpo_best_match=0.90)
    c = classify(v, a)
    assert "BS1" in [x.code for x in c.criteria if x.applies and x.met]
    assert c.tier == rules.VUS, f"a 2% allele classified {c.tier}: {c.rule_path}"
    assert "conflicting" in c.rule_path


def test_the_audit_trail_never_shows_evidence_the_verdict_ignored():
    """The trail is built from every met criterion. If the tier came from one side only while
    the other held Strong evidence, the rule label and the criteria table beside it disagree."""
    tier, path = _combine("PVS1", "BS2", "BP6")
    assert "PVS1" in path                       # the trail still names it...
    assert "LB-1" not in path and "Likely Benign" not in path   # ...and the verdict reflects it
    assert tier == rules.VUS


# ---------------------------------------------------------------------------------------
# Calibration: the gate must NOT fire. These are what break if it is widened to
# "any met criterion on both sides".
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("codes,expected", [
    # "rare" next to "predicted benign" is not a contradiction — a rare variant can be benign.
    (("PM2", "BP4", "BP6"), rules.LIKELY_BENIGN),
    (("PM2", "BS2", "BP4"), rules.LIKELY_BENIGN),
    (("PM4", "BS1", "BS2"), rules.BENIGN),
    # A single Supporting benign line must not veto a call built on two Strong criteria.
    (("PS1", "PS3", "BP4"), rules.PATHOGENIC),
])
def test_weak_opposing_evidence_does_not_force_a_conflict(codes, expected):
    tier, path = _combine(*codes)
    assert tier == expected and "conflicting" not in path


@pytest.mark.parametrize("codes,expected", [
    (("BA1",), rules.BENIGN),
    (("BS1", "BS2"), rules.BENIGN),
    (("BS1", "BP4"), rules.LIKELY_BENIGN),
    (("PVS1", "PM2"), rules.LIKELY_PATHOGENIC),
    (("PS1", "PS3"), rules.PATHOGENIC),
    (("PVS1", "PS1"), rules.PATHOGENIC),
])
def test_single_sided_evidence_is_untouched(codes, expected):
    """No opposing evidence at all — the gate must be invisible."""
    tier, path = _combine(*codes)
    assert tier == expected and "conflicting" not in path


def test_the_gate_is_symmetric():
    """Same two criteria, mirrored strengths: neither side gets to win by rule-shape accident."""
    assert _combine("PVS1", "BS1")[0] == rules.VUS
    assert _combine("PS1", "BA1")[0] == rules.VUS


def test_committed_exomes_flag_their_mixed_evidence_calls():
    """Live regression: 64 of 645 classifications on the committed exomes carry evidence on
    both sides and ZERO were flagged. The ones holding Strong evidence must now be."""
    import glob

    from vcf2report.pipeline import run_pipeline
    unflagged_strong = []
    for f in sorted(glob.glob("data/example/*.annotated.vcf.gz")):
        for c in run_pipeline(f).classifications:
            met = [(x.code, x.applied_strength or x.default_strength)
                   for x in c.criteria if x.applies and x.met]
            path = [s for code, s in met if code not in rules.BENIGN_CODES]
            ben = [s for code, s in met if code in rules.BENIGN_CODES]
            if not (path and ben):
                continue
            losing_strong = (c.tier in (rules.BENIGN, rules.LIKELY_BENIGN)
                             and any(s in ("very_strong", "strong") for s in path)) or \
                            (c.tier in (rules.PATHOGENIC, rules.LIKELY_PATHOGENIC)
                             and any(s in ("stand_alone", "strong") for s in ben))
            if losing_strong:
                unflagged_strong.append((c.variant.gene, c.tier, c.rule_path))
    assert not unflagged_strong, f"decisive evidence discarded: {unflagged_strong}"
