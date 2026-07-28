"""Coverage of the ACMG/AMP criterion set, and the doc that describes it.

The report claims a "fully auditable ACMG trail — all 28 criteria shown", and the laudo's own
limitations section names PP1 and BS4 as reported N/A. Both were false: only 22 criteria were
registered, so six (PP1, BS3, BS4, BP2, BP3, BP5) never appeared in any trail at all.

These tests make the claim enforceable — the registry must hold exactly the 28 criteria of Richards
et al. 2015, the engine must emit all of them in canonical order, and docs/ACMG_CRITERIA.md must
agree with the code about how each one is decided.
"""
import re
from pathlib import Path

import pytest

from vcf2report.acmg.criteria import all_criteria
from vcf2report.acmg.engine import evaluate_criteria
from vcf2report.models import Annotation, Variant

# Richards et al., Genet Med 2015 — the complete set, in the canonical reporting order.
ACMG_28 = [
    "PVS1", "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1", "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
]
DOC = Path(__file__).resolve().parent.parent / "docs" / "ACMG_CRITERIA.md"


def _v():
    return Variant(chrom="1", pos=100, ref="C", alt="T", gene="TESTG",
                   consequence="missense_variant", hgvs_p="p.Arg123Cys")


def _results():
    return evaluate_criteria(_v(), Annotation())


def _kind(cr) -> str:
    return "N/A" if not cr.applies else cr.adjudicated_by


def test_every_acmg_criterion_is_registered():
    missing = [c for c in ACMG_28 if c not in all_criteria()]
    assert not missing, f"criteria claimed by the report but not implemented: {missing}"


def test_no_criteria_beyond_the_acmg_set():
    extra = sorted(set(all_criteria()) - set(ACMG_28))
    assert not extra, f"unknown criteria in the registry: {extra}"


def test_engine_emits_all_28_in_canonical_order():
    assert [c.code for c in _results()] == ACMG_28


def test_every_criterion_carries_a_reason():
    # A trail entry with no reasoning is unauditable — the whole point of showing not-met criteria.
    for cr in _results():
        assert (cr.reasoning or "").strip(), f"{cr.code} has no reasoning"


def test_na_criteria_are_exactly_the_single_proband_undecidables():
    # These need trio / segregation / phasing data a single proband VCF cannot supply. The laudo's
    # limitations section names them, so the set must match what the engine actually reports.
    na = {c.code for c in _results() if not c.applies}
    assert na == {"PS2", "PM3", "PM6", "PP1", "BP2", "BS4"}


def test_model_adjudicated_criteria_are_the_judgement_calls():
    # Functional-assay literature (PS3/BS3), case-control data (PS4), a repeat/domain track (BP3),
    # and whole-case context (BP5). Everything else is decided from data.
    model = {c.code for c in _results() if c.applies and c.adjudicated_by == "model"}
    assert model == {"PS3", "PS4", "BS3", "BP3", "BP5"}


def test_model_and_na_criteria_never_default_to_met():
    for cr in _results():
        if _kind(cr) != "engine":
            assert cr.met is False, f"{cr.code} defaulted to met without engine evidence"


# --- the doc must not drift from the code -----------------------------------
def _doc_rows() -> dict[str, str]:
    """code -> 'decided by', parsed from the tables in docs/ACMG_CRITERIA.md."""
    rows = {}
    for m in re.finditer(r"^\|\s*`(\w+)`\s*\|[^|]*\|\s*([^|]+?)\s*\|", DOC.read_text(), re.M):
        rows[m.group(1)] = m.group(2)
    return rows


def test_doc_lists_every_criterion():
    assert set(_doc_rows()) == set(ACMG_28)


@pytest.mark.parametrize("cr", _results(), ids=lambda c: c.code)
def test_doc_agrees_with_the_code_on_how_each_is_decided(cr):
    assert _doc_rows()[cr.code] == _kind(cr), (
        f"docs/ACMG_CRITERIA.md says {cr.code} is '{_doc_rows()[cr.code]}', code says '{_kind(cr)}'"
    )


def test_doc_headline_counts_match_the_code():
    counts = {}
    for cr in _results():
        counts[_kind(cr)] = counts.get(_kind(cr), 0) + 1
    expected = f"**{counts['engine']} engine · {counts['model']} model · {counts['N/A']} N/A.**"
    assert expected in DOC.read_text(), f"doc headline should read: {expected}"


def test_every_registered_criterion_is_on_a_known_side():
    """`evaluate_criteria` auto-appends any newly registered code; `rules.BENIGN_CODES` is a
    hardcoded set. Nothing ties the two, and `side_of` treats an unlisted code as PATHOGENIC —
    so a future benign criterion would be scored as pathogenic evidence, silently. Two such
    criteria would combine to `Pathogenic [PATH-2 (>=2 Strong)]`. Pin the invariant."""
    from vcf2report.acmg.rules import BENIGN_CODES, side_of
    registered = set(all_criteria())
    for code in registered:
        expected = "benign" if code.startswith("B") else "pathogenic"
        assert side_of(code) == expected, (
            f"{code} is scored as {side_of(code)} evidence; add it to rules.BENIGN_CODES")
    # ...and nothing stale in the other direction (reserved simulation codes excepted).
    from vcf2report.acmg import rules
    stale = (BENIGN_CODES - registered) - {rules.HYPOTHETICAL_BENIGN}
    assert not stale, f"BENIGN_CODES names criteria that are not registered: {sorted(stale)}"


def test_reserved_simulation_codes_are_not_real_criteria():
    """explore.missing_evidence feeds these to rules.combine; if an evaluator ever claimed one
    of the names, the simulation would double-count a real criterion."""
    from vcf2report.acmg import rules
    for code in (rules.HYPOTHETICAL_PATHOGENIC, rules.HYPOTHETICAL_BENIGN):
        assert code not in all_criteria()
    assert rules.side_of(rules.HYPOTHETICAL_BENIGN) == "benign"
    assert rules.side_of(rules.HYPOTHETICAL_PATHOGENIC) == "pathogenic"
