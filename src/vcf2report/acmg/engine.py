"""ACMG classification engine.

Evaluates every registered criterion against a variant + its annotation, then
applies the combining rules to reach a 5-tier call. The output is the full
criterion trail (not just the tier) so the report can show *why*.
"""
from __future__ import annotations

from ..models import Annotation, Classification, CriterionResult, Variant
from . import criteria as _criteria
from . import rules as _rules


def evaluate_criteria(variant: Variant, annotation: Annotation) -> list[CriterionResult]:
    """Run all criterion evaluators; return results in a stable, readable order."""
    order = [
        "PVS1", "PS1", "PS2", "PS3", "PS4",
        "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
        "PP2", "PP3", "PP4", "PP5",
        "BA1", "BS1", "BS2", "BP4", "BP6", "BP7",
    ]
    fns = _criteria.all_criteria()
    results: list[CriterionResult] = []
    for code in order:
        fn = fns.get(code)
        if fn:
            results.append(fn(variant, annotation))
    # append any criteria not in the explicit order (future-proofing)
    for code, fn in fns.items():
        if code not in order:
            results.append(fn(variant, annotation))
    return results


def classify(variant: Variant, annotation: Annotation) -> Classification:
    """Classify a single variant into an auditable :class:`Classification`."""
    results = evaluate_criteria(variant, annotation)
    tier, rule_path = _rules.combine(results)
    return Classification(
        variant=variant,
        annotation=annotation,
        criteria=results,
        tier=tier,
        rule_path=rule_path,
    )
