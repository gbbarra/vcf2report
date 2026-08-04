"""The per-variant ACMG rationale is 97.8% of the report. What it repeats, and what it omits.

Measured on SYN-016 (1,178 classifications): the report was 6.76 MB, of which the rationale
section was 6.61 MB. Ten criteria produced byte-identical rows 1,178 times each (1.42 MB); seven
more fired for no variant at all but carried a per-variant evidence dict, so PM1 alone spent
356 KB printing `window=7, cutoff=3, enrichment_cutoff=2.0` for variants it never assessed
because the ClinVar residue index was absent.

The section is titled "auditable", so the tests here are mostly about what may NOT disappear.
"""

from __future__ import annotations

import pytest

from vcf2report.models import (
    Annotation,
    Classification,
    CriterionResult,
    QCSummary,
    Variant,
)
from vcf2report.report.assemble import ReportModel
from vcf2report.report.render import (
    HOIST_MIN_VARIANTS,
    _uninformative,
    render_markdown,
    run_constant_criteria,
)


N = HOIST_MIN_VARIANTS  # enough variants for the hoist to engage


def _cr(
    code, met=False, evidence=None, reasoning="because", citation=None, applies=True
):
    return CriterionResult(
        code,
        f"{code} name",
        "supporting",
        applies=applies,
        met=met,
        applied_strength="supporting" if met else None,
        evidence=evidence or {},
        citation=citation or [],
        reasoning=reasoning,
    )


def _cls(gene, criteria, tier="Uncertain Significance (VUS)"):
    return Classification(
        variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene=gene, hgvs_c="c.1A>T"),
        annotation=Annotation(),
        criteria=criteria,
        tier=tier,
        rule_path="x => y",
    )


def _report(classifications):
    return ReportModel(
        sample_id="S", hpo_terms=[], qc=QCSummary(), classifications=classifications
    )


# ------------------------------------------------------------------ what gets hoisted


def test_a_criterion_that_never_fires_is_stated_once():
    r = _report(
        [
            _cls(f"G{i}", [_cr("PS2", reasoning="Requires parental (trio) data")])
            for i in range(N)
        ]
    )
    assert run_constant_criteria(r.classifications) == {
        "PS2": "Requires parental (trio) data"
    }


def test_a_criterion_that_fires_even_once_is_never_hoisted():
    """The hoist must key on the RUN, not on the common case. One variant where PP3 fires
    makes PP3 a real per-variant line for all of them — otherwise the one row that carried
    the classification would vanish into a header note."""
    crits = [[_cr("PP3", met=(i == 3), evidence={"am": 0.9})] for i in range(N)]
    r = _report([_cls(f"G{i}", c) for i, c in enumerate(crits)])
    assert "PP3" not in run_constant_criteria(r.classifications)


def test_a_never_fired_criterion_with_evidence_is_still_hoisted():
    """This is the case the first implementation missed. PM1 records hotspot parameters for
    every variant, so it is not "uninformative" row-by-row — yet it fired for none of 1,178
    and the reason was the same absent residue index each time. 356 KB of parameters."""
    r = _report(
        [
            _cls(
                f"G{i}",
                [
                    _cr(
                        "PM1",
                        evidence={"window": 7, "cutoff": 3},
                        reasoning="ClinVar residue index unavailable",
                    )
                ],
            )
            for i in range(N)
        ]
    )
    assert run_constant_criteria(r.classifications) == {
        "PM1": "ClinVar residue index unavailable"
    }


def test_mixed_reasons_report_the_dominant_one_and_admit_the_rest():
    """PM1 had 12 distinct reasons on SYN-016 (mostly "not a missense variant"). Printing only
    the commonest would state something false about the other 131 variants."""
    crits = [
        [
            _cr(
                "PM1",
                reasoning="residue index unavailable" if i < N - 1 else "not missense",
            )
        ]
        for i in range(N)
    ]
    got = run_constant_criteria([_cls(f"G{i}", c) for i, c in enumerate(crits)])["PM1"]
    assert "residue index unavailable" in got
    assert f"for {N - 1} of {N} variants" in got and "1 other reason" in got


def test_no_classifications_hoists_nothing():
    assert run_constant_criteria([]) == {}


# --------------------------------------------------------------- what may not disappear


def test_a_row_that_did_not_fire_but_recorded_evidence_is_kept():
    """The line proving the engine LOOKED. PM2 with gnomad_af=0.0714 did not fire, and that
    is precisely the audit trail for why a variant is not rare."""
    assert not _uninformative(_cr("PM2", evidence={"gnomad_af": 0.0714}))
    assert _uninformative(_cr("PS2", evidence={}, citation=[]))
    assert not _uninformative(_cr("PP5", met=True))
    assert not _uninformative(_cr("BP6", evidence={}, citation=["VCV000123"]))


def test_the_rendered_report_declares_what_it_omitted():
    """Auditability means the omission is stated. A shorter report that silently drops rows
    is worse than a long one."""
    md = render_markdown(
        _report(
            [
                _cls(
                    f"G{i}",
                    [
                        _cr("PS2", reasoning="Requires parental (trio) data"),
                        _cr("PP3", met=True, evidence={"am": 0.9}),
                    ],
                )
                for i in range(N)
            ]
        )
    )
    assert "fired for **no** variant" in md
    assert "PS2" in md and "Requires parental (trio) data" in md
    assert "_results.json" in md, "must point at where the full 28-criterion set lives"
    assert "| **PP3** |" in md, "a criterion that fired was dropped from the table"
    assert md.count("| **PS2** |") == 0, (
        "the hoisted criterion is still printed per variant"
    )


def test_per_variant_omissions_are_counted_not_silent():
    """A criterion uninformative for THIS variant but not for the run gets no header note,
    so the block itself must say how many rows it left out."""
    md = render_markdown(
        _report(
            [
                _cls(
                    "A",
                    [
                        _cr("BP7", evidence={}, reasoning="not synonymous"),
                        _cr("PP3", met=True, evidence={"am": 0.9}),
                    ],
                )
            ]
            + [
                _cls(
                    f"G{i}",
                    [
                        _cr("BP7", met=True, evidence={"x": 1}),
                        _cr("PP3", met=True, evidence={"am": 0.8}),
                    ],
                )
                for i in range(N)
            ]
        )
    )
    assert "further criterion" in md, "the block dropped a row without saying so"


# --------------------------------------------------------- the two renderers must agree


def test_both_renderers_apply_the_same_trim():
    """Found the hard way: the trim was written into the builtin renderer while the report is
    actually produced by the Jinja template, so the first measurement showed a 2-byte saving
    on a 6.9 MB file. Whether a laudo is 4 MB or 7 MB must not depend on whether jinja2
    happens to be importable."""
    from vcf2report.report import render as R

    r = _report(
        [
            _cls(
                f"G{i}",
                [
                    _cr("PS2", reasoning="Requires parental (trio) data"),
                    _cr(
                        "PM1",
                        evidence={"window": 7},
                        reasoning="residue index unavailable",
                    ),
                    _cr("PP3", met=True, evidence={"am": 0.9}),
                ],
            )
            for i in range(N)
        ]
    )
    jinja_md = render_markdown(r)
    builtin_md = R._render_markdown_builtin(r)
    for md, who in ((jinja_md, "jinja"), (builtin_md, "builtin")):
        assert md.count("| **PS2** |") == 0, f"{who} still prints the hoisted PS2"
        assert md.count("| **PM1** |") == 0, f"{who} still prints the hoisted PM1"
        assert md.count("| **PP3** |") == N, f"{who} lost a criterion that fired"
        assert "fired for **no** variant" in md, f"{who} does not declare the omission"


def test_a_small_report_is_never_trimmed():
    """The hoist removes REPETITION. On a one-variant report "fired for no variant" is
    trivially true of nearly every criterion, and applying it would empty the table and call
    that a summary — which is how the first version broke test_no_rendered_report_prints_None.
    Below the threshold the section is tens of KB anyway; print it whole.
    """
    small = [
        _cls(
            "A",
            [
                _cr("PS2", reasoning="Requires parental (trio) data"),
                _cr("PM2", evidence={"gnomad_af": None}),
            ],
        )
    ]
    assert run_constant_criteria(small) == {}
    md = render_markdown(_report(small))
    assert "| **PS2** |" in md and "| **PM2** |" in md
    assert "fired for **no** variant" not in md
