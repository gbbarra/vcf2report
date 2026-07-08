"""Property-based / fuzz tests (hypothesis) for the parser and the ACMG engine.

These assert invariants over *generated* inputs: the parser never crashes on
arbitrary bytes, and the combining rules always return a valid tier and behave
monotonically (more pathogenic evidence never makes a call more benign, and
vice-versa).
"""
import os
import tempfile

from hypothesis import given, settings
from hypothesis import strategies as st

from vcf2report.acmg import rules
from vcf2report.acmg.rules import (BENIGN, LIKELY_BENIGN, LIKELY_PATHOGENIC,
                                   PATHOGENIC, VUS)
from vcf2report.models import CriterionResult
from vcf2report.vcf.parse import parse_vcf

_HEADER = "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
_RANK = {PATHOGENIC: 0, LIKELY_PATHOGENIC: 1, VUS: 2, LIKELY_BENIGN: 3, BENIGN: 4}

_PATHO = [("PVS1", "very_strong"), ("PS1", "strong"), ("PS3", "strong"),
          ("PM1", "moderate"), ("PM2", "moderate"), ("PM4", "moderate"),
          ("PP2", "supporting"), ("PP3", "supporting"), ("PP4", "supporting"), ("PP5", "supporting")]
_BENIGN = [("BA1", "stand_alone"), ("BS1", "strong"), ("BS2", "strong"),
           ("BP4", "supporting"), ("BP7", "supporting")]
_ALL = _PATHO + _BENIGN


def _mk(code, strength):
    return CriterionResult(code, code, strength, applies=True, met=True, applied_strength=strength)


def _write(text: str) -> str:
    fd, path = tempfile.mkstemp(suffix=".vcf")
    os.write(fd, text.encode("utf-8", "ignore"))
    os.close(fd)
    return path


# --- parser fuzz -------------------------------------------------------------
@settings(max_examples=250, deadline=None)
@given(st.text(max_size=400))
def test_parser_never_crashes_on_arbitrary_text(text):
    path = _write(text)
    try:
        variants, build, header = parse_vcf(path)
        assert isinstance(variants, list)
        assert build in (None, "GRCh38", "GRCh37")
        for v in variants:                       # any emitted variant is well-formed
            assert v.pos > 0 and v.ref and v.alt
    finally:
        os.unlink(path)


@settings(max_examples=150, deadline=None)
@given(
    chrom=st.sampled_from(["1", "2", "7", "X", "chr17"]),
    pos=st.integers(min_value=1, max_value=10 ** 8),
    ref=st.sampled_from(["A", "C", "G", "T"]),
    alt=st.sampled_from(["A", "C", "G", "T"]),
)
def test_wellformed_record_always_parses(chrom, pos, ref, alt):
    line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1"
    path = _write(_HEADER + line + "\n")
    try:
        variants, _, _ = parse_vcf(path)
        assert len(variants) == 1
        c = chrom[3:] if chrom.lower().startswith("chr") else chrom
        assert variants[0].key == f"{c}-{pos}-{ref}-{alt}"
    finally:
        os.unlink(path)


# --- ACMG combining invariants ----------------------------------------------
@settings(max_examples=400, deadline=None)
@given(st.lists(st.sampled_from(_ALL), unique_by=lambda x: x[0], max_size=12))
def test_combine_always_returns_valid_tier(codes):
    tier, path = rules.combine([_mk(c, s) for c, s in codes])
    assert tier in _RANK
    assert isinstance(path, str) and path


@settings(max_examples=400, deadline=None)
@given(st.lists(st.sampled_from(_PATHO), unique_by=lambda x: x[0], max_size=8),
       st.sampled_from(_PATHO))
def test_more_pathogenic_evidence_never_more_benign(base, extra):
    codes = {c for c, _ in base}
    r0 = _RANK[rules.combine([_mk(c, s) for c, s in base])[0]]
    if extra[0] not in codes:
        r1 = _RANK[rules.combine([_mk(c, s) for c, s in base] + [_mk(*extra)])[0]]
        assert r1 <= r0        # rank 0 = Pathogenic; adding pathogenic can't raise it


@settings(max_examples=400, deadline=None)
@given(st.lists(st.sampled_from(_BENIGN), unique_by=lambda x: x[0], max_size=4),
       st.sampled_from(_BENIGN))
def test_more_benign_evidence_never_more_pathogenic(base, extra):
    codes = {c for c, _ in base}
    r0 = _RANK[rules.combine([_mk(c, s) for c, s in base])[0]]
    if extra[0] not in codes:
        r1 = _RANK[rules.combine([_mk(c, s) for c, s in base] + [_mk(*extra)])[0]]
        assert r1 >= r0        # rank 4 = Benign; adding benign can't lower it
