"""Integrity tests: exhaustive ACMG combining matrix, determinism, robustness.

These lock the classification core (Richards 2015 Table 5) against regressions
and prove the parser degrades gracefully on malformed/edge input.
"""
import gzip

from vcf2report import config
from vcf2report.acmg import rules
from vcf2report.acmg.rules import (BENIGN, LIKELY_BENIGN, LIKELY_PATHOGENIC,
                                   PATHOGENIC, VUS)
from vcf2report.models import CriterionResult
from vcf2report.pipeline import run_pipeline
from vcf2report.vcf.parse import parse_vcf


def _mk(code: str, strength: str) -> CriterionResult:
    return CriterionResult(code, code, strength, applies=True, met=True,
                           applied_strength=strength)


# (list of (code, strength), expected tier) covering every Richards Table 5 rule.
_MATRIX = [
    # Pathogenic
    ([("PVS1", "very_strong"), ("PS1", "strong")], PATHOGENIC),                 # PATH-1
    ([("PVS1", "very_strong"), ("PM2", "moderate"), ("PM4", "moderate")], PATHOGENIC),
    ([("PVS1", "very_strong"), ("PM2", "moderate"), ("PP3", "supporting")], PATHOGENIC),
    ([("PVS1", "very_strong"), ("PP3", "supporting"), ("PP4", "supporting")], PATHOGENIC),
    ([("PS1", "strong"), ("PS3", "strong")], PATHOGENIC),                       # PATH-2
    ([("PS1", "strong"), ("PM1", "moderate"), ("PM2", "moderate"), ("PM4", "moderate")], PATHOGENIC),  # PATH-3
    ([("PS1", "strong"), ("PM2", "moderate"), ("PP3", "supporting"),
      ("PP4", "supporting"), ("PM4", "moderate")], PATHOGENIC),
    # Likely Pathogenic
    ([("PVS1", "very_strong"), ("PM2", "moderate")], LIKELY_PATHOGENIC),        # LP-1
    ([("PS1", "strong"), ("PM2", "moderate")], LIKELY_PATHOGENIC),              # LP-2
    ([("PS1", "strong"), ("PP3", "supporting"), ("PP4", "supporting")], LIKELY_PATHOGENIC),  # LP-3
    ([("PM1", "moderate"), ("PM2", "moderate"), ("PM4", "moderate")], LIKELY_PATHOGENIC),    # LP-4
    ([("PM2", "moderate"), ("PM4", "moderate"), ("PP3", "supporting"), ("PP4", "supporting")], LIKELY_PATHOGENIC),  # LP-5
    ([("PM2", "moderate"), ("PP2", "supporting"), ("PP3", "supporting"),
      ("PP4", "supporting"), ("PP5", "supporting")], LIKELY_PATHOGENIC),        # LP-6
    # Benign
    ([("BA1", "stand_alone")], BENIGN),                                         # BEN-1
    ([("BS1", "strong"), ("BS2", "strong")], BENIGN),                           # BEN-2
    ([("BS1", "strong"), ("BP4", "supporting")], LIKELY_BENIGN),               # LB-1
    ([("BP4", "supporting"), ("BP7", "supporting")], LIKELY_BENIGN),           # LB-2
    # VUS
    ([("PM2", "moderate")], VUS),                                               # insufficient
    # conflicting: PVS1+PM2 fires LP-1 AND BA1 fires BEN-1 -> VUS
    ([("PVS1", "very_strong"), ("PM2", "moderate"), ("BA1", "stand_alone")], VUS),
]


def test_acmg_combining_matrix():
    for codes, expected in _MATRIX:
        crits = [_mk(c, s) for c, s in codes]
        tier, path = rules.combine(crits)
        assert tier == expected, f"{[c for c, _ in codes]} -> {tier} (want {expected}); {path}"


def test_pipeline_is_deterministic():
    hpo = ["HP:0001250", "HP:0001263", "HP:0002133"]
    r1 = run_pipeline(config.SAMPLE_VCF, hpo_terms=hpo)
    r2 = run_pipeline(config.SAMPLE_VCF, hpo_terms=hpo)
    t1 = [(c.variant.gene, c.tier, tuple(c.met_codes)) for c in r1.classifications]
    t2 = [(c.variant.gene, c.tier, tuple(c.met_codes)) for c in r2.classifications]
    assert t1 == t2


MALFORMED = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
short\tline
2\t.\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
2\t100\t.\tA\t\t.\tPASS\t.\tGT\t0/1
2\t200\t.\tA\tG\t.\tPASS\tGENE=X;CSQ=missense_variant\tGT:DP:GQ:AD\t0/1:30:99:15,15
2\t300\t.\tC\tT\t.\tPASS\tGENE=Y;CSQ=missense_variant
"""


def test_parser_survives_malformed_input(tmp_path):
    p = tmp_path / "bad.vcf"
    p.write_text(MALFORMED)
    variants, build, _ = parse_vcf(p)
    keys = {v.key for v in variants}
    assert "2-200-A-G" in keys       # the one well-formed genotyped record
    assert "2-300-C-T" in keys        # valid record without a sample column
    assert build == "GRCh38"
    # malformed rows (short line, POS='.', empty ALT) are skipped, not crashing.
    assert all(v.pos > 0 for v in variants)


def test_empty_vcf(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    variants, _, _ = parse_vcf(p)
    assert variants == []


def test_gzipped_vcf_roundtrips(tmp_path):
    raw = (config.SAMPLE_VCF).read_text()
    gz = tmp_path / "s.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(raw)
    plain, _, _ = parse_vcf(config.SAMPLE_VCF)
    gzv, _, _ = parse_vcf(gz)
    assert [v.key for v in gzv] == [v.key for v in plain]
