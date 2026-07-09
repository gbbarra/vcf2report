"""End-to-end regression guard over the committed example exomes.

Runs the *actual* VCFs shipped in ``data/example/`` (real 1000G DRAGEN
backgrounds with real ClinVar pathogenics spiked in) through the full pipeline
and asserts, for every case, that:

* both spiked variants are surfaced with their expected ACMG tier, and
* the phenotype-matched DEE gene routes to the **primary** findings while the
  unrelated ACMG-SF gene routes to **secondary** (the split must not leak).

This ties the repo's demo inputs to their expected outputs, so a regression in
parsing, annotation, classification, or the primary/secondary gating fails a
test instead of silently changing a report. It exercises the same files a user
gets from ``git pull`` — no network, no external annotation.
"""
from pathlib import Path

import pytest

from vcf2report.cli import read_hpo_file
from vcf2report.pipeline import run_pipeline
from vcf2report.report.assemble import split_findings

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "data" / "example"
HPO = ROOT / "data" / "synthetic"

# case -> {gene: (expected_tier, expected_bucket)}
CASES = {
    "SYN-001": {"SCN1A": ("Pathogenic", "primary"), "RB1": ("Pathogenic", "secondary")},
    "SYN-002": {"KCNQ2": ("Pathogenic", "primary"), "APC": ("Likely Pathogenic", "secondary")},
    "SYN-003": {"SCN2A": ("Pathogenic", "primary"), "STK11": ("Pathogenic", "secondary")},
    "SYN-004": {"STXBP1": ("Pathogenic", "primary"), "WT1": ("Pathogenic", "secondary")},
    "SYN-005": {"SLC2A1": ("Pathogenic", "primary"), "FBN1": ("Pathogenic", "secondary")},
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_example_case_tiers_and_split(case):
    vcf = EXAMPLE / f"{case}.synthetic.vcf.gz"
    hpo = HPO / f"{case}.hpo.txt"
    assert vcf.exists(), f"missing committed example VCF {vcf}"

    report = run_pipeline(vcf, hpo_terms=read_hpo_file(hpo), sample_id=case)
    by_gene = {c.variant.gene: c for c in report.classifications}

    primary, secondary, _other = split_findings(report.classifications)
    bucket = {}
    for c in primary:
        bucket[c.variant.gene] = "primary"
    for c in secondary:
        bucket[c.variant.gene] = "secondary"

    for gene, (tier, want_bucket) in CASES[case].items():
        assert gene in by_gene, f"{case}: spiked {gene} not classified"
        assert by_gene[gene].tier == tier, (
            f"{case}: {gene} tier {by_gene[gene].tier!r} != expected {tier!r}"
        )
        assert bucket.get(gene) == want_bucket, (
            f"{case}: {gene} routed to {bucket.get(gene)!r}, expected {want_bucket!r}"
        )
