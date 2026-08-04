"""Inputs at the edges: nothing, one variant, the wrong build, a phenotype that matches nothing.

The happy path has been exercised to exhaustion. These had never been run — and two of the five
produced a conclusion that a clinician would read as a negative RESULT when no result existed.

The conclusion is the bottom-line-up-front line. Everything here is about what it may and may
not claim.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from vcf2report import config

REPO = config.REPO_ROOT
HDR = "##fileformat=VCFv4.2\n"
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
# BBS2 c.472-2A>G — the planted diagnosis in SYN-073, with a real SnpEff ANN.
BBS2_ROW = (
    "chr16\t56510923\t.\tT\tC\t500\tPASS\tANN=C|splice_acceptor_variant|HIGH|BBS2|BBS2|"
    "transcript|NM_031885.5|protein_coding|3/16|c.472-2A>G||||||\tGT:DP:GQ\t1/1:40:60\n"
)


def _run(vcf: Path, out: Path, hpo: Path | None = None, sample: str = "T"):
    argv = [
        sys.executable,
        "-m",
        "vcf2report.cli",
        str(vcf),
        "--sample-id",
        sample,
        "--out",
        str(out),
    ]
    if hpo:
        argv += ["--hpo", str(hpo)]
    p = subprocess.run(argv, cwd=REPO, capture_output=True, text=True, timeout=900)
    assert p.returncode == 0, p.stdout + p.stderr
    return (out / f"{sample}_report.md").read_text()


def _conclusion(md: str) -> str:
    return md.split("## Conclusion", 1)[1].split("## Quality", 1)[0]


# ------------------------------------------------------------------ nothing to analyse


def test_an_empty_vcf_does_not_report_a_negative_finding(tmp_path):
    """ "No Pathogenic / Likely Pathogenic finding" is a RESULT.

    With zero variants there is no result — only an absence of input, which is the one thing
    this report may never dress up as evidence. The QC section did carry a warning, but the
    conclusion, read first, announced a negative.
    """
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(HDR + "##reference=GRCh38\n" + COLS)
    concl = _conclusion(_run(vcf, tmp_path / "out"))

    assert "No analysis was performed" in concl
    assert "NOT a negative result" in concl
    assert "No Pathogenic / Likely Pathogenic finding" not in concl


def test_an_empty_vcf_still_warns_in_the_qc_section(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(HDR + "##reference=GRCh38\n" + COLS)
    md = _run(vcf, tmp_path / "out")
    assert "contained no variant records" in md


def test_a_single_variant_vcf_classifies_it(tmp_path):
    """One record must not trip the funnel, the QC percentages, or the quality panel."""
    vcf = tmp_path / "one.vcf"
    vcf.write_text(HDR + "##reference=GRCh38\n##contig=<ID=chr16>\n" + COLS + BBS2_ROW)
    md = _run(vcf, tmp_path / "out", hpo=config.DATA_DIR / "example/SYN-073.hpo.txt")
    assert "BBS2" in md
    # Ti/Tv is undefined on one SNV; the panel must omit it, not print a fabricated ratio.
    assert "Ti/Tv" not in md


# ------------------------------------------------------------------------- wrong build


@pytest.mark.parametrize(
    "ref", ["hs37d5.fa", "Homo_sapiens_assembly19.fasta", "hg19.fa"]
)
def test_a_grch37_vcf_is_detected_and_coordinate_lookups_are_skipped(tmp_path, ref):
    """The dangerous failure is silent: GRCh37 coordinates looked up against GRCh38 stores
    match the wrong loci and every frequency is wrong without anything looking wrong.

    hs37d5 and assembly19 are what legacy exomes actually declare, and both were missing from
    the detector until #24 — this exercises that fix end to end.
    """
    vcf = tmp_path / "b37.vcf"
    vcf.write_text(
        HDR
        + f"##reference=file:///refs/{ref}\n"
        + COLS
        + "16\t56544835\t.\tT\tC\t500\tPASS\t.\tGT:DP:GQ\t1/1:40:60\n"
    )
    md = _run(vcf, tmp_path / "out")

    assert "**Genome build:** GRCh37" in md, (
        "the report claimed the assumed build, not the real one"
    )
    assert "GRCh37, expected GRCh38" in md
    assert "SKIPPED" in md, (
        "coordinate annotation must be skipped, not run against the wrong build"
    )


# ---------------------------------------------------------------------------- phenotype


def test_a_phenotype_that_matches_nothing_does_not_hide_a_pathogenic_call(tmp_path):
    """Measured on SYN-073 with two unrelated HPO terms: the conclusion opened with "No
    Pathogenic / Likely Pathogenic finding" while the report below listed a Pathogenic BBS2.

    The qualifier "in phenotype-matched genes" is accurate and easy to miss — and an incomplete
    referral phenotype is the norm, not the exception, so this is the common way a real finding
    disappears from the line a clinician reads first.
    """
    hpo = tmp_path / "nomatch.hpo.txt"
    hpo.write_text("HP:0031936\nHP:0032101\n")  # deliberately unrelated to BBS2
    md = _run(
        config.DATA_DIR / "example/SYN-073.BBS2.annotated.vcf.gz",
        tmp_path / "out",
        hpo=hpo,
        sample="NOMATCH",
    )
    concl = _conclusion(md)

    assert "OUTSIDE the phenotype match" in concl, (
        "a Pathogenic call outside the phenotype match is absent from the conclusion"
    )
    assert "BBS2" in concl
    assert "referral phenotype is incomplete" in concl


def test_a_recessive_carrier_is_not_listed_as_a_hidden_pathogenic_finding(tmp_path):
    """SYN-073 also carries a lone het P/LP in SKIC2, an AR-only gene. That is carrier status
    and already has its own bucket; surfacing it as a missed diagnosis would undo the triage."""
    hpo = tmp_path / "nomatch.hpo.txt"
    hpo.write_text("HP:0031936\nHP:0032101\n")
    concl = _conclusion(
        _run(
            config.DATA_DIR / "example/SYN-073.BBS2.annotated.vcf.gz",
            tmp_path / "out",
            hpo=hpo,
            sample="NOMATCH",
        )
    )
    outside = [ln for ln in concl.splitlines() if "OUTSIDE the phenotype match" in ln]
    assert outside and "SKIC2" not in outside[0], (
        "the carrier leaked into the P/LP headline"
    )


def test_no_phenotype_at_all_says_so_and_still_surfaces_the_finding(tmp_path):
    """The no-HPO path was already right; pin it, because the fix above touches the same
    branch and this is the behaviour it must not regress."""
    concl = _conclusion(
        _run(
            config.DATA_DIR / "example/SYN-073.BBS2.annotated.vcf.gz",
            tmp_path / "out",
            sample="NOHPO",
        )
    )
    assert "No phenotype was available" in concl
    assert "BBS2" in concl and "phenotype overlap NOT assessed" in concl
