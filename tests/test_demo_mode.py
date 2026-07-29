"""Demo mode must relax the store gate for the repository's own fixtures and for NOTHING else,
and it must never produce an unstamped laudo.

The two properties under test are asymmetric on purpose. Letting a fixture through the gate is
a convenience; letting a *real* exome through, or letting a fixture laudo pass for a patient
result, is the failure the gate exists to prevent. So every test here is written from the
dangerous direction: what would have to break for a demo laudo to look real.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from vcf2report import config, demo, stores
from vcf2report.models import QCSummary
from vcf2report.report.assemble import build_report, summarize
from vcf2report.report.render import _render_markdown_builtin, render_markdown

REPO = config.REPO_ROOT
EXAMPLE = config.DATA_DIR / "example"
DEMO_VCF = EXAMPLE / "SYN-073.BBS2.annotated.vcf.gz"


# --------------------------------------------------------------------------- recognition

def test_committed_example_is_recognised_as_a_fixture():
    assert DEMO_VCF.exists(), "the committed demo fixture moved — update this test and the skill"
    assert demo.is_demo_vcf(DEMO_VCF)


def test_relative_path_to_the_same_file_is_recognised(monkeypatch):
    """The skill and the docs use repo-relative paths; resolution must not depend on cwd."""
    monkeypatch.chdir(REPO)
    assert demo.is_demo_vcf("data/example/SYN-073.BBS2.annotated.vcf.gz")


def test_a_path_outside_data_example_is_not_a_fixture(tmp_path):
    other = tmp_path / "patient.vcf.gz"
    other.write_bytes(b"")
    assert not demo.is_demo_vcf(other)


def test_traversal_out_of_data_example_is_not_a_fixture(tmp_path):
    """`data/example/../../patient.vcf` must not read as a fixture — hence the resolve()."""
    sneaky = EXAMPLE / ".." / ".." / "patient.vcf"
    assert not demo.is_demo_vcf(sneaky)


def test_a_non_vcf_inside_data_example_is_not_a_fixture():
    """The HPO term files live in the same directory and are not VCFs."""
    hpo = EXAMPLE / "SYN-073.hpo.txt"
    assert hpo.exists()
    assert not demo.is_demo_vcf(hpo)


def test_missing_and_none_paths_are_safe():
    assert not demo.is_demo_vcf(None)
    assert not demo.is_demo_vcf("")
    assert not demo.is_demo_vcf("/nonexistent/nowhere.vcf.gz")


# ------------------------------------------------------------------- the gate exemption

def test_demo_not_requested_leaves_the_gate_untouched():
    d = demo.decide(DEMO_VCF, demo_requested=False)
    assert not d.active and not d.refused


def test_demo_requested_on_a_fixture_is_active():
    d = demo.decide(DEMO_VCF, demo_requested=True)
    assert d.active and not d.refused
    assert "data/example" in d.vcf


def test_demo_requested_on_a_real_vcf_is_REFUSED_not_honoured(tmp_path):
    """The core safety property: --demo can never be the reason a real exome ran blind."""
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    d = demo.decide(patient, demo_requested=True)
    assert d.refused and not d.active
    assert "not a store override" in d.reason


def test_gate_refusal_reports_not_ready(tmp_path):
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    g = demo.gate(patient, demo_requested=True, measure=False)
    assert g["ready"] is False
    assert g["mode"] == "refused"
    assert g["exempted"] == []


def test_gate_exempts_blocking_stores_for_a_fixture(monkeypatch):
    monkeypatch.setattr(stores, "gate", lambda measure=True: {
        "ready": False, "blocking": ["gnomad", "clinvar"], "stale": [], "health": {}})
    g = demo.gate(DEMO_VCF, demo_requested=True, measure=False)
    assert g["ready"] is True
    assert g["mode"] == "demo"
    assert g["exempted"] == ["gnomad", "clinvar"]


def test_gate_still_blocks_a_real_vcf_when_stores_block(monkeypatch, tmp_path):
    monkeypatch.setattr(stores, "gate", lambda measure=True: {
        "ready": False, "blocking": ["gnomad"], "stale": [], "health": {}})
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    g = demo.gate(patient, demo_requested=False, measure=False)
    assert g["ready"] is False
    assert g["mode"] == "full"


def test_a_fixture_on_a_healthy_machine_is_still_flagged_demo(monkeypatch):
    """Nothing was exempted, but the input is still a fixture — the laudo must say so."""
    monkeypatch.setattr(stores, "gate", lambda measure=True: {
        "ready": True, "blocking": [], "stale": [], "health": {}})
    g = demo.gate(DEMO_VCF, demo_requested=True, measure=False)
    assert g["ready"] is True
    assert g["mode"] == "demo"
    assert g["exempted"] == []


# ------------------------------------------------------------------------- the stamp

def test_provenance_is_empty_for_a_real_vcf(tmp_path):
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    assert demo.provenance(patient) == {}


def test_provenance_stamps_a_fixture_with_NO_flag_set(monkeypatch):
    """The stamp must not depend on --demo or VCF2REPORT_DEMO.

    A forgotten flag has to fail toward "stamped", never toward "looks like a patient result".
    """
    monkeypatch.delenv("VCF2REPORT_DEMO", raising=False)
    prov = demo.provenance(DEMO_VCF)
    assert prov["mode"] == "demo"
    assert prov["vcf"].startswith("data/example/")


def test_stamp_names_the_criteria_each_absent_store_costs():
    prov = {"mode": "demo", "vcf": "data/example/x.vcf.gz",
            "stores_absent": ["gnomad", "clinvar"],
            "criteria_degraded": demo.degraded_criteria(["gnomad", "clinvar"])}
    line = demo.stamp_line(prov)
    # A vague "demo data" note would leave the reader guessing which rows to distrust.
    for code in ("PM2", "BA1", "BS1", "PS1", "PM5", "PP5"):
        assert code in line
    # AlphaMissense was present, so its criteria must NOT be listed as degraded.
    assert "BP4" not in line


def test_stamp_line_is_empty_without_provenance():
    assert demo.stamp_line({}) == ""


def test_degraded_criteria_dedupes_and_is_stable():
    assert demo.degraded_criteria([]) == []
    assert demo.degraded_criteria(["gnomad", "gnomad"]) == ["PM2", "BA1", "BS1", "BS2"]


# --------------------------------------------------------- the stamp reaches the reader

def _demo_report():
    qc = QCSummary(total_variants=1, build="GRCh38")
    prov = {"mode": "demo", "vcf": "data/example/SYN-073.BBS2.annotated.vcf.gz",
            "reason": "fixture", "stores_absent": ["gnomad"],
            "criteria_degraded": ["PM2", "BA1", "BS1", "BS2"]}
    return build_report("SYN-073", ["HP:0000041"], qc, [], provenance=prov)


def test_conclusion_leads_with_the_demo_stamp():
    """A reader who reads one bullet must already know this is not a patient."""
    lines = summarize(_demo_report())
    assert lines, "conclusion must not be empty"
    assert "DEMONSTRATION RUN" in lines[0]


def test_methods_records_the_data_mode():
    r = _demo_report()
    assert "DEMONSTRATION" in r.methods["data_mode"]
    assert r.methods["criteria_not_backed_by_full_stores"] == ["PM2", "BA1", "BS1", "BS2"]


def test_methods_is_untouched_for_a_normal_run():
    qc = QCSummary(total_variants=1, build="GRCh38")
    r = build_report("REAL-001", [], qc, [])
    assert "data_mode" not in r.methods
    assert r.provenance == {}


def test_provenance_survives_into_the_json():
    assert _demo_report().to_dict()["provenance"]["mode"] == "demo"


@pytest.mark.parametrize("render", [render_markdown, _render_markdown_builtin],
                         ids=["jinja-template", "builtin"])
def test_both_renderers_carry_the_demo_banner(render):
    """The repo ships templates/report.md.j2 AND a dependency-free fallback; a laudo rendered
    by either must be unmistakable. A banner in only one is how a fixture escapes."""
    md = render(_demo_report())
    assert "DEMONSTRATION RUN" in md
    assert "SYN-073.BBS2.annotated.vcf.gz" in md
    # It must appear ABOVE the findings, not buried in Methods at the bottom.
    assert md.index("DEMONSTRATION RUN") < md.index("Conclusion")


@pytest.mark.parametrize("render", [render_markdown, _render_markdown_builtin],
                         ids=["jinja-template", "builtin"])
def test_neither_renderer_stamps_a_normal_run(render):
    qc = QCSummary(total_variants=1, build="GRCh38")
    md = render(build_report("REAL-001", [], qc, []))
    assert "DEMONSTRATION" not in md


# ------------------------------------------------------------------- the CLI / gate script

def _run(*argv):
    return subprocess.run([sys.executable, *argv], cwd=REPO,
                          capture_output=True, text=True, timeout=600)


def test_gate_script_exits_zero_for_a_fixture_and_says_DEMO():
    p = _run("scripts/check_stores.py", "--gate", "--demo",
             "data/example/SYN-073.BBS2.annotated.vcf.gz")
    assert p.returncode == 0, p.stdout + p.stderr
    assert "DEMO MODE" in p.stdout
    # The workflow keys its narration off this phrase; keep them in sync.
    assert "stamped" in p.stdout.lower()


def test_gate_script_refuses_a_non_fixture_and_exits_nonzero(tmp_path):
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    p = _run("scripts/check_stores.py", "--gate", "--demo", str(patient))
    assert p.returncode != 0, p.stdout
    assert "REFUSED" in p.stdout


def test_cli_refuses_demo_on_a_non_fixture(tmp_path):
    patient = tmp_path / "patient.vcf.gz"
    patient.write_bytes(b"")
    p = _run("-m", "vcf2report.cli", str(patient), "--demo")
    assert p.returncode == 2, p.stdout + p.stderr
    assert "REFUSED" in p.stderr


def test_the_skill_documents_demo_mode():
    """The guided flow is what a user actually drives; if SKILL.md does not describe the
    exemption, the flow will keep dead-ending at the gate on a store-less machine."""
    skill = (REPO / ".claude/skills/vcf2report/SKILL.md").read_text()
    assert "demo mode" in skill.lower()
    assert "data/example/" in skill
    # It must document the refusal, not just the convenience.
    assert "refuse" in skill.lower()


def test_the_artifact_template_has_a_demo_banner_placeholder():
    html = (REPO / ".claude/skills/vcf2report/references/report_template.html").read_text()
    assert "{{DEMO_BANNER" in html
    assert ".demo{" in html, "the banner needs a style, or it renders as unstyled text"
