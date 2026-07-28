"""The MCP surface must be able to explore a finished run, not just produce one.

Claude Desktop drives the engine exclusively through these tools. Until now `run_report`
wrote only the Markdown laudo — `write_explore` was called from the CLI alone — and there
were no explore tools, so the conversational follow-up the persisted run exists for was
simply unavailable on that surface: a Desktop user could generate a report and then had no
way to ask why a criterion fired.

These tests are skipped when the MCP SDK is absent (it is an optional extra, and CI installs
the light dependency set), so they never make the suite depend on it.
"""
import pytest

pytest.importorskip("mcp", reason="MCP SDK is an optional extra")

from vcf2report import mcp_server as M  # noqa: E402


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    out = tmp_path_factory.mktemp("mcp")
    return M.run_report("data/example/SYN-051.PIGA.annotated.vcf.gz",
                        hpo_terms=["HP:0001250"], sample_id="T-PIGA", out_dir=str(out))


def test_run_report_writes_and_returns_the_queryable_json(run):
    # The Desktop path produced a laudo with no companion JSON, so nothing downstream could
    # read the run back. The path must also be RETURNED, or the client cannot chain to it.
    from pathlib import Path
    assert "results_json" in run, "run_report must tell the client where the data is"
    assert Path(run["results_json"]).is_file()


def test_explore_case_opens_the_run(run):
    o = M.explore_case(run["results_json"])
    assert o["sample_id"] == "T-PIGA"
    assert o["bucket_counts"]["primary"] >= 1
    assert isinstance(o["conclusion"], list) and o["conclusion"]


def test_explore_gene_gives_the_digest(run):
    g = M.explore_gene(run["results_json"], "PIGA")
    assert g["gene"] == "PIGA"
    assert g["variants"] and g["variants"][0]["hgvs_p"] == "p.Ser330Asn"


def test_explore_gene_with_a_criterion_gives_the_audit_trail(run):
    r = M.explore_gene(run["results_json"], "PIGA", "pm2")   # case-insensitive
    assert r["criterion"] == "PM2"
    basis = r["basis"][0]
    assert basis["code"] == "PM2" and basis["met"] is True
    assert basis["reasoning"], "an audit answer with no reasoning is not an audit answer"


def test_explore_missing_evidence_names_an_actionable_next_step(run):
    cases = M.explore_missing_evidence(run["results_json"], "PIGA")["cases"]
    up = [s for s in cases[0]["would_change_with"] if s["direction"] == "up"]
    assert up, "a VUS should report what would raise it"
    needs = " ".join(c["needs"] for c in up[0]["candidates"])
    assert "parental" in needs, "the answer must name a concrete next step, not a bare strength"


def test_explore_evidence_sources_views(run):
    rj = run["results_json"]
    assert M.explore_evidence_sources(rj, "clinvar")["view"] == "clinvar"
    assert M.explore_evidence_sources(rj, "missense")["view"] == "missense"
    # the audit view must surface MORE than the met-only view
    met_only = len(M.explore_evidence_sources(rj, "missense")["findings"])
    audit = len(M.explore_evidence_sources(rj, "missense", include_not_met=True)["findings"])
    assert audit >= met_only


def test_explore_evidence_sources_rejects_an_unknown_view(run):
    assert "error" in M.explore_evidence_sources(run["results_json"], "nonsense")


def test_setup_doc_lists_every_registered_mcp_tool():
    """SETUP.md's "Which tools the MCP server exposes" named 9 of 17 — the four explore_* tools
    SKILL.md stage 8 depends on were all missing. Assert BOTH directions so the list cannot go
    stale in either: nothing registered is unlisted, nothing listed is unregistered."""
    import pathlib
    import re

    server = pathlib.Path("src/vcf2report/mcp_server.py").read_text()
    registered = set(re.findall(r"@mcp\.tool\(\)\s*\ndef\s+(\w+)", server))
    assert registered, "no @mcp.tool() functions found — has the decorator changed?"

    doc = pathlib.Path("docs/SETUP.md").read_text()
    section = doc.split("## Which tools the MCP server exposes", 1)[1].split("\n## ", 1)[0]
    listed = set(re.findall(r"`(\w+)`", section))

    assert not registered - listed, f"registered but undocumented: {sorted(registered - listed)}"
    assert not (listed - registered) - {"tests/test_mcp_explore.py"}, \
        f"documented but not registered: {sorted(listed - registered)}"


def test_the_server_class_import_survives_both_sdk_generations():
    """The SDK removed `mcp.server.fastmcp.FastMCP` in v2 in favour of `mcp.server.MCPServer`.
    pyproject declares `mcp[cli]>=1.2`, so CI resolved v2 the day it shipped and the module
    raised SystemExit at import — an upstream rename taking the whole MCP surface down. The
    decorator API is identical across both, so the import accepts either."""
    import mcp.server

    available = [n for n in ("FastMCP", "MCPServer")
                 if hasattr(mcp.server, n)
                 or (n == "FastMCP" and _has_fastmcp_module())]
    assert available, "neither SDK generation's server class is importable"
    assert M.mcp is not None and hasattr(M.mcp, "tool")


def _has_fastmcp_module():
    try:
        import mcp.server.fastmcp  # noqa: F401
        return True
    except ImportError:
        return False
