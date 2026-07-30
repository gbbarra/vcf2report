"""Documentation-vs-reality guards.

A docs audit found 40+ claims the code did not support: a README quickstart block whose
pasted output no longer matched, a demo script naming a gene absent from the sample VCF, an
inverted description of the network default, a licence stated two incompatible ways, and
several counts that had drifted by 2x.

`docs/ACMG_CRITERIA.md` was the ONLY document with zero drift — because
`tests/test_acmg_coverage.py` parses its tables and asserts them against the registry. Where a
test existed, nothing drifted; where none existed, everything did. This file generalises that
pattern to the claims most likely to rot: things a reader types, paths they follow, and output
they are told to expect.

Deliberately NOT pinned here: the suite's own test count (circular, churns every PR) and any
figure that needs the 200-exome benchmark or a multi-GB store (unavailable in CI).
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = sorted((ROOT / "docs").glob("*.md")) + [ROOT / "README.md",
                                               ROOT / ".claude/skills/vcf2report/SKILL.md"]

# Paths a doc may reference that are generated, downloaded, or built on demand.
_GENERATED = (
    "data/gnomad/gnomad_parquet", "data/alphamissense", "data/clinvar/clinvar_parquet",
    "data/out/", "data/real/", "data/tools/", "synthetic_exomes/", "data/benchmark",
    "data/cache/", "data/clinvar/clinvar.tsv", "data/gnomad/gnomad.tsv",
    # Operator-supplied, never shipped: this project redistributes no population cohort.
    "data/local_cohort/",
)


_PACKAGE_SOURCE = "\n".join(
    p.read_text() for p in (ROOT / "src/vcf2report").rglob("*.py"))


def _fenced_bash(text: str) -> list[str]:
    return [b for b in re.findall(r"```(?:bash|sh|console)\n(.*?)```", text, re.S)]


# --------------------------------------------------------------------------------------
# What a reader types
# --------------------------------------------------------------------------------------

def test_every_documented_flag_exists_in_the_script_that_takes_it():
    """This project once shipped `--vcf` / `--out-dir` in the README when the real interface
    was positional + `--out`. Parse every documented invocation and check its flags."""
    missing = []
    for doc in DOCS:
        for block in _fenced_bash(doc.read_text()):
            for line in block.splitlines():
                line = line.strip()
                m = re.match(r"(?:python3?|bash)?\s*(scripts/[\w./-]+\.py)\s+(.*)", line)
                if not m:
                    continue
                script, rest = ROOT / m.group(1), m.group(2)
                if not script.exists():
                    continue
                # Several scripts are thin wrappers over a module (run_headless.py delegates to
                # vcf2report.cli), so the flag literals live in the package, not the script.
                src = script.read_text() + _PACKAGE_SOURCE
                for flag in re.findall(r"(?<![\w-])(--[a-z][\w-]*)", rest):
                    if f'"{flag}"' not in src and f"'{flag}'" not in src:
                        missing.append(f"{doc.relative_to(ROOT)}: {m.group(1)} {flag}")
    assert not missing, "documented flags that the script does not define:\n" + "\n".join(missing)


def test_shell_scripts_are_invoked_with_bash_and_python_scripts_with_python():
    """`bash scripts/build_v2_biallelic.py` is a syntax error — the file is Python and is not
    executable. Cheap to get wrong, cheap to catch."""
    wrong = []
    for doc in DOCS:
        for block in _fenced_bash(doc.read_text()):
            for line in block.splitlines():
                line = line.strip()
                if m := re.match(r"bash\s+(scripts/[\w./-]+)", line):
                    if m.group(1).endswith(".py"):
                        wrong.append(f"{doc.relative_to(ROOT)}: bash {m.group(1)}")
                if m := re.match(r"python3?\s+(scripts/[\w./-]+)", line):
                    if m.group(1).endswith(".sh"):
                        wrong.append(f"{doc.relative_to(ROOT)}: python {m.group(1)}")
    assert not wrong, "wrong interpreter for a documented script:\n" + "\n".join(wrong)


def test_every_repo_path_named_in_the_docs_exists():
    """Catches a doc pointing at a script or example that was renamed or never existed."""
    absent = []
    pattern = re.compile(r"`((?:src|scripts|tests|templates|docs|data)/[\w./-]+)`")
    for doc in DOCS:
        for path in set(pattern.findall(doc.read_text())):
            if path.startswith(_GENERATED) or path.rstrip("/").endswith(("*", "…")):
                continue
            if not (ROOT / path).exists():
                absent.append(f"{doc.relative_to(ROOT)}: {path}")
    assert not absent, "paths named in docs that do not exist:\n" + "\n".join(absent)


# --------------------------------------------------------------------------------------
# What a reader is told to expect
# --------------------------------------------------------------------------------------

def test_readme_quickstart_block_matches_the_real_demo_output():
    """The README's pasted terminal transcript drifted: it showed KCNQ2 as VUS long after the
    residue index started firing PM1 and lifting it to Likely Pathogenic. A pasted transcript
    nobody re-runs is this project's most reliable source of documentation drift."""
    from vcf2report import config
    from vcf2report.pipeline import run_pipeline

    hpo = [ln.split()[0] for ln in config.SAMPLE_HPO.read_text().splitlines()
           if ln.strip() and not ln.startswith("#")]
    real = {c.variant.gene: c.tier for c in run_pipeline(config.SAMPLE_VCF, hpo_terms=hpo).classifications}

    readme = (ROOT / "README.md").read_text()
    blocks = [b for b in re.findall(r"```[a-z]*\n(.*?)```", readme, re.S)
              if "candidates classified:" in b]
    assert blocks, "README.md no longer quotes a demo output block"
    block = blocks[0]

    quoted = dict(re.findall(r"^\s*-\s+(\w+):\s+(.+?)\s*$", block, re.M))
    assert quoted, "could not find the quoted `- GENE: Tier` lines in README.md"
    assert quoted == real, (
        f"README quickstart block is stale.\n  quoted: {quoted}\n  actual: {real}\n"
        "Re-run `python3 scripts/run_headless.py` and paste the real output.")
    assert f"candidates classified: {len(real)}" in block


def test_demo_script_only_names_genes_that_are_in_the_sample_vcf():
    """DEMO_SCRIPT.md told the presenter to point at `LDLR p.Arg350Ter` as the ACMG SF
    secondary finding. LDLR appears zero times in the sample VCF — the demo failed on stage."""
    from vcf2report import config

    vcf = config.SAMPLE_VCF.read_text()
    demo = (ROOT / "docs/DEMO_SCRIPT.md").read_text()
    section = demo.split("money shots", 1)[1].split("## ", 1)[0]
    # Gene-shaped tokens in bold, e.g. **RB1 p.Arg320Ter → ...**
    named = {m for m in re.findall(r"\*\*([A-Z][A-Z0-9]{2,9})\s+p\.", section)}
    assert named, "no gene names found in the demo's money-shots section"
    absent = sorted(g for g in named if g not in vcf)
    assert not absent, f"DEMO_SCRIPT names genes absent from the sample VCF: {absent}"


# --------------------------------------------------------------------------------------
# Claims that have already drifted once
# --------------------------------------------------------------------------------------

def test_network_is_opt_in_and_the_docs_say_so(monkeypatch):
    """ARCHITECTURE.md described network as opt-OUT (`OFFLINE=1` to disable) when it is
    opt-IN. The code was safe; the doc misdescribed the privacy model to exactly the reader
    assessing data governance."""
    from vcf2report import config

    for var in ("VCF2REPORT_ALLOW_NETWORK", "OFFLINE"):
        monkeypatch.delenv(var, raising=False)
    assert config.allow_network() is False, "the default run must not reach the network"
    monkeypatch.setenv("VCF2REPORT_ALLOW_NETWORK", "1")
    assert config.allow_network() is True

    arch = (ROOT / "docs/ARCHITECTURE.md").read_text()
    assert "VCF2REPORT_ALLOW_NETWORK" in arch
    assert "opt-in" in arch.lower()


def test_alphamissense_licence_is_stated_identically_everywhere():
    """The predictions were relicensed CC BY 4.0 in March 2024; CC BY-NC-SA 4.0 belongs to a
    different artifact (the AlphaMissense Database). README stated both, contradicting itself
    on a redistribution claim — the one kind of drift that is not merely embarrassing."""
    offenders = []
    for f in DOCS + [ROOT / "scripts/fetch_alphamissense.sh"]:
        text = f.read_text()
        for line in text.splitlines():
            if "CC BY-NC-SA" in line and "AlphaMissense" in text:
                # Permitted only where the line explicitly says it is the superseded licence.
                if not re.search(r"supersed|older|Database|Zenodo|relicens", line):
                    offenders.append(f"{f.relative_to(ROOT)}: {line.strip()[:90]}")
    assert not offenders, ("AlphaMissense predictions are CC BY 4.0; these lines still assert "
                           "the superseded NC-SA without saying so:\n" + "\n".join(offenders))


def test_phenotype_routing_uses_the_average_not_the_strongest_match():
    """HOW_IT_WORKS claimed routing used the single strongest match at >=0.50. It uses the
    best-match average at HPO_RELATED_MIN (0.60) — the switch off the max is what cut the
    decoy false-match rate from 62% to 22%, so a silent revert would undo a measured result."""
    from vcf2report import config
    from vcf2report.models import Annotation, Classification, Variant
    from vcf2report.report.assemble import split_findings

    assert config.HPO_RELATED_MIN == 0.6

    def cls(avg, best):
        return Classification(
            variant=Variant(chrom="1", pos=1, ref="A", alt="T", gene="G"),
            annotation=Annotation(hpo_match_score=avg, hpo_best_match=best),
            criteria=[], tier="Pathogenic", rule_path="")

    strong_max_low_avg = cls(0.55, 1.0)
    primary, _s, other = split_findings([strong_max_low_avg])
    assert strong_max_low_avg in other and strong_max_low_avg not in primary

    good_avg = cls(0.65, 0.65)
    primary, _s, _o = split_findings([good_avg])
    assert good_avg in primary


def test_documented_env_vars_are_read_and_read_env_vars_are_documented():
    """Four store-path overrides were the documented way to put a multi-GB store on another
    disk, and none of them appeared in any doc."""
    sources = "\n".join(p.read_text() for p in (ROOT / "src/vcf2report").rglob("*.py"))
    read = {v for v in re.findall(r"VCF2REPORT_[A-Z0-9_]+", sources)}
    documented = "\n".join(f.read_text() for f in DOCS) + (ROOT / ".env.example").read_text()
    # NO_CYVCF2 is a test-harness escape hatch, not user-facing configuration.
    undocumented = sorted(v for v in read - {"VCF2REPORT_NO_CYVCF2"} if v not in documented)
    assert not undocumented, f"env vars the code reads but no doc mentions: {undocumented}"


@pytest.mark.parametrize("doc", ["README.md", "docs/CONCORDANCE.md"])
def test_published_concordance_headline_matches_the_frozen_panel(doc):
    """The panel's numbers are quoted in two places and have drifted before (the docs published
    60% / LoF 84% against a code that produced neither)."""
    panel = ROOT / "data/concordance/ground_truth.tsv"
    if not panel.exists():
        pytest.skip("concordance panel not frozen")
    from vcf2report.concordance import evaluate_panel, load_panel

    metrics = evaluate_panel(load_panel()).metrics
    text = (ROOT / doc).read_text()
    assert metrics["gross_discordances"] == 0
    for pct in (round(metrics["pathogenic_sensitivity"] * 100),
                round(metrics["lof_pathogenic_sensitivity"] * 100)):
        assert f"{pct}%" in text, f"{doc} does not quote {pct}% from the live panel"


# --------------------------------------------------------------------- invocable entry point

def _slash_command_claims() -> list[Path]:
    """Every repo file that tells a reader to type a `/name` slash command."""
    candidates = list((ROOT / "docs").glob("*.md")) + [ROOT / "README.md", ROOT / "vcf2report.md"]
    candidates += list((ROOT / ".claude/skills").rglob("SKILL.md"))
    return [p for p in candidates if p.exists()]


def test_every_advertised_slash_command_actually_exists():
    """Four files told the reader to type `/vcf2report` while `.claude/commands/` did not exist,
    so the app answered "Unknown command" — a skill under `.claude/skills/` is model-invoked and
    does NOT register a typed slash command by itself. Any `/name` the repo advertises must have a
    backing `.claude/commands/name.md`.
    """
    # Built-ins the harness provides; the repo may reference them without shipping a file.
    builtin = {"clear", "help", "config", "init", "review", "compact", "model", "cost",
               "doctor", "login", "logout", "memory", "vim", "terminal-setup", "loop"}
    commands = {p.stem for p in (ROOT / ".claude/commands").glob("*.md")}
    missing: dict[str, list[str]] = {}
    for doc in _slash_command_claims():
        for name in re.findall(r"`/([a-z][a-z0-9-]{2,})`", doc.read_text()):
            if name in builtin or name in commands:
                continue
            missing.setdefault(name, []).append(str(doc.relative_to(ROOT)))
    assert not missing, (
        "slash commands advertised with no .claude/commands/<name>.md to back them "
        f"(the app answers 'Unknown command'): {missing}")


def test_the_vcf2report_command_delegates_to_the_skill():
    """The command must POINT AT the skill, not restate the flow. Two copies of an 8-stage
    clinical procedure is exactly the duplication that drifted everywhere else in this repo."""
    cmd = ROOT / ".claude/commands/vcf2report.md"
    assert cmd.exists(), "the /vcf2report slash command is the documented entry point"
    text = cmd.read_text()
    assert "vcf2report" in text and "Skill" in text, "must invoke the skill"
    # A copy of the flow would have to re-list the stages; a pointer does not.
    assert "SKILL.md" in text, "must name the skill as the source of truth"
    assert len(text.splitlines()) < 60, "this should be a thin pointer, not a second copy"


# ---------------------------------------------------------------- README discoverability

def test_readme_surfaces_the_phenotype_algorithm_before_the_fold():
    """A reader (or an agent) skimming the README must learn HOW phenotype→gene works
    without scrolling past the feature list.

    This is not cosmetic. The mechanism — Lin/IC similarity over the HPO `is_a` graph — was
    documented only in a bullet ~100 lines down, so a session that had the repo but skimmed
    the top proposed building "Exomiser-lite" as if phenotype scoring were missing. The
    capability existed; the README did not say so where anyone would look.
    """
    text = (ROOT / "README.md").read_text()
    fold = text.index("## What it does")
    head = text[:fold].lower()
    assert "hpo" in head
    assert "lin" in head and "information-content" in head, (
        "the README's opening does not name the phenotype-matching algorithm")
    assert "is_a" in head, "the opening does not say the match is over the HPO graph"


def test_readme_internal_anchors_resolve_to_real_headings():
    """A `(#anchor)` link that matches no heading is a dead link that renders as one."""
    text = (ROOT / "README.md").read_text()
    headings = set()
    for line in text.splitlines():
        if line.startswith("#"):
            title = line.lstrip("#").strip()
            slug = re.sub(r"[^a-z0-9 -]", "", title.lower()).replace(" ", "-")
            headings.add(slug)
    broken = [a for a in re.findall(r"\]\(#([a-z0-9-]+)\)", text) if a not in headings]
    assert not broken, f"README links to anchors with no matching heading: {broken}"
