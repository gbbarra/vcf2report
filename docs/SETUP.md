# Setup — run vcf2report from Claude

There are two Claude surfaces. **Claude Code is the quickest** (a one-step skill);
**Claude Desktop** adds a natural-language chat over the MCP server.

## Option A — Claude Code (recommended, one step)

Install the guided `analyze-vcf` skill once — it works in any session and
bootstraps everything (clone, install, run, render) itself:

```bash
mkdir -p ~/.claude/skills/analyze-vcf && curl -fsSL \
  https://raw.githubusercontent.com/gbbarra/vcf2report/main/.claude/skills/analyze-vcf/SKILL.md \
  -o ~/.claude/skills/analyze-vcf/SKILL.md
```

Restart Claude Code, then say *"analyze this VCF: /path/to/exome.vcf"* (or `/analyze-vcf`).
That's the whole setup — see [../vcf2report.md](../vcf2report.md). Everything below is
only for the Claude Desktop / MCP path.

## Option B — Claude Desktop (natural-language chat via MCP)

vcf2report runs **locally**. For the chat surface use **Claude Desktop** (not
claude.ai web — a browser can't read local files, run the annotators, or keep the
databases local). The MCP server runs as a local subprocess and is the bridge
between Claude and the tools/databases.

**Privacy — safe by default.** No patient data ever leaves the machine unless you
opt in: outbound gnomAD/NCBI lookups are **off** by default and only run when you
set `VCF2REPORT_ALLOW_NETWORK=1`. The VCF file itself is never transmitted; when
network is enabled, only *individual variant coordinates* (chrom-pos-ref-alt) are
sent to those public APIs. For fully local operation, annotate with the local
toolchain (below) and leave network off (or set `OFFLINE=1` to hard-guarantee it).

## 1. Install

```bash
git clone <repo> && cd vcf2report
python -m pip install -e ".[mcp]"        # engine + MCP SDK (jinja2 included)
# annotation toolchain (all MIT/permissive), via bioconda:
conda install -c bioconda bcftools snpeff vcfanno htslib
```

## 2. One-time data download (for real exomes)

```bash
scripts/setup_data.sh ./annotation_data                 # full
scripts/setup_data.sh ./annotation_data --panel genes.bed   # small subset (laptop/demo)
```

Downloads ClinVar, the SnpEff GRCh38 DB, and points you at gnomAD/ABraOM. Then
edit `scripts/vcfanno.conf.toml` so its `file` paths point at `./annotation_data`.
(You can skip this entirely to run the **offline synthetic demo** — the bundled
`data/` already has everything for that.)

## 3. Register the MCP server in Claude Desktop

Config file: macOS `~/Library/Application Support/Claude/claude_desktop_config.json`,
Windows `%APPDATA%\Claude\claude_desktop_config.json`. Copy the `mcpServers` block
from `claude_desktop_config.example.json`, fixing the absolute paths, then fully
quit and reopen Claude Desktop. `command` must be the interpreter where you ran
`pip install` (use the venv's python if you used one).

Verify: in a new chat, ask Claude to call **`data_status`** — it reports which
tools/data are ready.

## 4. Add the Skill

The Claude Desktop clinical SOP lives in `.claude/skills/` (`vcf2report-orchestrator`,
`acmg-classify`, `variant-report`) — these drive the flow through the MCP tools.
(The separate `analyze-vcf` skill is the terminal-based harness for Option A / Claude
Code and needs no MCP server.) Add the Desktop skills via the **Skills** capability
(Settings → Capabilities in the Claude apps), or run from Claude Code where
`.claude/skills/` is discovered automatically. If Skills aren't available in your
setup, the MCP tool descriptions still drive the flow.

## 5. Analysis flow (what the bench scientist does)

1. Open Claude Desktop and say: *"Analyze this VCF for a patient with seizures and
   developmental delay: /path/to/exome.vcf"* (or point to a GA4GH phenopacket —
   convert it with `scripts/phenopacket_to_inputs.py`).
2. The orchestrator calls **`annotate_and_report(vcf, hpo_terms, reference)`**:
   - already annotated (ANN/CSQ or gnomAD in INFO) → classify directly;
   - raw VCF + tools + a GRCh38 reference → annotate locally (SnpEff + vcfanno)
     first.
3. Get the draft report inline (tiers, auditable ACMG per variant, ABraOM
   callout, per-stage timings).
4. Drill into any candidate with the live `clinvar_lookup` / `gnomad_frequency`
   tools for up-to-the-minute confirmation.
5. Review, correct, and sign out in the lab's template.

## Which tools the MCP server exposes

`data_status`, `annotate_and_report`, `run_report`, `parse_vcf`,
`classify_variant`, `gnomad_frequency`, `clinvar_lookup`, `abraom_frequency`,
`hpo_phenotype_match`. Bulk annotation is local (fast, private); the live
gnomAD/ClinVar tools are for last-mile freshness on the shortlist — see
[ANNOTATION.md](ANNOTATION.md) for why local scales and live doesn't.
