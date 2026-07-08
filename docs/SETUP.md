# Setup & Claude Desktop integration

vcf2report runs **locally** so patient VCFs never leave the machine. The right
surface is **Claude Desktop** (not claude.ai web — a browser can't read local
files, run the annotators, or keep the databases local). The MCP server runs as a
local subprocess and is the bridge between Claude and the tools/databases.

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

The clinical SOP lives in `.claude/skills/` (`vcf2report-orchestrator`,
`acmg-classify`, `variant-report`). Add it to Claude via the **Skills** capability
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
