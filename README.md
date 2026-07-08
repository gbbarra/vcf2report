# vcf2report

**Turn a raw exome VCF into an auditable, clinical-grade ACMG variant report —
driven from Claude Desktop, no terminal required.**

Built with Claude — Life Sciences Hackathon.

Exome interpretation today is slow and manual: an analyst hand-filters thousands
of variants, cross-references several databases, and matches each candidate to
the patient's phenotype one by one. vcf2report turns that into minutes: you point
Claude Desktop at a VCF and its phenotype, and it runs the whole pipeline — parse,
QC, annotate, prioritize, classify, report — and hands back a **draft for expert
review**.

## What makes it different

- **Auditable ACMG, not a black box.** 20 of the 28 ACMG/AMP criteria are
  evaluated (the rest — PS2/PM3/PM6/PP1/BS4 and other trio/segregation/phasing
  rules — require data a single-proband VCF can't provide, and are reported N/A)
  each with the concrete evidence value, the source (DB + accession + date),
  and a one-line reason. Deterministic criteria are decided by the engine; genuine
  judgment criteria are tagged for expert/model adjudication. You see *why* a
  variant is Pathogenic, and the exact combining-rule path to the tier.
- **Brazilian population frequencies (ABraOM).** Alongside gnomAD, vcf2report
  checks the ABraOM (SABE) Brazilian cohort. A variant absent from gnomAD but
  common in admixed Brazilians is correctly down-weighted — the report names the
  spurious candidates a gnomAD-only pipeline would have kept.
- **Open and composable.** Plain Python engine + local datasets; swap in your
  lab's internal variant DB, thresholds, and report template.

## Architecture

```
Bench scientist (Claude Desktop)
   │  natural language + VCF path
   ▼
Agent Skills (the clinical SOP)  ──►  vcf2report MCP server (thin tool adapter)
                                          ▼
                                     vcf2report Python package
                                     parse ▸ qc ▸ annotate ▸ filter ▸ ACMG ▸ report
                                          ├─ real APIs: gnomAD, ClinVar, HPO
                                          └─ local: ClinVar slice, gnomAD cache,
                                             ABraOM, HPO, gene constraint, in-silico
```

The MCP server is a thin adapter; all logic lives in the importable `vcf2report`
package so it runs and is tested **headless**, without Claude. See
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Quick start (headless)

```bash
pip install -e .                     # core deps only (jinja2)
python scripts/run_headless.py       # runs on the bundled sample VCF
# -> writes data/out/sample_exome_report.md
```

Run on your own VCF + phenotype:

```bash
python scripts/run_headless.py path/to/exome.vcf --hpo hpo_terms.txt --stdout
```

## Use it from Claude Desktop

1. `pip install -e ".[mcp]"` to install the MCP SDK.
2. Copy the `mcpServers` block from
   [`claude_desktop_config.example.json`](claude_desktop_config.example.json)
   into your Claude Desktop config, fixing the absolute paths.
3. Restart Claude Desktop and say: *"Analyze this VCF for a patient with seizures
   and developmental delay: /path/to/exome.vcf"*.

The Agent Skills in `.claude/skills/` encode the SOP; the MCP tools (`parse_vcf`,
`gnomad_frequency`, `clinvar_lookup`, `abraom_frequency`, `hpo_phenotype_match`,
`classify_variant`, `run_report`) do the work.

## Demo

The bundled synthetic sample (`data/sample/sample_exome.vcf`, phenotype: seizures
+ developmental delay) yields:

| Gene | Variant | ACMG tier | Bucket | Why |
|---|---|---|---|---|
| SCN1A | p.Arg612Ter | **Pathogenic** | Primary | PVS1 + PM2 + PP4 + PP5 (matches the seizure phenotype) |
| KCNQ2 | p.Arg213Trp | **VUS** | Primary | PM2 + PP3 + PP4 + PP5 (ClinVar P is supporting, not strong) |
| CACNA1A | p.Ser34Pro | **VUS** | Primary | PM2 only |
| LDLR | p.Arg350Ter | **Likely Pathogenic** | **Secondary (ACMG SF)** | PVS1 + PM2 — actionable incidental finding on an SF v3.2 gene |
| PAX6 | p.Arg68Ter | **Likely Pathogenic** | Other | PVS1 + PM2 — incidental, but PAX6 is not on the ACMG SF list |
| OBSCN | p.Val100Ile | *dropped* | — | common in ABraOM, absent in gnomAD |

The report separates **primary** (phenotype-related) from **secondary** findings
gated on the real **ACMG SF v3.2** gene list — so LDLR (familial hypercholesterolemia,
actionable) is a reportable secondary finding while PAX6 (aniridia, not on the SF list)
is routed to "other". Each call carries its full sourced ACMG derivation.

See [docs/DEMO_SCRIPT.md](docs/DEMO_SCRIPT.md) for the walkthrough.

## Testing

```bash
pip install -e ".[dev]" && pytest
```

## Disclaimer

vcf2report is a **draft-generation aid, not a diagnostic device**. All output must
be reviewed and signed out by a qualified professional. See
[docs/DISCLAIMERS.md](docs/DISCLAIMERS.md).
