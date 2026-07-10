# Architecture

## Principle

The **MCP server is a thin adapter**. All real logic lives in the importable
`vcf2report` package, so the pipeline runs and is unit-tested headless — no Claude,
no network required. Skills orchestrate; MCP tools execute; the package does the
work.

```
src/vcf2report/
├── config.py        paths, dataset locations, API endpoints, thresholds
├── models.py        Variant, Annotation, CriterionResult, Classification, QCSummary
├── pipeline.py      the one orchestration entry point (parse→qc→annotate→filter→acmg→report)
├── cli.py           headless CLI (also the `vcf2report` console script)
├── mcp_server.py    FastMCP tool wrappers for Claude Desktop
├── vcf/             parse (pure-python + cyvcf2 fallback), qc, filter (tiering)
├── annotate/        gnomad, clinvar, abraom, hpo, alphamissense, extra (constraint/in-silico), cache
├── acmg/            criteria (28 evaluators), rules (Richards 2015 Table 5), engine
├── report/          assemble (ReportModel), render (Jinja2 + built-in fallback)
└── concordance.py   ClinVar-vs-engine validation panel (see docs/CONCORDANCE.md)
```

Two thin adapters drive the same package: the Claude Code **`analyze-vcf` skill**
(a terminal harness — clone/install/run/render, no MCP needed) and the Claude
Desktop **MCP server** (natural-language chat).

## Data flow

1. **parse** — read VCF (pure-Python reader by default; `cyvcf2` if installed),
   split multiallelics, detect genome build.
2. **qc** — drop by FILTER, DP, GQ, het allele balance; record reasons.
3. **annotate** — merge gnomAD (popmax AF + filtering AF, homozygotes), ClinVar
   (significance, accession), ABraOM (Brazilian AF), gene constraint (pLI/LOEUF),
   in-silico (**AlphaMissense** at ClinGen-calibrated PP3/BP4 strength, with
   REVEL/CADD as fallback), HPO phenotype match. Every field records its source.
   (AlphaMissense is looked up lazily — only for the surviving candidates.)
4. **filter** — funnel: rarity (gnomAD **and** ABraOM) → coding/splice impact →
   phenotype ranking. ClinVar P/LP bypass the funnel. Records ABraOM-specific drops.
5. **acmg** — evaluate 20 of the 28 criteria (the rest need trio/segregation/
   phasing data unavailable from a single proband and are reported N/A), apply
   the combining rules → 5-tier call
   with the full auditable trail.
6. **report** — assemble a `ReportModel`, render Markdown (native in Claude Desktop).

## Real APIs + local fallback

Each annotator resolves: on-disk cache (`data/cache/`) → live API (unless
`OFFLINE=1`) → bundled local snapshot. The live gnomAD (GraphQL) and ClinVar
(E-utilities) clients use only the standard library (`urllib`), so real calls need
no extra dependency. Set `OFFLINE=1` to force cache/local-only.
`scripts/warm_cache.py` pre-fills the cache so a demo is network-independent while
the code stays capable of real calls.

| Source | Live | Local fallback |
|---|---|---|
| gnomAD | GraphQL (`gnomad_r4`) | `data/gnomad/gnomad_cache.json` |
| ClinVar | NCBI E-utilities | `data/clinvar/clinvar_grch38_slice.tsv` |
| ABraOM | — (static dataset) | `data/abraom/abraom_sabe.tsv` |
| HPO | ontology.jax.org | `data/hpo/genes_to_phenotype.tsv` |

## Auditability

`CriterionResult` carries `code, met, applied_strength, evidence, citation,
reasoning, confidence, adjudicated_by`. `adjudicated_by` separates **engine**
(deterministic) from **model** (judgment) criteria — the report shows both, so a
reviewer can see exactly what the machine decided vs. what needs human/model
judgment. The combining rules (`acmg/rules.py`) are pure boolean logic over
strength counts, matching Richards 2015 Table 5, and return the rule label used
verbatim in the report.

## Known limitations / production notes

- **Genome build**: everything assumes GRCh38 (to match gnomAD r4 / ClinVar
  GRCh38). The pipeline warns on a mismatched or undeclared build.
- **Normalization**: production should run `bcftools norm` against a GRCh38 FASTA;
  the demo ships pre-normalized to avoid the ~3 GB reference.
- **Consequence annotation**: expected in the VCF INFO (VEP/SnpEff) or via ClinVar;
  the sample carries `CSQ`.
- **Single proband**: PS2/PM3/PM6/PP1/BS4 are reported N/A without trio/segregation.
