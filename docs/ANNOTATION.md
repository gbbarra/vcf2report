# Annotating a real exome for vcf2report

vcf2report classifies and reports; it does **not** re-implement variant
annotation. A raw VCF from a caller (GATK, DeepVariant) has no gene / consequence
/ HGVS / population frequency — annotate it first, then feed the annotated VCF in.

## Recommended toolchain (all open-source, MIT/permissive, fully local)

| Step | Tool | License | Adds |
|---|---|---|---|
| Normalize | **bcftools norm** | MIT | split multiallelics, left-align indels |
| Consequence + HGVS | **SnpEff** | MIT | `INFO/ANN` (gene, consequence, HGVS c./p.) |
| Population + clinical | **vcfanno** | MIT | gnomAD AF, ClinVar, REVEL/CADD, ABraOM from local files |

(Ensembl **VEP** — Apache-2.0 — is a fine alternative to SnpEff; vcf2report reads
its `CSQ` too. Avoid ANNOVAR: not open-source.)

**Why local, not live APIs:** a live gnomAD/ClinVar call per variant is
rate-limited and network-bound — infeasible for the ~20–100k variants in an
exome (hours). Annotating against local files is O(1) per variant (minutes total).
It is also **private**: local annotation sends nothing off-machine. The live MCP
tools (`gnomad_frequency`, `clinvar_lookup`) are for interactive drill-down on the
final shortlist, and they only run when you opt in with `VCF2REPORT_ALLOW_NETWORK=1`
— at which point only that variant's coordinates are sent. By default the whole
pipeline is offline.

## One command

```bash
scripts/annotate_vcf.sh raw.vcf.gz GRCh38.fa out.annotated.vcf.gz
python scripts/run_headless.py out.annotated.vcf.gz --hpo patient_hpo_terms.txt
```

`scripts/vcfanno.conf.toml` maps each data file's fields to the INFO names
vcf2report expects (`config.INFO_ALIASES`): `gnomad_AF`, `CLNSIG`, `CLNREVSTAT`,
`CLNDN`, `ABraOM_AF`, `REVEL`, `CADD_PHRED`. If your annotation uses different
INFO keys, add them to `INFO_ALIASES` — no code change needed elsewhere.

## Data files (GRCh38, free)

- **gnomAD v4 sites VCF** — https://gnomad.broadinstitute.org/downloads
- **ClinVar weekly VCF** — https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
- **ABraOM (SABE)** — http://abraom.ib.usp.br/ (convert to a bgzipped, tabix VCF)
- **REVEL / CADD** (optional) — from their sites or dbNSFP; build a small VCF.

## Expected timing (real exome, local)

parse (cyvcf2) seconds · **SnpEff ~1–3 min · vcfanno ~1–2 min (the bottleneck)** ·
vcf2report filter + ACMG + report seconds. End-to-end ≈ a few minutes. Use
`--timing` on `run_headless.py` to measure each stage on your machine.
