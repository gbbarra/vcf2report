# Running a phenotype-linked case (GA4GH Phenopacket)

Public patient VCFs paired with phenotype are rare (privacy). The practical source
of paired **HPO terms + causative variants** is the **GA4GH / Monarch Phenopacket
Store** (thousands of curated, published cases). A phenopacket carries the
phenotype and the variant coordinates/HGVS — but **not** the molecular
consequence, which is exactly what annotation adds.

## Flow

```bash
# 1. Phenopacket -> a VCF + an HPO-terms file
python scripts/phenopacket_to_inputs.py case.json --out-prefix data/out/case

# 2. Annotate the VCF (consequence + gnomAD/ClinVar) — see ANNOTATION.md
scripts/annotate_vcf.sh data/out/case.vcf GRCh38.fa data/out/case.annotated.vcf.gz

# 3. Classify + report
python scripts/run_headless.py data/out/case.annotated.vcf.gz --hpo data/out/case.hpo.txt
```

## Why step 2 matters (the annotation dependency)

Run on the **raw** phenopacket VCF (no consequence), the bundled SCN1A/Dravet
example (`data/sample/example_phenopacket.json`) classifies the causal variant as
**VUS**: it is retained (ClinVar P/LP bypasses the impact filter) and earns
PM2 + PP4 + PP5, but **PVS1 cannot fire without a molecular consequence**. After
annotation adds `stop_gained`, the same variant becomes **Pathogenic**
(PVS1 + PM2 + PP4 + PP5). This is locked by
`tests/test_phenopacket.py::test_phenopacket_end_to_end_and_annotation_dependency`.

Takeaway: phenopackets give you a real, phenotype-linked case to drive the whole
pipeline; annotation (a local, one-time-setup step) is what completes the
consequence-based ACMG criteria.
