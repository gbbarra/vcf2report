# Synthetic case definitions (5 diverse DEE exomes)

Per-sample **spike targets** and **HPO terms** for the synthetic-case pipeline
(`scripts/spike_pathogenic.py` + `make_synthetic_exomes.sh`). Each case is a real,
healthy 1000G exome (real diverse background) with two **real ClinVar pathogenic
variants** spiked at their true GRCh38 coordinates: one **primary** developmental
& epileptic encephalopathy (DEE) gene that matches the seizure phenotype, and one
**secondary** ACMG SF v3.2 gene that is an unrelated, actionable incidental finding.

The gene choices are **validated against vcf2report's engine + the full real HPO
ontology + the real gnomAD v2.1.1 constraint table**, so each routes correctly:

* **Primary** genes are LoF-intolerant (→ PVS1 on a spiked LoF) *and* share HPO
  terms with the seizure phenotype (→ PP4), so they land Pathogenic → **primary**.
* **Secondary** genes are LoF-intolerant (→ reach Pathogenic/LP via PVS1) *and*
  have **zero** HPO overlap with the seizure phenotype, so they route to
  **secondary** instead of leaking into primary.

HPO is deliberately **seizure-specific** (Seizure / Status epilepticus / Epileptic
spasms) — the generic "global developmental delay" term (HP:0001263) makes several
SF cancer genes match the phenotype under the full HPO and mis-route.

| Sample | Super-pop (1000G id) | Primary (DEE) | LOEUF | Secondary (ACMG SF) | LOEUF |
|---|---|---|---|---|---|
| SYN-001 | EUR / CEU (NA12878) | SCN1A | 0.07 | RB1 (retinoblastoma) | 0.13 |
| SYN-002 | AFR / YRI (NA19240) | KCNQ2 | 0.16 | APC (FAP) | 0.16 |
| SYN-003 | EAS / JPT (NA18939) | SCN2A | 0.13 | SMAD4 (JP/HHT) | 0.22 |
| SYN-004 | SAS / GIH (NA20845) | STXBP1 | 0.09 | WT1 (Wilms/Denys-Drash) | 0.25 |
| SYN-005 | AMR / CLM (HG01112) | SLC2A1 | 0.24 | FBN1 (Marfan) | 0.05 |

NA12878 is confirmed in 1000G + GIAB. Verify the other DRAGEN ids exist and swap
freely — the super-population labels are the intent, not a guarantee.

All secondary genes are on the ACMG SF v3.2 list and LoF-intolerant; all avoid the
seizure HPO terms. If you swap a gene, re-check both properties (there's a helper
in the repo history: match against the seizure HPO set + `gene_constraint`).

**Not clinical data.** Synthetic, de-identified; every output is marked
demonstration-only in-header.
