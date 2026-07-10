# Variant Interpretation Report — sample_exome

> **DRAFT — for expert review. Not for clinical use.** Auto-generated candidate
> interpretation to be verified and signed out by a qualified professional.

- **Genome build:** GRCh38
- **Pipeline:** vcf2report v0.1.0
- **Generated:** 2026-07-09T21:22:36+00:00
- **Patient HPO terms:** HP:0001250, HP:0002133, HP:0011097

## Quality control & filtering funnel

- Total variants: **11**
- PASS filter: **10**
- After QC (DP/GQ/AB): **9**
- After rarity: **7**
- After coding/splice impact: **5**
- **Candidates classified: 5**

### Brazilian-frequency filtering (ABraOM)

Spurious candidates a gnomAD-only pipeline would have kept, removed using ABraOM (SABE) local frequencies:
- OBSCN (p.Val100Ile): gnomAD AF=0.000000 but ABraOM AF=0.0320 — common in Brazilians, dropped


## Primary (diagnostic) findings

_Variants in genes overlapping the patient's phenotype._
| Gene | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | ABraOM AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|
| SCN1A | c.1834C>T p.Arg612Ter | het | stop_gained | Pathogenic | 0.000000 | 0.000000 | 1.0 | **Pathogenic** |
| KCNQ2 | c.637C>T p.Arg213Trp | het | missense_variant | Pathogenic | 0.000000 | 0.000000 | 1.0 | **Uncertain Significance (VUS)** |
| CACNA1A | c.100T>C p.Ser34Pro | het | missense_variant | — | 0.000000 | 0.000000 | 0.667 | **Uncertain Significance (VUS)** |


## Secondary findings (ACMG SF v3.2)

_P/LP variants in ACMG SF v3.2 genes, unrelated to the indication — reportable actionable secondary findings, subject to the patient's opt-in policy._
| Gene | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | ABraOM AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|
| RB1 | c.958C>T p.Arg320Ter | het | stop_gained | — | 0.000000 | 0.000000 | 0.0 | **Likely Pathogenic** |

## Other candidates

_Incidental P/LP not on the ACMG SF list, plus phenotype-unrelated uncertain/benign candidates. Not routinely reported._
| Gene | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | ABraOM AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|
| PAX6 | c.202C>T p.Arg68Ter | het | stop_gained | — | 0.000000 | 0.000000 | 0.0 | **Likely Pathogenic** |


## Per-variant ACMG rationale (auditable)

### SCN1A — p.Arg612Ter → Pathogenic

**Rule path:** `PVS1 + PM2 + PP4 + PP5 => Pathogenic [PATH-1 (PVS1 + strong/moderate/supporting)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True | gnomAD v2.1.1 LoF constraint (local) | engine | stop_gained is loss-of-function and SCN1A is LoF-intolerant |
| **PS1** | — | strong | hgvs_p=p.Arg612Ter | — | model | Requires a residue-level cross-match to a distinct ClinVar pathogenic variant — model adjudication (own ClinVar record is PP5) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, abraom_af=0.0 | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg612Ter | — | model | Domain/hotspot membership requires curated annotation — model adjudication |
| **PM2** | ✅ met | moderate | gnomad_af=0.0, abraom_af=0.0, ceiling=0.0001, moi=AD | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD popmax AF=0.000000, ABraOM AF=0.000000 — both at/under 0.0001 (SCN1A is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg612Ter | — | model | Requires residue-level ClinVar cross-check — model adjudication |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=SCN1A | — | model | Gene-level missense constraint requires curated metric — model adjudication |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | ✅ met | supporting | hpo_match_score=1.0, matched_terms=['HP:0001250 (Seizure)', 'HP:0002133 (Status epilepticus)', 'HP:0011097 (Epileptic spasm)'], cutoff=0.6 | HPO genes_to_phenotype (local) | engine | phenotype match 1.00 (terms: HP:0001250 (Seizure), HP:0002133 (Status epilepticus), HP:0011097 (Epileptic spasm)) |
| **PP5** | ✅ met | supporting | clinvar=Pathogenic, review_status=criteria_provided,_multiple_submitters,_no_conflicts | VCV000012345 | engine | ClinVar Pathogenic (criteria_provided,_multiple_submitters,_no_conflicts) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (SCN1A is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | — | engine | in-silico predictors not benign / unavailable |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### PAX6 — p.Arg68Ter → Likely Pathogenic

**Rule path:** `PVS1 + PM2 => Likely Pathogenic [LP-1 (PVS1 + 1 Moderate)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True | gnomAD v2.1.1 LoF constraint (local) | engine | stop_gained is loss-of-function and PAX6 is LoF-intolerant |
| **PS1** | — | strong | hgvs_p=p.Arg68Ter | — | model | Requires a residue-level cross-match to a distinct ClinVar pathogenic variant — model adjudication (own ClinVar record is PP5) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, abraom_af=0.0 | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg68Ter | — | model | Domain/hotspot membership requires curated annotation — model adjudication |
| **PM2** | ✅ met | moderate | gnomad_af=0.0, abraom_af=0.0, ceiling=0.0001, moi=AD | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD popmax AF=0.000000, ABraOM AF=0.000000 — both at/under 0.0001 (PAX6 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg68Ter | — | model | Requires residue-level ClinVar cross-check — model adjudication |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=PAX6 | — | model | Gene-level missense constraint requires curated metric — model adjudication |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | — | supporting | hpo_match_score=0.0, matched_terms=[], cutoff=0.6 | HPO genes_to_phenotype (local) | engine | phenotype match 0.00 below 0.6 |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (PAX6 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | — | engine | in-silico predictors not benign / unavailable |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### RB1 — p.Arg320Ter → Likely Pathogenic

**Rule path:** `PVS1 + PM2 => Likely Pathogenic [LP-1 (PVS1 + 1 Moderate)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True | gnomAD v2.1.1 LoF constraint (local) | engine | stop_gained is loss-of-function and RB1 is LoF-intolerant |
| **PS1** | — | strong | hgvs_p=p.Arg320Ter | — | model | Requires a residue-level cross-match to a distinct ClinVar pathogenic variant — model adjudication (own ClinVar record is PP5) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, abraom_af=0.0 | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg320Ter | — | model | Domain/hotspot membership requires curated annotation — model adjudication |
| **PM2** | ✅ met | moderate | gnomad_af=0.0, abraom_af=0.0, ceiling=0.0001, moi=AD | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD popmax AF=0.000000, ABraOM AF=0.000000 — both at/under 0.0001 (RB1 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg320Ter | — | model | Requires residue-level ClinVar cross-check — model adjudication |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=RB1 | — | model | Gene-level missense constraint requires curated metric — model adjudication |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | — | supporting | hpo_match_score=0.0, matched_terms=[], cutoff=0.6 | HPO genes_to_phenotype (local) | engine | phenotype match 0.00 below 0.6 |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (RB1 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | — | engine | in-silico predictors not benign / unavailable |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### KCNQ2 — p.Arg213Trp → Uncertain Significance (VUS)

**Rule path:** `PM2 + PP3 + PP4 + PP5 => criteria insufficient for a benign or pathogenic call => VUS`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | — | very_strong | consequence=missense_variant, gene_lof_intolerant=True | gnomAD v2.1.1 LoF constraint (local) | engine | missense_variant is not a qualifying null variant in a LoF-intolerant gene |
| **PS1** | — | strong | hgvs_p=p.Arg213Trp | — | model | Requires a residue-level cross-match to a distinct ClinVar pathogenic variant — model adjudication (own ClinVar record is PP5) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, abraom_af=0.0 | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=missense_variant, hgvs_p=p.Arg213Trp | — | model | Domain/hotspot membership requires curated annotation — model adjudication |
| **PM2** | ✅ met | moderate | gnomad_af=0.0, abraom_af=0.0, ceiling=0.0001, moi=AD | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD popmax AF=0.000000, ABraOM AF=0.000000 — both at/under 0.0001 (KCNQ2 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=missense_variant | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg213Trp | — | model | Requires residue-level ClinVar cross-check — model adjudication |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=missense_variant, gene=KCNQ2 | — | model | Gene-level missense constraint requires curated metric — model adjudication |
| **PP3** | ✅ met | supporting | revel=0.92, cadd_phred=32.0, revel_cutoff=0.7, cadd_cutoff=20.0 | REVEL/CADD (local) | engine | REVEL=0.92, CADD=32.0 above deleterious cutoffs |
| **PP4** | ✅ met | supporting | hpo_match_score=1.0, matched_terms=['HP:0001250 (Seizure)', 'HP:0002133 (Status epilepticus)', 'HP:0011097 (Epileptic spasm)'], cutoff=0.6 | HPO genes_to_phenotype (local) | engine | phenotype match 1.00 (terms: HP:0001250 (Seizure), HP:0002133 (Status epilepticus), HP:0011097 (Epileptic spasm)) |
| **PP5** | ✅ met | supporting | clinvar=Pathogenic, review_status=criteria_provided,_single_submitter | VCV000067890 | engine | ClinVar Pathogenic (criteria_provided,_single_submitter) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (KCNQ2 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BP4** | — | supporting | revel=0.92, cadd_phred=32.0, revel_cutoff=0.15, cadd_cutoff=10.0 | — | engine | in-silico predictors not benign / unavailable |
| **BP7** | — | supporting | consequence=missense_variant | — | engine | not a synonymous variant |

### CACNA1A — p.Ser34Pro → Uncertain Significance (VUS)

**Rule path:** `PM2 + PP4 => criteria insufficient for a benign or pathogenic call => VUS`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | — | very_strong | consequence=missense_variant, gene_lof_intolerant=True | gnomAD v2.1.1 LoF constraint (local) | engine | missense_variant is not a qualifying null variant in a LoF-intolerant gene |
| **PS1** | — | strong | hgvs_p=p.Ser34Pro | — | model | Requires a residue-level cross-match to a distinct ClinVar pathogenic variant — model adjudication (own ClinVar record is PP5) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, abraom_af=0.0 | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=missense_variant, hgvs_p=p.Ser34Pro | — | model | Domain/hotspot membership requires curated annotation — model adjudication |
| **PM2** | ✅ met | moderate | gnomad_af=0.0, abraom_af=0.0, ceiling=0.0001, moi=AD | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD popmax AF=0.000000, ABraOM AF=0.000000 — both at/under 0.0001 (CACNA1A is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=missense_variant | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Ser34Pro | — | model | Requires residue-level ClinVar cross-check — model adjudication |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=missense_variant, gene=CACNA1A | — | model | Gene-level missense constraint requires curated metric — model adjudication |
| **PP3** | — | supporting | revel=0.45, cadd_phred=18.0, revel_cutoff=0.7, cadd_cutoff=20.0 | REVEL/CADD (local) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | ✅ met | supporting | hpo_match_score=0.667, matched_terms=['HP:0001250 (Seizure)', 'HP:0002133 (Status epilepticus)'], cutoff=0.6 | HPO genes_to_phenotype (local) | engine | phenotype match 0.67 (terms: HP:0001250 (Seizure), HP:0002133 (Status epilepticus)) |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/ABraOM popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); ABraOM SABE (not observed) | engine | gnomAD/ABraOM popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (CACNA1A is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BP4** | — | supporting | revel=0.45, cadd_phred=18.0, revel_cutoff=0.15, cadd_cutoff=10.0 | — | engine | in-silico predictors not benign / unavailable |
| **BP7** | — | supporting | consequence=missense_variant | — | engine | not a synonymous variant |

## Methods
- **genome_build:** GRCh38
- **qc_thresholds:** {'min_DP': 10, 'min_GQ': 20}
- **rarity_cutoff_popmax_af:** 0.005
- **ba1_cutoff:** 0.05
- **databases:** ['ClinVar', 'gnomAD r4', 'ABraOM (SABE)', 'HPO', 'gnomAD constraint']
- **standards:** ['ACMG/AMP variant classification (Richards et al., Genet Med 2015)', 'ClinGen SVI criteria refinements', 'ACMG secondary-findings list (SF v3.2, Miller et al. 2023)', 'HGVS nomenclature', 'GA4GH Phenopackets (phenotype exchange)']

## Performance (this run)
- **parse:** 0.0004 s- **qc:** 0.0 s- **annotate:** 0.2367 s- **filter:** 0.0 s- **alphamissense:** 0.001 s- **classify:** 0.0002 s- **total:** 0.2383 s- **variants per:** 46.2
## Limitations & disclaimers

- Single-proband analysis: criteria requiring parental/segregation/phasing data (PS2, PM3, PM6, PP1, BS4) are reported as N/A.
- Judgment criteria (PS3, PS4, PM1, PM5, PP2) are surfaced for expert/model adjudication and default to not-met unless explicitly supported.
- Population and clinical databases are versioned snapshots; re-check before sign-out.
- **This is a draft-generation aid, not a diagnostic device.**