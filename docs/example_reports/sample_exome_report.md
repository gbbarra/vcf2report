# Variant Interpretation Report — sample_exome

> **DRAFT — for expert review. Not for clinical use.** Auto-generated candidate
> interpretation to be verified and signed out by a qualified professional.

- **Genome build:** GRCh38
- **Pipeline:** vcf2report v0.1.0
- **Generated:** 2026-07-28T20:37:14+00:00
- **Patient HPO terms:** HP:0001250, HP:0002133, HP:0011097

## Conclusion (draft interpretation)

- Likely explanatory finding for the clinical indication: **SCN1A — Pathogenic; KCNQ2 — Likely Pathogenic** (in a gene overlapping the patient's phenotype) — confirm and review.
- Reportable **secondary finding** (ACMG SF v3.2 — actionable, subject to the patient's opt-in policy): RB1 — Likely Pathogenic.
- Additional **Pathogenic / Likely Pathogenic** variant(s) not matching the stated phenotype and not on the ACMG SF actionable list: PAX6 — Likely Pathogenic. Clinical relevance to the indication is uncertain — review in context.
- Sequencing depth at called sites is adequate (median 44x); note that a variant-only VCF conveys no breadth of coverage, so poorly-covered regions cannot be assessed from this input.
- Single-proband analysis: de novo / segregation / phasing criteria (PS2, PM3, PM6, PP1, BS4) are N/A — parental or trio testing could upgrade candidates or resolve VUS.
- **Recommended next steps:** expert review and sign-out; orthogonal confirmation (e.g. Sanger) of any reported P/LP variant; segregation / functional evidence to resolve variants of uncertain significance.

## Quality control & filtering funnel

- Total variants: **11**
- PASS filter: **10**
- After QC (DP/GQ/AB): **9**
- After rarity: **8**
- After coding/splice impact: **6**
- **Candidates classified: 6**

## Sequencing quality (estimated from variant sites)

- **Assay (by variant count):** small / demo VCF (11 variants)
- **Depth at called sites:** 41.3x mean / 44x median — 90.9% ≥10x, 90.9% ≥20x
- **Genotype quality:** median 99, 100.0% ≥20
- **Ti/Tv (SNVs):** 4.5 (11 SNVs)
- **Indel:SNV ratio:** 0.0 (0 indels)
- **Multiallelic sites:** 0.0% (0 / 11)
- **Het allele balance:** 100.0% balanced
- **FILTER = PASS:** 90.9%
- _Depth is measured only at called variant sites — a proxy for sequencing quality, not genome-wide breadth of coverage (a variants-only VCF cannot give breadth; a gVCF or BAM would)._
- _Ti/Tv = 4.5 (expected ~3.0 for exome, ~2.0-2.1 for whole genome; a much lower value suggests false-positive calls)._
- _VCF is not dbSNP-annotated (few/no rsIDs in the ID column) — novelty rate not computed (annotate with dbSNP IDs to enable it)._


## Primary (diagnostic) findings

_Variants in genes overlapping the patient's phenotype._
| Gene | Transcript | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | local cohort AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|---|
| SCN1A | — | c.1834C>T p.Arg612Ter | het | stop_gained | Pathogenic | 0.000000 | n/a | 1.0 | **Pathogenic** |
| KCNQ2 | — | c.637C>T p.Arg213Trp | het | missense_variant | Pathogenic | 0.000000 | n/a | 1.0 | **Likely Pathogenic** |
| CACNA1A | — | c.100T>C p.Ser34Pro | het | missense_variant | — | 0.000000 | n/a | 0.93 | **Uncertain Significance (VUS)** |


## Secondary findings (ACMG SF v3.2)

_P/LP variants in ACMG SF v3.2 genes, unrelated to the indication — reportable actionable secondary findings, subject to the patient's opt-in policy._
| Gene | Transcript | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | local cohort AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|---|
| RB1 | — | c.958C>T p.Arg320Ter | het | stop_gained | — | 0.000000 | n/a | 0.198 | **Likely Pathogenic** |

## Other candidates

_Incidental P/LP not on the ACMG SF list, plus phenotype-unrelated uncertain/benign candidates. Not routinely reported._
| Gene | Transcript | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | local cohort AF | HPO | ACMG |
|---|---|---|---|---|---|---|---|---|---|
| PAX6 | — | c.202C>T p.Arg68Ter | het | stop_gained | — | 0.000000 | n/a | 0.199 | **Likely Pathogenic** |
| OBSCN | — | c.298G>A p.Val100Ile | het | missense_variant | — | 0.000000 | n/a | 0.187 | **Uncertain Significance (VUS)** |


## Per-variant ACMG rationale (auditable)

### SCN1A — c.1834C>T (p.Arg612Ter) → Pathogenic

**Rule path:** `PVS1 + PM2 + PP4 + PP5 => Pathogenic [PATH-1 (PVS1 + strong/moderate/supporting)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True, lof_mechanism_basis=gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease), gene_moi=AD, exon=None, pvs1_strength=very_strong | ClinGen Dosage Sensitivity (HI=3, local); gnomAD v2.1.1 constraint (local) | engine | stop_gained is loss-of-function and LoF is a disease mechanism in SCN1A (gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease)) |
| **PS1** | — | strong | hgvs_p=p.Arg612Ter, ps1_match=None, own_clinvar_plp=True | — | engine | no distinct ClinVar pathogenic variant with the same amino-acid change |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg612Ter, hotspot_residues=0, hotspot_changes=0, enrichment=0.0, gene_baseline=None, window=None, cutoff=3, enrichment_cutoff=2.0 | — | engine | stop_gained is not a missense variant |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.0001, moi=AD, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.0001 (SCN1A is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg612Ter, pm5_match=None, own_clinvar_plp=True, pm5_strength=None | — | engine | variant's own ClinVar assertion is pathogenic — captured by PP5; PM5 withheld to avoid double-counting the same ClinVar evidence |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=SCN1A, mis_z=5.221, mis_z_cutoff=3.09, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | ✅ met | supporting | hpo_match_score=1.0, matched_terms=['HP:0001250→HP:0001250 (Seizure, 1.00)', 'HP:0002133→HP:0002133 (Status epilepticus, 1.00)', 'HP:0011097→HP:0011097 (Epileptic spasm, 1.00)'], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 1.00 (terms: HP:0001250→HP:0001250 (Seizure, 1.00), HP:0002133→HP:0002133 (Status epilepticus, 1.00), HP:0011097→HP:0011097 (Epileptic spasm, 1.00)) |
| **PP5** | ✅ met | supporting | clinvar=Pathogenic, review_status=criteria_provided,_multiple_submitters,_no_conflicts | VCV000012345 | engine | ClinVar Pathogenic (criteria_provided,_multiple_submitters,_no_conflicts) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (SCN1A is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=stop_gained, gene=SCN1A, gene_lof_intolerant=True, oe_mis_upper=0.59, oe_mis_cutoff=1.0, missense_tolerant=False, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=stop_gained, hgvs_p=p.Arg612Ter | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | in-silico (none) | engine | in-silico predictors not benign / unavailable |
| **BP5** | — | supporting | gene=SCN1A | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=Pathogenic, review_status=criteria_provided,_multiple_submitters,_no_conflicts | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### KCNQ2 — c.637C>T (p.Arg213Trp) → Likely Pathogenic

**Rule path:** `PM1 + PM2 + PP3 + PP4 + PP5 => Likely Pathogenic [LP-6 (1 Moderate + >=4 Supporting)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | — | very_strong | consequence=missense_variant, gene_lof_intolerant=True, lof_mechanism_basis=gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease), gene_moi=AD, exon=None, pvs1_strength=None | — | engine | missense_variant is not a qualifying null variant |
| **PS1** | — | strong | hgvs_p=p.Arg213Trp, ps1_match={'alt_aa': 'Trp', 'ref_aa': 'Arg', 'stars': 2, 'genomic_key': '20-63444712-G-A', 'accession': 'VCV000021795'}, own_clinvar_plp=True | — | engine | variant's own ClinVar assertion is pathogenic — captured by PP5; PS1 withheld to avoid double-counting the same ClinVar evidence |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | ✅ met | moderate | consequence=missense_variant, hgvs_p=p.Arg213Trp, hotspot_residues=14, hotspot_changes=29, enrichment=3.22, gene_baseline=0.2895, window=7, cutoff=3, enrichment_cutoff=2.0 | ClinVar residue index (local) | engine | 14 distinct pathogenic-missense residues within ±7 aa (29 pathogenic changes), 3.22× denser than KCNQ2's own baseline — mutational hotspot by ClinVar density |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.0001, moi=AD, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.0001 (KCNQ2 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=missense_variant | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg213Trp, pm5_match=None, own_clinvar_plp=True, pm5_strength=None | — | engine | variant's own ClinVar assertion is pathogenic — captured by PP5; PM5 withheld to avoid double-counting the same ClinVar evidence |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=missense_variant, gene=KCNQ2, mis_z=4.041, mis_z_cutoff=3.09, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | regional hotspot evidence applies (PM1, Moderate) — PP2 stands down so the same missense-intolerance signal is not counted at two granularities |
| **PP3** | ✅ met | supporting | revel=0.92, cadd_phred=32.0, revel_cutoff=0.7, cadd_cutoff=20.0 | REVEL/CADD (local) | engine | REVEL=0.92, CADD=32.0 above deleterious cutoffs |
| **PP4** | ✅ met | supporting | hpo_match_score=1.0, matched_terms=['HP:0001250→HP:0001250 (Seizure, 1.00)', 'HP:0002133→HP:0002133 (Status epilepticus, 1.00)', 'HP:0011097→HP:0011097 (Epileptic spasm, 1.00)'], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 1.00 (terms: HP:0001250→HP:0001250 (Seizure, 1.00), HP:0002133→HP:0002133 (Status epilepticus, 1.00), HP:0011097→HP:0011097 (Epileptic spasm, 1.00)) |
| **PP5** | ✅ met | supporting | clinvar=Pathogenic, review_status=criteria_provided,_single_submitter | VCV000067890 | engine | ClinVar Pathogenic (criteria_provided,_single_submitter) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (KCNQ2 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=missense_variant, gene=KCNQ2, gene_lof_intolerant=True, oe_mis_upper=0.57, oe_mis_cutoff=1.0, missense_tolerant=False, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | KCNQ2 does not tolerate missense (not depleted-free) — BP1 not supported |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=missense_variant, hgvs_p=p.Arg213Trp | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | — | supporting | revel=0.92, cadd_phred=32.0, revel_cutoff=0.15, cadd_cutoff=10.0 | REVEL/CADD (local) | engine | in-silico predictors not benign / unavailable |
| **BP5** | — | supporting | gene=KCNQ2 | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=Pathogenic, review_status=criteria_provided,_single_submitter | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=missense_variant | — | engine | not a synonymous variant |

### PAX6 — c.202C>T (p.Arg68Ter) → Likely Pathogenic

**Rule path:** `PVS1 + PM2 => Likely Pathogenic [LP-1 (PVS1 + Supporting)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True, lof_mechanism_basis=gene is LoF-intolerant (population constraint), gene_moi=AD, exon=None, pvs1_strength=very_strong | gnomAD v2.1.1 constraint (local) | engine | stop_gained is loss-of-function and LoF is a disease mechanism in PAX6 (gene is LoF-intolerant (population constraint)) |
| **PS1** | — | strong | hgvs_p=p.Arg68Ter, ps1_match=None, own_clinvar_plp=False | — | engine | ClinVar residue index unavailable — PS1 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg68Ter, hotspot_residues=0, hotspot_changes=0, enrichment=0.0, gene_baseline=None, window=None, cutoff=3, enrichment_cutoff=2.0 | — | engine | stop_gained is not a missense variant |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.0001, moi=AD, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.0001 (PAX6 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg68Ter, pm5_match=None, own_clinvar_plp=False, pm5_strength=None | — | engine | ClinVar residue index unavailable — PM5 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=PAX6, mis_z=2.817, mis_z_cutoff=3.09, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | — | supporting | hpo_match_score=0.199, matched_terms=['HP:0001250→HP:0100022 (Abnormality of movement, 0.31)'], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 0.20 below 0.6 |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (PAX6 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=stop_gained, gene=PAX6, gene_lof_intolerant=True, oe_mis_upper=0.573, oe_mis_cutoff=1.0, missense_tolerant=False, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=stop_gained, hgvs_p=p.Arg68Ter | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | in-silico (none) | engine | in-silico predictors not benign / unavailable |
| **BP5** | — | supporting | gene=PAX6 | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### RB1 — c.958C>T (p.Arg320Ter) → Likely Pathogenic

**Rule path:** `PVS1 + PM2 => Likely Pathogenic [LP-1 (PVS1 + Supporting)]`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | ✅ met | very_strong | consequence=stop_gained, gene_lof_intolerant=True, lof_mechanism_basis=gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease), gene_moi=AD, exon=None, pvs1_strength=very_strong | ClinGen Dosage Sensitivity (HI=3, local); gnomAD v2.1.1 constraint (local) | engine | stop_gained is loss-of-function and LoF is a disease mechanism in RB1 (gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease)) |
| **PS1** | — | strong | hgvs_p=p.Arg320Ter, ps1_match=None, own_clinvar_plp=False | — | engine | no distinct ClinVar pathogenic variant with the same amino-acid change |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=stop_gained, hgvs_p=p.Arg320Ter, hotspot_residues=0, hotspot_changes=0, enrichment=0.0, gene_baseline=None, window=None, cutoff=3, enrichment_cutoff=2.0 | — | engine | stop_gained is not a missense variant |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.0001, moi=AD, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.0001 (RB1 is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=stop_gained | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Arg320Ter, pm5_match=None, own_clinvar_plp=False, pm5_strength=None | — | engine | no other ClinVar pathogenic missense at this residue |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=stop_gained, gene=RB1, mis_z=2.673, mis_z_cutoff=3.09, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **PP3** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.7, cadd_cutoff=20.0 | in-silico (none) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | — | supporting | hpo_match_score=0.198, matched_terms=['HP:0001250→HP:0001249 (Intellectual disability, 0.31)'], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 0.20 below 0.6 |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (RB1 is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=stop_gained, gene=RB1, gene_lof_intolerant=True, oe_mis_upper=0.719, oe_mis_cutoff=1.0, missense_tolerant=False, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | stop_gained is not a missense variant |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=stop_gained, hgvs_p=p.Arg320Ter | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | — | supporting | revel=None, cadd_phred=None, revel_cutoff=0.15, cadd_cutoff=10.0 | in-silico (none) | engine | in-silico predictors not benign / unavailable |
| **BP5** | — | supporting | gene=RB1 | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=stop_gained | — | engine | not a synonymous variant |

### CACNA1A — c.100T>C (p.Ser34Pro) → Uncertain Significance (VUS)

**Rule path:** `PM2 + PP2 + PP4 => criteria insufficient for a benign or pathogenic call => VUS`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | — | very_strong | consequence=missense_variant, gene_lof_intolerant=True, lof_mechanism_basis=gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease), gene_moi=AD, exon=None, pvs1_strength=None | — | engine | missense_variant is not a qualifying null variant |
| **PS1** | — | strong | hgvs_p=p.Ser34Pro, ps1_match=None, own_clinvar_plp=False | — | engine | ClinVar residue index unavailable — PS1 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=missense_variant, hgvs_p=p.Ser34Pro, hotspot_residues=0, hotspot_changes=0, enrichment=0.0, gene_baseline=0.0, window=7, cutoff=3, enrichment_cutoff=2.0 | — | engine | ClinVar residue index unavailable — PM1 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.0001, moi=AD, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.0001 (CACNA1A is AD) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=missense_variant | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Ser34Pro, pm5_match=None, own_clinvar_plp=False, pm5_strength=None | — | engine | ClinVar residue index unavailable — PM5 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | ✅ met | supporting | consequence=missense_variant, gene=CACNA1A, mis_z=5.785, mis_z_cutoff=3.09, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | missense in CACNA1A, a missense-constrained gene (gnomAD mis_z=5.79 ≥ 3.09) |
| **PP3** | — | supporting | revel=0.45, cadd_phred=18.0, revel_cutoff=0.7, cadd_cutoff=20.0 | REVEL/CADD (local) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | ✅ met | supporting | hpo_match_score=0.93, matched_terms=['HP:0001250→HP:0001250 (Seizure, 1.00)', 'HP:0002133→HP:0002133 (Status epilepticus, 1.00)', 'HP:0011097→HP:0032794 (Myoclonic seizure, 0.79)'], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 0.93 (terms: HP:0001250→HP:0001250 (Seizure, 1.00), HP:0002133→HP:0002133 (Status epilepticus, 1.00), HP:0011097→HP:0032794 (Myoclonic seizure, 0.79)) |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.001, moi=AD, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.001 BS1 cutoff (CACNA1A is AD) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=missense_variant, gene=CACNA1A, gene_lof_intolerant=True, oe_mis_upper=0.609, oe_mis_cutoff=1.0, missense_tolerant=False, missense_constrained=True | gnomAD v2.1.1 constraint (local) | engine | CACNA1A does not tolerate missense (not depleted-free) — BP1 not supported |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=missense_variant, hgvs_p=p.Ser34Pro | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | — | supporting | revel=0.45, cadd_phred=18.0, revel_cutoff=0.15, cadd_cutoff=10.0 | REVEL/CADD (local) | engine | in-silico predictors not benign / unavailable |
| **BP5** | — | supporting | gene=CACNA1A | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=missense_variant | — | engine | not a synonymous variant |

### OBSCN — c.298G>A (p.Val100Ile) → Uncertain Significance (VUS)

**Rule path:** `PM2 + BP4 => criteria insufficient for a benign or pathogenic call => VUS`

| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |
|---|---|---|---|---|---|---|
| **PVS1** | — | very_strong | consequence=missense_variant, gene_lof_intolerant=False, lof_mechanism_basis=gene has an established autosomal-recessive phenotype (HPO), gene_moi=AR, exon=None, pvs1_strength=None | — | engine | missense_variant is not a qualifying null variant |
| **PS1** | — | strong | hgvs_p=p.Val100Ile, ps1_match=None, own_clinvar_plp=False | — | engine | ClinVar residue index unavailable — PS1 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PS2** | N/A | strong | — | — | engine | Requires parental (trio) data — not available from a single proband VCF |
| **PS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **PS4** | — | strong | gnomad_af=0.0, local_cohort_af=None | — | model | Needs case-control data; population absence alone is captured by PM2 |
| **PM1** | — | moderate | consequence=missense_variant, hgvs_p=p.Val100Ile, hotspot_residues=0, hotspot_changes=0, enrichment=0.0, gene_baseline=0.0, window=7, cutoff=3, enrichment_cutoff=2.0 | — | engine | ClinVar residue index unavailable — PM1 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PM2** | ✅ met | supporting | gnomad_af=0.0, local_cohort_af=None, local_cohort_checked=False, ceiling=0.001, moi=AR, strength_model=richards | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD popmax AF=0.000000, local cohort not checked — gnomAD at/under 0.001 (OBSCN is AR) |
| **PM3** | N/A | moderate | — | — | engine | Requires phasing / a second variant — not determinable from this VCF alone |
| **PM4** | — | moderate | consequence=missense_variant | — | engine | no protein-length-changing consequence |
| **PM5** | — | moderate | hgvs_p=p.Val100Ile, pm5_match=None, own_clinvar_plp=False, pm5_strength=None | — | engine | ClinVar residue index unavailable — PM5 not assessed (build it: scripts/fetch_clinvar_residue.py) |
| **PM6** | N/A | moderate | — | — | engine | Requires parental data — not available from a single proband VCF |
| **PP1** | N/A | supporting | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **PP2** | — | supporting | consequence=missense_variant, gene=OBSCN, mis_z=-1.145, mis_z_cutoff=3.09, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | OBSCN missense z-score -1.15 is below the 3.09 missense-constraint threshold |
| **PP3** | — | supporting | revel=0.08, cadd_phred=4.5, revel_cutoff=0.7, cadd_cutoff=20.0 | REVEL/CADD (local) | engine | in-silico predictors below deleterious cutoffs / unavailable |
| **PP4** | — | supporting | hpo_match_score=0.187, matched_terms=[], cutoff=0.6 | HPO ontology-aware (Lin/IC, local) | engine | phenotype match 0.19 below 0.6 |
| **PP5** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar pathogenic assertion (or 0-star) |
| **BA1** | — | stand_alone | af=0.0, cutoff=0.05, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 at or below 0.05 |
| **BS1** | — | strong | af=0.0, cutoff=0.01, moi=AR, basis=gnomAD/local cohort popmax AF (no faf95 available) | gnomAD gnomad_r4 (local snapshot); local cohort (not in local table) | engine | gnomAD/local cohort popmax AF (no faf95 available) = 0.0000 under the 0.01 BS1 cutoff (OBSCN is AR) |
| **BS2** | — | strong | gnomad_homozygotes=0, cutoff=2 | gnomAD gnomad_r4 (local snapshot) | engine | 0 homozygotes (below 2) |
| **BS3** | — | strong | — | — | model | Requires literature review of functional assays — left for expert/model adjudication |
| **BS4** | N/A | strong | — | — | engine | Requires genotypes for affected relatives — not available from a single proband VCF |
| **BP1** | — | supporting | consequence=missense_variant, gene=OBSCN, gene_lof_intolerant=False, oe_mis_upper=1.068, oe_mis_cutoff=1.0, missense_tolerant=True, missense_constrained=False | gnomAD v2.1.1 constraint (local) | engine | OBSCN is not LoF-intolerant — truncating-only mechanism not established |
| **BP2** | N/A | supporting | — | — | engine | Requires phasing / parental data to establish trans or cis — not available from a single proband VCF |
| **BP3** | — | supporting | consequence=missense_variant, hgvs_p=p.Val100Ile | — | model | Requires a repeat/domain annotation the engine does not carry — model adjudication |
| **BP4** | ✅ met | supporting | revel=0.08, cadd_phred=4.5, revel_cutoff=0.15, cadd_cutoff=10.0 | REVEL/CADD (local) | engine | REVEL=0.08, CADD=4.5 below benign cutoffs |
| **BP5** | — | supporting | gene=OBSCN | — | model | Requires whole-case context (another finding explaining the phenotype) plus clinical judgement — left for expert/model adjudication |
| **BP6** | — | supporting | clinvar=None, review_status=None | — | engine | no reviewed ClinVar benign assertion (or 0-star) |
| **BP7** | — | supporting | consequence=missense_variant | — | engine | not a synonymous variant |

## Methods
- **genome_build:** GRCh38
- **qc_thresholds:** {'min_DP': 10, 'min_GQ': 20}
- **rarity_cutoff_popmax_af:** 0.005
- **ba1_cutoff:** 0.05
- **databases:** ['ClinVar', 'gnomAD r4', 'local cohort', 'HPO', 'gnomAD constraint']
- **standards:** ['ACMG/AMP variant classification (Richards et al., Genet Med 2015)', 'ClinGen SVI criteria refinements', 'ACMG secondary-findings list (SF v3.2, Miller et al. 2023)', 'HGVS nomenclature', 'GA4GH Phenopackets (phenotype exchange)']

## Performance (this run)
- **parse:** 0.0004 s
- **qc:** 0.0 s
- **gnomad prime:** 0.0 s
- **clinvar prime:** 0.0971 s
- **annotate:** 0.3893 s
- **filter:** 0.0001 s
- **alphamissense:** 0.0006 s
- **clinvar residue:** 0.007 s
- **classify:** 0.1943 s
- **total:** 0.6888 s
- **variants per:** 16.0

## Limitations & disclaimers

- Single-proband analysis: criteria requiring parental/segregation/phasing data (PS2, PM3, PM6, PP1, BP2, BS4) are reported as N/A.
- Judgment criteria (PS3, PS4, BS3, BP3, BP5) are surfaced for expert/model adjudication and default to not-met unless explicitly supported.
- Population and clinical databases are versioned snapshots; re-check before sign-out.
- **This is a draft-generation aid, not a diagnostic device.**