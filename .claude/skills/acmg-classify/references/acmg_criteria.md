# ACMG/AMP criteria & combining rules (Richards et al., Genet Med 2015)

Reference for grounding classification. Strengths: Very Strong (PVS), Strong
(PS/BS), Moderate (PM), Supporting (PP/BP), Stand-alone (BA).

## Pathogenic criteria
- **PVS1** (very strong) — null variant (nonsense, frameshift, canonical ±1/2
  splice, initiation codon, single/multi-exon deletion) in a gene where loss of
  function is a known disease mechanism.
- **PS1** (strong) — same amino-acid change as an established pathogenic variant.
- **PS2** (strong) — de novo (maternity & paternity confirmed) in a patient with
  no family history. *Needs a trio.*
- **PS3** (strong) — well-established functional studies show a damaging effect.
- **PS4** (strong) — prevalence in affecteds significantly increased vs controls.
- **PM1** (moderate) — mutational hotspot / critical, well-established functional
  domain without benign variation.
- **PM2** (moderate) — absent (or at extremely low frequency) from population
  databases. *vcf2report requires absence in BOTH gnomAD and ABraOM.*
- **PM3** (moderate) — for recessive disorders, detected in trans with a
  pathogenic variant. *Needs phasing.*
- **PM4** (moderate) — protein length change (in-frame indel / stop-loss).
- **PM5** (moderate) — novel missense at a residue where a different pathogenic
  missense has been seen.
- **PM6** (moderate) — assumed de novo without confirmation of parentage.
- **PP1** (supporting) — cosegregation with disease. *Needs a family.*
- **PP2** (supporting) — missense in a gene with low benign-missense rate where
  missense is a common mechanism.
- **PP3** (supporting) — multiple computational lines support a deleterious effect.
- **PP4** (supporting) — patient phenotype/family history highly specific for a
  disease with a single genetic aetiology.
- **PP5** (supporting) — reputable source reports pathogenic (deprecated).

## Benign criteria
- **BA1** (stand-alone) — allele frequency > 5% in a population database.
- **BS1** (strong) — allele frequency greater than expected for the disorder.
- **BS2** (strong) — observed in a healthy adult (homozygous/hemizygous/het as
  appropriate) for a fully-penetrant early-onset disorder.
- **BS3** (strong) — well-established functional studies show no damaging effect.
- **BS4** (strong) — lack of segregation in affected family members. *Needs family.*
- **BP1** (supporting) — missense in a gene where only LoF causes disease.
- **BP2** (supporting) — observed in trans/cis inconsistent with the mechanism.
- **BP3** (supporting) — in-frame indel in a repeat region.
- **BP4** (supporting) — multiple computational lines suggest no impact.
- **BP5** (supporting) — found in a case with an alternate molecular cause.
- **BP6** (supporting) — reputable source reports benign (deprecated).
- **BP7** (supporting) — synonymous with no predicted splice impact.

## Combining rules (Table 5)
**Pathogenic** if:
1. 1 PVS1 + (≥1 PS OR ≥2 PM OR 1 PM+1 PP OR ≥2 PP)
2. ≥2 PS
3. 1 PS + (≥3 PM OR 2 PM+≥2 PP OR 1 PM+≥4 PP)

**Likely Pathogenic** if:
1. 1 PVS1 + 1 PM
2. 1 PS + 1–2 PM
3. 1 PS + ≥2 PP
4. ≥3 PM
5. 2 PM + ≥2 PP
6. 1 PM + ≥4 PP

**Benign** if: BA1 (stand-alone) OR ≥2 BS.

**Likely Benign** if: 1 BS + 1 BP OR ≥2 BP.

**Uncertain Significance (VUS)** if criteria are insufficient OR pathogenic and
benign evidence contradict.
