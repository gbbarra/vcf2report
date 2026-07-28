# The 28 ACMG/AMP criteria — what this engine decides, and from what

Every criterion in Richards et al. (Genet Med 2015) is evaluated and shown on every variant, in the
order below. Nothing is silently omitted: a criterion the engine cannot decide says so, with the
reason, rather than disappearing from the trail.

Each falls into exactly one class:

| class | meaning |
|---|---|
| **engine** | decided deterministically from local data (thresholds, lookups). Reproducible and unit-tested. |
| **model** | genuine clinical/literature judgement. The evidence is surfaced and the criterion is tagged `adjudicated_by="model"`, defaulting to not-met, so a human or model decides it transparently instead of the engine pretending it is a fact. |
| **N/A** | undecidable from a single-proband VCF (trio, segregation, phasing). Reported with `applies=False` and an explicit reason. |

**17 engine · 5 model · 6 N/A.**

> This table is generated from the criterion registry and asserted against it by
> `tests/test_acmg_coverage.py`, so it cannot drift from the code.

## Pathogenic

| code | strength | decided by | criterion |
|---|---|---|---|
| `PVS1` | very strong | engine | Null variant in a gene where LoF is a known disease mechanism |
| `PS1` | strong | engine | Same amino-acid change as a DIFFERENT established pathogenic variant |
| `PS2` | strong | N/A | De novo (confirmed) in a patient |
| `PS3` | strong | model | Well-established functional studies show a damaging effect |
| `PS4` | strong | model | Prevalence in affected significantly increased vs controls |
| `PM1` | moderate | engine | Located in a mutational hotspot / critical functional domain |
| `PM2` | moderate | engine | Absent or ultra-rare in population databases (gnomAD + local cohort) |
| `PM3` | moderate | N/A | Detected in trans with a pathogenic variant (recessive) |
| `PM4` | moderate | engine | Protein length change (in-frame indel / stop-loss) in a non-repeat region |
| `PM5` | moderate | engine | Novel missense at a residue where a different pathogenic missense is known |
| `PM6` | moderate | N/A | Assumed de novo (parentage not confirmed) |
| `PP1` | supporting | N/A | Co-segregation with disease in multiple affected family members |
| `PP2` | supporting | engine | Missense in a gene with a low rate of benign missense, where missense is a mechanism |
| `PP3` | supporting | engine | Multiple in-silico lines of evidence support a deleterious effect |
| `PP4` | supporting | engine | Patient phenotype highly specific for the gene (HPO match) |
| `PP5` | supporting | engine | Reputable source (ClinVar) classifies the variant as pathogenic |

## Benign

| code | strength | decided by | criterion |
|---|---|---|---|
| `BA1` | stand alone | engine | Allele frequency > 5% in a population database |
| `BS1` | strong | engine | Allele frequency greater than expected for the disorder |
| `BS2` | strong | engine | Observed in healthy adult homozygotes |
| `BS3` | strong | model | Well-established functional studies show NO damaging effect |
| `BS4` | strong | N/A | Lack of segregation in affected members of a family |
| `BP1` | supporting | engine | Missense in a gene where primarily truncating variants cause disease |
| `BP2` | supporting | N/A | Observed in trans with a pathogenic variant (dominant gene), or in cis with one |
| `BP3` | supporting | model | In-frame indel in a repetitive region without a known function |
| `BP4` | supporting | engine | Multiple in-silico lines of evidence suggest no impact |
| `BP5` | supporting | model | Variant found in a case with an alternate molecular basis for disease |
| `BP6` | supporting | engine | Reputable source (ClinVar) classifies the variant as benign |
| `BP7` | supporting | engine | Synonymous variant with no predicted splice impact |

## What drives each engine criterion

| data source | criteria |
|---|---|
| gnomAD allele frequency (+ local cohort) | `PM2`, `BA1`, `BS1`, `BS2` |
| gnomAD per-gene LoF constraint (pLI / LOEUF) + ClinGen dosage + HPO inheritance | `PVS1` |
| gnomAD per-gene **missense** constraint (`mis_z`, `oe_mis_upper`) | `PP2`, `BP1` |
| ClinVar — the variant's **own** assertion | `PP5`, `BP6` |
| ClinVar **residue index** (same residue / neighbourhood) | `PS1`, `PM5`, `PM1` |
| AlphaMissense (ClinGen-calibrated) / REVEL / CADD | `PP3`, `BP4`, `BP7` |
| HPO phenotype similarity | `PP4` |
| Consequence + exon rank from the annotator | `PVS1` strength, `PM4`, `BP7` |

## Each line of evidence is counted once

ACMG's combining rules assume independent evidence. Three double-counts have been found and closed;
`tests/test_criteria_invariants.py` asserts the relationships across value ranges so they cannot
reappear:

- **`PP5` ⊥ `PS1`/`PM5`** — a variant ClinVar already calls P/LP is covered by PP5, so the residue
  criteria stand down rather than counting the same ClinVar record twice.
- **`PS1` › `PM5` › `PM1`** — same residue + same change / same residue + different change /
  neighbourhood only. Exactly one fires, so ACMG's rule against combining PM1 with PM5 holds by
  construction.
- **`PP2` ⊥ `PM1`** — the same missense-intolerance claim at gene-wide vs. regional granularity; the
  more specific PM1 wins.
- **`PP2` ⊥ `BP1`**, **`PP3` ⊥ `BP4`**, **`PP5` ⊥ `BP6`**, **`PVS1` ⊥ `PM4`**, and at most one of
  **`PM2`/`BS1`/`BA1`**.

An annotation with nothing looked up fires **no criterion at all** — absence of data is never read as
evidence of absence.

## Known approximations

Stated here rather than buried, because each is a place where the engine is doing something the
guideline describes more loosely:

- **`PM1`'s "without benign variation"** has no residue-level benign index; it is approximated at
  gene level (genes that tolerate missense are excluded).
- **`PM1`'s hotspot thresholds** (±7 aa, ≥3 pathogenic residues, ≥2× the gene's own baseline density)
  are an empirical calibration against over-call, not a published cut. Raw density alone fired on 61%
  of FBN1 and 49% of SCN1A novel residues, because ClinVar density tracks how thoroughly a gene has
  been *studied*.
- **`BP1`'s "primarily truncating"** has no curated gene list; it is proxied by LoF-intolerant *and*
  missense-tolerant constraint.
- **`PVS1`'s mechanism gate** accepts ClinGen HI=3, population constraint, *or* an established
  recessive phenotype — constraint alone is structurally blind to recessive disease.
- **`PP5`/`BP6`** were deprecated by the ClinGen SVI; they are retained as transparent, ≥1★-gated
  Supporting lines rather than removed.

## Combining

Two models, toggled by `VCF2REPORT_ACMG_MODEL`: **Richards 2015 Table 5** (the conservative default)
and the **ClinGen/Tavtigian 2020 points** system. Both emit the rule path that produced the tier, so
the classification is auditable end to end.
