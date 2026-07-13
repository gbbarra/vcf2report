# Benchmark — does vcf2report flag real pathogenic variants?

A validation harness that answers the only question that matters for a variant-interpretation
tool: **when a real pathogenic variant is present, does the tool bring it to the clinician's
attention?** It uses real, published cases — not hand-invented variants.

## Data (all openly licensed)

| piece | source | license |
|---|---|---|
| variants + phenotypes | **GA4GH Phenopacket Store** (10,377 curated cases from 959 publications; GRCh38 coordinates, HGVS, HPO terms, ACMG P/LP) | BSD-3 |
| healthy background exome | **GIAB NA12878 / HG001** (NIST) | public domain |
| ClinVar status | ClinVar GRCh38 VCF (NCBI) | public domain |
| phenotype ontology | HPO `genes_to_phenotype` | open |

Each synthetic case = the NA12878 exome + **one** real pathogenic variant spiked in (with a
constructed SnpEff `ANN`) + that case's HPO terms. Half the cases are fed a ClinVar
Pathogenic assertion ("with ClinVar"), half are not ("without" → the engine must classify
from mechanism: PVS1/PM2/PP3/PP4). All variants are absent from the local gnomAD store (they
are rare disease alleles) — verified: on the 12-case set, 11/12 are absent from **all** of
gnomAD and 0/12 are genomes-only.

## The metric: *surfaced to the clinician*, not the raw ACMG tier

vcf2report classifies **independently** and conservatively: it does not echo ClinVar into the
ACMG tier (ClinGen SVI deprecated PP5/BP6 for circularity), and with an exomes+MANE gnomAD
store it will not assert a false absence, so **PM2 does not fire for a variant merely absent
from the store**. So the raw tier is deliberately strict. The clinically-meaningful metric is
whether the variant is **brought to attention** by *any* of:

- the engine's ACMG tier is Pathogenic / Likely Pathogenic, **or**
- it is flagged as a ≥2-star **ClinVar Pathogenic** (surfaced in the conclusion), **or**
- it is a **phenotype-matched** primary candidate (HPO overlap ≥ threshold, not benign).

## Adversarial review (23 agents, 4 lenses)

An adversarial review of the 12-case run found **0 code bugs** and confirmed the 0/12
"auto-Pathogenic" rate is **correct-conservative behavior**, not a defect, plus one
clinical-safety gap and a set of tunable defaults:

- **FIXED — ClinVar surface (clinical-safety).** A ≥2-star ClinVar Pathogenic variant the
  engine tiered VUS was being reported as "no Pathogenic finding". The conclusion now flags
  it explicitly ("⚠️ Classified Pathogenic in ClinVar … DO NOT dismiss") without touching the
  ACMG math (no circularity). This alone lifted *surfaced* from **0/12 → 9/12**.
- **Recommended (lab's call, not auto-applied):**
  - *Region-aware PM2 / a joint store.* The exomes+MANE store is `mode=partial`, so PM2 can
    never fire on genuine absence — the tool's primary use case. Building the **joint**
    (exomes+genomes) preset over the MANE panel and asserting absence per covered BED
    interval would let PM2 fire soundly (empirically 0 false positives on the benchmark).
    Flipping the exomes-only store to `mode=full` is **not** safe (off-panel false absence).
  - *ClinGen points model as default.* `VCF2REPORT_ACMG_MODEL=clingen` recovers ~4/12 (PVS1
    + phenotype/ClinVar reaches LP), but makes bare PVS1-alone → LP — a real over-call
    tradeoff, so it should ship with an SVI-style "needs corroboration" guard. Kept as an
    opt-in for now; Richards remains the conservative default.

## Results

**12-case set (mixed with/without ClinVar):** found 12/12 · engine P/LP 0/12 · **surfaced
9/12**. The 3 residual: 2 LoF with no phenotype match (would need region-aware PM2) + 1
non-coding-RNA gene (out of the protein-coding ACMG scope).

**284-case set:** _(running — numbers filled on completion)_

## Honest limitations

- Single-proband: PS2/PM3/PM6/PP1/BS4 (de novo, in-trans, segregation) are N/A.
- The exomes+MANE store keeps PM2 conservative (see above) — a deliberate no-false-absence
  choice, quantified here rather than hidden.
- The spike-in uses a constructed `ANN`; a real SnpEff/VEP annotation may add exon rank
  (PVS1 strength modulation) that this harness omits.

## Reproduce

The harness (survey phenopackets → curate a balanced set → spike into NA12878 → run + score)
is deterministic given the three public inputs above; see `scripts/` and this doc.
