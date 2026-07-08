---
name: acmg-classify
description: >-
  Apply auditable ACMG/AMP variant classification. Use when interpreting a genetic
  variant's pathogenicity, evaluating ACMG criteria (PVS1, PM2, PP3, BA1, ...), or
  explaining why a variant is Pathogenic / Likely Pathogenic / VUS / Likely Benign /
  Benign. Pairs with the vcf2report engine, which decides deterministic criteria.
---

# Auditable ACMG/AMP classification

Classify variants transparently — a glass box, not a score. The `vcf2report`
engine (`classify_variant` MCP tool) decides the **deterministic** criteria from
data; you **explain** them and **adjudicate** the judgment criteria.

## Division of labour
- **Engine decides (adjudicated_by=engine):** PVS1 (LoF mechanics + gene
  constraint), PM2 (gnomAD **and** ABraOM absence), PM4, PP3/BP4 (in-silico
  cutoffs, mutually exclusive), PP4 (HPO match), PP5 (reviewed ClinVar P/LP as a
  *supporting* reputable-source line), BA1/BS1/BS2 (frequency), BP7.
- **You adjudicate (adjudicated_by=model):** PS1 (same aa change as a *different*
  established pathogenic variant — not the variant's own ClinVar record), PM1
  (hotspot/domain), PS3 (functional studies), PP2 (gene missense constraint), PM5
  (residue-level cross-check), PS4 (case-control). Only mark these met with an
  explicit, cited reason. Default to not-met when evidence is absent.
- **N/A for a single proband:** PS2, PM3, PM6, PP1, BS4 — never assert these
  without trio / segregation / phasing data.

## Rules
1. Ground every criterion in `references/acmg_criteria.md` (definitions + the
   Richards 2015 combining rules). Do not rely on memory for thresholds.
2. Each criterion needs: state (met / not / N/A), the concrete evidence value,
   the source (DB + accession + date), and a one-line reason.
3. Apply the combining rules exactly to reach the 5-tier call; if pathogenic and
   benign evidence both fire, report **VUS (conflicting)**.
4. If your adjudication of a judgment criterion changes the tier, say so and show
   the before/after.
