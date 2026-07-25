# Disclaimers

**vcf2report is a research and draft-generation aid. It is NOT a medical device
and NOT for clinical diagnosis.**

- Every report is a **draft** that must be reviewed, corrected, and signed out by a
  qualified professional (clinical geneticist / molecular pathologist).
- The ACMG classification automates deterministic criteria and surfaces judgment
  criteria for human/model adjudication. It does not replace expert curation.
- Population and clinical databases are **versioned snapshots**; re-check current
  ClinVar/gnomAD before any clinical decision.
- Single-proband analysis cannot assess criteria that require trio, segregation,
  or phasing data (PS2, PM3, PM6, PP1, BP2, BS4); these are reported as N/A.
- All 28 ACMG/AMP criteria are evaluated and shown on every variant — 17 decided by
  the engine, 5 left for expert/model adjudication, 6 reported N/A. See
  [ACMG_CRITERIA.md](ACMG_CRITERIA.md) for the per-criterion breakdown, the data
  behind each engine decision, and the documented approximations.
- The bundled sample data is **synthetic and de-identified** — not real patient
  data. Gene/variant coordinates in the sample are illustrative.
- No patient-identifying information should be committed to this repository or sent
  to any external service. Handle real VCFs under your institution's data-governance
  and privacy (LGPD/HIPAA) policies.

By using vcf2report you accept that the authors provide no warranty and assume no
liability for clinical use.
