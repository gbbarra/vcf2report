# Report style guide

Adjust the generated report to match the lab's sign-out format. The Jinja2
template at `templates/report.md.j2` is the source of truth — edit it to change
section order, headings, or wording for the whole pipeline.

## Conventions
- Language: match the lab's report language (pt-BR or en). The template ships in
  English; translate section headings as needed.
- Variant nomenclature: HGVS with the transcript, e.g. `NM_006920.6:c.1834C>T
  (p.Arg612Ter)`. Always pair c. and p.
- Tier wording: Pathogenic / Likely Pathogenic / Uncertain Significance (VUS) /
  Likely Benign / Benign.
- Always keep: the DRAFT banner, the per-variant criterion table, the methods and
  limitations sections.
- Never remove the disclaimer that this is a draft-generation aid, not a
  diagnostic device.
