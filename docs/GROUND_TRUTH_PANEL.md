# Ground-truth concordance panel (build locally)

Measures whether vcf2report's ACMG engine agrees with **expert-reviewed ClinVar**
on real variants — using only *orthogonal* evidence (frequency, consequence/LoF,
constraint), with **ClinVar withheld from the engine** (PP5 off) so the check is
honest, not circular. The test (`tests/test_ground_truth_concordance.py`) runs
fully offline against a frozen panel; this doc builds that panel.

Everything here is **local** — no live API. The one flaky dependency (remote
gnomAD) is replaced by your local gnomAD sites VCF (the same file vcfanno uses).

## 1. Distill a compact ClinVar TSV

Stream the public ClinVar VCF once and keep only expert-reviewed SNVs with a
definite classification (`allele_id chrom pos ref alt gene clnsig clnrevstat mc`):

```bash
# NCBI weekly VCF (or the Broad Funcotator mirror used in make_synthetic_exomes.sh)
CLINVAR=clinvar_GRCh38.vcf            # uncompressed, or `zcat clinvar.vcf.gz`
awk -F'\t' '
/^#/ {next}
{ info=$8
  if (info !~ /CLNVC=single_nucleotide_variant/ || info !~ /MC=/) next
  match(info,/CLNSIG=[^;]*/); s=substr(info,RSTART+7,RLENGTH-7)
  if (s !~ /^(Pathogenic|Likely_pathogenic|Benign|Likely_benign|Pathogenic\/Likely_pathogenic|Benign\/Likely_benign)$/) next
  match(info,/CLNREVSTAT=[^;]*/); rev=substr(info,RSTART+11,RLENGTH-11)
  if (rev !~ /criteria_provided|reviewed_by_expert_panel|practice_guideline/) next
  match(info,/GENEINFO=[^;]*/); gi=substr(info,RSTART+9,RLENGTH-9); split(gi,g,":")
  match(info,/MC=[^;]*/); mc=substr(info,RSTART+3,RLENGTH-3)
  match(info,/ALLELEID=[^;]*/); aid=substr(info,RSTART+9,RLENGTH-9)
  print aid"\t"$1"\t"$2"\t"$4"\t"$5"\t"g[1]"\t"s"\t"rev"\t"mc }' "$CLINVAR" > clinvar_compact.tsv
```

## 2. Freeze the panel with LOCAL gnomAD (no live API)

Point `--gnomad-vcf` at your local gnomAD v4.1 sites VCF (bgzipped + tabixed —
the same one referenced by `scripts/vcfanno.conf.toml`). Real AF + `faf95` are
read locally, per variant, no network:

```bash
python scripts/build_truth_panel.py clinvar_compact.tsv \
    --gnomad-vcf gnomad.genomes.v4.1.sites.GRCh38.vcf.gz
# -> data/truth/clinvar_panel.json  (~100 variants: Pathogenic LoF, Pathogenic
#    missense, Benign — balanced, deterministic, capped per gene)
```

> Omitting `--gnomad-vcf` falls back to the live remote path
> (`VCF2REPORT_ALLOW_NETWORK=1`). Prefer the local file — the remote API is slow
> and flaky in locked-down networks.

## 3. Run the concordance check

```bash
python -m pytest tests/test_ground_truth_concordance.py -v -s
```

It asserts the safety-critical properties and prints the confusion matrix:

* **zero directional contradictions** — no ClinVar-Pathogenic called Benign/LB,
  no ClinVar-Benign called Pathogenic/LP (VUS is a safe non-contradiction).
* **BA1-band benign recall** — benign at AF ≥ 5% must land Benign/LB.

Without ClinVar as evidence the engine is deliberately conservative: pathogenic
missense with no in-silico/functional support lands **VUS**, not Pathogenic —
that's ACMG behaving correctly, not a miss. The matrix makes this explicit.

The frozen `data/truth/clinvar_panel.json` is small and safe to commit (public,
de-identified). Rebuild it when you refresh ClinVar or gnomAD.
