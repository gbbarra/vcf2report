# gnomAD frequency store — data provenance & license

The DuckDB/Parquet store used by vcf2report for population allele frequencies is
**derived from gnomAD v4.1**. The store shipped/auto-detected today is built from the
**exomes** release, MANE-sliced. Building from the **joint** (exomes+genomes) release
— `build_gnomad_parquet.py --preset joint` (the script default) — is recommended and
removes the exome-only false-absence limitation documented in `docs/BENCHMARK.md`
(a variant present in genomes but not exomes would otherwise be treated as absent).

## Source

- **Dataset:** Genome Aggregation Database (gnomAD) v4.1 (exomes for the current store;
  joint exomes+genomes for a `--preset joint` rebuild).
- **Publisher:** Broad Institute — <https://gnomad.broadinstitute.org>
- **Retrieved from:** `gs://gcp-public-data--gnomad/release/4.1/vcf/{exomes,joint}/`

## What is included

Only **allele-frequency fields** are extracted (one row per gnomAD variant):
`chrom, pos, ref, alt, filter, af, af_grpmax, ac, an, nhomalt, faf95, grpmax_pop`, and
per-ancestry AFs (`af_afr, af_amr, af_asj, af_eas, af_fin, af_mid, af_nfe, af_sas,
af_remaining`, plus `af_ami` in the joint release). No restrictively-licensed
annotations (e.g. SpliceAI) are included.

## License

gnomAD data is released under the **Open Data Commons Open Database License (ODbL) v1.0**
— <https://opendatacommons.org/licenses/odbl/1-0/>. Redistribution (including this
derived Parquet) is permitted, **provided**:

1. **Attribution** — credit gnomAD v4.1 (this notice).
2. **Share-Alike** — any redistributed derivative database stays under ODbL-1.0.

Citation: Chen, S., Francioli, L.C., Goodrich, J.K. et al. *A genomic mutational
constraint map using variation in 76,156 human genomes.* Nature 625, 92–100 (2024).

## Reproduce it yourself

This store is **not committed** (~1 GB). Rebuild it from the public bucket with
`scripts/build_gnomad_parquet.py`, or download the published, checksummed copy with
`scripts/fetch_gnomad_parquet.sh`. Either lands it at `data/gnomad/gnomad_parquet/`,
where vcf2report auto-detects it.
