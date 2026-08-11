# Whole-genome (WGS) on a small-RAM machine — the chunked path

Running the standard pipeline (`annotate_vcf.sh` + `run_headless.py`) directly on a
whole-genome VCF (~5M variants) OOMs a low-RAM (e.g. 8 GB) machine, in two places:

1. **SnpEff heap.** `annotate_vcf.sh` used a fixed `java -Xmx8g` — larger than total
   RAM on an 8 GB box. It is now overridable: **`SNPEFF_XMX`** (default `8g`). `SNPEFF_XMX=3g`
   runs comfortably (peak RSS ~1.5 GB measured, even on chr1).
2. **In-memory variant list.** `parse_vcf()` loads every variant into one Python list;
   ~5M `Variant` objects exhaust RAM regardless of how efficient the DuckDB/Parquet
   stores are (those scan columnar and are not the bottleneck).

Both are avoided by processing **one chromosome at a time** and stitching the results
back together. ACMG classification is per-variant with no cross-chromosome context, so
this is exact — only the sequencing-quality panel needs count-weighted recombination.

## One command

```bash
scripts/run_wgs_chunked.sh RAW.wgs.vcf.gz GRCh38.fa patient_hpo.txt OUTDIR
# env: SNPEFF_XMX (default 3g), CHROMS (override the contig list)
```

Per chromosome it runs `bcftools view -r` → `annotate_vcf.sh` (norm + chr↔Ensembl
rename + SnpEff at `SNPEFF_XMX` + vcfanno) → `run_headless.py`, then
`merge_wgs_results.py` stitches the per-chromosome `results.json` into
`OUTDIR/merged.results.json`. Render the laudo from that exactly like any single run.

## What it does, and why each step

- **Scope to the 25 primary contigs** (`chr1..22, X, Y, M`). They hold ~98.6% of a
  typical WGS VCF's variants; the ~195 alt/decoy/unplaced/HLA contigs are not worth the
  per-chunk overhead. Naming (`chr1` vs `1`) is auto-detected from the VCF.
- **Annotate per chromosome** at a reduced heap — peak RSS ~1.5 GB.
- **Classify per chromosome** — peak RSS ~515 MB on the largest chromosome. Safe because
  ACMG is local to each variant.
- **Merge** (`scripts/merge_wgs_results.py`, logic in `vcf2report/report/merge.py`):
  - QC funnel counters **summed**; finding lists (`classifications`, every bucket,
    `clinvar_do_not_dismiss`, `conclusion`, `qc_rescued`) **concatenated**.
  - Sequencing-quality panel: counts summed; means/medians/percentages **weighted by each
    chunk's variant/site count** — never a naive mean of per-chromosome means. Ratios that
    derive exactly from counts (het:hom, indel:SNV, %multiallelic) are recomputed exactly.

## Manual / partial runs

`CHROMS="chr21 chr22" scripts/run_wgs_chunked.sh ...` restricts the run. To merge an
arbitrary set of per-chromosome results yourself:

```bash
python scripts/merge_wgs_results.py OUTDIR/chr*_results.json --out OUTDIR/merged.results.json
```

## Roadmap

Making `parse_vcf()` a streaming generator (the `iter_variants()` iterator already exists)
would remove the in-memory-list ceiling and let a single non-chunked run scale, at which
point chunking would remain useful only for the annotation heap and for parallelism.
