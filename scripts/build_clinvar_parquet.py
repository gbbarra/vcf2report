#!/usr/bin/env python3
"""Convert the local ClinVar tabix TSV into a locus-partitioned Parquet store for the
batch DuckDB annotate join (see docs/DATA_ARCHITECTURE.md).

Reads the store built by scripts/build_clinvar_local.py (chrom pos ref alt significance
review_status accession condition), adds a precomputed review_stars column (0-4, the
≥2★ safety-flag gate), chr-prefixes the contig to match the gnomAD/query key, and writes
one partition per chromosome. This is the WEEKLY-refresh store — rebuilding it touches no
other source (gnomAD/AlphaMissense stay frozen).

    python3 scripts/build_clinvar_parquet.py [SRC.tsv.gz] [OUT_DIR]

Defaults: SRC = config.CLINVAR_TABIX, OUT = config.CLINVAR_PARQUET.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

try:
    import duckdb
except ImportError:
    raise SystemExit("Needs duckdb: pip install duckdb")

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else config.CLINVAR_TABIX
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else config.CLINVAR_PARQUET

if not SRC.exists():
    raise SystemExit(f"ClinVar TSV not found: {SRC}\n"
                     f"Build it first: python3 scripts/build_clinvar_local.py <clinvar.vcf.gz>")

TMP = OUT.with_name(OUT.name + ".building")
OUT.parent.mkdir(parents=True, exist_ok=True)

# review_stars mirrors report.assemble.clinvar_stars (kept simple + stable). review_status
# is already space-normalised by build_clinvar_local.py.
con = duckdb.connect()
con.execute(f"""
COPY (
  SELECT CASE WHEN chrom LIKE 'chr%' THEN chrom ELSE 'chr' || chrom END AS chrom,
         pos, ref, alt, significance, review_status,
         CASE
           WHEN review_status ILIKE '%practice guideline%'                                   THEN 4
           WHEN review_status ILIKE '%reviewed by expert panel%'                             THEN 3
           WHEN review_status ILIKE '%multiple submitters%' AND review_status ILIKE '%no conflict%' THEN 2
           WHEN review_status ILIKE 'criteria provided%' OR review_status ILIKE '%single submitter%'
                OR review_status ILIKE '%conflicting%'                                        THEN 1
           ELSE 0
         END AS review_stars,
         accession, condition
  FROM read_csv('{SRC.as_posix()}', delim='\t', header=false,
                auto_detect=false, quote='', escape='', strict_mode=false,
                columns={{'chrom':'VARCHAR','pos':'INTEGER','ref':'VARCHAR','alt':'VARCHAR',
                          'significance':'VARCHAR','review_status':'VARCHAR',
                          'accession':'VARCHAR','condition':'VARCHAR'}})
) TO '{TMP.as_posix()}' (FORMAT PARQUET, PARTITION_BY (chrom), OVERWRITE_OR_IGNORE);
""")
n = con.execute(f"SELECT count(*) FROM read_parquet('{TMP.as_posix()}/**/*.parquet')").fetchone()[0]
con.close()

import shutil  # noqa: E402
if OUT.exists():
    shutil.rmtree(OUT)
TMP.rename(OUT)
print(f"[clinvar-parquet] wrote {n:,} variants (review_stars precomputed) -> {OUT}")
from vcf2report import stores  # noqa: E402
stores.write_manifest("clinvar", path=str(OUT))
print("[clinvar-parquet] _manifest.json stamped — verify with scripts/check_stores.py")
