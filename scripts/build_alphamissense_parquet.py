#!/usr/bin/env python3
"""Convert the AlphaMissense hg38 tabix TSV into a locus-partitioned Parquet store for
the batch DuckDB annotate join (see docs/DATA_ARCHITECTURE.md).

Pre-aggregates to ONE row per (chrom,pos,ref,alt) = the MAX am_pathogenicity across
transcripts and its class, so the LEFT JOIN in the annotate stage stays 1:1 (no row
fan-out) and the per-candidate tabix loop disappears. Missense-SNV-only: most gnomAD/
ClinVar loci correctly LEFT-JOIN to NULL (AlphaMissense absence is never asserted).

    python3 scripts/build_alphamissense_parquet.py [SRC.tsv.gz] [OUT_DIR]

Defaults: SRC = config.ALPHAMISSENSE_LOCAL, OUT = config.ALPHAMISSENSE_PARQUET.

⚠️  LICENSE: the AlphaMissense file header declares CC BY-NC-SA 4.0 (non-commercial,
share-alike). Do NOT commit or redistribute the derived Parquet; it is git-ignored.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

try:
    import duckdb
except ImportError:
    raise SystemExit("Needs duckdb: pip install duckdb")

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else config.ALPHAMISSENSE_LOCAL
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else config.ALPHAMISSENSE_PARQUET

if not SRC.exists():
    raise SystemExit(f"AlphaMissense TSV not found: {SRC}\n"
                     f"Fetch it first: bash scripts/fetch_alphamissense.sh")

# Build into a temp dir then swap, so a reader never sees a half-written store.
TMP = OUT.with_name(OUT.name + ".building")
OUT.parent.mkdir(parents=True, exist_ok=True)

con = duckdb.connect()
# comment='#' drops the 3 licence lines AND the #CHROM header; we supply the schema.
con.execute(f"""
COPY (
  SELECT chrom, pos, ref, alt,
         MAX(am_pathogenicity)                    AS am_pathogenicity,
         arg_max(am_class, am_pathogenicity)      AS am_class
  FROM read_csv('{SRC.as_posix()}', delim='\t', header=false, comment='#',
                auto_detect=false, quote='', escape='', strict_mode=false,
                columns={{'chrom':'VARCHAR','pos':'INTEGER','ref':'VARCHAR','alt':'VARCHAR',
                          'genome':'VARCHAR','uniprot_id':'VARCHAR','transcript_id':'VARCHAR',
                          'protein_variant':'VARCHAR','am_pathogenicity':'DOUBLE','am_class':'VARCHAR'}})
  GROUP BY chrom, pos, ref, alt
) TO '{TMP.as_posix()}' (FORMAT PARQUET, PARTITION_BY (chrom), OVERWRITE_OR_IGNORE);
""")
n = con.execute(f"SELECT count(*) FROM read_parquet('{TMP.as_posix()}/**/*.parquet')").fetchone()[0]
con.close()

import shutil  # noqa: E402
if OUT.exists():
    shutil.rmtree(OUT)
TMP.rename(OUT)
print(f"[am-parquet] wrote {n:,} loci (MAX-per-locus) -> {OUT}")
from vcf2report import stores  # noqa: E402
stores.write_manifest("alphamissense", path=str(OUT))
print("[am-parquet] _manifest.json stamped — verify with scripts/check_stores.py")
