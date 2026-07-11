"""gnomAD frequencies from a DuckDB/Parquet store — batch, offline, instant.

This is the genomic-lakehouse approach: instead of ~11k per-variant tabix lookups
(latency-bound), load all the post-QC variants into DuckDB once and answer them with a
single vectorised LEFT JOIN against a columnar gnomAD Parquet — partition pruning +
row-group skipping make a whole exome's frequencies come back in ~seconds.

Usage in the pipeline: ``prime(variants)`` runs the one join and fills a per-key cache;
``gnomad.lookup`` then reads that cache first (before the local tabix / remote path).
Absent config, no duckdb, or no parquet -> a no-op / None, so everything falls through
unchanged (behaviour-preserving).

Point at it with ``VCF2REPORT_GNOMAD_PARQUET`` — a single ``.parquet`` file or a
Hive-partitioned dir (``chrom=chrN/*.parquet``). Expected columns: chrom, pos, ref, alt,
af, af_grpmax, ac, an, nhomalt, faf95 (optional), grpmax_pop (optional).

Safety: like the full-mode tabix, absence (af 0.0) is asserted ONLY for a contig the
parquet actually covers; a variant on an uncovered contig is left unprimed so the
caller falls back — never a fabricated absence.
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from .. import config
from ..models import Variant

_duckdb = None
_duckdb_tried = False
_primed: dict[str, dict] = {}
_lock = threading.Lock()


def _get_duckdb():
    global _duckdb, _duckdb_tried
    if not _duckdb_tried:
        _duckdb_tried = True
        try:
            import duckdb  # noqa: F401
            _duckdb = duckdb
        except Exception:
            _duckdb = None
    return _duckdb


def available() -> bool:
    p = config.GNOMAD_PARQUET
    return bool(p) and _get_duckdb() is not None and Path(p).exists()


def _norm_chrom(chrom: str) -> str:
    """Match the parquet's 'chr'-prefixed contig naming."""
    c = str(chrom)
    return c if c.lower().startswith("chr") else f"chr{c}"


def _source_expr() -> str:
    """A read_parquet(...) expression: a single file, or a dir's real parquet files
    (excluding macOS '._' AppleDouble sidecars that break the glob on exFAT)."""
    p = Path(config.GNOMAD_PARQUET)
    if p.is_dir():
        files = [str(f) for f in sorted(p.rglob("*.parquet")) if not f.name.startswith("._")]
        lst = ", ".join("'" + f.replace("'", "''") + "'" for f in files)
        return f"read_parquet([{lst}], hive_partitioning=true, union_by_name=true)"
    return f"read_parquet('{str(p)}')"


def prime(variants) -> int:
    """Batch-resolve every variant's gnomAD frequency in one DuckDB join and cache it.
    Returns the number of variants primed (0 if the feature is off)."""
    if not available() or not variants:
        return 0
    import csv
    import os
    import tempfile
    duckdb = _get_duckdb()
    con = duckdb.connect()
    tmp = None
    try:
        src = _source_expr()
        # Schema-aware: a lakehouse parquet may lack faf95/grpmax_pop; a from-scratch
        # build (build_gnomad_parquet.py) carries them. Select NULL for missing columns.
        schema = {r[0].lower() for r in con.execute(f"DESCRIBE SELECT * FROM {src}").fetchall()}
        col = lambda name: f"g.{name}" if name in schema else "NULL"
        covered = {r[0] for r in con.execute(f"SELECT DISTINCT chrom FROM {src}").fetchall()}
        # Bulk-load the query variants via a temp TSV — executemany is per-row and would
        # take minutes on a whole exome; read_csv loads ~24k rows in ~1 s.
        fd, tmp = tempfile.mkstemp(suffix=".tsv")
        with os.fdopen(fd, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            for v in variants:
                w.writerow([_norm_chrom(v.chrom), v.pos, v.ref, v.alt, v.key])
        con.execute(
            "CREATE TABLE q AS SELECT column0 AS chrom, CAST(column1 AS INTEGER) AS pos, "
            "column2 AS ref, column3 AS alt, column4 AS key "
            "FROM read_csv(?, header=false, delim='\t', all_varchar=true)", [tmp])
        rows = con.execute(f"""
            SELECT q.chrom, q.key, {col('af')}, {col('af_grpmax')}, {col('ac')},
                   {col('an')}, {col('nhomalt')},
                   TRY_CAST({col('faf95')} AS DOUBLE), {col('grpmax_pop')}
            FROM q LEFT JOIN {src} g
              ON g.chrom = q.chrom AND g.pos = q.pos AND g.ref = q.ref AND g.alt = q.alt
        """).fetchall()
    except Exception:
        con.close()
        if tmp and os.path.exists(tmp):
            os.remove(tmp)
        return 0
    con.close()
    if tmp and os.path.exists(tmp):
        os.remove(tmp)

    n = 0
    with _lock:
        for chrom, key, af, af_grpmax, ac, an, nhomalt, faf95, pop in rows:
            matched = not (af is None and af_grpmax is None and ac is None)
            if matched:
                _primed[key] = {"af": af_grpmax if af_grpmax is not None else af,
                                "ac": ac, "an": an, "hom": nhomalt, "faf95": faf95, "pop": pop}
                n += 1
            elif chrom in covered:
                # covered contig, no match -> genuine absence.
                _primed[key] = {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
                n += 1
            # uncovered contig -> leave unprimed so gnomad.lookup falls through.
    return n


def get(key: str) -> Optional[dict]:
    """A primed frequency dict for a variant key, or None (fall through)."""
    return _primed.get(key)


def _reset_for_tests() -> None:
    global _duckdb, _duckdb_tried
    with _lock:
        _primed.clear()
        _duckdb, _duckdb_tried = None, False
