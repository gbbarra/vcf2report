"""AlphaMissense from a DuckDB/Parquet store — batch, offline, instant.

Mirrors :mod:`gnomad_parquet`: :func:`prime` resolves every candidate's AlphaMissense
score in ONE chr-pruned DuckDB LEFT JOIN against the locus-partitioned Parquet
(``build_alphamissense_parquet.py``, MAX-per-locus so the join stays 1:1) and returns
``{key: {am_pathogenicity, am_class} or None}``. ``alphamissense.prime`` fills its cache
from this and falls back to the per-variant tabix loop when the store / duckdb is absent —
behaviour-preserving. AlphaMissense absence is never asserted: a locus with no missense
score simply maps to ``None`` (exactly as the tabix path's ``_best`` returns None).
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from .. import config

_duckdb = None
_duckdb_tried = False
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
    p = config.ALPHAMISSENSE_PARQUET
    return bool(p) and _get_duckdb() is not None and Path(p).exists()


def _norm_chrom(chrom: str) -> str:
    c = str(chrom)
    return c if c.lower().startswith("chr") else f"chr{c}"


def _q(path: str) -> str:
    return path.replace("'", "''")


def _source_expr() -> str:
    """read_parquet(...) over the store's real parquet files (excluding macOS '._' sidecars)."""
    p = Path(config.ALPHAMISSENSE_PARQUET)
    if p.is_dir():
        files = [str(f) for f in sorted(p.rglob("*.parquet")) if not f.name.startswith("._")]
        lst = ", ".join("'" + _q(f) + "'" for f in files)
        return f"read_parquet([{lst}], hive_partitioning=true, union_by_name=true)"
    return f"read_parquet('{_q(str(p))}')"


def prime(variants) -> Optional[dict]:
    """One chr-pruned LEFT JOIN → ``{key: {am_pathogenicity, am_class} or None}`` for every
    variant. None on any failure so the caller falls back to the tabix loop."""
    if not available() or not variants:
        return None
    import csv
    import os
    import tempfile
    duckdb = _get_duckdb()
    con = duckdb.connect()
    tmp = None
    try:
        src = _source_expr()
        fd, tmp = tempfile.mkstemp(suffix=".tsv")
        chroms = set()
        with os.fdopen(fd, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            for v in variants:
                c = _norm_chrom(v.chrom)
                chroms.add(c)
                w.writerow([c, v.pos, v.ref.upper(), v.alt.upper(), v.key])
        con.execute(
            "CREATE TABLE q AS SELECT column0 AS chrom, CAST(column1 AS INTEGER) AS pos, "
            "column2 AS ref, column3 AS alt, column4 AS key "
            "FROM read_csv(?, header=false, delim='\t', all_varchar=true)", [tmp])
        # Prune the store to the candidate chroms (partition pruning -> sub-second).
        chlist = ", ".join("'" + _q(c) + "'" for c in sorted(chroms))
        rows = con.execute(f"""
            SELECT q.key, a.am_pathogenicity, a.am_class
            FROM q LEFT JOIN (SELECT * FROM {src} WHERE chrom IN ({chlist})) a
              ON a.chrom = q.chrom AND a.pos = q.pos
             AND upper(a.ref) = q.ref AND upper(a.alt) = q.alt
        """).fetchall()
    except Exception:
        con.close()
        if tmp and os.path.exists(tmp):
            os.remove(tmp)
        return None
    con.close()
    if tmp and os.path.exists(tmp):
        os.remove(tmp)
    out: dict = {}
    for key, am_path, am_class in rows:
        out[key] = ({"am_pathogenicity": float(am_path), "am_class": am_class}
                    if am_path is not None else None)
    return out


def _reset_for_tests() -> None:
    global _duckdb, _duckdb_tried
    _duckdb, _duckdb_tried = None, False
