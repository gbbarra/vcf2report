"""ClinVar from a DuckDB/Parquet store — batch, offline, instant.

Mirrors :mod:`gnomad_parquet`: :func:`prime` resolves every post-QC variant's ClinVar
record in ONE chr-pruned DuckDB LEFT JOIN against the locus-partitioned Parquet
(``build_clinvar_parquet.py``, weekly), replacing ClinVar's per-variant tabix loop —
the one source that previously had no batch path. ``clinvar.lookup`` reads :func:`get`
FIRST (before the tabix / live / slice fallback), so a store miss falls through exactly
like a tabix miss (never a fabricated 'no record'). Absent store / duckdb -> a no-op.

The returned fields (significance, review_status, accession, condition) are the SAME
values the tabix ``_tabix_lookup`` serves — both derive from ``clinvar_grch38.tsv.gz`` —
so the classification is byte-identical to the tabix path.
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

from .. import config

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
    p = config.CLINVAR_PARQUET
    return bool(p) and _get_duckdb() is not None and Path(p).exists()


def _norm_chrom(chrom: str) -> str:
    c = str(chrom)
    return c if c.lower().startswith("chr") else f"chr{c}"


def _q(path: str) -> str:
    return path.replace("'", "''")


def _source_expr() -> str:
    p = Path(config.CLINVAR_PARQUET)
    if p.is_dir():
        files = [str(f) for f in sorted(p.rglob("*.parquet")) if not f.name.startswith("._")]
        lst = ", ".join("'" + _q(f) + "'" for f in files)
        return f"read_parquet([{lst}], hive_partitioning=true, union_by_name=true)"
    return f"read_parquet('{_q(str(p))}')"


def prime(variants) -> int:
    """One chr-pruned LEFT JOIN → cache the ClinVar record for every MATCHED variant.
    A miss is left unprimed so :func:`get` returns None and ``clinvar.lookup`` falls
    through to tabix/live/slice. Returns the number of records cached (0 if off)."""
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
        chlist = ", ".join("'" + _q(c) + "'" for c in sorted(chroms))
        rows = con.execute(f"""
            SELECT q.key, c.significance, c.review_status, c.accession, c.condition
            FROM q JOIN (SELECT * FROM {src} WHERE chrom IN ({chlist})) c
              ON c.chrom = q.chrom AND c.pos = q.pos
             AND upper(c.ref) = q.ref AND upper(c.alt) = q.alt
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
        for key, sig, rev, acc, cond in rows:
            if not sig:
                continue
            _primed[key] = {"significance": sig, "review_status": rev or None,
                            "accession": acc or None, "condition": cond or None,
                            "date": "local snapshot"}
            n += 1
    return n


def get(key: str) -> Optional[dict]:
    """A primed ClinVar record for a variant key, or None (fall through to tabix/live/slice)."""
    return _primed.get(key)


def _reset_for_tests() -> None:
    global _duckdb, _duckdb_tried
    with _lock:
        _primed.clear()
        _duckdb, _duckdb_tried = None, False
