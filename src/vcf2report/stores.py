"""Health, integrity, and freshness of the local annotation Parquet stores.

One place to answer: are gnomAD / AlphaMissense / ClinVar present, intact, complete, and
current? :func:`store_health` measures each store directly (size, row count, chromosomes,
schema, readability) and cross-checks a ``_manifest.json`` sidecar (build date, source, the
row count at build time) to detect truncation/corruption and to decide, by each source's
update cadence (ClinVar weekly; gnomAD v4.1 / AlphaMissense frozen), whether a refresh is due.

MCP-free so both the MCP ``data_status`` tool and ``scripts/check_stores.py`` can call it.
``write_manifest`` stamps a store after a build (or a one-off ``scripts/stamp_store_manifest.py``).
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from . import config

MANIFEST = "_manifest.json"
# Completeness floor: every store must cover the 24 primary contigs; extra contigs
# (ClinVar's MT / alt scaffolds, etc.) are a bonus, never a requirement.
CORE_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

_duckdb = None
_duckdb_tried = False


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


def _registry() -> dict:
    """name -> spec, with paths resolved at call time (env overrides stay live)."""
    return {
        "gnomad": {
            "path": config._resolve_gnomad_parquet(),
            "key_columns": {"chrom", "pos", "ref", "alt", "af"},
            "cadence": "frozen", "stale_after_days": None,
            "source": {"name": "gnomAD", "release": "v4.1",
                       "url": "gs://gcp-public-data--gnomad", "note": "frozen release"},
            "enables": "PM2 / BA1 / BS1 (population frequency)",
        },
        "alphamissense": {
            "path": config.ALPHAMISSENSE_PARQUET,
            "key_columns": {"chrom", "pos", "ref", "alt", "am_pathogenicity", "am_class"},
            "cadence": "frozen", "stale_after_days": None,
            "source": {"name": "AlphaMissense hg38", "release": "2023",
                       "license": "CC BY-NC-SA 4.0", "note": "frozen release"},
            "enables": "PP3 / BP4 (missense pathogenicity)",
        },
        "clinvar": {
            "path": config.CLINVAR_PARQUET,
            "key_columns": {"chrom", "pos", "ref", "alt", "significance",
                            "review_status", "review_stars"},
            "cadence": "weekly", "stale_after_days": 14,
            "source": {"name": "ClinVar GRCh38", "release": "weekly",
                       "url": "ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"},
            "enables": "PS1 / PM5 / PP5 / BP6 + the >=2-star ClinVar safety flag",
        },
    }


def _parquet_files(path: Path) -> list[str]:
    if path.is_dir():
        return [str(f) for f in sorted(path.rglob("*.parquet")) if not f.name.startswith("._")]
    return [str(path)] if path.exists() else []


def _dir_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    total = 0
    for f in path.rglob("*"):
        if f.is_file() and not f.name.startswith("._"):
            try:
                total += f.stat().st_size
            except OSError:
                pass
    return total


def _human(n: int) -> str:
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.0f} {unit}" if unit == "B" else f"{x:.1f} {unit}"
        x /= 1024
    return f"{x:.1f} TB"


def _measure(path: Path) -> dict:
    """Row count, chromosomes, and schema via DuckDB — and thus whether the store READS
    (a corrupt/truncated Parquet raises here). Returns readable=False + error on failure."""
    duckdb = _get_duckdb()
    if duckdb is None:
        return {"readable": None, "error": "duckdb not installed", "rows": None,
                "chroms": [], "schema": []}
    files = _parquet_files(path)
    if not files:
        return {"readable": False, "error": "no parquet files", "rows": None,
                "chroms": [], "schema": []}
    q = ", ".join("'" + f.replace("'", "''") + "'" for f in files)
    src = f"read_parquet([{q}], hive_partitioning=true, union_by_name=true)"
    con = duckdb.connect()
    try:
        schema = [r[0].lower() for r in con.execute(f"DESCRIBE SELECT * FROM {src}").fetchall()]
        rows = con.execute(f"SELECT count(*) FROM {src}").fetchone()[0]
        chroms = [r[0] for r in con.execute(
            f"SELECT DISTINCT chrom FROM {src}").fetchall()] if "chrom" in schema else []
    except Exception as exc:
        con.close()
        return {"readable": False, "error": str(exc)[:200], "rows": None,
                "chroms": [], "schema": []}
    con.close()
    return {"readable": True, "error": None, "rows": rows,
            "chroms": sorted(chroms), "schema": schema, "n_files": len(files)}


def _read_manifest(path: Path) -> Optional[dict]:
    mf = (path / MANIFEST) if path.is_dir() else path.with_name(path.name + MANIFEST)
    try:
        if mf.exists():
            return json.loads(mf.read_text())
    except Exception:
        return None
    return None


def _age_days(built_at: Optional[str]) -> Optional[float]:
    if not built_at:
        return None
    try:
        dt = datetime.fromisoformat(built_at.replace("Z", "+00:00"))
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return (datetime.now(timezone.utc) - dt).total_seconds() / 86400.0
    except Exception:
        return None


def store_health(name: Optional[str] = None, measure: bool = True) -> dict:
    """Health of one store (``name``) or all three. Each entry reports present, size, rows,
    chromosomes, schema_ok, readable (integrity), the manifest (build date + source), age,
    completeness, and an update recommendation by cadence. ``measure=False`` skips the DuckDB
    row-count/integrity read (presence + size + manifest only) for a fast dependency check."""
    reg = _registry()
    names = [name] if name else list(reg)
    out = {}
    for nm in names:
        spec = reg[nm]
        path = Path(spec["path"]) if spec["path"] else None
        present = bool(path and path.exists())
        e: dict = {"present": present, "path": str(path) if path else None,
                   "enables": spec["enables"], "cadence": spec["cadence"]}
        if not present:
            e.update(status="missing", update_recommended=True,
                     reason="store absent — build it (see docs/DATA_ARCHITECTURE.md)")
            out[nm] = e
            continue
        e["size_bytes"] = _dir_size(path)
        e["size"] = _human(e["size_bytes"])
        mf = _read_manifest(path)
        e["manifest"] = mf
        built_at = (mf or {}).get("built_at")
        e["built_at"] = built_at
        e["source"] = (mf or {}).get("source", spec["source"])
        age = _age_days(built_at)
        e["age_days"] = round(age, 1) if age is not None else None

        schema_ok = readable = complete = None
        if measure:
            m = _measure(path)
            e["rows"] = m["rows"]
            e["chroms_present"] = len(m["chroms"])
            e["readable"] = m["readable"]
            if m["error"]:
                e["error"] = m["error"]
            schema_ok = bool(m["schema"]) and spec["key_columns"].issubset(set(m["schema"]))
            e["schema_ok"] = schema_ok
            missing_chroms = sorted(CORE_CHROMS - set(m["chroms"]))
            e["missing_core_chroms"] = missing_chroms
            readable = m["readable"]
            # Completeness: manifest row count must match (detects truncation), core contigs present.
            expect_rows = (mf or {}).get("rows")
            rows_match = (expect_rows is None) or (m["rows"] == expect_rows)
            e["rows_expected"] = expect_rows
            complete = bool(readable and schema_ok and not missing_chroms and rows_match)
            e["complete"] = complete
            if expect_rows is not None and m["rows"] != expect_rows:
                e["rows_mismatch"] = {"expected": expect_rows, "actual": m["rows"]}

        # Update recommendation by cadence.
        stale = False
        reason = "frozen release — no routine update"
        if spec["cadence"] == "weekly":
            if age is None:
                reason = "no build date recorded — stamp/rebuild to track freshness"
            elif age > spec["stale_after_days"]:
                stale = True
                reason = (f"built {int(age)} d ago; ClinVar releases weekly — "
                          f"rebuild (older than {spec['stale_after_days']} d)")
            else:
                reason = f"built {int(age)} d ago — within the weekly window"
        e["update_recommended"] = bool(stale) or (measure and complete is False) or \
            (readable is False)
        e["reason"] = reason
        # Overall status.
        if measure and readable is False:
            e["status"] = "corrupt"
        elif measure and complete is False:
            e["status"] = "incomplete"
        elif stale:
            e["status"] = "stale"
        else:
            e["status"] = "ok"
        out[nm] = e
    return out


# The Parquet stores an analysis needs before it may run. Missing / corrupt / incomplete
# among these BLOCK the run; merely stale (an old but intact ClinVar) only warns.
REQUIRED = ("gnomad", "alphamissense", "clinvar")
_BLOCKING = {"missing", "corrupt", "incomplete"}


def gate(required=REQUIRED, measure: bool = True) -> dict:
    """Analysis-readiness gate for the guided flow's Stage 1.

    Returns ``{ready, blocking, stale, health}``: ``ready`` is False (do NOT run the analysis)
    when any ``required`` store is missing / corrupt / incomplete; a stale store is listed in
    ``stale`` but does not block (it is intact, just past its refresh window)."""
    health = store_health(measure=measure)
    blocking = [n for n in required if health.get(n, {}).get("status") in _BLOCKING]
    stale = [n for n, e in health.items() if e.get("status") == "stale"]
    return {"ready": not blocking, "blocking": blocking, "stale": stale, "health": health}


def write_manifest(name: str, source: Optional[dict] = None,
                   built_at: Optional[str] = None, path: Optional[str] = None) -> dict:
    """Measure a built store and write its ``_manifest.json`` (build date, source, row count,
    chromosomes, schema). Called by the build scripts (which pass their output ``path``) and by
    scripts/stamp_store_manifest.py (which stamps the registry path)."""
    spec = _registry()[name]
    path = Path(path) if path else (Path(spec["path"]) if spec["path"] else None)
    if not path or not path.exists():
        raise SystemExit(f"{name}: store not found at {path}")
    m = _measure(path)
    if not m["readable"]:
        raise SystemExit(f"{name}: store unreadable — {m.get('error')}")
    manifest = {
        "store": name,
        "format": "parquet",
        "built_at": built_at or datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": source or spec["source"],
        "rows": m["rows"],
        "chroms": m["chroms"],
        "schema": m["schema"],
        "update": {"cadence": spec["cadence"], "stale_after_days": spec["stale_after_days"]},
    }
    mf = (path / MANIFEST) if path.is_dir() else path.with_name(path.name + MANIFEST)
    mf.write_text(json.dumps(manifest, indent=2))
    return manifest


def _reset_for_tests() -> None:
    global _duckdb, _duckdb_tried
    _duckdb, _duckdb_tried = None, False
