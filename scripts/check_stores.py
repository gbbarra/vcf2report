#!/usr/bin/env python3
"""Health check for the annotation Parquet stores (gnomAD · AlphaMissense · ClinVar).

Reports, per store: presence, size, row count, integrity (reads cleanly = not corrupt),
completeness (all 24 core contigs present AND the row count matches the build manifest), the
build date + source version, and — by each source's cadence — whether a refresh is due
(ClinVar weekly; gnomAD v4.1 / AlphaMissense frozen).

    python3 scripts/check_stores.py             # human report
    python3 scripts/check_stores.py --json        # machine-readable (monitoring)
    python3 scripts/check_stores.py --quick        # presence + size + freshness (no row scan)
    python3 scripts/check_stores.py clinvar         # one store

Exit code is non-zero if any store is missing / corrupt / incomplete / stale, so it doubles
as a cron / CI monitoring probe.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import stores  # noqa: E402

_ICON = {"ok": "OK ", "stale": "STALE", "incomplete": "BAD", "corrupt": "BAD", "missing": "—"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("store", nargs="?", choices=["gnomad", "alphamissense", "clinvar"])
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--quick", action="store_true", help="skip the row-count / integrity scan")
    a = ap.parse_args()

    health = stores.store_health(name=a.store, measure=not a.quick)

    if a.json:
        print(json.dumps(health, indent=2, default=str))
    else:
        print(f"\n  vcf2report data stores{'  (quick)' if a.quick else ''}")
        print("  " + "-" * 70)
        for name, e in health.items():
            tag = _ICON.get(e["status"], "?")
            line = f"  [{tag:5s}] {name:13s} {e.get('size', '—'):>9}"
            if e.get("rows") is not None:
                line += f"  {e['rows']:>12,} rows  {e.get('chroms_present', '?')} chr"
            print(line)
            if not e["present"]:
                print(f"          {e['reason']}")
                continue
            det = []
            if e.get("readable") is not None:
                det.append(f"integrity={'ok' if e['readable'] else 'FAIL'}")
            if e.get("complete") is not None:
                det.append(f"complete={'yes' if e['complete'] else 'NO'}")
            if e.get("built_at"):
                age = f" ({int(e['age_days'])}d ago)" if e.get("age_days") is not None else ""
                det.append(f"built={e['built_at'][:10]}{age}")
            src = (e.get("source") or {})
            if src.get("release"):
                det.append(f"{src.get('name', '')} {src['release']}")
            if det:
                print(f"          {' · '.join(det)}")
            print(f"          → {e['reason']}")
            if e.get("missing_core_chroms"):
                print(f"          ⚠ missing core contigs: {e['missing_core_chroms']}")
            if e.get("rows_mismatch"):
                print(f"          ⚠ row-count mismatch (corrupt/partial): {e['rows_mismatch']}")
        print()

    bad = {n: e["status"] for n, e in health.items() if e["status"] != "ok"}
    if bad and not a.json:
        print(f"  ⚠ action needed: {bad}\n")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
