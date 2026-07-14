#!/usr/bin/env python3
"""Write / refresh a store's ``_manifest.json`` (build date, source, row count, schema) for an
already-built store, so the health check (scripts/check_stores.py) can verify completeness and
freshness. The build scripts stamp automatically; run this to stamp stores built earlier.

    python3 scripts/stamp_store_manifest.py [gnomad|alphamissense|clinvar|all]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import stores  # noqa: E402

if __name__ == "__main__":
    targets = sys.argv[1:] or ["all"]
    names = list(stores._registry()) if "all" in targets else targets
    rc = 0
    for n in names:
        try:
            mf = stores.write_manifest(n)
            print(f"[{n}] stamped — {mf['rows']:,} rows, built_at {mf['built_at']}")
        except SystemExit as exc:
            print(f"[{n}] skipped — {exc}")
            rc = 1
    sys.exit(rc)
