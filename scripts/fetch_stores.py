#!/usr/bin/env python3
"""Provision the three annotation Parquet stores using only Python — no `gh`, no `zstd` binary.

``scripts/setup_stores.sh`` already does this, but it hard-requires the GitHub CLI and the
zstd command-line tool, and exits if either is missing. Neither is present in a bare container,
which is where provisioning matters most: it is the difference between the store gate blocking
every real analysis and the pipeline actually running.

Neither requirement was necessary. The gnomAD and ClinVar stores live on this repository's
*public* releases, so their assets download unauthenticated over plain HTTPS; and zstd
decompression is available as the pure-Python `zstandard` wheel. Measured on a cold container:

    gnomAD v4.1     1.4 GB    21 s     69,898,057 rows
    AlphaMissense   0.64 GB    6 s     71,034,269 rows   (+ local build)
    ClinVar          50 MB    1.4 s     4,195,020 rows

    python3 scripts/fetch_stores.py               # all three, skipping any already present
    python3 scripts/fetch_stores.py gnomad clinvar
    python3 scripts/fetch_stores.py --force clinvar

AlphaMissense is fetched from DeepMind under **CC BY 4.0** and built locally; this project does
not redistribute it. gnomAD is ODbL-1.0 (attribution + share-alike); ClinVar is public domain.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import subprocess
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

RELEASES = "https://api.github.com/repos/gbbarra/vcf2report/releases/tags/{tag}"
ALPHAMISSENSE_SRC = (
    "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
)

STORES = {
    "gnomad": {
        "dest": ROOT / "data/gnomad/gnomad_parquet",
        "extract_to": ROOT / "data/gnomad",
        "tag": "gnomad-parquet-v4.1",
        "licence": "gnomAD v4.1 — ODbL-1.0 (attribution + share-alike)",
    },
    "clinvar": {
        "dest": ROOT / "data/clinvar/clinvar_parquet",
        "extract_to": ROOT / "data/clinvar",
        "tag": "clinvar-parquet-latest",
        "licence": "ClinVar (NCBI) — public domain",
    },
    "alphamissense": {
        "dest": ROOT / "data/alphamissense/am_parquet",
        "extract_to": ROOT / "data/alphamissense",
        "tag": None,  # built locally from the DeepMind source
        "licence": "AlphaMissense hg38 — CC BY 4.0, fetched from DeepMind, NOT redistributed",
    },
}


def _get(url: str, dest: Path) -> int:
    """Stream a URL to disk. Chunked so a 1.4 GB asset never sits in memory."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=1800) as r, open(dest, "wb") as fh:
        while chunk := r.read(1 << 22):
            fh.write(chunk)
    return dest.stat().st_size


def _asset_url(tag: str) -> str:
    """The public browser_download_url — no token, no `gh`. The releases are public."""
    with urllib.request.urlopen(RELEASES.format(tag=tag), timeout=120) as r:
        rel = json.load(r)
    for a in rel.get("assets", []):
        if a["name"].endswith(".tar.zst"):
            return a["browser_download_url"]
    raise SystemExit(f"no .tar.zst asset on release {tag!r}")


def _unzstd_tar(blob_path: Path, into: Path) -> None:
    """Decompress+untar, preferring the zstd binary and falling back to the wheel."""
    into.mkdir(parents=True, exist_ok=True)
    if shutil.which("zstd"):
        subprocess.run(
            f"zstd -dc {blob_path} | tar -x -C {into}", shell=True, check=True
        )
        return
    try:
        import zstandard
    except ImportError:
        raise SystemExit(
            "neither the `zstd` binary nor the `zstandard` package is available.\n"
            "  pip install zstandard      (or: apt install zstd / brew install zstd)"
        )
    with (
        open(blob_path, "rb") as fh,
        zstandard.ZstdDecompressor().stream_reader(fh) as sr,
    ):
        with tarfile.open(fileobj=sr, mode="r|") as tf:
            tf.extractall(into)


def fetch_release_store(name: str, force: bool) -> bool:
    spec = STORES[name]
    if spec["dest"].exists() and not force:
        print(
            f"  {name}: already present at {spec['dest'].relative_to(ROOT)} — skipping"
        )
        return False
    if force and spec["dest"].exists():
        shutil.rmtree(spec["dest"])
    print(f"  {name}: {spec['licence']}")
    url = _asset_url(spec["tag"])
    tmp = ROOT / f".{name}.tar.zst"
    try:
        t0 = time.perf_counter()
        size = _get(url, tmp)
        print(f"    downloaded {size / 1e6:.0f} MB in {time.perf_counter() - t0:.0f}s")
        _unzstd_tar(tmp, spec["extract_to"])
    finally:
        tmp.unlink(missing_ok=True)
    return True


def fetch_alphamissense(force: bool) -> bool:
    spec = STORES["alphamissense"]
    if spec["dest"].exists() and not force:
        print(f"  alphamissense: already present — skipping")
        return False
    print(f"  alphamissense: {spec['licence']}")
    src = spec["extract_to"] / "AlphaMissense_hg38.tsv.gz"
    if not src.exists():
        t0 = time.perf_counter()
        size = _get(ALPHAMISSENSE_SRC, src)
        print(f"    downloaded {size / 1e9:.2f} GB in {time.perf_counter() - t0:.0f}s")
    print("    building the Parquet store (this is the slow step)…")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/build_alphamissense_parquet.py"),
            str(src),
        ],
        check=True,
    )
    return True


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # No argparse `choices` here: with nargs="*" it also validates the default, which turns
    # "run with no arguments" — the common case — into a usage error.
    ap.add_argument(
        "stores",
        nargs="*",
        metavar="STORE",
        help=f"which stores to provision: {', '.join(STORES)} (default: all)",
    )
    ap.add_argument(
        "--force", action="store_true", help="re-fetch even if already present"
    )
    a = ap.parse_args()
    wanted = a.stores or list(STORES)
    if unknown := [s for s in wanted if s not in STORES]:
        raise SystemExit(f"unknown store(s) {unknown}; choose from {list(STORES)}")

    print(
        "\nProvisioning annotation stores. Nothing here needs `gh` or the zstd binary.\n"
    )
    for name in wanted:
        if name == "alphamissense":
            fetch_alphamissense(a.force)
        else:
            fetch_release_store(name, a.force)

    print("\n== store gate ==")
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_stores.py"), "--gate"]
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
