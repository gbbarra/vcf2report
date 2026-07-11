"""Local gnomAD frequency client (reduced tabix table).

Reads a small, offline, tabix-indexed TSV built by ``scripts/build_gnomad_local.py``
that carries only the fields the ACMG engine cites — grpmax AF, AC/AN, homozygote
count, and the ClinGen filtering AF (faf95) — one row per gnomAD variant. It is the
same reduction ``gnomad_remote`` computes per record, so local and remote agree; the
local table just makes it offline and instant (no 150 GB full download, no network).

Absent file or missing pysam -> returns None so ``gnomad.lookup`` falls back to the
remote/live/bundled path exactly as before (behaviour-preserving).

Schema (tab-separated, bgzipped, ``tabix -s1 -b2 -e2``; ``#`` header skipped):
    #chrom  pos  ref  alt  af  ac  an  hom  faf95  pop
``chrom`` has no ``chr`` prefix (matches :pyattr:`Variant.key`); ``pop`` may be empty.

Miss semantics (the safety invariant — a table must NEVER fabricate an absence for a
variant that is actually present in gnomAD):
  * **partial** (panel / per-VCF): a value is returned ONLY on an exact match; any
    miss -> None (the caller falls back). It never asserts absence.
  * **full**: asserts absence (af 0.0) on a miss ONLY for a contig the build actually
    covered (recorded in the sidecar ``contigs``). A query on an uncovered contig
    (chrM/MT, a chromosome whose stream failed, an alt/decoy) returns None -> fallback.
"""
from __future__ import annotations

import threading
from typing import Optional

from .. import config
from ..models import Variant

_tabix = None
_tabix_tried = False
_pysam = None
_pysam_tried = False
_mode: Optional[str] = None                 # "full" | "partial"
_contigs: frozenset = frozenset()           # covered contigs (no 'chr') for a full table
_lock = threading.Lock()


def _get_pysam():
    global _pysam, _pysam_tried
    if not _pysam_tried:
        _pysam_tried = True
        try:
            import pysam  # noqa: F401
            _pysam = pysam
        except Exception:
            _pysam = None
    return _pysam


def _norm_chrom(chrom: str) -> str:
    chrom = str(chrom)
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _read_meta(fp) -> tuple[str, frozenset]:
    """(mode, covered-contigs) from the ``<table>.meta`` sidecar. Safe defaults —
    ``partial`` and an empty contig set — so a missing/corrupt sidecar never lets the
    client assert an absence."""
    mode, contigs = "partial", frozenset()
    try:
        import json
        meta = fp.with_name(fp.name + ".meta")
        if meta.exists():
            d = json.loads(meta.read_text())
            m = (d.get("mode") or "").lower()
            if m in ("full", "partial"):
                mode = m
            cs = d.get("contigs")
            if isinstance(cs, list):
                contigs = frozenset(_norm_chrom(str(c)) for c in cs)
    except Exception:
        pass
    return mode, contigs


def _open():
    """Lazily open the tabix file once; None if unavailable (missing file/pysam)."""
    global _tabix, _tabix_tried, _mode, _contigs
    if _tabix_tried:
        return _tabix
    with _lock:
        if _tabix_tried:
            return _tabix
        pysam = _get_pysam()
        fp = config.GNOMAD_LOCAL_TABIX
        tabix, mode, contigs = None, "partial", frozenset()
        if pysam is not None and fp.exists():
            try:
                tabix = pysam.TabixFile(str(fp))
                mode, contigs = _read_meta(fp)
            except Exception:
                tabix = None
        # Publish the fully-initialised state, THEN flip the tried flag — so a
        # concurrent reader on the unlocked fast path never sees tried=True while
        # _tabix / _mode / _contigs are still half-assigned (DCL-safe).
        _tabix, _mode, _contigs = tabix, mode, contigs
        _tabix_tried = True
        return _tabix


def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _to_int(s: str) -> Optional[int]:
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


def query(variant: Variant) -> Optional[dict]:
    """Return {'af','ac','an','hom','faf95','pop'} for an exact CHROM-POS-REF-ALT
    match, ``{...af:0.0...}`` for a genuine absence a full table can vouch for, or
    None when no local answer is available (caller falls back)."""
    tabix = _open()
    if tabix is None:
        return None
    chrom = _norm_chrom(variant.chrom)
    rows: list = []
    for name in (chrom, "chr" + chrom):   # tolerate either contig naming in the file
        try:
            rows = list(tabix.fetch(name, variant.pos - 1, variant.pos))
        except (ValueError, KeyError):
            rows = []
        if rows:
            break
    ref, alt = variant.ref.upper(), variant.alt.upper()   # alleles are case-insensitive
    for row in rows:
        f = row.split("\t")
        if len(f) < 4:
            continue
        if _to_int(f[1]) != variant.pos or f[2].upper() != ref or f[3].upper() != alt:
            continue
        return {"af": _to_float(f[4]) if len(f) > 4 else None,
                "ac": _to_int(f[5]) if len(f) > 5 else None,
                "an": _to_int(f[6]) if len(f) > 6 else None,
                "hom": _to_int(f[7]) if len(f) > 7 else None,
                "faf95": _to_float(f[8]) if len(f) > 8 else None,
                "pop": (f[9].strip() or None) if len(f) > 9 else None}
    # No exact allele match. A partial table cannot assert absence at all; a full
    # table may — but ONLY for a contig it actually covered (else chrM/MT, an
    # alt/decoy, or a chromosome whose build stream failed would be fabricated as 0.0).
    if _mode == "full" and chrom in _contigs:
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
    return None


def _reset_for_tests() -> None:
    global _tabix, _tabix_tried, _mode, _contigs
    with _lock:
        _tabix, _tabix_tried, _mode, _contigs = None, False, None, frozenset()
