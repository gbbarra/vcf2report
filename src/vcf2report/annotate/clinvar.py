"""ClinVar clinical-significance client.

Resolution order: on-disk cache -> live NCBI E-utilities (unless OFFLINE) ->
bundled ClinVar slice (TSV keyed by CHROM-POS-REF-ALT). The local slice is the
authoritative, offline, deterministic coordinate lookup; the live path is a
best-effort enrichment that only accepts an exact location match. Returns
significance, review status, accession, condition, and date.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from .. import config
from ..models import Variant
from . import cache

_SOURCE = "clinvar"
_local: Optional[dict] = None
_COLUMNS = ["key", "significance", "review_status", "accession", "condition", "date"]


def _load_local() -> dict:
    global _local
    if _local is None:
        _local = {}
        fp = config.CLINVAR_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                row = dict(zip(_COLUMNS, parts))
                if row.get("key"):
                    _local[row["key"]] = row
    return _local


def _chrpos_field() -> str:
    return "chrpos37" if config.GENOME_BUILD == "GRCh37" else "chrpos38"


def _ncbi_params() -> dict:
    return {"db": "clinvar", "retmode": "json", "tool": "vcf2report",
            "email": config.NCBI_EMAIL, "api_key": config.NCBI_API_KEY}


def _loc_matches(variant: Variant, docsum: dict) -> bool:
    """True if the ClinVar record's location matches this variant's chr/pos(/alt)."""
    chrom = variant.chrom[3:] if variant.chrom.lower().startswith("chr") else variant.chrom
    for vset in docsum.get("variation_set", []) or []:
        for loc in vset.get("variation_loc", []) or []:
            if (loc.get("assembly_name") or "").upper() != config.GENOME_BUILD.upper():
                continue
            if str(loc.get("chr")) != str(chrom):
                continue
            if str(loc.get("start")) != str(variant.pos):
                continue
            alt = loc.get("alt")
            if alt and str(alt).upper() != variant.alt.upper():
                continue  # position matched but a different alt allele
            return True
    return False


def _extract(docsum: dict) -> dict:
    """Pull significance/review/condition from new or legacy esummary shapes."""
    germ = docsum.get("germline_classification") or {}
    legacy = docsum.get("clinical_significance") or {}
    significance = germ.get("description") or legacy.get("description")
    review = germ.get("review_status") or legacy.get("review_status")
    date = germ.get("last_evaluated") or legacy.get("last_evaluated")
    traits = germ.get("trait_set") or docsum.get("trait_set") or []
    condition = "; ".join(
        t.get("trait_name", "") for t in traits if t.get("trait_name")
    ) or None
    accession = docsum.get("accession") or docsum.get("obj_type")
    return {"significance": significance, "review_status": review,
            "accession": accession, "condition": condition, "date": date}


def _live(variant: Variant) -> Optional[dict]:
    """Best-effort ClinVar lookup via NCBI E-utilities (esearch -> esummary).

    Coordinate-driven and conservative: only returns a record whose location
    (assembly/chr/pos and, when present, alt) matches the query variant. Any
    ambiguity, transport error, or no-match returns None so the caller falls
    back to the authoritative local ClinVar slice.
    """
    from . import _http

    chrom = variant.chrom[3:] if variant.chrom.lower().startswith("chr") else variant.chrom
    # NCBI: 3 req/s anonymous, 10 req/s with an API key.
    interval = 0.11 if config.NCBI_API_KEY else 0.34

    _http.throttle("ncbi", interval)
    search = _http.get_json(
        f"{config.NCBI_EUTILS}/esearch.fcgi",
        {**_ncbi_params(), "retmax": 20,
         "term": f"{variant.pos}[{_chrpos_field()}] AND {chrom}[chr]"},
    )
    if not search:
        return None
    ids = (((search.get("esearchresult") or {}).get("idlist")) or [])
    if not ids:
        return {"significance": None, "review_status": None, "accession": None,
                "condition": None, "date": None}  # searched, genuinely absent

    _http.throttle("ncbi", interval)
    summary = _http.get_json(
        f"{config.NCBI_EUTILS}/esummary.fcgi",
        {**_ncbi_params(), "id": ",".join(ids)},
    )
    if not summary:
        return None
    result = summary.get("result") or {}
    for uid in result.get("uids", []) or []:
        docsum = result.get(uid) or {}
        if _loc_matches(variant, docsum):
            return _extract(docsum)
    return None  # candidates found but none matched exactly -> fall back


def lookup(variant: Variant) -> dict:
    cached = cache.get(_SOURCE, variant.key)
    if cached is not None:
        return {**cached, "_source": "ClinVar (cache)"}

    if not config.offline():
        live = _live(variant)
        if live is not None:
            cache.put(_SOURCE, variant.key, live)
            return {**live, "_source": "ClinVar (live E-utilities)"}

    local = _load_local().get(variant.key)
    if local is not None:
        date = local.get("date", "")
        return {**local, "_source": f"ClinVar slice ({date})"}

    return {"significance": None, "review_status": None, "accession": None,
            "condition": None, "date": None, "_source": "ClinVar (no record)"}
