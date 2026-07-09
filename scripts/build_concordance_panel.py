#!/usr/bin/env python3
"""Build (freeze) the ClinVar-vs-engine concordance panel — run ONCE with network.

Harvests real, expert-classified variants from ClinVar for a curated gene panel
(coordinates always from ClinVar — never fabricated), then freezes each variant's
gnomAD grpmax frequency via remote tabix. Writes two files under
``data/concordance/``:

    ground_truth.tsv    key  gene  consequence  hgvs_p  clinvar_sig  review  acc  condition
    gnomad_frozen.json  { key: {af, ac, an, hom, pop, faf95, release} }

After this runs once (~1-2 min), ``tests/test_concordance_panel.py`` and
``scripts/run_concordance.py`` execute fully offline against the frozen files.

Idempotent / resumable: existing entries are kept; only missing gnomAD records
are fetched. Re-run any time to top up the panel.

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_concordance_panel.py \
        [--per-group 50] [--genes SCN1A,KCNQ2,...]

Network egress here is to PUBLIC ClinVar / gnomAD records only (no patient data),
so it is safe to run; the pipeline's patient-facing egress gate is separate.
"""
from __future__ import annotations

import argparse
import json
import signal
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import config  # noqa: E402
from vcf2report.annotate import _http, clinvar, gnomad, gnomad_remote  # noqa: E402
from vcf2report.models import Variant  # noqa: E402

# Curated gene panel: dominant LoF-driven disease genes + ACMG SF genes that also
# carry well-classified common benign variation — a mix that exercises both
# deterministic axes (LoF pathogenicity and common-benign frequency).
DEFAULT_GENES = [
    "SCN1A", "SCN2A", "KCNQ2", "STXBP1", "CACNA1A", "SLC2A1", "PAX6",
    "RB1", "APC", "MUTYH", "TP53", "STK11", "PTEN", "VHL", "NF1",
    "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2",
    "FBN1", "LDLR", "MYH7", "MYBPC3", "KCNQ1", "KCNH2", "SCN5A",
    "CFTR", "HFE", "GJB2", "ATP7B", "GAA",
]

# ClinVar molecular-consequence text -> our consequence vocabulary.
_CONSEQUENCE_MAP = {
    "nonsense": "stop_gained",
    "frameshift variant": "frameshift_variant",
    "frameshift": "frameshift_variant",
    "splice donor variant": "splice_donor_variant",
    "splice acceptor variant": "splice_acceptor_variant",
    "missense variant": "missense_variant",
    "synonymous variant": "synonymous_variant",
    "stop lost": "stop_lost",
    "initiator codon variant": "start_lost",
    "start lost": "start_lost",
    "inframe deletion": "inframe_deletion",
    "inframe insertion": "inframe_insertion",
}

# When ClinVar lists several consequences for one variant, keep the most severe
# so a co-listed milder term (e.g. "5 prime UTR variant") never masks the LoF
# axis the panel reports on. Most severe first.
_CONSEQUENCE_SEVERITY = [
    "splice_donor_variant", "splice_acceptor_variant", "stop_gained",
    "frameshift_variant", "start_lost", "stop_lost",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "synonymous_variant",
]

_GT_HEADER = ("# Concordance-panel ground truth (real ClinVar, GRCh38). "
              "Columns: key gene consequence hgvs_p clinvar_significance "
              "review_status accession condition")


def _online() -> bool:
    return config.allow_network()


def _interval() -> float:
    # NCBI: 3 req/s anonymous, 10 req/s with an API key.
    return 0.11 if config.NCBI_API_KEY else 0.34


def _sig_term(group: str) -> str:
    if group == "pathogenic":
        return ('("pathogenic"[Clinical significance] OR '
                '"likely pathogenic"[Clinical significance])')
    return ('("benign"[Clinical significance] OR '
            '"likely benign"[Clinical significance])')


def _reviewed(review: str | None) -> bool:
    """>=1-star ClinVar assertion (same gate PP5 uses)."""
    r = (review or "").lower().replace("_", " ").strip()
    return (r.startswith("criteria provided")
            or "reviewed by expert" in r
            or "practice guideline" in r)


def _group_of(significance: str | None) -> str | None:
    """Collapse a ClinVar significance to 'pathogenic' / 'benign' / None."""
    sig = (significance or "").strip().lower()
    if sig.startswith("pathogenic") or sig.startswith("likely pathogenic"):
        return "pathogenic"
    if sig.startswith("benign") or sig.startswith("likely benign"):
        return "benign"
    return None


def _consequence_from(docsum: dict) -> str:
    # Real ClinVar esummary carries the SO terms under molecular_consequence_list
    # (a list); the singular key is kept only as a defensive fallback.
    mc = docsum.get("molecular_consequence_list") or docsum.get("molecular_consequence")
    if isinstance(mc, list):
        mapped = {_CONSEQUENCE_MAP.get(str(term).strip().lower()) for term in mc}
        mapped.discard(None)
        for cons in _CONSEQUENCE_SEVERITY:
            if cons in mapped:
                return cons
    # Fall back to the protein change encoded in the record title.
    title = (docsum.get("title") or "").lower()
    if "ter)" in title or "*)" in title or "=stop" in title:
        return "stop_gained"
    if "fs" in title and "p." in title:
        return "frameshift_variant"
    if "p.(=)" in title or "p.=" in title:
        return "synonymous_variant"
    return ""


def _hgvs_p_from(docsum: dict) -> str:
    title = docsum.get("title") or ""
    if "(p." in title and title.endswith(")"):
        return "p." + title.rsplit("(p.", 1)[1][:-1]
    return ""


def _grch38_snv(variant_set: list) -> tuple[str, str, str] | None:
    """Return (key, ref, alt) for a simple GRCh38 SNV, else None.

    v1 restricts the panel to single-nucleotide variants so the frozen gnomAD
    tabix match is unambiguous (no left-alignment / representation mismatch).
    """
    for vset in variant_set or []:
        spdi = (vset.get("canonical_spdi") or "").split(":")
        spdi_ref = spdi[-2] if len(spdi) >= 4 else None
        spdi_alt = spdi[-1] if len(spdi) >= 4 else None
        for loc in vset.get("variation_loc", []) or []:
            if (loc.get("assembly_name") or "").upper() != "GRCH38":
                continue
            chrom = str(loc.get("chr") or "").replace("chr", "")
            start = loc.get("start")
            ref = (loc.get("ref") or spdi_ref or "").upper()
            alt = (loc.get("alt") or spdi_alt or "").upper()
            if not (chrom and start and ref and alt):
                continue
            if len(ref) != 1 or len(alt) != 1 or ref == alt:
                continue  # SNV only for v1
            if ref not in "ACGT" or alt not in "ACGT":
                continue
            return f"{chrom}-{start}-{ref}-{alt}", ref, alt
    return None


def _esearch_ids(params: dict, term: str, retmax: int) -> list[str]:
    _http.throttle("ncbi", _interval())
    search = _http.get_json(
        f"{config.NCBI_EUTILS}/esearch.fcgi",
        {**params, "retmax": retmax, "term": term})
    return (((search or {}).get("esearchresult") or {}).get("idlist")) or []


def _harvest_gene(gene: str, group: str, want: int, seen: set[str]) -> list[dict]:
    """Fetch up to ``want`` reviewed GRCh38 SNVs of ``group`` for ``gene``.

    The engine only accepts a variant whose *extracted* significance actually
    collapses to the requested group and that carries a >=1-star review — so a
    mis-tagged search hit can never land in the wrong bucket. If the field-qualified
    significance term returns nothing (Entrez field quirks), it retries with a plain
    keyword term and leans entirely on that post-filter.
    """
    params = {"db": "clinvar", "retmode": "json", "tool": "vcf2report",
              "email": config.NCBI_EMAIL, "api_key": config.NCBI_API_KEY}
    retmax = max(want * 6, 30)
    ids = _esearch_ids(params, f"{gene}[gene] AND {_sig_term(group)}", retmax)
    if not ids:
        ids = _esearch_ids(params, f"{gene}[gene] AND {group}", retmax)
    if not ids:
        return []

    _http.throttle("ncbi", _interval())
    summary = _http.get_json(
        f"{config.NCBI_EUTILS}/esummary.fcgi", {**params, "id": ",".join(ids)})
    result = (summary or {}).get("result") or {}

    rows: list[dict] = []
    for uid in result.get("uids", []) or []:
        if len(rows) >= want:
            break
        docsum = result.get(uid) or {}
        extracted = clinvar._extract(docsum)
        if _group_of(extracted.get("significance")) != group:
            continue  # wrong significance bucket -> never mis-file it
        if not _reviewed(extracted.get("review_status")):
            continue
        snv = _grch38_snv(docsum.get("variation_set", []))
        if not snv:
            continue
        key, _ref, _alt = snv
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "key": key,
            "gene": gene,
            "consequence": _consequence_from(docsum),
            "hgvs_p": _hgvs_p_from(docsum),
            "clinvar_significance": extracted["significance"],
            "review_status": extracted.get("review_status") or "",
            "accession": extracted.get("accession") or "",
            "condition": (extracted.get("condition") or "").replace("\t", " "),
        })
    return rows


class _FreezeTimeout(Exception):
    pass


def _call_bounded(seconds: int, fn, *args):
    """Run ``fn(*args)`` with a hard SIGALRM deadline (main thread only).

    The remote-tabix path opens gnomAD VCFs over HTTPS via htslib, which carries
    no read timeout — a slow/stalled GCS read would otherwise hang the whole
    freeze. On the deadline we return None so the caller falls back cleanly.
    """
    def _handler(signum, frame):
        raise _FreezeTimeout()

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        return fn(*args)
    except _FreezeTimeout:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _freeze_gnomad(key: str) -> dict | None:
    chrom, pos, ref, alt = key.split("-")
    v = Variant(chrom=chrom, pos=int(pos), ref=ref, alt=alt)

    # GraphQL first: urllib with a 15s timeout + retry/backoff (annotate._http),
    # so it never hangs and never raises. A genuine "variant not found" returns a
    # real absent record (af 0.0), which is correct for a rare pathogenic.
    _http.throttle("gnomad", 0.5)
    live = gnomad._live(v)
    if live is not None and live.get("af") is not None:
        live = dict(live)
        live.setdefault("faf95", None)
        live["release"] = config.GNOMAD_DATASET
        return live

    # Fallback: remote tabix (adds the filtering AF, faf95), bounded so a stalled
    # GCS read can't hang the build.
    remote = _call_bounded(20, gnomad_remote.query, v)
    if remote is not None:
        remote = dict(remote)
        remote["release"] = gnomad_remote.RELEASE
        return remote
    return None


def _load_existing_truth(fp: Path) -> tuple[list[dict], set[str]]:
    rows: list[dict] = []
    seen: set[str] = set()
    if not fp.exists():
        return rows, seen
    cols = ["key", "gene", "consequence", "hgvs_p",
            "clinvar_significance", "review_status", "accession", "condition"]
    for line in fp.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        row = dict(zip(cols, line.split("\t")))
        if row.get("key"):
            rows.append(row)
            seen.add(row["key"])
    return rows, seen


def _write_truth(fp: Path, rows: list[dict]) -> None:
    cols = ["key", "gene", "consequence", "hgvs_p",
            "clinvar_significance", "review_status", "accession", "condition"]
    lines = [_GT_HEADER]
    for r in sorted(rows, key=lambda x: x["key"]):
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    fp.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Freeze the ClinVar-vs-engine concordance panel.")
    ap.add_argument("--per-group", type=int, default=50,
                    help="target pathogenic AND benign variants (each) across the panel")
    ap.add_argument("--per-gene", type=int, default=4,
                    help="max variants to take from one gene per significance group")
    ap.add_argument("--genes", default="",
                    help="comma-separated gene override (default: curated panel)")
    ap.add_argument("--refresh", action="store_true",
                    help="re-fetch gnomAD even for already-frozen variants (apply client fixes)")
    args = ap.parse_args()

    if not _online():
        print("ERROR: network is disabled. Re-run with "
              "VCF2REPORT_ALLOW_NETWORK=1 (and OFFLINE unset).", file=sys.stderr)
        return 2

    genes = [g.strip().upper() for g in args.genes.split(",") if g.strip()] or DEFAULT_GENES
    out_dir = config.DATA_DIR / "concordance"
    out_dir.mkdir(parents=True, exist_ok=True)
    truth_fp = out_dir / "ground_truth.tsv"
    frozen_fp = out_dir / "gnomad_frozen.json"

    truth_rows, seen = _load_existing_truth(truth_fp)
    frozen: dict = json.loads(frozen_fp.read_text()) if frozen_fp.exists() else {}
    counts = {"pathogenic": 0, "benign": 0}
    for r in truth_rows:
        sig = (r.get("clinvar_significance") or "").lower()
        if sig.startswith("pathogenic") or sig.startswith("likely pathogenic"):
            counts["pathogenic"] += 1
        elif sig.startswith("benign") or sig.startswith("likely benign"):
            counts["benign"] += 1

    t0 = time.perf_counter()
    print(f"Harvesting ClinVar across {len(genes)} genes "
          f"(target {args.per_group}/group; have "
          f"{counts['pathogenic']} path / {counts['benign']} benign)...")

    for group in ("pathogenic", "benign"):
        for gene in genes:
            if counts[group] >= args.per_group:
                break
            want = min(args.per_gene, args.per_group - counts[group])
            try:
                rows = _harvest_gene(gene, group, want, seen)
            except Exception as exc:  # noqa: BLE001 — one bad gene must not abort the build
                print(f"  ! {gene} ({group}): harvest failed: {exc}", file=sys.stderr)
                continue
            for row in rows:
                truth_rows.append(row)
                counts[group] += 1
            if rows:
                print(f"  {gene} ({group}): +{len(rows)} "
                      f"[{counts['pathogenic']}p/{counts['benign']}b]")

    _write_truth(truth_fp, truth_rows)

    # Freeze gnomAD for every key not already frozen (or all, with --refresh).
    to_freeze = [r["key"] for r in truth_rows if args.refresh or r["key"] not in frozen]
    print(f"Freezing gnomAD for {len(to_freeze)} new variants "
          f"({len(frozen)} already cached)...")
    for i, key in enumerate(to_freeze, 1):
        try:
            rec = _freeze_gnomad(key)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! {key}: gnomAD lookup failed: {exc}", file=sys.stderr)
            rec = None
        if rec is not None:
            frozen[key] = {k: rec.get(k) for k in
                           ("af", "ac", "an", "hom", "pop", "faf95", "release")}
        if i % 10 == 0 or i == len(to_freeze):
            frozen_fp.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
            print(f"  frozen {i}/{len(to_freeze)}")

    frozen_fp.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
    dt = time.perf_counter() - t0
    print(f"\nDone in {dt:.0f}s: {len(truth_rows)} truth variants "
          f"({counts['pathogenic']} path / {counts['benign']} benign), "
          f"{len(frozen)} gnomAD records frozen.")
    print(f"  {truth_fp}")
    print(f"  {frozen_fp}")
    print("Now run:  python scripts/run_concordance.py   (offline)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
