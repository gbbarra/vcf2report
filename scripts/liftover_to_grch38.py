#!/usr/bin/env python3
"""Lift a GRCh37/hg19 VCF over to GRCh38/hg38 so the pipeline's GRCh38-only
coordinate lookups (gnomAD r4, ClinVar) can run.

Coordinate conversion uses ``pyliftover`` (a pure-Python UCSC liftOver port).
The output carries build provenance in its header (``##reference=GRCh38`` +
``##liftover=...``) so :func:`vcf2report.vcf.parse.detect_build` reports GRCh38.

    python3 scripts/liftover_to_grch38.py IN.vcf[.gz] OUT.vcf [--chain PATH]

Chain resolution order:
  1. ``--chain PATH`` if given;
  2. else ``<repo>/data/liftover/hg19ToHg38.over.chain.gz``;
  3. else, only when ``VCF2REPORT_ALLOW_NETWORK`` is truthy (1/true/yes),
     download that UCSC chain into the default path.

Notes / caveats:
  * Output is written UNSORTED — liftOver can reorder coordinates. Run
    ``bcftools sort`` (and tabix) afterwards if you need a sorted/indexed file;
    vcf2report's own reader does not require sorting.
  * ``pyliftover`` is 0-based; UCSC chains name contigs ``chr1``..``chrX/Y/M``.
    We query ``chr<core>`` at ``POS-1`` and emit ``tpos0+1``.
  * On a '-' strand hit only SNVs are safe (reverse-complemented); indels,
    MNVs and symbolic ALTs on '-' are dropped and counted (skipped_strand).
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CHAIN = REPO_ROOT / "data" / "liftover" / "hg19ToHg38.over.chain.gz"
CHAIN_URL = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"
             "hg19ToHg38.over.chain.gz")

_ACGT = frozenset("ACGT")
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

_PROVENANCE = (
    "##reference=GRCh38\n"
    "##liftover=hg19ToHg38 via pyliftover\n"
)


def _truthy(val: str | None) -> bool:
    return (val or "").strip().lower() in {"1", "true", "yes"}


def _open_text(path: Path):
    # gzip by extension; both paths yield universal-newline text streams.
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def resolve_chain(chain_arg: str | None) -> Path | None:
    """Return a usable chain-file path, or None (after printing why) if none."""
    if chain_arg:
        p = Path(chain_arg).expanduser()
        if not p.exists():
            print(f"ERROR: --chain path does not exist: {p}", file=sys.stderr)
            return None
        return p

    if DEFAULT_CHAIN.exists():
        return DEFAULT_CHAIN

    if _truthy(os.environ.get("VCF2REPORT_ALLOW_NETWORK")):
        print(f"Chain absent; downloading {CHAIN_URL}\n  -> {DEFAULT_CHAIN}",
              file=sys.stderr)
        try:
            DEFAULT_CHAIN.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(CHAIN_URL) as resp, \
                    open(DEFAULT_CHAIN, "wb") as out:
                shutil.copyfileobj(resp, out)
        except Exception as exc:  # network / IO — leave no partial file behind
            if DEFAULT_CHAIN.exists():
                try:
                    DEFAULT_CHAIN.unlink()
                except OSError:
                    pass
            print(f"ERROR: chain download failed: {exc}", file=sys.stderr)
            return None
        return DEFAULT_CHAIN

    print(
        f"ERROR: no liftover chain found at {DEFAULT_CHAIN} and no --chain given.\n"
        f"  Provide one with --chain PATH, or set VCF2REPORT_ALLOW_NETWORK=1 to\n"
        f"  auto-download it from UCSC, or fetch it manually:\n"
        f"    curl -L -o {DEFAULT_CHAIN} {CHAIN_URL}",
        file=sys.stderr,
    )
    return None


def _match_style(tchrom: str, has_chr: bool) -> str:
    """Render ``tchrom`` (always 'chrN' from the chain) in the input's style."""
    core = tchrom[3:] if tchrom[:3].lower() == "chr" else tchrom
    return f"chr{core}" if has_chr else core


def _lift_record(fields: list[str], lo) -> list[str] | str | None:
    """Lift one split data record.

    Returns the new field list on success, or a sentinel string reason
    ('unmapped' / 'strand' / 'malformed') on drop.
    """
    if len(fields) < 8:
        return "malformed"

    raw_chrom = fields[0]
    try:
        pos = int(fields[1])
    except ValueError:
        return "malformed"

    has_chr = raw_chrom[:3].lower() == "chr"
    core = raw_chrom[3:] if has_chr else raw_chrom
    query_chrom = "chr" + core

    # pyliftover: 0-based query, returns None or a list of hit tuples.
    hits = lo.convert_coordinate(query_chrom, pos - 1)
    if not hits:
        return "unmapped"

    tchrom, tpos0, strand, _ = hits[0]
    ref, alt = fields[3], fields[4]

    if strand == "-":
        is_snv = (len(ref) == 1 and len(alt) == 1
                  and ref.upper() in _ACGT and alt.upper() in _ACGT)
        if not is_snv:
            return "strand"  # indel/MNV/symbolic on '-' is not safe to flip
        ref = _COMP[ref.upper()]
        alt = _COMP[alt.upper()]
    # strand == "+" (or anything else): REF/ALT unchanged.

    out = list(fields)
    out[0] = _match_style(tchrom, has_chr)
    out[1] = str(tpos0 + 1)  # 0-based -> 1-based
    out[3] = ref
    out[4] = alt
    return out


def _write_header_line(line: str, fout, state: dict) -> None:
    """Emit a header line, injecting provenance once near the top (before #CHROM)."""
    if state["injected"]:
        fout.write(line)
        return
    # Keep ##fileformat as the very first line (VCF spec); provenance right after.
    if state["first"] and line.lower().startswith("##fileformat"):
        fout.write(line)
        fout.write(_PROVENANCE)
    else:
        fout.write(_PROVENANCE)
        fout.write(line)
    state["injected"] = True
    state["first"] = False


def liftover(in_path: Path, out_path: Path, lo) -> dict:
    counts = {"total": 0, "lifted": 0, "unmapped": 0, "strand": 0, "malformed": 0}
    state = {"first": True, "injected": False}

    with _open_text(in_path) as fin, open(out_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                _write_header_line(line, fout, state)
                state["first"] = False
                continue

            # Data record.
            counts["total"] += 1
            row = line.rstrip("\r\n")
            if not row:
                counts["malformed"] += 1
                continue
            fields = row.split("\t")

            try:
                result = _lift_record(fields, lo)
            except Exception:
                # Never crash on a single bad record — drop and keep going.
                counts["malformed"] += 1
                continue

            if result == "malformed":
                counts["malformed"] += 1
            elif result == "unmapped":
                counts["unmapped"] += 1
            elif result == "strand":
                counts["strand"] += 1
            else:
                fout.write("\t".join(result) + "\n")
                counts["lifted"] += 1

        # Degenerate input with no #CHROM line: still stamp provenance.
        if not state["injected"]:
            fout.write(_PROVENANCE)

    return counts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Lift a GRCh37/hg19 VCF over to GRCh38/hg38 (pyliftover).")
    ap.add_argument("input", help="input VCF (plain or .gz)")
    ap.add_argument("output", help="output VCF (plain text)")
    ap.add_argument("--chain", help="UCSC hg19->hg38 chain (.over.chain[.gz])")
    args = ap.parse_args(argv)

    try:
        from pyliftover import LiftOver
    except ImportError:
        print("ERROR: pyliftover is not installed. Run: pip install pyliftover",
              file=sys.stderr)
        return 2

    in_path = Path(args.input).expanduser()
    if not in_path.exists():
        print(f"ERROR: input VCF not found: {in_path}", file=sys.stderr)
        return 2

    chain = resolve_chain(args.chain)
    if chain is None:
        return 2

    try:
        lo = LiftOver(str(chain))
    except Exception as exc:
        print(f"ERROR: could not load chain {chain}: {exc}", file=sys.stderr)
        return 2

    out_path = Path(args.output).expanduser()
    counts = liftover(in_path, out_path, lo)

    print(
        f"[liftover] records={counts['total']} lifted={counts['lifted']} "
        f"skipped_unmapped={counts['unmapped']} skipped_strand={counts['strand']} "
        f"malformed={counts['malformed']}",
        file=sys.stderr,
    )
    print(
        "[liftover] WARNING: output is UNSORTED (liftOver can reorder positions). "
        "If you need a sorted/tabix'd file, run e.g. "
        "`bcftools sort -o sorted.vcf OUT.vcf`. vcf2report's own reader does not "
        "require sorting.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
