"""Headless command-line entry point.

Runs the full VCF -> report pipeline without Claude, for fast iteration and CI.
This is the same pipeline the MCP server exposes to Claude Desktop.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import config
from .pipeline import run_pipeline
from .report.render import render_markdown, write_report


def read_hpo_file(path: str | Path) -> list[str]:
    """Resolve patient HPO terms from either a file (one HP:id per line) OR an inline
    list of HP: ids (comma/space/semicolon-separated). Inline is accepted because
    ``--hpo HP:0001250,HP:0001263`` is a natural thing to type; treating it strictly as
    a path would silently drop the terms and disable PP4 — dangerous for a clinical tool.
    """
    p = Path(path)
    if p.exists():
        terms: list[str] = []
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            terms.append(line.split()[0])  # first token = HP:xxxxxxx
        return terms
    # Not a file -> parse HP: ids straight out of the argument (e.g. "HP:1,HP:2").
    import re
    return [t.upper() for t in re.findall(r"HP:\d+", str(path), flags=re.IGNORECASE)]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="vcf2report",
        description="Turn an exome VCF into an auditable ACMG variant report.",
    )
    parser.add_argument("vcf", nargs="?", default=str(config.SAMPLE_VCF),
                        help="Path to the input VCF (defaults to the bundled sample).")
    parser.add_argument("--hpo", default=None,
                        help="Patient HPO terms: a file (one HP:id per line) OR an inline "
                             "list like HP:0001250,HP:0001263 (defaults to the sample).")
    parser.add_argument("--sample-id", default=None, help="Override the sample identifier.")
    parser.add_argument("--out", default=None, help="Output directory for the report.")
    parser.add_argument("--stdout", action="store_true",
                        help="Print the Markdown report to stdout instead of writing a file.")
    parser.add_argument("--timing", action="store_true",
                        help="Print per-stage timings (parse/QC/annotate/filter/classify).")
    args = parser.parse_args(argv)

    hpo_path = args.hpo
    if hpo_path is None and args.vcf == str(config.SAMPLE_VCF):
        hpo_path = str(config.SAMPLE_HPO)
    hpo_terms = read_hpo_file(hpo_path) if hpo_path else []

    report = run_pipeline(args.vcf, hpo_terms=hpo_terms, sample_id=args.sample_id)

    def _print_timing() -> None:
        if not report.timings:
            return
        print("  timing (s):")
        for k, v in report.timings.items():
            print(f"    {k}: {v}")

    if args.stdout:
        sys.stdout.write(render_markdown(report))
        if args.timing:
            _print_timing()
        return 0

    out_dir = Path(args.out) if args.out else config.OUTPUT_DIR
    fp = write_report(report, out_dir)
    from .report.explore import write_explore
    jp = write_explore(report, str(fp).replace("_report.md", "_results.json"))
    print(f"Report written to {fp}")
    print(f"  explorable data: {jp}")
    print(f"  candidates classified: {report.qc.candidates}")
    for c in report.classifications:
        print(f"  - {c.variant.gene}: {c.tier}")
    if args.timing:
        _print_timing()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
