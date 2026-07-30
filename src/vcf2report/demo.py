"""Demo mode — run the guided flow end-to-end on the repo's own committed example VCFs
without the full Parquet stores, and make every resulting laudo say so unmissably.

Why this exists: the Stage-1 store gate is deliberately hard, because a laudo produced with
gnomAD absent has PM2/BA1/BS1 blind and over-calls. But the gate cannot tell a real patient
exome from `data/example/SYN-073.BBS2.annotated.vcf.gz` — a file committed to this repository
purely so the pipeline can be demonstrated. That left the demo cases unrunnable through the
guided flow on any machine without ~1 GB of stores (a fresh container, a phone-driven cloud
session), which is the one place a demo is most useful.

Two properties make the exemption safe, and both are enforced here rather than by convention:

1. **It can only ever apply to a committed example file.** ``--demo`` is not a gate override:
   pointed at anything outside ``data/example/`` it is REFUSED, not honoured. So the flag can
   never be the reason a real patient VCF was analysed with blind stores.
2. **The exemption is stamped, not silent.** Whenever demo mode is what allowed a run to
   proceed, the provenance travels with the report model into the conclusion, the Methods
   block, the JSON, and both renderers. A demo laudo cannot be mistaken for a real one, and
   the stamp is not suppressible from the CLI.

The failure mode is asymmetric by design: a genuine analysis mislabelled "demo" merely
under-claims, while a demo laudo passing as real is the dangerous direction — so the check
errs toward stamping.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from . import config, stores

# The one directory whose VCFs are demonstration fixtures. These files are committed to the
# repository; they are synthetic exomes with a known planted diagnosis (see
# docs/SYNTHETIC_CASES.md), never patient data.
DEMO_VCF_DIR = config.DATA_DIR / "example"

_VCF_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz")


def _is_vcf(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(s) for s in _VCF_SUFFIXES)


def is_demo_vcf(vcf: str | Path | None) -> bool:
    """True when ``vcf`` is one of the repository's committed example VCFs.

    Resolved before comparison so a relative path, a symlink, or ``../`` cannot smuggle a
    file in from outside ``data/example/``.
    """
    if not vcf:
        return False
    try:
        p = Path(vcf).resolve()
        demo_dir = DEMO_VCF_DIR.resolve()
    except OSError:
        return False
    if not _is_vcf(p):
        return False
    return p.parent == demo_dir


def requested() -> bool:
    """Whether demo mode was asked for via the environment (``VCF2REPORT_DEMO=1``).

    The environment is the channel that reaches every surface — CLI, headless script, MCP
    server — without each having to thread a flag down to the pipeline.
    """
    return os.environ.get("VCF2REPORT_DEMO", "").strip().lower() in {"1", "true", "yes", "on"}


@dataclass
class DemoDecision:
    """Outcome of asking "may this run proceed under the demo exemption?"."""

    active: bool = False           # demo mode is in force for this run
    refused: bool = False          # --demo was asked for but must NOT be honoured
    reason: str = ""               # human-readable, shown in the progress surface
    vcf: str = ""                  # the example file, repo-relative, for the stamp

    def to_dict(self) -> dict:
        return {"active": self.active, "refused": self.refused,
                "reason": self.reason, "vcf": self.vcf}


def decide(vcf: str | Path | None, demo_requested: bool | None = None) -> DemoDecision:
    """Resolve whether demo mode applies, WITHOUT consulting the stores.

    ``demo_requested`` defaults to the environment. A request for a file outside
    ``data/example/`` comes back ``refused`` — the caller must treat that as a hard stop, not
    fall through to a normal run, so an operator who believed they were in demo mode is told
    plainly that they are not.
    """
    if demo_requested is None:
        demo_requested = requested()
    if not demo_requested:
        return DemoDecision(reason="demo mode not requested")
    if not is_demo_vcf(vcf):
        return DemoDecision(
            refused=True,
            reason=(f"demo mode REFUSED — {vcf!r} is not one of this repository's committed "
                    f"example VCFs ({DEMO_VCF_DIR}). Demo mode is not a store override: it "
                    "exists only to demonstrate the pipeline on its own fixtures. Build the "
                    "Parquet stores to analyse this file."),
        )
    try:
        rel = str(Path(vcf).resolve().relative_to(config.REPO_ROOT.resolve()))
    except (ValueError, OSError):
        rel = str(vcf)
    return DemoDecision(
        active=True,
        reason=(f"demo mode ACTIVE — {rel} is a committed synthetic example with a known "
                "planted diagnosis. The store gate is relaxed for it and the laudo is "
                "stamped as a demonstration."),
        vcf=rel,
    )


def gate(vcf: str | Path | None = None, demo_requested: bool | None = None,
         measure: bool = True) -> dict:
    """The guided flow's Stage-1 store gate, with the demo exemption applied.

    Returns the fields of :func:`stores.gate` plus:

    ``mode``      ``"full"`` (all stores intact) | ``"demo"`` (exemption carried the run)
    ``demo``      the :class:`DemoDecision`, as a dict
    ``exempted``  stores that were blocking and were let through by the exemption

    ``ready`` is True under the exemption *even though stores block*, which is the whole
    point — but ``mode`` is then ``"demo"``, and the pipeline stamps the report from the same
    decision, so nothing downstream can quietly treat it as a full run.
    """
    g = stores.gate(measure=measure)
    d = decide(vcf, demo_requested)
    g["demo"] = d.to_dict()
    g["exempted"] = []
    g["mode"] = "full"

    if d.refused:
        g["ready"] = False
        g["mode"] = "refused"
        return g
    if not d.active or not g["blocking"]:
        # Nothing to exempt: either demo wasn't asked for, or every store is already intact.
        # A demo run against a fully-built machine still gets mode="demo" so the laudo is
        # stamped — the file is a fixture regardless of how healthy the stores are.
        if d.active:
            g["mode"] = "demo"
        return g

    g["exempted"] = list(g["blocking"])
    g["ready"] = True
    g["mode"] = "demo"
    return g


def provenance(vcf: str | Path | None = None, measure: bool = False) -> dict:
    """The stamp the report carries. Empty dict for an ordinary run on a real VCF.

    Deliberately does NOT consult ``VCF2REPORT_DEMO`` or any flag: a file in ``data/example/``
    is a demonstration fixture however it was invoked, so the stamp is unconditional. Making it
    flag-driven would mean a forgotten ``--demo`` — on the CLI, or on an MCP surface with no
    environment to set — produced an *unstamped* fixture laudo, which is exactly the dangerous
    direction. The flag governs only whether the store gate relaxes.

    Records which stores were actually absent, so the laudo names the specific criteria the
    reader must not trust rather than a vague "demo data" note. ``measure=False`` by default:
    presence is all the stamp needs, and a row-count scan would cost seconds per run.
    """
    if not is_demo_vcf(vcf):
        return {}
    try:
        rel = str(Path(vcf).resolve().relative_to(config.REPO_ROOT.resolve()))
    except (ValueError, OSError):
        rel = str(vcf)
    health = stores.store_health(measure=measure)
    absent = [n for n in stores.REQUIRED if health.get(n, {}).get("status") in stores._BLOCKING]
    return {
        "mode": "demo",
        "vcf": rel,
        "reason": (f"{rel} is a synthetic exome committed to this repository with a known "
                   "planted diagnosis (docs/SYNTHETIC_CASES.md) — a demonstration fixture, "
                   "not patient data."),
        "stores_absent": absent,
        "criteria_degraded": degraded_criteria(absent),
    }


# What each absent store costs, in criteria. Named explicitly so the demo stamp tells the
# reader which lines of the trail to distrust instead of leaving them to guess.
_STORE_COST = {
    "gnomad": ("PM2", "BA1", "BS1", "BS2"),
    "alphamissense": ("PP3", "BP4"),
    "clinvar": ("PS1", "PM1", "PM5", "PP5", "BP6"),
}


def degraded_criteria(absent) -> list[str]:
    out: list[str] = []
    for name in absent or ():
        for code in _STORE_COST.get(name, ()):
            if code not in out:
                out.append(code)
    return out


def stamp_line(prov: dict) -> str:
    """The one-sentence banner every renderer puts at the top of a demo laudo."""
    if not prov:
        return ""
    absent = prov.get("stores_absent") or []
    degraded = prov.get("criteria_degraded") or []
    tail = (
        f" The full Parquet store(s) {', '.join(absent)} are NOT present, so "
        f"{', '.join(degraded)} rest on the repository's frozen demo slices rather than the "
        "complete databases — treat every one of them as unverified."
        if absent else
        " Stores are present, but the input is still a fixture, not a patient."
    )
    return (
        "**DEMONSTRATION RUN — NOT A PATIENT RESULT.** This laudo was produced from "
        f"`{prov.get('vcf', 'a committed example VCF')}`, a synthetic exome committed to this "
        "repository with a known planted diagnosis." + tail
    )
