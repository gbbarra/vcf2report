"""Ground-truth concordance: our ACMG engine vs expert-reviewed ClinVar.

Loads the frozen panel (``data/truth/clinvar_panel.json`` — real ClinVar SNVs
with real gnomAD frequencies; build it with ``scripts/build_truth_panel.py``,
see docs/GROUND_TRUTH_PANEL.md) and classifies every variant with **ClinVar
withheld from the engine** (no ``clinvar_significance`` → PP5 cannot fire), so
the comparison is honest, not circular — it measures whether our *orthogonal*
evidence (frequency, consequence/LoF, constraint) agrees with expert curation.

Asserts the safety-critical property and the deterministic-recall floors:

* **zero directional contradictions** — no ClinVar-Pathogenic called Benign/
  Likely Benign, and no ClinVar-Benign called Pathogenic/Likely Pathogenic.
  (VUS is a safe non-contradiction: honest uncertainty, not a wrong call.)
* **BA1-band benign recall** — a benign variant at gnomAD AF ≥ 5% must reach
  Benign/Likely Benign (BA1 is deterministic).

The panel is optional: if it hasn't been built on this machine the test skips,
so CI stays green while the concordance check is available to run locally.
"""
import json
from pathlib import Path

import pytest

from vcf2report.acmg.engine import classify
from vcf2report.annotate import extra
from vcf2report.config import AF_BA1
from vcf2report.models import Annotation, Variant

PANEL = Path(__file__).resolve().parent.parent / "data" / "truth" / "clinvar_panel.json"
PATHOGENIC = {"Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"}
BENIGN = {"Benign", "Benign/Likely_benign", "Likely_benign"}
PATH_TIERS = {"Pathogenic", "Likely Pathogenic"}
BENIGN_TIERS = {"Benign", "Likely Benign"}
MIN_PANEL = 40  # below this the metrics aren't meaningful → skip


def _load():
    if not PANEL.exists():
        pytest.skip(f"{PANEL} not built — run scripts/build_truth_panel.py (see "
                    "docs/GROUND_TRUTH_PANEL.md)")
    data = json.loads(PANEL.read_text())
    variants = data.get("variants", [])
    if len(variants) < MIN_PANEL:
        pytest.skip(f"panel has {len(variants)} < {MIN_PANEL} variants — rebuild it fuller")
    return variants


def _classify_without_clinvar(v: dict) -> str:
    """Classify a panel row with ClinVar withheld (PP5 off) — orthogonal evidence only."""
    variant = Variant(chrom=v["chrom"], pos=v["pos"], ref=v["ref"], alt=v["alt"],
                      gene=v["gene"], consequence=v["consequence"])
    # gene_lof_intolerant comes from the committed local constraint table (drives PVS1);
    # clinvar_significance is deliberately left None so PP5 never fires.
    con = extra.gene_constraint(v["gene"])
    ann = Annotation(gnomad_af=v.get("gnomad_af"), gnomad_faf95=v.get("gnomad_faf95"),
                     gene_lof_intolerant=con.get("lof_intolerant"))
    return classify(variant, ann).tier


def test_no_directional_contradiction():
    variants = _load()
    bad = []
    for v in variants:
        tier = _classify_without_clinvar(v)
        truth = v["clinvar_sig"]
        if truth in PATHOGENIC and tier in BENIGN_TIERS:
            bad.append((v["gene"], v["chrom"], v["pos"], truth, tier))
        if truth in BENIGN and tier in PATH_TIERS:
            bad.append((v["gene"], v["chrom"], v["pos"], truth, tier))
    assert not bad, f"directional contradictions (truth vs call): {bad}"


def test_ba1_band_benign_recall():
    """Every benign variant common enough for BA1 (AF ≥ 5%) must land Benign/LB."""
    variants = _load()
    common_benign = [v for v in variants
                     if v["clinvar_sig"] in BENIGN and (v.get("gnomad_af") or 0) >= AF_BA1]
    if not common_benign:
        pytest.skip("no benign variants at AF ≥ 5% in the panel")
    misses = [(v["gene"], v.get("gnomad_af"), _classify_without_clinvar(v))
              for v in common_benign
              if _classify_without_clinvar(v) not in BENIGN_TIERS]
    assert not misses, f"BA1-band benign not called benign: {misses}"


def test_report_concordance_summary(capsys):
    """Print the full confusion matrix (informational — always passes if loaded)."""
    variants = _load()
    from collections import Counter
    matrix = Counter()
    for v in variants:
        truth = "Pathogenic" if v["clinvar_sig"] in PATHOGENIC else "Benign"
        tier = _classify_without_clinvar(v)
        side = ("Pathogenic" if tier in PATH_TIERS else
                "Benign" if tier in BENIGN_TIERS else "VUS")
        matrix[(truth, side)] += 1
    with capsys.disabled():
        print("\nGround-truth concordance (ClinVar withheld / PP5 off):")
        for truth in ("Pathogenic", "Benign"):
            row = {s: matrix[(truth, s)] for s in ("Pathogenic", "Benign", "VUS")}
            print(f"  ClinVar {truth:11s} -> {row}")
