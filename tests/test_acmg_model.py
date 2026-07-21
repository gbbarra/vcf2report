"""Configurable ACMG combining model: Richards Table 5 vs ClinGen/Tavtigian points."""
import pytest

from vcf2report import config
from vcf2report.acmg import rules
from vcf2report.acmg.engine import classify
from vcf2report.acmg.rules import BENIGN, LIKELY_BENIGN, LIKELY_PATHOGENIC, PATHOGENIC, VUS
from vcf2report.models import Annotation, CriterionResult, Variant


def _cr(code, strength):
    return CriterionResult(code=code, name="", default_strength=strength, applies=True,
                           met=True, applied_strength=strength)


# ---------------------------------------------------------------------------
# Model toggle
# ---------------------------------------------------------------------------
def test_acmg_model_toggle(monkeypatch):
    # The Table-5-vs-points combiner and the PM2 strength are SEPARATE knobs now.
    monkeypatch.delenv("VCF2REPORT_PM2_STRENGTH", raising=False)
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "richards")
    assert config.acmg_model() == "richards"
    assert config.pm2_strength() == "supporting"   # default: ClinGen SVI 2020, decoupled from the combiner
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    assert config.acmg_model() == "clingen" and config.pm2_strength() == "supporting"
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "points")
    assert config.acmg_model() == "clingen"
    # PM2 strength override restores the legacy Richards-2015 Moderate.
    monkeypatch.setenv("VCF2REPORT_PM2_STRENGTH", "moderate")
    assert config.pm2_strength() == "moderate"


# ---------------------------------------------------------------------------
# Points math (Tavtigian 2020): VS8 S4 M2 P1; BA-8 BS-4 BP-1; P>=10 LP6-9 VUS0-5 LB-1..-6 B<=-7
# ---------------------------------------------------------------------------
def test_points_thresholds():
    P = rules._combine_points
    assert P([_cr("PVS1", "very_strong")])[0] == LIKELY_PATHOGENIC            # 8
    assert P([_cr("PVS1", "very_strong"), _cr("PM2", "moderate")])[0] == PATHOGENIC   # 10
    assert P([_cr("PM2", "supporting"), _cr("PP3", "strong")])[0] == VUS      # 5 (one short of LP)
    assert P([_cr("PM2", "supporting"), _cr("PP3", "strong"),
              _cr("PP4", "supporting")])[0] == LIKELY_PATHOGENIC              # 6
    assert P([_cr("BA1", "stand_alone")])[0] == BENIGN                        # -8
    assert P([_cr("BS1", "strong")])[0] == LIKELY_BENIGN                      # -4
    assert P([_cr("BS1", "strong"), _cr("BS2", "strong")])[0] == BENIGN       # -8
    # conflicting nets out
    assert P([_cr("PVS1", "very_strong"), _cr("BS1", "strong")])[0] == VUS    # +4


# ---------------------------------------------------------------------------
# End-to-end under the ClinGen model
# ---------------------------------------------------------------------------
def _ann(**kw):
    base = dict(gene_lof_intolerant=False, abraom_af=0.0, source={})
    base.update(kw)
    return Annotation(**base)


def test_clingen_lof_stays_lp(monkeypatch):
    # PVS1(8) alone is already LP under points, so PM2->Supporting does not hurt LoF.
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="G", consequence="stop_gained")
    c = classify(v, _ann(gene_lof_intolerant=True, gnomad_af=0.0, gnomad_faf95=0.0))
    assert "PM2" in c.met_codes and c.tier == LIKELY_PATHOGENIC


def test_clingen_missense_strong_am_is_vus(monkeypatch):
    # PM2_Supporting(1) + PP3_Strong(4) = 5 points -> VUS (the missense-recovery loss).
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="G", consequence="missense_variant")
    c = classify(v, _ann(gnomad_af=0.0, gnomad_faf95=0.0, am_pathogenicity=0.999))
    assert c.tier == VUS


def test_richards_missense_pm2_strength(monkeypatch):
    # Richards Table-5 combiner, rare missense + strong AlphaMissense = PM2 + PP3_Strong.
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "richards")
    v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="G", consequence="missense_variant")
    a = _ann(gnomad_af=0.0, gnomad_faf95=0.0, am_pathogenicity=0.999)
    # Default PM2 Supporting: one point short of LP -> held at VUS (the probable-pathogenic VUS).
    monkeypatch.delenv("VCF2REPORT_PM2_STRENGTH", raising=False)
    assert classify(v, a).tier == VUS
    # Legacy PM2 Moderate override -> LP-2.
    monkeypatch.setenv("VCF2REPORT_PM2_STRENGTH", "moderate")
    assert classify(v, a).tier == LIKELY_PATHOGENIC


# ---------------------------------------------------------------------------
# Safety invariant holds under BOTH models on the frozen panel
# ---------------------------------------------------------------------------
def test_panel_zero_gross_both_models(monkeypatch):
    from vcf2report import concordance
    if not concordance.FROZEN_GNOMAD.exists() or not concordance.GROUND_TRUTH.exists():
        pytest.skip("panel not frozen")
    entries = concordance.load_panel()
    if not entries:
        pytest.skip("empty panel")
    for model in ("richards", "clingen"):
        monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", model)
        res = concordance.evaluate_panel(entries)
        assert res.metrics["gross_discordances"] == 0, f"{model}: {res.to_markdown()}"
        assert res.metrics["benign_precision"] >= 0.99
