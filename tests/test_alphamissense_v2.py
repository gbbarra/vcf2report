"""v2: AlphaMissense-driven, ClinGen-calibrated PP3/BP4 + the client.

Locks the calibrated evidence-strength mapping, the missense-recovery it enables
(PM2 + PP3_Strong -> Likely Pathogenic), the backward-compatible REVEL/CADD
fallback, and — critically — that a *miscalibrated* high AlphaMissense score on a
rare benign variant is turned into a detectable gross discordance rather than
silently shipped.
"""
import pytest

from vcf2report import config
from vcf2report.acmg import criteria
from vcf2report.acmg.engine import classify
from vcf2report.acmg.rules import LIKELY_PATHOGENIC, VUS
from vcf2report.annotate import alphamissense
from vcf2report.models import Annotation, Variant


def _ann(**kw) -> Annotation:
    base = dict(gene_lof_intolerant=False, abraom_af=0.0,
                source={"alphamissense": "test", "insilico": "test"})
    base.update(kw)
    return Annotation(**base)


# ---------------------------------------------------------------------------
# config: calibrated strength mapping
# ---------------------------------------------------------------------------
def test_am_pp3_strength_bands():
    assert config.am_pp3_strength(0.999) == "strong"
    assert config.am_pp3_strength(config.AM_PP3_STRONG) == "strong"
    assert config.am_pp3_strength(0.95) == "moderate"
    assert config.am_pp3_strength(0.60) == "supporting"
    assert config.am_pp3_strength(0.45) is None       # ambiguous band
    assert config.am_pp3_strength(None) is None


def test_am_bp4_strength_bands():
    assert config.am_bp4_strength(0.10) == "supporting"
    assert config.am_bp4_strength(config.AM_BP4_SUPPORTING) == "supporting"
    assert config.am_bp4_strength(0.50) is None       # not benign enough
    assert config.am_bp4_strength(None) is None


# ---------------------------------------------------------------------------
# criteria: PP3/BP4 read AlphaMissense at variable strength
# ---------------------------------------------------------------------------
def _v(cons="missense_variant", gene="TESTG"):
    return Variant(chrom="1", pos=100, ref="A", alt="T", gene=gene, consequence=cons)


def test_pp3_uses_alphamissense_strength():
    r = criteria.pp3(_v(), _ann(am_pathogenicity=0.995))
    assert r.met and r.applied_strength == "strong"
    assert "AlphaMissense" in r.reasoning
    r2 = criteria.pp3(_v(), _ann(am_pathogenicity=0.93))
    assert r2.met and r2.applied_strength == "moderate"
    r3 = criteria.pp3(_v(), _ann(am_pathogenicity=0.20))
    assert not r3.met            # benign-leaning score -> PP3 does not fire


def test_bp4_uses_alphamissense():
    r = criteria.bp4(_v(), _ann(am_pathogenicity=0.05))
    assert r.met and r.applied_strength == "supporting"
    r2 = criteria.bp4(_v(), _ann(am_pathogenicity=0.95))
    assert not r2.met            # pathogenic score -> BP4 does not fire


def test_pp3_falls_back_to_revel_when_no_alphamissense():
    # am_pathogenicity None -> REVEL/CADD path (supporting), unchanged v1 behaviour.
    r = criteria.pp3(_v(), _ann(am_pathogenicity=None, revel=0.9))
    assert r.met and r.applied_strength == "supporting"
    assert "REVEL" in r.reasoning


@pytest.mark.parametrize("am", [0.0, 0.20, 0.34, 0.45, 0.564, 0.70, 0.90, 0.99, 1.0])
def test_pp3_and_bp4_never_both_fire(am):
    """A single AlphaMissense score can support at most one direction."""
    a = _ann(am_pathogenicity=am)
    assert not (criteria.pp3(_v(), a).met and criteria.bp4(_v(), a).met)


def test_ambiguous_band_fires_neither():
    a = _ann(am_pathogenicity=0.45)  # 0.34 < am < 0.564 -> AlphaMissense 'ambiguous'
    assert not criteria.pp3(_v(), a).met
    assert not criteria.bp4(_v(), a).met


# ---------------------------------------------------------------------------
# The recovery: a rare missense with a strong AlphaMissense score reaches LP
# ---------------------------------------------------------------------------
def test_rare_missense_strong_am_triaged_not_lp(monkeypatch):
    v = _v()  # missense, non-LoF gene
    a = _ann(gnomad_af=0.0, gnomad_faf95=0.0, am_pathogenicity=0.999)
    c = classify(v, a)
    assert "PM2" in c.met_codes and "PP3" in c.met_codes
    # Default (PM2 Supporting, ClinGen SVI 2020): PP3_Strong + PM2_Supporting is one point short of
    # LP -> VUS. Surfaced by the probable-pathogenic-VUS triage, not called LP on in-silico alone.
    assert c.tier == VUS
    # Legacy Richards-2015 Moderate override still recovers it to LP-2.
    monkeypatch.setenv("VCF2REPORT_PM2_STRENGTH", "moderate")
    assert classify(v, a).tier == LIKELY_PATHOGENIC


def test_rare_missense_moderate_am_stays_vus():
    # PP3_Moderate + PM2 = two Moderate -> not enough for LP under Richards.
    c = classify(_v(), _ann(gnomad_af=0.0, gnomad_faf95=0.0, am_pathogenicity=0.93))
    assert c.tier == VUS


# ---------------------------------------------------------------------------
# Adversarial: a miscalibrated high AM on a rare *benign* variant would flip it
# to LP — the engine must at least produce that (so the panel can catch it).
# ---------------------------------------------------------------------------
def test_high_am_on_rare_variant_pm2_strength(monkeypatch):
    from vcf2report.concordance import collapse_engine_tier, PATH
    a = _ann(gnomad_af=0.0, gnomad_faf95=0.0, am_pathogenicity=0.9995)
    # Default (PM2 Supporting): rare + a single in-silico predictor no longer reaches LP, so this is
    # NOT a potential gross discordance — the engine holds it at VUS.
    assert collapse_engine_tier(classify(_v(), a).tier) != PATH
    # Under the Moderate override it reaches LP/PATH — the failure mode the panel guards against
    # (a miscalibrated high AM on a ClinVar-benign variant would then be a gross discordance).
    monkeypatch.setenv("VCF2REPORT_PM2_STRENGTH", "moderate")
    assert collapse_engine_tier(classify(_v(), a).tier) == PATH


# ---------------------------------------------------------------------------
# AlphaMissense client: multi-transcript max + REF/ALT matching + no file
# ---------------------------------------------------------------------------
def test_alphamissense_best_takes_max_matching_transcript():
    v = Variant(chrom="17", pos=7670000, ref="C", alt="T")
    rows = [
        "chr17\t7670000\tC\tT\thg38\tP04637\tENST1\tp.A\t0.71\tlikely_pathogenic",
        "chr17\t7670000\tC\tT\thg38\tP04637\tENST2\tp.A\t0.98\tlikely_pathogenic",
        "chr17\t7670000\tC\tG\thg38\tP04637\tENST3\tp.A\t0.99\tlikely_pathogenic",  # wrong ALT
    ]
    best = alphamissense._best(rows, v)
    assert best["am_pathogenicity"] == pytest.approx(0.98)
    assert best["am_class"] == "likely_pathogenic"


def test_alphamissense_best_none_when_no_match():
    v = Variant(chrom="17", pos=7670000, ref="C", alt="T")
    assert alphamissense._best(["chr17\t7670000\tG\tA\th\tu\tt\tp\t0.9\tlikely_pathogenic"], v) is None


def test_alphamissense_lookup_no_local_file(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ALPHAMISSENSE_LOCAL", tmp_path / "missing.tsv.gz")
    monkeypatch.setattr(alphamissense, "_tabix", None)
    monkeypatch.setattr(alphamissense, "_tabix_tried", False)
    r = alphamissense.lookup(Variant(chrom="1", pos=1, ref="A", alt="T"))
    assert r["am_pathogenicity"] is None and "unavailable" in r["_source"]
