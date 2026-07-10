"""Concordance panel: harness math, engine path, builder parsing, and the
forever-offline safety guard on the real frozen panel.

The harness logic is proven on controlled synthetic inputs (so the confusion
matrix / metrics / gross-discordance detection are locked regardless of bundled
data), and the real frozen panel — once built — is guarded by a hard invariant:
the engine must never flip a ClinVar-benign variant to Pathogenic or vice-versa.
"""
import importlib.util
import types
from pathlib import Path

import pytest

from vcf2report import concordance
from vcf2report.acmg.rules import (
    BENIGN, LIKELY_BENIGN, LIKELY_PATHOGENIC, PATHOGENIC, VUS,
)
from vcf2report.concordance import (
    BEN, PATH, UNCERTAIN, PanelEntry, collapse_clinvar, collapse_engine_tier,
)
from vcf2report.models import Variant

_SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load_builder():
    fp = _SCRIPTS / "build_concordance_panel.py"
    spec = importlib.util.spec_from_file_location("build_concordance_panel", fp)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _entry(key, gene, consequence, truth_sig):
    chrom, pos, ref, alt = key.split("-")
    return PanelEntry(
        variant=Variant(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                        gene=gene, consequence=consequence),
        truth_significance=truth_sig,
        truth_class=collapse_clinvar(truth_sig),
    )


# ---------------------------------------------------------------------------
# Collapsing helpers
# ---------------------------------------------------------------------------
def test_collapse_engine_tier():
    assert collapse_engine_tier(PATHOGENIC) == PATH
    assert collapse_engine_tier(LIKELY_PATHOGENIC) == PATH
    assert collapse_engine_tier(VUS) == UNCERTAIN
    assert collapse_engine_tier(LIKELY_BENIGN) == BEN
    assert collapse_engine_tier(BENIGN) == BEN


def test_collapse_clinvar():
    assert collapse_clinvar("Pathogenic") == PATH
    assert collapse_clinvar("Likely pathogenic") == PATH
    assert collapse_clinvar("Pathogenic/Likely pathogenic") == PATH
    assert collapse_clinvar("Benign") == BEN
    assert collapse_clinvar("Likely benign") == BEN
    assert collapse_clinvar("Conflicting classifications of pathogenicity") == UNCERTAIN
    assert collapse_clinvar("Uncertain significance") == UNCERTAIN
    assert collapse_clinvar("") is None
    assert collapse_clinvar(None) is None


# ---------------------------------------------------------------------------
# Confusion matrix + metrics (engine mocked → pure harness math)
# ---------------------------------------------------------------------------
def test_matrix_and_metrics(monkeypatch):
    entries = [
        _entry("1-1-A-T", "G1", "missense_variant", "Pathogenic"),        # e1 -> PATH
        _entry("1-2-A-T", "G2", "missense_variant", "Pathogenic"),        # e2 -> VUS
        _entry("1-3-A-T", "G3", "stop_gained", "Pathogenic"),             # e3 -> PATH (LoF)
        _entry("1-4-A-T", "G4", "missense_variant", "Benign"),            # e4 -> BEN
        _entry("1-5-A-T", "G5", "missense_variant", "Benign"),            # e5 -> VUS
        _entry("1-6-A-T", "G6", "missense_variant", "Benign"),            # e6 -> PATH (GROSS)
    ]
    tiers = {
        "1-1-A-T": LIKELY_PATHOGENIC, "1-2-A-T": VUS, "1-3-A-T": PATHOGENIC,
        "1-4-A-T": BENIGN, "1-5-A-T": VUS, "1-6-A-T": PATHOGENIC,
    }
    monkeypatch.setattr(concordance, "classify_entry", lambda e, withhold_clinvar=True:
                        types.SimpleNamespace(
                            tier=tiers[e.variant.key], rule_path="mock", met_codes=[]))

    res = concordance.evaluate_panel(entries)
    m = res.metrics
    assert m["n"] == 6 and m["n_pathogenic"] == 3 and m["n_benign"] == 3
    assert m["n_lof_pathogenic"] == 1
    assert res.matrix[PATH] == {PATH: 2, UNCERTAIN: 1, BEN: 0}
    assert res.matrix[BEN] == {PATH: 1, UNCERTAIN: 1, BEN: 1}
    assert m["concordance"] == pytest.approx(0.5)
    assert m["pathogenic_sensitivity"] == pytest.approx(2 / 3, abs=1e-4)
    assert m["lof_pathogenic_sensitivity"] == pytest.approx(1.0)
    assert m["benign_agreement"] == pytest.approx(1 / 3, abs=1e-4)
    assert m["benign_specificity"] == pytest.approx(2 / 3, abs=1e-4)
    assert m["gross_discordances"] == 1
    # Decisiveness / precision framing (4 decisive: e1,e3,e4,e6).
    assert m["decisiveness"] == pytest.approx(4 / 6, abs=1e-4)
    assert m["concordance_when_decisive"] == pytest.approx(3 / 4)  # e6 is the miss
    assert m["pathogenic_precision"] == pytest.approx(2 / 3, abs=1e-4)  # engine PATH: e1,e3,e6
    assert m["benign_precision"] == pytest.approx(1.0)  # engine BEN: e4


def test_gross_discordance_surfaced(monkeypatch):
    entries = [_entry("1-6-A-T", "G6", "missense_variant", "Benign")]
    monkeypatch.setattr(concordance, "classify_entry", lambda e, withhold_clinvar=True:
                        types.SimpleNamespace(tier=PATHOGENIC, rule_path="x", met_codes=[]))
    res = concordance.evaluate_panel(entries)
    gross = res.gross_discordances
    assert len(gross) == 1 and gross[0].truth_class == BEN and gross[0].engine_class == PATH
    assert "Gross discordances" in res.to_markdown()


def test_gross_discordance_pathogenic_called_benign(monkeypatch):
    """The most dangerous direction: engine calls a ClinVar-pathogenic variant Benign."""
    entries = [_entry("1-7-A-T", "G7", "missense_variant", "Pathogenic")]
    monkeypatch.setattr(concordance, "classify_entry", lambda e, withhold_clinvar=True:
                        types.SimpleNamespace(tier=BENIGN, rule_path="x", met_codes=[]))
    res = concordance.evaluate_panel(entries)
    assert res.matrix[PATH][BEN] == 1
    assert res.metrics["gross_discordances"] == 1
    gross = res.gross_discordances
    assert len(gross) == 1 and gross[0].truth_class == PATH and gross[0].engine_class == BEN
    assert "Gross discordances" in res.to_markdown()


# ---------------------------------------------------------------------------
# The real engine path on controlled inputs (constraint / ABraOM patched)
# ---------------------------------------------------------------------------
def _patch_local(monkeypatch, lof_intolerant, abraom_af=0.0):
    monkeypatch.setattr(concordance.extra, "gene_constraint",
                        lambda gene: {"lof_intolerant": lof_intolerant, "_source": "test"})
    monkeypatch.setattr(concordance.abraom, "lookup",
                        lambda v: {"af": abraom_af, "_source": "test"})


def test_engine_calls_rare_lof_pathogenic(monkeypatch):
    _patch_local(monkeypatch, lof_intolerant=True)
    e = _entry("2-100-C-T", "TESTLOF", "stop_gained", "Pathogenic")
    e.frozen_gnomad = {"af": 0.0, "faf95": 0.0, "ac": 0, "an": 152000, "hom": 0, "pop": None}
    c = concordance.classify_entry(e, withhold_clinvar=True)
    assert collapse_engine_tier(c.tier) == PATH
    assert "PVS1" in c.met_codes and "PM2" in c.met_codes


def test_engine_calls_common_benign(monkeypatch):
    _patch_local(monkeypatch, lof_intolerant=False)
    e = _entry("3-200-A-G", "TESTBEN", "missense_variant", "Benign")
    e.frozen_gnomad = {"af": 0.12, "faf95": 0.11, "ac": 18000, "an": 150000,
                       "hom": 50, "pop": "nfe"}
    c = concordance.classify_entry(e, withhold_clinvar=True)
    assert collapse_engine_tier(c.tier) == BEN
    assert "BA1" in c.met_codes


def test_clinvar_withheld_suppresses_pp5(monkeypatch):
    _patch_local(monkeypatch, lof_intolerant=True)
    e = _entry("2-166003360-C-T", "SCN1A", "stop_gained", "Pathogenic")
    e.review_status = "criteria provided, multiple submitters, no conflicts"
    e.accession = "VCV000012345"
    e.frozen_gnomad = {"af": 0.0, "faf95": 0.0, "ac": 0, "an": 152000, "hom": 0, "pop": None}

    withheld = concordance.classify_entry(e, withhold_clinvar=True)
    assert "PP5" not in withheld.met_codes

    included = concordance.classify_entry(e, withhold_clinvar=False)
    assert "PP5" in included.met_codes  # wiring sanity: ClinVar does reach the engine


def test_frozen_alphamissense_recovers_missense(monkeypatch):
    _patch_local(monkeypatch, lof_intolerant=False)
    e = _entry("5-100-A-T", "TESTG", "missense_variant", "Pathogenic")
    e.frozen_gnomad = {"af": 0.0, "faf95": 0.0, "ac": 0, "an": 152000, "hom": 0, "pop": None}
    e.frozen_alphamissense = {"am_pathogenicity": 0.999, "am_class": "likely_pathogenic"}
    c = concordance.classify_entry(e, withhold_clinvar=True)
    assert "PP3" in c.met_codes and "PM2" in c.met_codes
    assert collapse_engine_tier(c.tier) == PATH   # VUS -> PATH via calibrated PP3_Strong


def test_load_panel_attaches_alphamissense(tmp_path):
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text("# h\n2-166003360-C-T\tSCN1A\tmissense_variant\t\tPathogenic\tcrit\tVCV1\tx\n")
    frozen = tmp_path / "gnomad.json"
    frozen.write_text('{"2-166003360-C-T": {"af": 0.0}}')
    am = tmp_path / "am.json"
    am.write_text('{"2-166003360-C-T": {"am_pathogenicity": 0.97, "am_class": "likely_pathogenic"}}')
    entries = concordance.load_panel(gt, frozen, am)
    assert len(entries) == 1
    assert entries[0].frozen_alphamissense["am_pathogenicity"] == pytest.approx(0.97)


# ---------------------------------------------------------------------------
# load_panel: keeps only definite PATH/BEN truth, attaches frozen gnomAD
# ---------------------------------------------------------------------------
def test_load_panel_filters_and_attaches(tmp_path):
    gt = tmp_path / "ground_truth.tsv"
    gt.write_text(
        "# header\n"
        "2-166003360-C-T\tSCN1A\tstop_gained\tp.Arg612Ter\tPathogenic\tcriteria_provided\tVCV1\tDravet\n"
        "7-100-A-G\tCFTR\tmissense_variant\t\tBenign\tcriteria_provided\tVCV2\t-\n"
        "1-9-A-T\tXX\tmissense_variant\t\tConflicting classifications of pathogenicity\t-\tVCV3\t-\n"
    )
    frozen = tmp_path / "gnomad_frozen.json"
    frozen.write_text('{"2-166003360-C-T": {"af": 0.0, "faf95": 0.0}}')

    entries = concordance.load_panel(gt, frozen)
    keys = {e.variant.key for e in entries}
    assert keys == {"2-166003360-C-T", "7-100-A-G"}  # conflicting row dropped
    scn = next(e for e in entries if e.variant.key == "2-166003360-C-T")
    assert scn.truth_class == PATH and scn.frozen_gnomad == {"af": 0.0, "faf95": 0.0}
    cftr = next(e for e in entries if e.variant.key == "7-100-A-G")
    assert cftr.truth_class == BEN and cftr.frozen_gnomad == {}  # no frozen record


# ---------------------------------------------------------------------------
# Builder parsing (network mocked)
# ---------------------------------------------------------------------------
def test_builder_grch38_snv_and_consequence():
    build = _load_builder()
    assert build._reviewed("criteria provided, single submitter") is True
    assert build._reviewed("no assertion criteria provided") is False
    vset = [{"variation_loc": [{"assembly_name": "GRCh38", "chr": "2",
                                "start": "166003360", "ref": "C", "alt": "T"}]}]
    assert build._grch38_snv(vset) == ("2-166003360-C-T", "C", "T")
    # indel rejected in v1
    vset_indel = [{"variation_loc": [{"assembly_name": "GRCh38", "chr": "2",
                                      "start": "10", "ref": "CA", "alt": "C"}]}]
    assert build._grch38_snv(vset_indel) is None
    # Real ClinVar shape is molecular_consequence_list (a list).
    assert build._consequence_from({"molecular_consequence_list": ["nonsense"]}) == "stop_gained"
    # Most severe consequence wins when several are co-listed (LoF not masked).
    assert build._consequence_from(
        {"molecular_consequence_list": ["missense variant", "splice donor variant"]}
    ) == "splice_donor_variant"
    # The singular key is still honoured as a defensive fallback.
    assert build._consequence_from({"molecular_consequence": ["nonsense"]}) == "stop_gained"


def test_builder_harvest_gene_mocked(monkeypatch):
    build = _load_builder()
    esummary = {"result": {"uids": ["1"], "1": {
        "accession": "VCV000012345",
        "molecular_consequence_list": ["nonsense"],
        "title": "NM_006920.6(SCN1A):c.1834C>T (p.Arg612Ter)",
        "germline_classification": {
            "description": "Pathogenic",
            "review_status": "criteria provided, multiple submitters, no conflicts",
            "trait_set": [{"trait_name": "Dravet syndrome"}]},
        "variation_set": [{"variation_loc": [{
            "assembly_name": "GRCh38", "chr": "2", "start": "166003360",
            "ref": "C", "alt": "T"}]}],
    }}}
    monkeypatch.setattr(build._http, "throttle", lambda *a, **k: None)
    monkeypatch.setattr(build._http, "get_json", lambda url, params, **k: (
        {"esearchresult": {"idlist": ["1"]}} if "esearch" in url else esummary))

    rows = build._harvest_gene("SCN1A", "pathogenic", want=5, seen=set())
    assert len(rows) == 1
    r = rows[0]
    assert r["key"] == "2-166003360-C-T"
    assert r["consequence"] == "stop_gained"
    assert r["clinvar_significance"] == "Pathogenic"
    assert r["hgvs_p"] == "p.Arg612Ter"


# ---------------------------------------------------------------------------
# The forever-offline guard on the REAL frozen panel (skips until built)
# ---------------------------------------------------------------------------
def test_frozen_panel_safety_invariant():
    if not concordance.FROZEN_GNOMAD.exists() or not concordance.GROUND_TRUTH.exists():
        pytest.skip("panel not frozen — run scripts/build_concordance_panel.py")
    entries = concordance.load_panel()
    if not entries:
        pytest.skip("frozen panel has no PATH/BEN-labelled variants")
    res = concordance.evaluate_panel(entries)
    # Safety invariant: the engine never flips PATH<->BEN vs ClinVar.
    assert res.metrics["gross_discordances"] == 0, res.to_markdown()
    # A ClinVar-benign variant is (almost) never called Pathogenic. Only assert
    # when the panel actually has benign truth — else the metric is undefined
    # (_rate(0, 0) == 0.0) and would spuriously fail without any real violation.
    if res.metrics["n_benign"]:
        assert res.metrics["benign_specificity"] >= 0.9, res.to_markdown()
