"""Explore — persist a run as queryable JSON and answer follow-up questions off it, no re-run.

Locks the persisted schema (so a downstream conversation can rely on it) and the read helpers that
answer the conversational questions the module exists for: "show gene X", "why did Y get PM2",
"which findings rest on ClinVar", "summarise / open this case".
"""
import pytest

from vcf2report.models import Annotation, Classification, CriterionResult, QCSummary, Variant
from vcf2report.report.assemble import build_report
from vcf2report.report.explore import (BUCKETS, build_explore, criterion_basis, explain,
                                       findings_citing_clinvar, findings_for_gene, load_explore,
                                       missense_evidence, missing_evidence, overview,
                                       variants_in_bucket, write_explore)

_VUS = "Uncertain Significance (VUS)"


def _crit(code, *, met=True, strength="supporting", citation=None, evidence=None):
    return CriterionResult(code=code, name=f"{code} criterion", default_strength=strength,
                           applies=True, met=met, applied_strength=strength,
                           citation=citation or [], evidence=evidence or {})


def _c(gene, *, tier, hpo, zyg="het", gnomad=1e-6, clinvar=None, review=None, criteria=None):
    return Classification(
        variant=Variant(chrom="1", pos=100, ref="G", alt="A", gene=gene,
                        consequence="missense_variant", hgvs_p="p.Gly97Arg", zygosity=zyg),
        annotation=Annotation(hpo_match_score=hpo, gnomad_af=gnomad,
                              clinvar_significance=clinvar, clinvar_review_status=review),
        criteria=criteria or [], tier=tier, rule_path=f"{tier} path")


def _report():
    # SCN1A: phenotype-matched Pathogenic, with PP5 citing a ClinVar VCV accession + a met PM2.
    scn1a = _c("SCN1A", tier="Pathogenic", hpo=1.0, zyg="hom",
               criteria=[_crit("PM2", citation=["gnomAD v4 (local)"], evidence={"popmax_af": 0.0}),
                         _crit("PP5", citation=["VCV000012345"]),
                         _crit("BS1", met=False)])
    # FLAG: engine holds it at VUS, but ClinVar calls it Pathogenic with expert-panel (3-star)
    # review — the do-not-dismiss safety net. Unrelated phenotype so it routes to `other`.
    flag = _c("FLAG", tier=_VUS, hpo=0.0, clinvar="Pathogenic",
              review="reviewed by expert panel", criteria=[_crit("PM2")])
    return build_report("CASE-1", ["HP:0001250"], QCSummary(candidates=2), [scn1a, flag])


# --- write side: the persisted schema -------------------------------------------------
def test_build_explore_has_the_documented_shape():
    d = build_explore(_report())
    for k in ("sample_id", "build", "classifications", "conclusion", "buckets",
              "clinvar_do_not_dismiss"):
        assert k in d
    assert set(d["buckets"]) == set(BUCKETS)
    assert d["buckets"]["primary"] == ["SCN1A"]
    assert "FLAG" in d["buckets"]["other"]
    assert isinstance(d["conclusion"], list) and d["conclusion"]


def test_clinvar_do_not_dismiss_is_structured_not_a_repr_string():
    # Regression: the list held raw Classification objects, which json.dump(default=str) would have
    # flattened to a dataclass repr — unqueryable. It must be compact dicts.
    d = build_explore(_report())
    dnd = d["clinvar_do_not_dismiss"]
    assert len(dnd) == 1 and isinstance(dnd[0], dict)
    assert dnd[0]["gene"] == "FLAG"
    assert dnd[0]["clinvar_stars"] == 3           # "reviewed by expert panel"
    assert dnd[0]["engine_tier"] == _VUS


def test_write_and_load_roundtrip(tmp_path):
    p = tmp_path / "CASE-1_results.json"
    write_explore(_report(), p)
    d = load_explore(p)
    assert d["sample_id"] == "CASE-1"
    assert {c["variant"]["gene"] for c in d["classifications"]} == {"SCN1A", "FLAG"}


# --- read side: the conversational queries --------------------------------------------
def test_findings_for_gene_is_case_insensitive():
    d = build_explore(_report())
    hits = findings_for_gene(d, "scn1a")
    assert len(hits) == 1 and hits[0]["variant"]["gene"] == "SCN1A"


def test_variants_in_bucket_and_unknown_bucket_raises():
    d = build_explore(_report())
    assert [c["variant"]["gene"] for c in variants_in_bucket(d, "primary")] == ["SCN1A"]
    with pytest.raises(ValueError):
        variants_in_bucket(d, "not_a_bucket")


def test_criterion_basis_answers_why_gene_got_a_code():
    d = build_explore(_report())
    basis = criterion_basis(d, "SCN1A", "pm2")   # case-insensitive
    assert len(basis) == 1
    assert basis[0]["code"] == "PM2" and basis[0]["met"] is True
    assert basis[0]["citation"] == ["gnomAD v4 (local)"]


def test_findings_citing_clinvar_catches_pp5_via_vcv_accession():
    # PP5 cites "VCV000012345", not the literal word ClinVar — a text-only match would miss it.
    d = build_explore(_report())
    cv = findings_citing_clinvar(d)
    genes = {f["gene"] for f in cv}
    assert "SCN1A" in genes and "FLAG" not in genes
    codes = {cr["code"] for f in cv if f["gene"] == "SCN1A" for cr in f["criteria"]}
    assert "PP5" in codes


def test_explain_gives_a_gene_digest_with_bucket_membership():
    d = build_explore(_report())
    e = explain(d, "SCN1A")
    assert e["gene"] == "SCN1A" and "primary" in e["buckets"]
    assert e["variants"][0]["tier"] == "Pathogenic"
    assert set(e["variants"][0]["met_codes"]) == {"PM2", "PP5"}


def test_overview_counts_buckets_and_carries_the_conclusion():
    d = build_explore(_report())
    o = overview(d)
    assert o["sample_id"] == "CASE-1" and o["n_candidates"] == 2
    assert o["bucket_counts"]["primary"] == 1
    assert len(o["clinvar_do_not_dismiss"]) == 1
    assert isinstance(o["conclusion"], list) and o["conclusion"]


# --- gene/residue missense evidence (PP2/BP1/PS1/PM5) ---------------------------------
def _missense_report():
    # HOT: a novel missense elevated by residue + constraint evidence (PM5 + PP2), no ClinVar
    # record of its own. COLD: a missense where neither fired (the metric/index said no).
    hot = _c("HOT", tier="Likely Pathogenic", hpo=1.0,
             criteria=[_crit("PM5", strength="moderate", citation=["VCV000067890"]),
                       _crit("PP2"),
                       _crit("PS1", met=False, strength="strong"),
                       _crit("BP1", met=False)])
    cold = _c("COLD", tier=_VUS, hpo=0.9,
              criteria=[_crit("PP2", met=False), _crit("PM5", met=False), _crit("PM2")])
    return build_report("CASE-2", ["HP:0001250"], QCSummary(candidates=2), [hot, cold])


def test_missense_evidence_lists_only_fired_criteria_by_default():
    d = build_explore(_missense_report())
    ev = missense_evidence(d)
    assert [f["gene"] for f in ev] == ["HOT"]          # COLD has none met -> omitted
    codes = {cr["code"] for cr in ev[0]["criteria"]}
    assert codes == {"PM5", "PP2"}
    assert ev[0]["tier"] == "Likely Pathogenic"
    assert ev[0]["hgvs_p"] == "p.Gly97Arg"


def test_missense_evidence_all_criteria_shows_the_not_met_audit_view():
    d = build_explore(_missense_report())
    ev = missense_evidence(d, met_only=False)
    genes = {f["gene"] for f in ev}
    assert genes == {"HOT", "COLD"}                     # COLD now visible, with its not-met trail
    cold = next(f for f in ev if f["gene"] == "COLD")
    assert {cr["code"] for cr in cold["criteria"]} == {"PP2", "PM5"}
    assert all(cr["met"] is False for cr in cold["criteria"])


def test_ps1_pm5_count_as_clinvar_citing_criteria():
    # PS1/PM5 rest on a residue-level cross-match to a ClinVar record, so a "what rests on ClinVar?"
    # question must surface them alongside PP5/BP6.
    d = build_explore(_missense_report())
    cv = findings_citing_clinvar(d)
    hot = next(f for f in cv if f["gene"] == "HOT")
    assert "PM5" in {cr["code"] for cr in hot["criteria"]}


# --- "what would move this off VUS?" ---------------------------------------
def _na_crit(code, strength, reasoning="Requires parental (trio) data — not available from a "
                                       "single proband VCF"):
    return CriterionResult(code=code, name=f"{code} criterion", default_strength=strength,
                           applies=False, met=False, reasoning=reasoning)


def _gap_report():
    # A VUS resting on one Moderate + two Supporting: one more Moderate reaches Likely Pathogenic
    # (Richards LP-5). PM6 is carried as N/A so it can be named as the concrete way to get there.
    vus = _c("GAPG", tier=_VUS, hpo=0.9,
             criteria=[_crit("PM2", strength="moderate"), _crit("PP3"), _crit("PP4"),
                       _na_crit("PM6", "moderate"), _na_crit("PS2", "strong")])
    # A Pathogenic call that no single supporting/strong line can move.
    path = _c("SOLID", tier="Pathogenic", hpo=1.0,
              criteria=[_crit("PVS1", strength="very_strong"), _crit("PM2", strength="moderate"),
                        _crit("PP5")])
    return build_report("CASE-3", ["HP:0001250"], QCSummary(candidates=2), [vus, path])


def test_missing_evidence_names_the_weakest_addition_that_would_change_the_tier():
    d = build_explore(_gap_report())
    gap = missing_evidence(d, "GAPG")[0]
    assert gap["tier"] == _VUS
    up = [s for s in gap["would_change_with"] if s["direction"] == "up"]
    assert len(up) == 1, "should report exactly one (weakest) upgrade step"
    assert up[0]["strength"] == "moderate" and up[0]["side"] == "pathogenic"
    assert up[0]["would_become"] == "Likely Pathogenic"


def test_missing_evidence_names_the_concrete_next_step():
    # The answer must be actionable — "order a trio" — not just an abstract strength.
    d = build_explore(_gap_report())
    up = [s for s in missing_evidence(d, "GAPG")[0]["would_change_with"] if s["direction"] == "up"]
    codes = {c["code"] for c in up[0]["candidates"]}
    assert "PM6" in codes
    needs = next(c["needs"] for c in up[0]["candidates"] if c["code"] == "PM6")
    assert "parental" in needs


def test_pending_criteria_are_read_off_the_trail_not_a_hardcoded_list():
    """The "what's still available" set must come from the criteria themselves.

    A hardcoded map in the report layer goes stale the moment a criterion changes class — as it
    did when PM1 became engine-decided but was still advertised here as an orderable next step.
    Deriving it from `applies is False or adjudicated_by == "model"` means a criterion is listed
    if and only if it really is still open, whatever `acmg/criteria.py` does next.
    """
    from vcf2report.acmg.engine import evaluate_criteria
    from vcf2report.models import Variant as V

    var = V(chrom="1", pos=100, ref="G", alt="A", gene="REALG",
            consequence="missense_variant", hgvs_p="p.Gly97Arg", zygosity="het")
    real = Classification(variant=var, annotation=Annotation(),
                          criteria=evaluate_criteria(var, Annotation()),
                          tier=_VUS, rule_path="VUS path")
    d = build_explore(build_report("CASE-R", [], QCSummary(candidates=1), [real]))
    pending = {p["code"] for p in missing_evidence(d)[0]["pending_criteria"]}

    # exactly the single-proband-undecidable + judgement criteria, from the engine itself
    assert pending == {"PS2", "PM3", "PM6", "PP1", "BP2", "BS4",     # N/A
                       "PS3", "PS4", "BS3", "BP3", "BP5"}            # model-adjudicated
    # engine-decided criteria that merely did not fire are NOT offered as next steps
    assert not pending & {"PM1", "PS1", "PM5", "PP2", "BP1", "PVS1", "PM2", "PP3", "PP5"}


def test_missing_evidence_routes_the_benign_hypothetical_to_the_benign_side():
    # Regression: rules.combine decides the side by CODE membership, so a made-up benign code was
    # scored as PATHOGENIC — reporting that a benign criterion would make a VUS Likely Pathogenic.
    d = build_explore(_gap_report())
    for r in missing_evidence(d):
        for s in r["would_change_with"]:
            if s["side"] == "benign":
                assert s["direction"] == "down", f"benign line raised the tier: {s}"


def test_missing_evidence_skips_the_nonexistent_benign_moderate_bucket():
    # Richards Table 5 has no benign Moderate; escalating one would be a silent no-op step.
    d = build_explore(_gap_report())
    for r in missing_evidence(d):
        assert not [s for s in r["would_change_with"]
                    if s["side"] == "benign" and s["strength"] == "moderate"]


def test_missing_evidence_shows_a_robust_pathogenic_call_as_hard_to_move():
    d = build_explore(_gap_report())
    solid = next(r for r in missing_evidence(d) if r["gene"] == "SOLID")
    assert solid["tier"] == "Pathogenic"
    assert not [s for s in solid["would_change_with"] if s["direction"] == "up"]
    down = [s for s in solid["would_change_with"] if s["direction"] == "down"]
    assert down and down[0]["strength"] == "stand_alone"


def test_missing_evidence_covers_every_variant_when_no_gene_given():
    d = build_explore(_gap_report())
    assert {r["gene"] for r in missing_evidence(d)} == {"GAPG", "SOLID"}
