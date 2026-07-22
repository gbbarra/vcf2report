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
                                       overview, variants_in_bucket, write_explore)

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
