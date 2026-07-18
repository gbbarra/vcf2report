"""VUS triage — surface the phenotype-relevant, molecularly-suggestive VUS for expert review.

The engine correctly holds a variant with only Supporting evidence at Uncertain Significance;
this triage does NOT change that tier, it ranks which VUS are worth a human+model second look,
gated to the indication (phenotype-relevant) so a whole exome's incidental VUS are not noise.
"""
from vcf2report.models import Annotation, Classification, Variant
from vcf2report.report.vus_triage import probable_pathogenic_vus, exploration_prompt

_VUS = "Uncertain Significance (VUS)"


def _c(gene, *, tier=_VUS, am=None, hpo=0.9, clinvar=None, lof_intol=False, hp="p.Gly97Arg"):
    return Classification(
        variant=Variant(chrom="1", pos=1, ref="G", alt="A", gene=gene,
                        consequence="missense_variant", hgvs_p=hp),
        annotation=Annotation(am_pathogenicity=am, hpo_match_score=hpo,
                              clinvar_significance=clinvar, gene_lof_intolerant=lof_intol),
        criteria=[], tier=tier, rule_path="")


def test_strong_phenotype_matched_vus_is_prioritised():
    # RBSN-like: AM likely-path + phenotype 1.0 + ClinVar Pathogenic — held at VUS, but the top
    # thing to explore for this indication.
    r = probable_pathogenic_vus([_c("RBSN", am=0.968, hpo=1.0, clinvar="Pathogenic")])
    assert len(r) == 1 and r[0]["classification"].variant.gene == "RBSN"
    names = {s["signal"] for s in r[0]["signals"]}
    assert {"AlphaMissense likely-pathogenic", "phenotype match", "ClinVar Pathogenic assertion"} <= names


def test_does_not_change_tier():
    """Triage never promotes: the input stays VUS, it is only surfaced."""
    c = _c("RBSN", am=0.968, hpo=1.0, clinvar="Pathogenic")
    probable_pathogenic_vus([c])
    assert c.tier == _VUS


def test_phenotype_gate_excludes_incidental_vus():
    # A deleterious-looking VUS with NO phenotype overlap is an incidental finding, not a triage
    # priority for THIS indication — a whole exome carries dozens of these.
    assert probable_pathogenic_vus([_c("INCID", am=0.99, hpo=0.0, clinvar="Pathogenic")]) == []


def test_molecular_gate_excludes_bare_phenotype_overlap():
    # Phenotype overlap alone (no molecular signal) is not "probable-pathogenic".
    assert probable_pathogenic_vus([_c("BARE", am=None, hpo=1.0)]) == []


def test_only_vus_tier_considered():
    assert probable_pathogenic_vus([_c("PLP", tier="Pathogenic", am=0.99, hpo=1.0)]) == []


def test_ranked_by_evidence_weight():
    strong = _c("STRONG", am=0.97, hpo=1.0, clinvar="Pathogenic")     # AM3 + hpo2 + clinvar2 = 7
    weak = _c("WEAK", am=0.40, hpo=1.0)                               # AM-ambiguous1 + hpo2 = 3 -> gated out
    mid = _c("MID", am=0.60, hpo=1.0, lof_intol=True)                 # AM3 + hpo2 + lof1 = 6
    r = probable_pathogenic_vus([weak, mid, strong])
    genes = [e["classification"].variant.gene for e in r]
    assert genes == ["STRONG", "MID"]      # weak gated out (molecular < 2), strong ranks above mid


def test_conflicting_clinvar_is_a_signal():
    r = probable_pathogenic_vus([_c("CONF", am=0.6, hpo=1.0, clinvar="Conflicting interpretations")])
    assert r and any("conflicting" in s["signal"].lower() for s in r[0]["signals"])


def test_exploration_prompt_names_the_actions():
    r = probable_pathogenic_vus([_c("RBSN", am=0.968, hpo=1.0, clinvar="Pathogenic")])
    p = exploration_prompt(r[0])
    assert "RBSN" in p and "p.Gly97Arg" in p and "exploration" in p.lower()
