"""Annotation resolves from local datasets offline; end-to-end pipeline check."""

import os

from vcf2report import config
from vcf2report.annotate import annotate_variant
from vcf2report.models import Variant
from vcf2report.pipeline import run_pipeline


def test_offline_annotation_uses_local_snapshots(monkeypatch):
    monkeypatch.setenv("OFFLINE", "1")
    assert config.offline() is True
    v = Variant(chrom="2", pos=178562809, ref="G", alt="A", gene="TTN")
    a = annotate_variant(v, [])
    assert a.gnomad_af == 0.081  # from local snapshot
    # No cohort is shipped, so the local leg is UNKNOWN — never a checked zero.
    assert a.local_cohort_af is None
    assert "not in local table" in a.source["local_cohort"]
    assert "local" in a.source["gnomad"]


def test_clinvar_local_lookup():
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    a = annotate_variant(v, [])
    assert a.clinvar_significance == "Pathogenic"
    assert a.clinvar_accession == "VCV000012345"


def test_end_to_end_pipeline_tiers():
    hpo = ["HP:0001250", "HP:0002133", "HP:0011097"]
    report = run_pipeline(config.SAMPLE_VCF, hpo_terms=hpo)
    tiers = {c.variant.gene: c.tier for c in report.classifications}
    assert tiers["SCN1A"] == "Pathogenic"
    # PAX6: LoF in a LoF-intolerant gene, absent, no phenotype/ClinVar support
    # -> PVS1 + PM2 = Likely Pathogenic (an incidental finding, not over-called).
    assert tiers["PAX6"] == "Likely Pathogenic"
    # KCNQ2 p.Arg213Trp: the ClinVar assertion still contributes ONLY PP5 (supporting) — the
    # audit-#13 rule that a reputable-source call must not become a Strong line. PS1 is withheld
    # for exactly that reason. What carries it to Likely Pathogenic is INDEPENDENT evidence: PM1,
    # from 14 pathogenic-missense residues NEIGHBOURING Arg213 (the hotspot window excludes the
    # query residue, so it cannot restate the ClinVar record). Asserting the criteria, not just
    # the tier, keeps that distinction under test.
    kcnq2 = next(c for c in report.classifications if c.variant.gene == "KCNQ2")
    met = {x.code: x for x in kcnq2.criteria if x.met}
    assert met["PP5"].applied_strength == "supporting"
    assert "PS1" not in met, "ClinVar's own record must not also earn PS1"
    assert met["PM1"].applied_strength == "moderate"
    assert tiers["KCNQ2"] == "Likely Pathogenic"
    assert "VUS" in tiers["CACNA1A"]
    # TTN is dropped by gnomAD rarity alone (8.1%). OBSCN is NOT: it is absent from gnomAD
    # and this project ships no local cohort, so nothing establishes it as common and it
    # reaches the shortlist as a VUS. That is the honest result of having one frequency
    # source — a second cohort is the operator's to supply (annotate/local_cohort.py).
    assert "TTN" not in tiers
    assert tiers["OBSCN"] == "Uncertain Significance (VUS)"
    assert tiers["RB1"] == "Likely Pathogenic"  # incidental ACMG SF finding
    assert report.qc.candidates == 6
    assert report.qc.local_cohort_filtered == []
