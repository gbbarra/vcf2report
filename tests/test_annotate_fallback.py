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
    assert a.gnomad_af == 0.081          # from local snapshot
    assert a.abraom_af == 0.075          # from ABraOM local
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
    # KCNQ2: ClinVar P contributes only PP5 (supporting), so 1 PM + 3 PP -> VUS,
    # not the old PS1-driven Likely Pathogenic.
    assert "VUS" in tiers["KCNQ2"]
    assert "VUS" in tiers["CACNA1A"]
    # OBSCN dropped by ABraOM, TTN dropped by rarity -> not classified
    assert "OBSCN" not in tiers
    assert "TTN" not in tiers
    assert tiers["RB1"] == "Likely Pathogenic"   # incidental ACMG SF finding
    assert report.qc.candidates == 5
    assert any("OBSCN" in n for n in report.qc.abraom_filtered)
