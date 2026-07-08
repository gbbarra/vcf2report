"""T4: GA4GH Phenopacket -> pipeline inputs."""
from vcf2report import config
from vcf2report.phenopacket import load_phenopacket, write_inputs
from vcf2report.vcf.parse import parse_vcf

EXAMPLE = config.DATA_DIR / "sample" / "example_phenopacket.json"


def test_load_extracts_hpo_and_variant():
    data = load_phenopacket(EXAMPLE)
    assert data["subject_id"] == "DEMO-001"
    # excluded feature (Hydrocephalus) must be skipped
    assert data["hpo_terms"] == ["HP:0001250", "HP:0001263", "HP:0002133"]
    assert len(data["variants"]) == 1
    v = data["variants"][0]
    assert v["chrom"] == "2" and v["pos"] == 166003360 and v["ref"] == "C" and v["alt"] == "T"
    assert v["gene"] == "SCN1A"
    assert v["hgvs_p"] == "p.Arg612Ter"
    assert v["zygosity"] == "het"


def test_phenopacket_end_to_end_and_annotation_dependency(tmp_path):
    """A real-shaped phenopacket runs end-to-end; consequence (annotation) is the
    completing step that moves the causal variant from VUS to Pathogenic."""
    from vcf2report.phenopacket import load_phenopacket, write_inputs
    from vcf2report.pipeline import run_pipeline

    data = load_phenopacket(EXAMPLE)
    vcf, hpo = tmp_path / "c.vcf", tmp_path / "c.hpo.txt"
    write_inputs(data, vcf, hpo)

    # Raw phenopacket VCF has no molecular consequence. SCN1A is retained (ClinVar
    # P/LP bypasses the impact filter) but PVS1 can't fire -> VUS.
    raw = run_pipeline(vcf, hpo_terms=data["hpo_terms"])
    scn_raw = next(c for c in raw.classifications if c.variant.gene == "SCN1A")
    assert "VUS" in scn_raw.tier

    # After annotation (the consequence SnpEff/VEP would add) -> Pathogenic.
    annotated = tmp_path / "c.annotated.vcf"
    annotated.write_text(vcf.read_text().replace("GENE=SCN1A;", "GENE=SCN1A;CSQ=stop_gained;"))
    ann = run_pipeline(annotated, hpo_terms=data["hpo_terms"])
    scn_ann = next(c for c in ann.classifications if c.variant.gene == "SCN1A")
    assert scn_ann.tier == "Pathogenic"


def test_write_inputs_roundtrips_through_parser(tmp_path):
    data = load_phenopacket(EXAMPLE)
    vcf, hpo = tmp_path / "c.vcf", tmp_path / "c.hpo.txt"
    write_inputs(data, vcf, hpo)

    assert hpo.read_text().split() == ["HP:0001250", "HP:0001263", "HP:0002133"]
    variants, build, _ = parse_vcf(vcf)
    assert build == "GRCh38"
    assert len(variants) == 1
    v = variants[0]
    assert v.key == "2-166003360-C-T"
    assert v.gene == "SCN1A"
    assert v.hgvs_p == "p.Arg612Ter"
    assert v.zygosity == "het"
