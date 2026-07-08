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
