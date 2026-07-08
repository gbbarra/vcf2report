"""Regression tests for the robustness-audit fixes (pure parser)."""
from vcf2report.annotate import from_vcf
from vcf2report.models import Variant
from vcf2report.pipeline import run_pipeline
from vcf2report.vcf import annparse
from vcf2report.vcf.parse import parse_vcf

H = "##fileformat=VCFv4.2\n##reference=GRCh38\n"


def _w(tmp_path, body, name="v.vcf", header_extra=""):
    p = tmp_path / name
    p.write_text(H + header_extra + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + body)
    return p


# --- star / symbolic ALT (#1) ------------------------------------------------
def test_star_allele_skipped_not_annotated(tmp_path):
    body = ("1\t100\t.\tC\tT,*\t800\tPASS\tANN=T|missense_variant|MODERATE|BRCA1|"
            "G|transcript|TX|protein_coding|1/2|c.1C>T|p.Arg1Cys|||||\tGT\t0/2\n")
    variants, _, _ = parse_vcf(_w(tmp_path, body))
    keys = {v.key for v in variants}
    assert "1-100-C-*" not in keys          # star allele never emitted
    assert "1-100-C-T" in keys              # real allele kept + correctly annotated
    t = next(v for v in variants if v.key == "1-100-C-T")
    assert t.gene == "BRCA1" and t.consequence == "missense_variant"


def test_symbolic_alt_skipped(tmp_path):
    body = "1\t100\t.\tC\t<DEL>\t800\tPASS\t.\tGT\t0/1\n"
    variants, _, _ = parse_vcf(_w(tmp_path, body))
    assert variants == []


def test_multiallelic_no_annotation_match_returns_none(tmp_path):
    # ANN only carries the T allele; the G allele must NOT inherit T's gene.
    body = ("1\t100\t.\tC\tG,T\t800\tPASS\tANN=T|stop_gained|HIGH|GENET|x|transcript|"
            "TX|protein_coding|1/2|c.1C>T|p.X|||||\tGT\t1/2\n")
    variants, _, _ = parse_vcf(_w(tmp_path, body))
    g = next(v for v in variants if v.key == "1-100-C-G")
    t = next(v for v in variants if v.key == "1-100-C-T")
    assert t.gene == "GENET"
    assert g.gene is None and g.consequence is None   # not borrowed from T


# --- multi-sample proband selection (#2) ------------------------------------
_TRIO = ("##fileformat=VCFv4.2\n##reference=GRCh38\n"
         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\tMOTHER\tPROBAND\n"
         "2\t166003360\t.\tC\tT\t800\tPASS\tGENE=SCN1A;CSQ=stop_gained\tGT:DP:GQ\t0/0:40:99\t0/0:40:99\t0/1:40:99\n")


def test_multi_sample_selects_named_proband(tmp_path):
    p = tmp_path / "trio.vcf"
    p.write_text(_TRIO)
    # default (first column = FATHER, 0/0) -> non-carrier, zygosity None
    default = parse_vcf(p)[0][0]
    assert default.zygosity is None
    # naming the proband -> the real heterozygous call
    proband = parse_vcf(p, sample="PROBAND")[0][0]
    assert proband.zygosity == "het"


def test_pipeline_warns_on_multisample_without_selection(tmp_path):
    p = tmp_path / "trio.vcf"
    p.write_text(_TRIO)
    report = run_pipeline(p, hpo_terms=["HP:0001250"])
    assert any("Multi-sample" in w for w in report.qc.warnings)


# --- report build reflects the DETECTED build (found via real-VCF validation) --
def test_report_build_matches_detected_build(tmp_path):
    body = "2\t100\t.\tA\tG\t.\tPASS\tGENE=X;CSQ=missense_variant\tGT\t0/1\n"
    p37 = _w(tmp_path, body, name="b37.vcf", header_extra="##reference=hg19\n")
    r37 = run_pipeline(p37, hpo_terms=[])
    assert r37.build == "GRCh37"        # not the assumed GRCh38
    p_unknown = tmp_path / "u.vcf"
    p_unknown.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
                         "FILTER\tINFO\tFORMAT\tS\n" + body)
    assert run_pipeline(p_unknown, hpo_terms=[]).build == "unknown"


# --- VEP CSQ specifics (#3, #4, #10) ----------------------------------------
def test_vep_allele_num_matching():
    fmt = ["Allele", "Consequence", "SYMBOL", "HGVSc", "HGVSp", "ALLELE_NUM"]
    csq = "T|missense_variant|GENEA|c.1A>C|p.X|1,T|stop_gained|GENEB|c.1A>T|p.Y|2"
    r1 = annparse.parse_vep(csq, "?", fmt, ref="A", alt_index=0, n_alt=2)
    r2 = annparse.parse_vep(csq, "?", fmt, ref="A", alt_index=1, n_alt=2)
    assert r1["gene"] == "GENEA"        # ALLELE_NUM 1
    assert r2["gene"] == "GENEB"        # ALLELE_NUM 2


def test_vep_percent_decoding():
    fmt = ["Allele", "Consequence", "SYMBOL", "HGVSc", "HGVSp"]
    csq = "A|synonymous_variant|MLH1|c.655C%3DT|p.%3D"   # %3D -> '='
    r = annparse.parse_vep(csq, "A", fmt)
    assert r["hgvs_c"] == "c.655C=T" and r["hgvs_p"] == "p.="


def test_vep_mane_select_when_no_pick_or_canonical():
    fmt = ["Allele", "Consequence", "SYMBOL", "HGVSc", "HGVSp", "MANE_SELECT"]
    csq = "A|missense_variant|G|c.1|p.1|,A|missense_variant|G|c.2|p.2|NM_000.1"
    r = annparse.parse_vep(csq, "A", fmt)
    assert r["hgvs_c"] == "c.2"          # the MANE_SELECT transcript


# --- from_vcf INFO parsing (#5, #6, #7) -------------------------------------
def test_pick_no_broadcast_on_array_length_mismatch():
    v = Variant(chrom="1", pos=1, ref="A", alt="G",
                info={"gnomad_AF": "0.30,0.40"}, alt_index=2)   # only 2 values, idx 2
    assert from_vcf.extract(v).get("gnomad_af") is None
    v2 = Variant(chrom="1", pos=1, ref="A", alt="G",
                 info={"gnomad_AF": "0.05"}, alt_index=1)       # scalar broadcasts
    assert from_vcf.extract(v2)["gnomad_af"] == 0.05


def test_clnsig_with_internal_comma_not_indexed():
    v = Variant(chrom="1", pos=1, ref="A", alt="G",
                info={"CLNSIG": "Pathogenic,_low_penetrance"}, alt_index=1)
    assert from_vcf.extract(v)["clinvar_significance"] == "Pathogenic, low penetrance"


def test_revel_multi_transcript_takes_max():
    v = Variant(chrom="1", pos=1, ref="A", alt="G",
                info={"REVEL": "0.10;0.90;."}, alt_index=0)
    assert from_vcf.extract(v)["revel"] == 0.90
