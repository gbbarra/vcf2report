"""The benchmark's --withhold-clinvar re-classification (run_benchmark.withhold_own_clinvar_tier).

The hpo-spiked-exomes benchmark plants KNOWN ClinVar pathogenics, so PP5 (+ the ClinVar safety
flag) recover them no matter what the residue/constraint criteria do — masking whether PS1/PM5/PP2
can carry a variant alone. Nulling the plant's own ClinVar assertion turns it into the novel-variant
scenario those criteria exist for. This locks the re-classification logic: with ClinVar removed, a
variant at an established pathogenic residue is recovered to P/LP by PS1 (which the own-record guard
had suppressed while the ClinVar assertion was present).
"""
import importlib.util
from pathlib import Path

from vcf2report.acmg.engine import classify
from vcf2report.models import Annotation, Variant

_RB = Path(__file__).resolve().parent.parent / "scripts" / "run_benchmark.py"
_spec = importlib.util.spec_from_file_location("run_benchmark", _RB)
rb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rb)


def _kcnq2_arg213trp():
    # A missense at an established pathogenic residue in a missense-constrained gene, rare in
    # gnomAD, and itself catalogued by ClinVar as Pathogenic — the typical benchmark plant.
    v = Variant(chrom="20", pos=63446204, ref="G", alt="A", gene="KCNQ2",
                consequence="missense_variant", hgvs_p="p.Arg213Trp")
    ps1 = {"alt_aa": "Trp", "ref_aa": "Arg", "stars": 2,
           "genomic_key": "20-63444712-G-A", "accession": "VCV000021795"}
    a = Annotation(gnomad_af=0.0, gnomad_faf95=0.0,
                   clinvar_significance="Pathogenic",
                   clinvar_review_status="criteria provided, single submitter",
                   clinvar_accession="VCV1", clinvar_ps1=ps1, clinvar_residue_available=True,
                   gene_mis_z=5.2, gene_missense_constrained=True)  # PP2
    return classify(v, a)


def test_normal_run_is_pp5_driven_ps1_suppressed():
    # With the ClinVar assertion present: PP5 fires, PS1 is withheld (double-count guard).
    c = _kcnq2_arg213trp()
    codes = {x.code: x for x in c.criteria}
    assert codes["PP5"].met is True
    assert codes["PS1"].met is False


def test_withhold_flips_pp5_to_ps1():
    # The mechanism, tier-independent: nulling the own ClinVar assertion turns PP5 off and lets
    # PS1 fire from the residue index (no longer suppressed as a double-count).
    from dataclasses import replace
    c = _kcnq2_arg213trp()
    ann = replace(c.annotation, clinvar_significance=None, clinvar_review_status=None,
                  clinvar_accession=None, clinvar_condition=None)
    codes = {x.code: x for x in classify(c.variant, ann).criteria}
    assert codes["PP5"].met is False
    assert codes["PS1"].met is True


def test_withhold_recovers_plp_via_ps1():
    # Null the own ClinVar assertion -> PS1(strong) + PM2(supporting) + PP2(supporting) reaches
    # Likely Pathogenic on the residue/constraint/rarity evidence alone. This is what the mode measures.
    c = _kcnq2_arg213trp()
    withheld = rb.withhold_own_clinvar_tier(c)
    assert withheld.startswith(("Pathogenic", "Likely Pathogenic")), withheld


def test_withhold_does_not_mutate_original():
    # dataclasses.replace must leave the original classification's annotation untouched.
    c = _kcnq2_arg213trp()
    rb.withhold_own_clinvar_tier(c)
    assert c.annotation.clinvar_significance == "Pathogenic"
