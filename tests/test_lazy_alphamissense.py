"""Perf: AlphaMissense is queried lazily — only for the surviving candidates.

The score feeds PP3/BP4 at classification time and NEVER the rarity/impact filter,
so deferring the (~1 GB tabix) lookup from the whole post-QC set to just the
candidates cannot change any classification — it only avoids thousands of wasted
lookups on a real exome.
"""
from vcf2report import config, pipeline
from vcf2report.annotate import add_alphamissense, alphamissense, annotate_variant
from vcf2report.models import Annotation, Variant


def _v(**kw):
    base = dict(chrom="1", pos=100, ref="A", alt="T", gene="G", consequence="missense_variant")
    base.update(kw)
    return Variant(**base)


def test_annotate_variant_defers_client(monkeypatch):
    def boom(v):
        raise AssertionError("client must NOT be called when with_alphamissense=False")
    monkeypatch.setattr(alphamissense, "lookup", boom)
    a = annotate_variant(_v(), with_alphamissense=False)
    assert a.am_pathogenicity is None
    assert "deferred" in a.source["alphamissense"]


def test_deferred_still_reads_am_from_vcf_info(monkeypatch):
    def boom(v):
        raise AssertionError("client must NOT be called when INFO already has the score")
    monkeypatch.setattr(alphamissense, "lookup", boom)
    v = _v(info={"am_pathogenicity": "0.9", "am_class": "likely_pathogenic"})
    a = annotate_variant(v, with_alphamissense=False)
    assert a.am_pathogenicity == 0.9   # cheap INFO path is always honoured


def test_add_alphamissense_enriches_candidate(monkeypatch):
    monkeypatch.setattr(alphamissense, "lookup", lambda v: {
        "am_pathogenicity": 0.97, "am_class": "likely_pathogenic", "_source": "mock"})
    a = Annotation()
    add_alphamissense(_v(), a)
    assert a.am_pathogenicity == 0.97 and a.am_class == "likely_pathogenic"
    assert a.source["alphamissense"] == "mock"


def test_add_alphamissense_is_noop_when_present(monkeypatch):
    def boom(v):
        raise AssertionError("must NOT re-query when a score is already present")
    monkeypatch.setattr(alphamissense, "lookup", boom)
    a = Annotation(am_pathogenicity=0.5)
    add_alphamissense(_v(), a)
    assert a.am_pathogenicity == 0.5


def test_pipeline_queries_am_only_for_candidates(monkeypatch):
    calls = []
    monkeypatch.setattr(alphamissense, "lookup", lambda v: (
        calls.append(v.key), {"am_pathogenicity": None, "am_class": None, "_source": "mock"})[1])
    report = pipeline.run_pipeline(config.SAMPLE_VCF)
    n_candidates = len(report.classifications)
    assert n_candidates > 0
    # One lookup per candidate — NOT one per post-QC variant (that is the whole win).
    assert len(calls) == n_candidates
    assert n_candidates <= report.qc.after_qc <= report.qc.total_variants


def test_lazy_am_reaches_every_candidate(monkeypatch):
    """End-to-end: the deferred lookup actually populates each candidate's annotation
    (so a candidate is classified with the SAME score eager annotation would give it)."""
    monkeypatch.setattr(alphamissense, "lookup", lambda v: {
        "am_pathogenicity": 0.999, "am_class": "likely_pathogenic", "_source": "mock"})
    report = pipeline.run_pipeline(config.SAMPLE_VCF)
    assert report.classifications
    for c in report.classifications:
        assert c.annotation.am_pathogenicity == 0.999


def test_lazy_tiers_equal_eager_tiers(monkeypatch):
    """The refactor's central invariant: lazy classifications == eager classifications.

    Uses a PER-VARIANT score keyed on position, so an enrichment that routed the
    wrong variant's AlphaMissense score onto a candidate would break equivalence
    (a constant mock could not catch that).
    """
    from vcf2report.acmg.engine import classify
    from vcf2report.vcf.filter import filter_variants
    from vcf2report.vcf.parse import parse_vcf
    from vcf2report.vcf.qc import apply_qc

    def keyed(v):
        s = 0.999 if (v.pos % 2 == 0) else 0.05
        return {"am_pathogenicity": s, "_source": "mock",
                "am_class": "likely_pathogenic" if s > 0.5 else "likely_benign"}
    monkeypatch.setattr(alphamissense, "lookup", keyed)

    # Lazy = the production pipeline.
    lazy_report = pipeline.run_pipeline(config.SAMPLE_VCF)
    lazy = {c.variant.key: c.tier for c in lazy_report.classifications}

    # Eager reference = annotate the WHOLE post-QC set with the client, then classify.
    variants, _build, _h = parse_vcf(config.SAMPLE_VCF)
    kept, _ = apply_qc(variants)
    annotated = [(v, annotate_variant(v, [], build_trusted=True, with_alphamissense=True))
                 for v in kept]
    candidates, _ = filter_variants(annotated, max_af=config.AF_RECESSIVE_MAX)
    eager = {v.key: classify(v, a).tier for v, a in candidates}

    assert lazy == eager
    # Non-trivial: AlphaMissense actually drove a criterion on at least one candidate.
    assert any(("PP3" in c.met_codes) or ("BP4" in c.met_codes)
               for c in lazy_report.classifications)


def test_default_annotate_variant_queries_client(monkeypatch):
    """The default (with_alphamissense=True) path — used by the MCP 'annotate' tool
    and single-variant callers — must still query the client."""
    calls = []
    monkeypatch.setattr(alphamissense, "lookup", lambda v: (
        calls.append(v.key),
        {"am_pathogenicity": 0.88, "am_class": "ambiguous", "_source": "mock"})[1])
    v = _v()
    a = annotate_variant(v)   # defaults
    assert calls == [v.key]
    assert a.am_pathogenicity == 0.88
