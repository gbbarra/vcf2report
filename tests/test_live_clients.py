"""Live gnomAD / ClinVar clients: parsing, matching, and fallback invariants.

The HTTP layer (annotate._http) is mocked with realistic payloads — no network.
These tests lock the response parsing and the cache-first / OFFLINE / local-
fallback behaviour that the clients must preserve.
"""
from vcf2report.annotate import _http, clinvar, gnomad, gnomad_remote
from vcf2report.models import Variant

# Cache isolation + offline-by-default come from tests/conftest.py; the live
# tests below opt into "online" via _online() (network is opt-in by default).


def _online(mp):
    mp.setenv("OFFLINE", "")
    mp.setenv("VCF2REPORT_ALLOW_NETWORK", "1")
    # These tests exercise the GraphQL + local-fallback layer specifically, so
    # disable the remote-tabix path (covered by test_gnomad_remote.py) to keep
    # them hermetic — otherwise gnomad.lookup would try the real GCS bucket first.
    mp.setattr(gnomad_remote, "query", lambda variant: None)

# --- gnomAD -----------------------------------------------------------------
_GNOMAD_PAYLOAD = {
    "data": {"variant": {
        "variant_id": "2-166003360-C-T",
        "genome": None,
        "exome": {
            "ac": 3, "an": 152000, "homozygote_count": 1,
            "populations": [
                {"id": "nfe", "ac": 3, "an": 68000},
                {"id": "afr", "ac": 0, "an": 40000},
                {"id": "amr", "ac": 0, "an": 20000},
                {"id": "sas", "ac": 5, "an": 10000},     # excluded from popmax
                {"id": "nfe_bgr", "ac": 2, "an": 5000},  # sub-pop, excluded
            ],
        },
    }},
}


def test_gnomad_live_popmax(monkeypatch):
    _online(monkeypatch)
    monkeypatch.setattr(_http, "post_json", lambda *a, **k: _GNOMAD_PAYLOAD)
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T")
    r = gnomad.lookup(v)
    # sas is a full continental group INCLUDED in popmax; nfe_bgr (underscore) is
    # a sub-population and excluded.
    assert r["pop"] == "sas"
    assert abs(r["af"] - 5 / 10000) < 1e-12
    assert r["hom"] == 1
    assert "live" in r["_source"]


def test_gnomad_errors_block_falls_back(monkeypatch):
    """A 200 carrying a non-'not found' GraphQL error must fall back, not fake AF 0."""
    _online(monkeypatch)
    monkeypatch.setattr(_http, "post_json", lambda *a, **k: {
        "data": {"variant": None},
        "errors": [{"message": "Query timed out. Please try again."}],
    })
    v = Variant(chrom="2", pos=178562809, ref="G", alt="A")   # TTN, in local snapshot
    r = gnomad.lookup(v)
    assert r["af"] == 0.081
    assert "local snapshot" in r["_source"]


def test_gnomad_unknown_when_lookup_fails_and_no_local(monkeypatch):
    """Transport failure with no local record -> AF None (unknown), not absent."""
    _online(monkeypatch)
    monkeypatch.setattr(_http, "post_json", lambda *a, **k: None)
    r = gnomad.lookup(Variant(chrom="9", pos=999999, ref="A", alt="G"))
    assert r["af"] is None
    assert "unavailable" in r["_source"]


def test_gnomad_not_found_is_absent(monkeypatch):
    _online(monkeypatch)
    monkeypatch.setattr(_http, "post_json", lambda *a, **k: {"data": {"variant": None}})
    r = gnomad.lookup(Variant(chrom="1", pos=1, ref="A", alt="T"))
    assert r["af"] == 0.0


def test_gnomad_transport_error_falls_back_to_local(monkeypatch):
    _online(monkeypatch)
    monkeypatch.setattr(_http, "post_json", lambda *a, **k: None)  # network failure
    v = Variant(chrom="2", pos=178562809, ref="G", alt="A")       # TTN, in local snapshot
    r = gnomad.lookup(v)
    assert r["af"] == 0.081
    assert "local snapshot" in r["_source"]


def test_gnomad_offline_never_calls_network(monkeypatch):
    monkeypatch.setenv("OFFLINE", "1")

    def _boom(*a, **k):
        raise AssertionError("network must not be called in OFFLINE mode")

    monkeypatch.setattr(_http, "post_json", _boom)
    r = gnomad.lookup(Variant(chrom="2", pos=178562809, ref="G", alt="A"))
    assert r["af"] == 0.081


# --- ClinVar ----------------------------------------------------------------
_ESEARCH = {"esearchresult": {"idlist": ["12345"]}}
_ESUMMARY_MATCH = {"result": {"uids": ["12345"], "12345": {
    "accession": "VCV000012345",
    "germline_classification": {
        "description": "Pathogenic",
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "last_evaluated": "2024/11/01",
        "trait_set": [{"trait_name": "Dravet syndrome"}],
    },
    "variation_set": [{"variation_loc": [{
        "assembly_name": "GRCh38", "chr": "2", "start": "166003360",
        "stop": "166003360", "ref": "C", "alt": "T"}]}],
}}}


def _clinvar_router(match=True):
    summary = _ESUMMARY_MATCH if match else {"result": {"uids": ["12345"], "12345": {
        "accession": "VCV999", "variation_set": [{"variation_loc": [{
            "assembly_name": "GRCh38", "chr": "2", "start": "999999", "alt": "G"}]}]}}}

    def router(url, params, **k):
        return _ESEARCH if "esearch" in url else summary
    return router


def test_clinvar_live_exact_match(monkeypatch):
    _online(monkeypatch)
    monkeypatch.setattr(_http, "get_json", _clinvar_router(match=True))
    monkeypatch.setattr(_http, "throttle", lambda *a, **k: None)
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    r = clinvar.lookup(v)
    assert r["significance"] == "Pathogenic"
    assert r["accession"] == "VCV000012345"
    assert "Dravet" in (r["condition"] or "")


def test_clinvar_live_mismatch_falls_back_to_local(monkeypatch):
    _online(monkeypatch)
    # Live returns a record at a different position/allele -> must not be accepted;
    # falls back to the authoritative local slice (which has SCN1A Pathogenic).
    monkeypatch.setattr(_http, "get_json", _clinvar_router(match=False))
    monkeypatch.setattr(_http, "throttle", lambda *a, **k: None)
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    r = clinvar.lookup(v)
    assert r["significance"] == "Pathogenic"
    assert "slice" in r["_source"]


def test_clinvar_legacy_field_shape():
    docsum = {"accession": "VCV1", "clinical_significance": {
        "description": "Benign", "review_status": "single submitter"}}
    out = clinvar._extract(docsum)
    assert out["significance"] == "Benign"
    assert out["review_status"] == "single submitter"


def test_clinvar_empty_search_falls_back_to_local(monkeypatch):
    """esearch with no hits must fall back to the local slice, not cache a miss."""
    _online(monkeypatch)
    monkeypatch.setattr(_http, "throttle", lambda *a, **k: None)
    monkeypatch.setattr(_http, "get_json",
                        lambda url, params, **k: {"esearchresult": {"idlist": []}})
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    r = clinvar.lookup(v)
    assert r["significance"] == "Pathogenic"   # from local slice
    assert "slice" in r["_source"]


def test_clinvar_rejects_same_pos_unconfirmed_allele(monkeypatch):
    """Same position but no confirmable ref/alt must be REJECTED (not attached)."""
    _online(monkeypatch)
    monkeypatch.setattr(_http, "throttle", lambda *a, **k: None)
    summary = {"result": {"uids": ["1"], "1": {
        "accession": "VCV_OTHER",
        "germline_classification": {"description": "Benign"},
        # matching chr/pos but NO ref/alt and NO canonical_spdi -> cannot confirm
        "variation_set": [{"variation_loc": [{
            "assembly_name": "GRCh38", "chr": "2", "start": "166003360"}]}],
    }}}
    monkeypatch.setattr(_http, "get_json", lambda url, params, **k: (
        {"esearchresult": {"idlist": ["1"]}} if "esearch" in url else summary))
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    r = clinvar.lookup(v)
    # Rejected the unconfirmed record -> fell back to authoritative local slice.
    assert r["accession"] == "VCV000012345"
    assert "slice" in r["_source"]


def test_clinvar_matches_via_canonical_spdi(monkeypatch):
    """Allele can be positively confirmed from canonical_spdi when loc lacks alt."""
    _online(monkeypatch)
    monkeypatch.setattr(_http, "throttle", lambda *a, **k: None)
    summary = {"result": {"uids": ["1"], "1": {
        "accession": "VCV000012345",
        "germline_classification": {"description": "Pathogenic"},
        "variation_set": [{
            "canonical_spdi": "NC_000002.12:166003359:C:T",
            "variation_loc": [{"assembly_name": "GRCh38", "chr": "2",
                               "start": "166003360"}]}],
    }}}
    monkeypatch.setattr(_http, "get_json", lambda url, params, **k: (
        {"esearchresult": {"idlist": ["1"]}} if "esearch" in url else summary))
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A")
    r = clinvar.lookup(v)
    assert r["significance"] == "Pathogenic"
    assert "live" in r["_source"]
