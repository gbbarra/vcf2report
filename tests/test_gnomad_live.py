"""The gnomAD remote path, against the real gnomAD. No fixture, no fake, no mock.

Why this file exists as its own module, with real network:

`gnomad_remote` returned None for EVERY gnomAD v4.1 variant, because `rec.info.get(k)` raises
`ValueError("Invalid header")` for a key the header never declares — while returning None for a
key it declares and leaves empty. The frequency was read correctly and thrown away on the way
out. Nothing in the suite noticed, because every test of that module used a hand-written stand-in
for pysam, and a stand-in can only ever assert its author's model of the library.

So these tests range-query `storage.googleapis.com/gcp-public-data--gnomad` with real pysam and
assert the real numbers. They are slow (a few seconds per chromosome handle) and they need the
network — which is the point. In CI they run in the `tabix` job with
``VCF2REPORT_LIVE_GNOMAD=1``, and a network failure fails the job rather than skipping, because
a skipped network test is the artefact this file exists to stop shipping.

Locally they are opt-in (the env var), so a developer offline is not blocked.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from vcf2report.models import Variant

pytestmark = pytest.mark.skipif(
    os.environ.get("VCF2REPORT_LIVE_GNOMAD", "").strip() not in {"1", "true", "yes"},
    reason="live gnomAD test — set VCF2REPORT_LIVE_GNOMAD=1 (CI's tabix job does)",
)

pysam = pytest.importorskip("pysam", reason="the live path IS pysam")

# BBS2 c.472-2A>G, the planted diagnosis in data/example/SYN-073. The committed fixture baked
# gnomad_AF=8.99282e-07 / AN=1461862 / nhomalt=0 from the full store on another machine, so this
# doubles as an independent check that the fixture is faithful.
URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/"
    "gnomad.exomes.v4.1.sites.chr16.vcf.bgz"
)
BBS2 = Variant(chrom="chr16", pos=56510923, ref="T", alt="C")
BBS2_AF_GRPMAX = 8.9928198576672e-07


@pytest.fixture(scope="module", autouse=True)
def _network_or_fail():
    """Fail loudly, not skip, when the bucket is unreachable — probing THE SAME TRANSPORT.

    The first version of this used urllib. urllib rides Python's ssl module and the system CA
    store; pysam rides htslib's own libcurl, which on a GitHub runner could not find a CA
    bundle at all (``Libcurl reported error 77``). So the pre-check passed while every test
    failed — a guard that verified a different transport from the one under test, which is the
    same class of mistake as asserting against a hand-written stand-in.
    """
    try:
        vf = pysam.VariantFile(URL)
        next(iter(vf.fetch("chr16", BBS2.pos - 1, BBS2.pos)))
    except Exception as e:  # noqa: BLE001
        pytest.fail(
            f"gnomAD public bucket not readable BY PYSAM: {type(e).__name__}: {e}\n"
            "If this is 'Libcurl reported error 77', htslib cannot find a CA bundle — set "
            "CURL_CA_BUNDLE / SSL_CERT_FILE (the tabix CI job does)."
        )


def test_live_query_returns_the_real_frequency():
    """The end-to-end assertion the fake could not make: query() gives a number, not None."""
    from vcf2report.annotate import gnomad_remote

    got = gnomad_remote.query(BBS2)
    assert got is not None, (
        "query() returned None for a variant gnomAD publishes — the failure mode that shipped"
    )
    assert got["af"] == pytest.approx(BBS2_AF_GRPMAX, rel=1e-6)
    assert got["pop"] == "nfe"
    assert got["hom"] == 0


def test_the_live_value_matches_the_committed_fixture():
    """data/example/SYN-073 baked its gnomAD numbers from the full store. If the live release
    disagrees, either the fixture is stale or the reduction changed — both worth knowing."""
    from vcf2report.annotate import gnomad_remote

    assert gnomad_remote.query(BBS2)["af"] == pytest.approx(8.99282e-07, rel=1e-4)


def test_pysam_raises_on_an_undeclared_info_key():
    """Pin the library behaviour the whole bug rests on, against the real library and a real
    record — so the offline stand-in in test_annotate_absence.py can be checked against it
    rather than trusted."""
    rec = next(iter(pysam.VariantFile(URL).fetch("chr16", BBS2.pos - 1, BBS2.pos)))

    # declared and populated -> a value
    assert rec.info.get("AF_grpmax")[0] == pytest.approx(BBS2_AF_GRPMAX)
    # declared and empty on this record -> None, NOT a raise
    assert rec.info.get("fafmax_faf95_max") is None
    # never declared in this header -> raises. This asymmetry is the bug.
    for undeclared in ("faf95_grpmax", "faf95_max"):
        with pytest.raises(ValueError):
            rec.info.get(undeclared)


def test_the_offline_stand_in_matches_the_real_library():
    """The offline test uses a fake pysam record. This is what keeps the fake honest.

    A stand-in that drifts from the library is worse than no test: it stays green while the
    behaviour it models changes underneath. Here the same probes run against both, and any
    divergence fails.
    """
    # Load the sibling module BY PATH. `from tests.… import` needs `tests` to be an importable
    # package, which it is not — it has no __init__.py, and it only appeared to work locally
    # because pytest had inserted the repo root. In CI the whole live module failed at import,
    # before touching the network, so the job proved nothing about gnomAD at all.
    import importlib.util

    path = Path(__file__).with_name("test_annotate_absence.py")
    spec = importlib.util.spec_from_file_location("_absence_mod", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _FakeInfo = mod._FakeInfo

    real = next(
        iter(pysam.VariantFile(URL).fetch("chr16", BBS2.pos - 1, BBS2.pos))
    ).info
    fake = _FakeInfo(
        declared={"AF_grpmax": (BBS2_AF_GRPMAX,), "fafmax_faf95_max": None},
        undeclared={"faf95_grpmax", "faf95_max"},
    )

    def behaviour(info, key):
        try:
            return ("value", info.get(key))
        except Exception as e:  # noqa: BLE001
            return ("raises", type(e).__name__)

    for key in ("AF_grpmax", "fafmax_faf95_max", "faf95_grpmax", "faf95_max"):
        assert behaviour(real, key) == behaviour(fake, key), (
            f"the offline stand-in no longer matches real pysam for {key!r} — "
            "fix the stand-in, or the offline test is proving nothing"
        )
