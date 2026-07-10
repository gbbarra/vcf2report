"""Privacy: network egress is opt-in; no patient coordinates leave by default."""
from vcf2report import config
from vcf2report.annotate import _http, gnomad
from vcf2report.models import Variant


def test_network_off_by_default(monkeypatch):
    monkeypatch.delenv("OFFLINE", raising=False)
    monkeypatch.delenv("VCF2REPORT_ALLOW_NETWORK", raising=False)
    assert config.offline() is True           # safe default
    assert config.allow_network() is False

    def _boom(*a, **k):
        raise AssertionError("no network egress is allowed by default")

    monkeypatch.setattr(_http, "post_json", _boom)
    monkeypatch.setattr(_http, "get_json", _boom)
    # A variant absent from the local snapshot must NOT trigger a live call.
    r = gnomad.lookup(Variant(chrom="9", pos=999999, ref="A", alt="G"))
    assert "live" not in r["_source"]


def test_allow_network_enables_and_offline_overrides(monkeypatch):
    monkeypatch.delenv("OFFLINE", raising=False)
    monkeypatch.setenv("VCF2REPORT_ALLOW_NETWORK", "1")
    assert config.allow_network() is True
    monkeypatch.setenv("OFFLINE", "1")          # hard override wins
    assert config.allow_network() is False
