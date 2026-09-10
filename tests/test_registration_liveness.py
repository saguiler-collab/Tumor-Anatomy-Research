"""
A registration receipt must be checked against reality, not just against a hash.

This exists because of a specific, self-inflicted failure. On 2026-09-10 the project
recorded a receipt for https://osf.io/nj78m, reported REGISTERED, and stamped that into
every artefact a run produced. The registration had been deleted when its approval was
rejected — the URL returned HTTP 410 — and nothing noticed, because `status()` only ever
compared two hashes. The result was an automatic, written, false claim of
pre-registration in every result file. These tests make that state unreachable.
"""
from __future__ import annotations

import json

import pytest

from ivygap.anatomic import constraints as K
from ivygap.anatomic import registration as R


@pytest.fixture
def receipt(tmp_path, monkeypatch):
    """A well-formed receipt with the correct live hash, pointed at a given URL."""
    def _make(url: str):
        path = tmp_path / "REGISTRATION.json"
        path.write_text(json.dumps({
            "registry": "OSF", "url": url,
            "registered_utc": "2026-09-10T00:00:00Z",
            "constraint_freeze_hash": K.freeze_hash(),
        }))
        monkeypatch.setattr(R, "REGISTRATION_PATH", path)
        return path
    return _make


def _stub(monkeypatch, state, why="stubbed"):
    monkeypatch.setattr(R, "check_url_live", lambda url, timeout=15.0: (state, why))


def test_a_withdrawn_registration_is_never_reported_as_registered(receipt, monkeypatch):
    """The exact failure. 410 must not read as REGISTERED."""
    receipt("https://osf.io/nj78m")
    _stub(monkeypatch, "WITHDRAWN", "HTTP 410 — deleted")
    s = R.status(check_results=False)
    assert s.state == "WITHDRAWN"
    assert not s.is_registered
    assert "NOT registered" in s.verdict


def test_a_missing_registration_is_not_reported_as_registered(receipt, monkeypatch):
    receipt("https://osf.io/zzzzz")
    _stub(monkeypatch, "UNREACHABLE", "HTTP 404 — no such registration")
    s = R.status(check_results=False)
    assert s.state == "UNREACHABLE" and not s.is_registered


def test_a_pending_registration_is_its_own_state(receipt, monkeypatch):
    """
    Submitted but not public is neither REGISTERED nor UNREGISTERED. The timestamp is
    fixed; what is missing is the readability that makes it verifiable.
    """
    receipt("https://osf.io/dm2t8")
    _stub(monkeypatch, "PENDING", "HTTP 401 — not public yet")
    s = R.status(check_results=False)
    assert s.state == "PENDING"
    assert not s.is_registered, "PENDING must not satisfy is_registered"
    assert "not public" in s.verdict or "publicly" in s.verdict


def test_a_live_public_registration_still_passes(receipt, monkeypatch):
    receipt("https://osf.io/aaaaa")
    _stub(monkeypatch, "REGISTERED", "HTTP 200 — public")
    s = R.status(check_results=False)
    assert s.state == "REGISTERED" and s.is_registered


def test_an_unanswerable_liveness_check_does_not_downgrade(receipt, monkeypatch):
    """
    Offline is not evidence of withdrawal. A check that cannot run must leave the
    recorded state alone rather than silently invalidating a real registration.
    """
    receipt("https://osf.io/aaaaa")
    _stub(monkeypatch, None, "no network")
    s = R.status(check_results=False)
    assert s.state == "REGISTERED"


def test_hash_mismatch_still_outranks_liveness(receipt, monkeypatch, tmp_path):
    """A changed constraint file is fatal regardless of whether the URL resolves."""
    path = tmp_path / "REGISTRATION.json"
    path.write_text(json.dumps({
        "registry": "OSF", "url": "https://osf.io/aaaaa",
        "registered_utc": "2026-09-10T00:00:00Z",
        "constraint_freeze_hash": "0" * 64,
    }))
    monkeypatch.setattr(R, "REGISTRATION_PATH", path)
    _stub(monkeypatch, "REGISTERED", "HTTP 200")
    s = R.status(check_results=False)
    assert s.state == "HASH_MISMATCH"


def test_check_live_can_be_switched_off(receipt, monkeypatch):
    """Offline runs and tests must not require the network."""
    receipt("https://osf.io/nj78m")
    def boom(*a, **k):
        raise AssertionError("network must not be touched when check_live=False")
    monkeypatch.setattr(R, "check_url_live", boom)
    s = R.status(check_results=False, check_live=False)
    assert s.state == "REGISTERED"


@pytest.mark.parametrize("code,expect", [
    (200, "REGISTERED"), (401, "PENDING"), (403, "PENDING"),
    (404, "UNREACHABLE"), (410, "WITHDRAWN"),
])
def test_http_codes_map_to_the_right_states(code, expect):
    assert R.LIVENESS[code][0] == expect


def test_a_non_osf_url_is_not_guessed_at():
    state, why = R.check_url_live("https://example.com/whatever")
    assert state is None and "not an OSF" in why
