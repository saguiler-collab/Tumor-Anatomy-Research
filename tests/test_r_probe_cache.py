"""
The R package probe must be cached, and must not report a failed probe as an absence.

Two failures this pins, both observed:

  * An uncached probe spawned a fresh R process on every `check()`. A pytest run spent
    over twenty minutes almost entirely in repeated `requireNamespace` calls while the
    machine was busy with a re-measurement.
  * The old probe caught every exception and returned False, which `check()` reported as
    "R package X is not installed". Under load a probe can time out while the package is
    installed and usable. The method then runs its Python reimplementation and the
    recorded reason is false — and the leaderboard becomes a function of machine load.
"""
from __future__ import annotations

import subprocess

import pytest

from ivygap.deconv import r_bridge


@pytest.fixture(autouse=True)
def _clean_cache():
    r_bridge.clear_package_probe_cache()
    yield
    r_bridge.clear_package_probe_cache()


def test_probe_runs_once_per_package(monkeypatch):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="yes", stderr="")

    monkeypatch.setattr(r_bridge.shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.setattr(r_bridge.subprocess, "run", fake_run)

    assert r_bridge.r_package_available("MuSiC") is True
    for _ in range(20):
        assert r_bridge.r_package_available("MuSiC") is True
    assert len(calls) == 1, f"probed {len(calls)} times; the cache is not holding"

    r_bridge.r_package_available("SCDC")
    assert len(calls) == 2, "a different package must still be probed"


def test_absent_package_says_absent(monkeypatch):
    monkeypatch.setattr(r_bridge.shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.setattr(r_bridge.subprocess, "run",
                        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, "no", ""))
    ok, why = r_bridge.probe_r_package("NotARealPackage")
    assert ok is False
    assert "not installed" in why
    assert "PROBE FAILURE" not in why


def test_timeout_is_not_reported_as_absence(monkeypatch):
    """The failure that would make the leaderboard depend on machine load."""
    def boom(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, 120)

    monkeypatch.setattr(r_bridge.shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.setattr(r_bridge.subprocess, "run", boom)

    ok, why = r_bridge.probe_r_package("MuSiC")
    assert ok is False, "a method must not run when the probe could not confirm the package"
    assert "PROBE FAILURE" in why, (
        "a timed-out probe reported as an absence puts a false reason into "
        f"implementation_report.json; got: {why}")
    assert "not installed" not in why.replace("not evidence the package is absent", "")


def test_check_propagates_the_probe_reason(monkeypatch):
    def boom(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, 120)

    monkeypatch.setattr(r_bridge.shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.setattr(r_bridge.subprocess, "run", boom)

    method = next(iter(r_bridge.R_PACKAGES))
    avail = r_bridge.check(method, "gbmap")
    assert avail.available is False
    assert "PROBE FAILURE" in avail.reason, (
        "check() must carry the probe's reason through to the recorded availability, "
        f"got: {avail.reason}")


def test_an_unexpected_probe_result_is_not_silently_an_absence(monkeypatch):
    monkeypatch.setattr(r_bridge.shutil, "which", lambda _: "/usr/bin/Rscript")
    monkeypatch.setattr(r_bridge.subprocess, "run",
                        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "boom"))
    ok, why = r_bridge.probe_r_package("MuSiC")
    assert ok is False
    assert "PROBE FAILURE" in why
