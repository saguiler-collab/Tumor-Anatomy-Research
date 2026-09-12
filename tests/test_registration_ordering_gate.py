"""
The pre-registration ordering verdict must not be published when it is stale.

The failure this pins actually happened. `run_anatomic` called `registration.status()`
partway through the run, before most artefacts had been rewritten, so the mtimes it
compared against the registration were the PREVIOUS run's. The 2026-09-10 run published

    "REGISTERED, but 25 result file(s) are OLDER than the registration. ... Those files
     were produced before the constraints were registered and must be regenerated or
     labelled superseded."

into RESULTS.md, the primary results document, when the true count was zero. A protocol
violation asserted that did not happen is as damaging as one missed: a reader has no way
to tell which kind of error they are looking at.

Two defences, both tested here:
  * the status file is now written LAST, after every other artefact;
  * the summariser refuses to reprint an ordering verdict whose own file is older than the
    results it claims to have judged.
"""
from __future__ import annotations

import importlib.util
import json
import os
import time
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "summarize_results",
    Path(__file__).resolve().parent.parent / "scripts" / "summarize_results.py")
summarize = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summarize)


def _tree(tmp_path: Path, verdict: str, predating: list[str]) -> Path:
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    (anat / "registration_status.json").write_text(json.dumps({
        "state": "REGISTERED", "verdict": verdict,
        "results_predating_registration": predating,
        "freeze_hash": "a" * 64, "recorded_hash": "a" * 64,
    }))
    return anat


def test_a_status_file_older_than_its_results_is_flagged(tmp_path):
    anat = _tree(tmp_path, "REGISTERED, but 25 result file(s) are OLDER...", ["x"] * 25)
    assert summarize._ordering_verdict_is_stale(anat) == 0, "nothing newer yet"

    # A later run writes its artefacts without rewriting the status file.
    time.sleep(0.01)
    later = time.time() + 10
    for name in ("acs_leaderboard.csv", "acs_per_tumor.csv"):
        p = anat / name
        p.write_text("method,acs\nmusic,1.0\n")
        os.utime(p, (later, later))

    assert summarize._ordering_verdict_is_stale(anat) == 2


def test_the_stale_verdict_is_not_reprinted_without_a_warning(tmp_path):
    anat = _tree(tmp_path, "REGISTERED, but 25 result file(s) are OLDER than the "
                           "registration.", ["x"] * 25)
    later = time.time() + 10
    p = anat / "acs_leaderboard.csv"
    p.write_text("method,acs\nmusic,1.0\n")
    os.utime(p, (later, later))

    text = summarize.integrity_section(anat)
    assert "STALE" in text, (
        "a stale ordering verdict was reprinted with no warning; this is the exact text "
        f"that misreported a protocol violation:\n{text[:400]}")
    assert "must not be read" in text


def test_a_current_verdict_is_printed_clean(tmp_path):
    anat = _tree(tmp_path, "REGISTERED: ... every result file postdates it.", [])
    p = anat / "acs_leaderboard.csv"
    p.write_text("method,acs\nmusic,1.0\n")
    old = time.time() - 100
    os.utime(p, (old, old))

    text = summarize.integrity_section(anat)
    assert summarize._ordering_verdict_is_stale(anat) == 0
    assert "STALE" not in text
    assert "every result file postdates it" in text


def test_the_shipped_status_artefact_is_not_stale():
    """The repository's own results tree must not be carrying a stale verdict."""
    anat = Path(__file__).resolve().parent.parent / "results" / "anatomic"
    if not (anat / "registration_status.json").exists():
        pytest.skip("no results tree in this checkout")
    n = summarize._ordering_verdict_is_stale(anat)
    assert n == 0, (
        f"{n} result file(s) are newer than registration_status.json, so its ordering "
        "verdict judged an earlier state of the tree. Recompute it — do not edit it.")
