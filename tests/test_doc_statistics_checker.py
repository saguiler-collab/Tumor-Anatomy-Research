"""The statistics checker must FAIL when a statistic is wrong.

WHY THIS EXISTS. The first version of `scripts/check_doc_statistics.py` harvested every float
in the results tree and asked whether a quoted number appeared anywhere in that set. With
755,490 distinct values, essentially every 3- and 4-decimal number matched something, so it
reported CLEAN when a rho was replaced with an invented number AND when a sign was flipped. It
passed because it could not fail.

Two later bugs in the rewrite were also found only by running the controls, not by reading:

1. The "the current value is quoted alongside it" skip used absolute values, so when current
   and superseded differed ONLY in sign -- exactly OPEN_DEFECTS D21 -- the skip matched the
   flipped value and waved it through.
2. The sign matcher was written `"-|−" + digits`, and alternation binds looser than
   concatenation, so it meant "any hyphen OR the number" and fired on every line with a dash.

So the controls are the test. A checker is only worth its exit code if something makes it
non-zero.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "check_doc_statistics.py"
DOC = ROOT / "docs" / "ENDPOINT.md"


def _run() -> tuple[int, str]:
    r = subprocess.run([sys.executable, str(SCRIPT), "--strict"],
                       capture_output=True, text=True, cwd=ROOT, timeout=300)
    return r.returncode, r.stdout + r.stderr


@pytest.fixture
def doc_restored():
    """Edit ENDPOINT.md in place and always put it back byte-for-byte."""
    original = DOC.read_bytes()
    try:
        yield
    finally:
        DOC.write_bytes(original)


@pytest.mark.skipif(not SCRIPT.exists() or not DOC.exists(), reason="checker or doc absent")
def test_the_checker_is_clean_on_the_shipped_documents():
    code, out = _run()
    assert code == 0, f"the shipped documents do not pass the statistics check:\n{out[-2500:]}"


@pytest.mark.skipif(not SCRIPT.exists() or not DOC.exists(), reason="checker or doc absent")
def test_reverting_a_statistic_to_its_superseded_value_is_caught(doc_restored):
    """Staleness: the commonest failure, and what the log-`X` carryover actually was."""
    t = DOC.read_text()
    needle = "**Spearman rho = 0.6372** between"
    assert needle in t, "anchor moved; update this test rather than deleting it"
    DOC.write_text(t.replace(needle, "**Spearman rho = 0.7501** between", 1))
    code, out = _run()
    assert code != 0, "a statistic reverted to its superseded value was not caught"
    assert "0.7501" in out


@pytest.mark.skipif(not SCRIPT.exists() or not DOC.exists(), reason="checker or doc absent")
def test_a_sign_flip_is_caught_even_though_the_magnitude_is_correct(doc_restored):
    """The D21 shape. The magnitude is right, so only a sign-aware check can see it."""
    t = DOC.read_text()
    needle = "| **nothing** | **+0.0810** (p = 0.80)"
    assert needle in t, "anchor moved; update this test rather than deleting it"
    DOC.write_text(t.replace(needle, "| **nothing** | **-0.0810** (p = 0.80)", 1))
    code, out = _run()
    assert code != 0, (
        "a sign flip on a statistic whose magnitude is correct was not caught — this is "
        "exactly the defect the checker exists for")


@pytest.mark.skipif(not SCRIPT.exists(), reason="checker absent")
def test_every_canonical_statistic_resolves_to_an_artefact_field():
    """A statistic the checker cannot resolve is a statistic it is not checking."""
    code, out = _run()
    assert "ARTEFACT FIELD ABSENT" not in out, (
        "a canonical statistic has no artefact field, so it is silently unchecked:\n"
        + "\n".join(l for l in out.split("\n") if "ABSENT" in l))


@pytest.mark.skipif(not SCRIPT.exists(), reason="checker absent")
def test_the_sign_matcher_does_not_fire_on_every_dash():
    """Guards the alternation-precedence bug directly, not just through its symptom."""
    src = SCRIPT.read_text()
    bad = re.search(r'\(\s*"-\|−"\s+if', src)
    assert bad is None, (
        'the sign matcher is built as "-|−" + digits again; alternation binds looser than '
        'concatenation, so that pattern matches any hyphen. Wrap it: "(?:-|−)".')
