"""The timeout in `_run_bounded` must fire, and large output must not deadlock.

**Read D18's 2026-09-19 correction before trusting the framing here.** These tests were written
to demonstrate a fix for an observed 9-hour hang; checking them against the pre-fix code showed
the pre-fix code passes the timeout cases too. The file-backed capture is a genuine robustness
improvement -- the 2 MB test below is real, and a pipe with no reader does block at ~64 KB -- but
it is NOT established as the cause of the hang.

OPEN_DEFECTS D18. `_run_bounded` put the child in its own session and killed the GROUP on
timeout, which is correct and was still not enough: an R `parallel` socket cluster's workers
survived, re-parented to ppid 1, holding the inherited stdout pipe. The pipe never reached EOF,
so the parent blocked past its budget — a 2,400 s budget ran 93 minutes, and the project's own
synthetic validation gate sat for 8 h 54 min having written nothing.

These tests reproduce that shape with `sh` instead of R: a child that spawns a detached
grandchild which holds the inherited stdout open and sleeps far longer than the budget.
"""
from __future__ import annotations

import subprocess
import time

import pytest

from ivygap.deconv.r_bridge import _run_bounded


def test_a_quick_command_returns_its_output():
    """The ordinary path still works, including output larger than a pipe buffer."""
    got = _run_bounded(["sh", "-c", "echo hello; echo oops >&2"], budget=30)
    assert got.returncode == 0
    assert "hello" in got.stdout
    assert "oops" in got.stderr


def test_output_far_larger_than_a_pipe_buffer_does_not_deadlock():
    """A pipe would block the writer at ~64 KB with nobody draining. A file does not.

    This is the positive control for the fix: 2 MB through the same path.
    """
    got = _run_bounded(
        ["sh", "-c", "i=0; while [ $i -lt 2000 ]; do "
                     "printf '%s\\n' 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'; "
                     "i=$((i+1)); done"],
        budget=60)
    assert got.returncode == 0
    assert len(got.stdout) > 500_000, "large output must survive, not deadlock"


def test_timeout_fires_when_a_detached_grandchild_holds_the_output_handle():
    """Timeout must fire with a detached grandchild holding the output handle.

    **This is NOT a regression test for D18, and saying so would be false.** It was written as
    one, then checked against the pre-fix implementation, which *also* raises at 3.0 s on this
    input. `communicate(timeout=...)` fires correctly even when a grandchild keeps the pipe
    open, so the pipe-EOF story in D18's original wording does not explain the observed hang.
    See D18's 2026-09-19 correction.

    It is kept because it pins behaviour that must not regress, not because it demonstrates a
    fix.
    """
    start = time.monotonic()
    with pytest.raises(subprocess.TimeoutExpired):
        # `sleep 120` is detached into its own session, so killpg on the child's group need not
        # reach it -- which is precisely the situation the R cluster workers created.
        _run_bounded(["sh", "-c", "setsid sleep 120 & sleep 120"], budget=3)
    took = time.monotonic() - start
    # generous ceiling: the point is that it returns in seconds, not in two minutes
    assert took < 45, f"timeout took {took:.1f}s against a 3s budget — it did not fire"


def test_a_timeout_still_reports_what_the_child_managed_to_write():
    """OPEN_DEFECTS D19: a fallback whose reason is unrecorded cannot be acted on.

    The raised TimeoutExpired must carry the partial output, so "timed out" is distinguishable
    from "raised" without re-running anything.
    """
    with pytest.raises(subprocess.TimeoutExpired) as ei:
        _run_bounded(["sh", "-c", "echo progress-marker; sleep 120"], budget=3)
    text = ei.value.output or ""
    if isinstance(text, bytes):                            # pragma: no cover
        text = text.decode("utf-8", "replace")
    assert "progress-marker" in text, "partial output must survive the timeout"
