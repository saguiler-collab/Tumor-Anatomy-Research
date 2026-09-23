"""The reference list must describe itself accurately.

`docs/DATA_SOURCES.md` ends in the paper's reference list and states, in prose, how many
entries are still unverified. That count is hand-typed. A reference list whose own summary
disagrees with its contents is the kind of thing a reviewer finds in thirty seconds, and it
costs nothing to make it impossible: the summary names both a number and the entries, so
both are checkable against the markers themselves.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

DOC = Path(__file__).resolve().parent.parent / "docs" / "DATA_SOURCES.md"
MARKER = "[VERIFY]"


def _lines() -> list[str]:
    return DOC.read_text(errors="replace").split("\n")


def _summary_line(lines: list[str]) -> tuple[int, str]:
    hits = [(i, ln) for i, ln in enumerate(lines)
            if re.search(r"Outstanding:\s*\d+", ln)]
    assert len(hits) == 1, f"expected exactly one 'Outstanding: N' summary, found {len(hits)}"
    return hits[0]


@pytest.mark.skipif(not DOC.exists(), reason="DATA_SOURCES.md absent")
def test_the_outstanding_count_matches_the_markers():
    lines = _lines()
    idx, summary = _summary_line(lines)
    stated = int(re.search(r"Outstanding:\s*(\d+)", summary).group(1))
    # every marker EXCEPT the one inside the summary sentence itself
    actual = sum(MARKER in ln for i, ln in enumerate(lines) if i != idx)
    assert stated == actual, (
        f"{DOC.name} says 'Outstanding: {stated}' but carries {actual} {MARKER} "
        f"entries. Update the summary, or verify an entry and remove its marker.")


@pytest.mark.skipif(not DOC.exists(), reason="DATA_SOURCES.md absent")
def test_no_panel_method_is_left_unverified():
    """The summary's load-bearing clause is 'none of them a panel method'.

    That is the claim that matters: an unverified citation for a DATASET is untidy, an
    unverified citation for a method whose numbers are in the leaderboard is a hole in the
    paper. So it is tested, not asserted in prose.
    """
    lines = _lines()
    idx, _ = _summary_line(lines)
    panel = ("music", "scdc", "bisque", "epic", "quantiseq", "dwls", "bayesprism",
             "cibersortx", "nnls", "svr", "elastic net", "elastic_net")
    offenders = []
    for i, ln in enumerate(lines):
        if i == idx or MARKER not in ln:
            continue
        low = ln.lower()
        hit = [p for p in panel if re.search(rf"\b{re.escape(p)}\b", low)]
        if hit:
            offenders.append((i + 1, sorted(set(hit)), ln.strip()[:70]))
    assert not offenders, (
        "a panel method's citation is still unverified, which the summary denies: "
        + "; ".join(f"line {n} names {h}: {s}" for n, h, s in offenders))


@pytest.mark.skipif(not DOC.exists(), reason="DATA_SOURCES.md absent")
def test_the_summary_names_every_outstanding_entry():
    """It claims to list them. If it lists fewer than it counts, the list is the stale half."""
    lines = _lines()
    idx, summary = _summary_line(lines)
    # the sentence may wrap, so read to the end of the paragraph
    tail = []
    for ln in lines[idx:]:
        if not ln.strip():
            break
        tail.append(ln)
    named = " ".join(tail)
    stated = int(re.search(r"Outstanding:\s*(\d+)", summary).group(1))
    # count the em-dash-separated roll call, which is comma/'and'-separated
    after = named.split("—", 1)
    assert len(after) == 2, "the summary should name the outstanding entries after an em dash"
    roll = [t.strip(" .*") for t in re.split(r",|\band\b", after[1]) if t.strip(" .*")]
    assert len(roll) == stated, (
        f"the summary counts {stated} outstanding entries but names {len(roll)}: {roll}")
