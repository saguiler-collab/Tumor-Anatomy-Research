"""docs/REFERENCES.md must describe itself accurately.

The reference list is the last thing a reviewer checks and the first place a small error is
fatal to trust. Three things are checkable without judgement and therefore should be:
numbering that is contiguous and starts at 1, a self-reported unverified count that matches
the markers, and the guarantee that no method in the leaderboard panel rests on a citation
nobody has opened. The first draft of this file said "Six entries are unverified" and carried
eight, which is exactly the class of error this catches.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOC = ROOT / "docs" / "REFERENCES.md"
#: The published tools in the leaderboard panel, by reference number in this document.
PANEL_REFS = set(range(13, 21))


def _entries(text: str) -> list[tuple[int, str, str]]:
    """(number, citation, status) for every numbered row."""
    out = []
    for m in re.finditer(r"\*\*\[(\d+)\]\*\*\s*\|(.+?)\|\s*([^|]+?)\s*\|", text):
        out.append((int(m.group(1)), m.group(2).strip(), m.group(3).strip()))
    return out


@pytest.fixture(scope="module")
def doc() -> str:
    if not DOC.exists():
        pytest.skip("docs/REFERENCES.md absent")
    return DOC.read_text()


def test_numbering_is_contiguous_from_one(doc):
    nums = [n for n, _, _ in _entries(doc)]
    assert nums, "no numbered references parsed — the table format changed"
    assert nums == list(range(1, len(nums) + 1)), (
        f"reference numbers are not contiguous from 1: {nums}")


def test_the_stated_unverified_count_matches_the_markers(doc):
    actual = [n for n, _, s in _entries(doc) if "unverified" in s.lower()]
    stated = re.search(r"\*\*(\d+) unverified", doc)
    assert stated, "the Outstanding section no longer states a count"
    assert int(stated.group(1)) == len(actual), (
        f"REFERENCES.md says {stated.group(1)} unverified but carries {len(actual)}: {actual}")
    prose = re.search(r"(\d+) entries are unverified", doc)
    assert prose and int(prose.group(1)) == len(actual), (
        "the header count disagrees with the markers")


def test_no_panel_method_rests_on_an_unverified_citation(doc):
    """The load-bearing guarantee. A benchmark whose methods cite unchecked papers is weak."""
    bad = [n for n, _, s in _entries(doc)
           if n in PANEL_REFS and "unverified" in s.lower()]
    assert not bad, (
        f"reference(s) {bad} are panel methods AND unverified. Every method in the "
        f"leaderboard must cite a paper someone has actually opened.")


def test_every_entry_carries_a_status_and_a_year(doc):
    for n, cite, status in _entries(doc):
        assert status, f"[{n}] has no status"
        assert re.search(r"\b(19|20)\d{2}\b", cite), f"[{n}] has no year: {cite[:70]}"


def test_the_load_bearing_table_names_only_real_reference_numbers(doc):
    nums = {n for n, _, _ in _entries(doc)}
    seg = doc[doc.index("load-bearing"):] if "load-bearing" in doc else ""
    for m in re.finditer(r"\*\*\[(\d+)\]\*\*", seg):
        assert int(m.group(1)) in nums, (
            f"the load-bearing table cites [{m.group(1)}], which is not in the list")
