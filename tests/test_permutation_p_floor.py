"""A permutation p-value at its floor must be reported as `< x`, never as a point value.

WHY. The ACS leaderboard computes `null_p = (hits + 1) / (draws + 1)` over 10,000 draws, so
with zero hits it takes the value 1/10001 = 9.999e-05 and CANNOT go lower. **14 of 17 methods
sit exactly there.** Several documents printed that as "p = 0.0001", which reads as a value
measured to four decimal places. What it actually means is that no permutation out of 10,000
reached the observed ACS, so the honest statement is `p < 1e-4` — the test cannot resolve
finer without more draws. The 20,000-draw constraint checks have the same shape at 1/20000.

A reviewer who divides 1 by the draw count finds the floor in seconds, and a point estimate
sitting on it looks like overstatement of precision. This test keeps the convention.
"""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
LB = ROOT / "results" / "anatomic" / "acs_leaderboard.csv"
DOCS = ["README.md", "RESULTS.md", "docs/ENDPOINT.md", "docs/METHOD_LIMITATIONS.md",
        "docs/DATA_SOURCES.md", "docs/PAPER_OUTLINE.md", "docs/FIGURES.md"]


def _summarize():
    spec = importlib.util.spec_from_file_location(
        "summarize_results", ROOT / "scripts" / "summarize_results.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_formatter_marks_the_floor_and_leaves_real_values_alone():
    f = _summarize().fmt_perm_p
    assert f(1 / 10001, 10_000) == "< 1e-4"
    assert f(9.999000099990002e-05, 10_000) == "< 1e-4"
    assert f(0.00005, 20_000) == "< 5e-5"
    # a genuinely measured p is printed as itself
    assert f(0.3087691230876912, 10_000) == "0.3088"
    assert f(0.0143, 10_000) == "0.0143"
    assert f(None, 10_000) == "—"


@pytest.mark.skipif(not LB.exists(), reason="no leaderboard in this checkout")
def test_the_leaderboard_really_is_at_the_floor():
    """Guards the premise. If the null ever gains draws, this test says so."""
    lb = pd.read_csv(LB)
    floor = 1 / 10_001
    at = (lb["null_p"] - floor).abs() < 1e-12
    assert at.any(), (
        "no method sits at 1/10001 any more — the draw count or the p formula changed. "
        "Re-derive the floor before trusting fmt_perm_p's default.")
    assert lb.loc[at, "null_p"].min() >= floor - 1e-12


@pytest.mark.parametrize("rel", DOCS)
def test_no_document_prints_a_floor_p_as_a_point_value(rel):
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} absent")
    bad = []
    for ln, line in enumerate(p.read_text(errors="replace").split("\n"), 1):
        # "p = 0.0001" / "p = 0.00005" and the bare table cell "| 0.0001 |"
        if re.search(r"\bp\s*[=＝]\s*0\.0001\b", line) or \
           re.search(r"\bp\s*[=＝]\s*0\.00005\b", line) or \
           re.search(r"\|\s*0\.0001\s*\|", line):
            bad.append((ln, line.strip()[:100]))
    assert not bad, (
        f"{rel} prints a permutation p at its floor as a point value; use `< 1e-4` "
        f"(10,000 draws) or `< 5e-5` (20,000 draws):\n"
        + "\n".join(f"  line {n}: {s}" for n, s in bad))
