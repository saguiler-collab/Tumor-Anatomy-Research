"""A quantity that differs between cohorts must not be reported as one number.

FOUND 2026-09-28. EpiDISH ships 333 reference CpGs and complete-case filtering keeps a
different subset in each matrix: **258 in GBM, 255 in LGG**. The manuscript's Methods row
covers BOTH cohorts and printed 258; `docs/DATA_SOURCES.md` covered the same span and printed
255. Each was correct for one cohort and wrong for the other, and neither said which.

The error is easy to make and invisible once made, because both numbers are real and both
come from an artefact. What makes it wrong is scope: the sentence spans two cohorts and the
number spans one. This test pins the CpG case and guards the general habit -- a document may
state a cohort-varying count only if it names the cohort or gives both.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GBM = ROOT / "results" / "methylation_celltypes.json"
LGG = ROOT / "results" / "methylation_celltypes_lgg.json"
#: Documents that describe BOTH cohorts, so a bare single count is out of scope there.
SPANNING = ["docs/MANUSCRIPT.md", "docs/DATA_SOURCES.md"]
#: Documents scoped to one cohort, where a single count is correct.
SCOPED = ["docs/IMMUNE_ARM.md", "docs/PAPER_OUTLINE.md"]


@pytest.fixture(scope="module")
def counts() -> tuple[int, int]:
    if not (GBM.exists() and LGG.exists()):
        pytest.skip("methylation artefacts absent")
    g = json.loads(GBM.read_text())["n_reference_cpgs_used"]
    l = json.loads(LGG.read_text())["n_reference_cpgs_used"]
    return int(g), int(l)


def test_the_two_cohorts_really_do_differ(counts):
    """Guards the premise. If they ever match, a single number becomes correct again."""
    g, l = counts
    assert g != l, (
        f"GBM and LGG now use the same CpG count ({g}). A single number is then correct and "
        f"this test's reasoning no longer applies — relax it deliberately, do not delete it.")


@pytest.mark.parametrize("rel", SPANNING)
def test_a_cohort_spanning_document_gives_both_counts_or_names_the_cohort(rel, counts):
    g, l = counts
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} absent")
    bad = []
    for ln, line in enumerate(p.read_text(errors="replace").split("\n"), 1):
        if "333" not in line:
            continue
        has_g, has_l = str(g) in line, str(l) in line
        names_cohort = re.search(r"\bGBM\b|\bLGG\b|glioblastoma|lower-grade", line, re.I)
        # correct if it gives both, or gives one and says which cohort
        if (has_g and has_l) or ((has_g or has_l) and names_cohort):
            continue
        if has_g or has_l:
            bad.append((ln, line.strip()[:110]))
    assert not bad, (
        f"{rel} states a CpG count that differs by cohort (GBM {g}, LGG {l}) without giving "
        f"both or naming the cohort:\n"
        + "\n".join(f"  line {n}: {s}" for n, s in bad))


@pytest.mark.parametrize("rel", SCOPED)
def test_a_cohort_scoped_document_still_names_its_cohort(rel, counts):
    """A single count is fine here, but only because the cohort is stated nearby."""
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} absent")
    text = p.read_text(errors="replace")
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if "333" not in line or not any(str(c) in line for c in counts):
            continue
        window = " ".join(lines[max(0, i - 4):i + 5])
        assert re.search(r"\bGBM\b|\bLGG\b|glioblastoma|lower-grade", window, re.I), (
            f"{rel}:{i+1} gives a cohort-varying CpG count with no cohort named within four "
            f"lines: {line.strip()[:110]}")
