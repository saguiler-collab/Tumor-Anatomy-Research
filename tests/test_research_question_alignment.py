"""The manuscript must answer the question that was registered, and say which clause it is on.

WHY THIS EXISTS. The registered question in `Anatomy_Test.md` has two clauses:

    "Can a tumor's own anatomy stand in for ground truth when choosing a cell-type
     deconvolution method -- AND does a method that gets the anatomy right also get the
     biology right?"

They are different questions and this study answers them differently. Clause 1 is answered
`it detects, it does not rank`, with the ranking half INCONCLUSIVE at n = 12. Clause 2 is
answered `no`, categorically, by the pre-registered T > B criterion.

Before 2026-09-24 the manuscript answered clause 2 only through the clause-1 correlation --
its weakest evidence -- and titled itself with a claim about deconvolution's capability rather
than about the Anatomy Test. Both are supported, but a reader comparing registration against
report sees a study that asked one thing and reported another. The headline was NOT promoted
post hoc (the T > B criterion carries a timestamp 95 minutes ahead of its data), and the
manuscript has to show that rather than rely on the reader checking.

These tests fail if that alignment is lost again.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
MS = ROOT / "docs" / "MANUSCRIPT.md"
RQ = ROOT / "docs" / "RESEARCH_QUESTION.md"
AVB = ROOT / "results" / "anatomy_vs_biology.json"
PRESPEC = ROOT / "prespecified" / "immune_failure_factors.md"


@pytest.mark.skipif(not MS.exists(), reason="manuscript absent")
def test_the_manuscript_names_both_clauses_of_the_registered_question():
    t = MS.read_text()
    assert "registered question has two clauses" in t, (
        "the manuscript no longer tells the reader the registered question has two clauses; "
        "without it, Result 1 reads as the whole answer when it is half of one")
    for clause in ("stand in for ground truth when choosing",
                   "does a method that gets the anatomy right also get the biology right"):
        assert clause in t, f"the manuscript no longer quotes the registered clause: {clause!r}"


@pytest.mark.skipif(not MS.exists(), reason="manuscript absent")
def test_clause_two_is_answered_directly_not_only_by_the_correlation():
    t = MS.read_text()
    assert "## Result 1b" in t, (
        "Result 1b is gone. Clause 2 would then rest only on the ACS-versus-purity "
        "correlation, which is underpowered (n = 12, p = 0.80) and cannot carry the "
        "study's registered null.")


@pytest.mark.skipif(not (MS.exists() and AVB.exists()), reason="manuscript or artefact absent")
def test_result_1b_matches_its_artefact():
    avb = json.loads(AVB.read_text())
    t = MS.read_text()
    sec = t[t.index("## Result 1b"):t.index("## Result 2")]
    for cohort, b in avb["cohorts"].items():
        assert f"{b['n_getting_T_over_B_right']} of {b['n_methods']}" in sec, (
            f"{cohort}: the count of methods reproducing T > B does not match the artefact")
        assert b["top_method_by_acs"] in sec, (
            f"{cohort}: the top-ACS method named in the manuscript is not the artefact's")
        if b["top_method_frac_placing_B_over_T"] is not None:
            pct = f"{b['top_method_frac_placing_B_over_T']*100:.1f}%"
            assert pct in sec, f"{cohort}: the B-over-T fraction {pct} is not in Result 1b"


@pytest.mark.skipif(not MS.exists(), reason="manuscript absent")
def test_the_underpowered_trend_is_flagged_as_unclaimable():
    """The ACS-vs-lymphoid-failure correlation points the wrong way for the hypothesis.

    That is interesting and NOT significant at n = 12. A draft that states it without the
    caveat would be claiming the strongest available version of the study's own thesis on
    p = 0.14, which is exactly the error the project exists to document.
    """
    t = MS.read_text()
    assert "What must NOT be claimed" in t, "the unclaimable-trend warning is gone"
    sec = t[t.index("What must NOT be claimed"):]
    assert "not significant" in sec[:900].lower()


@pytest.mark.skipif(not MS.exists(), reason="manuscript absent")
def test_gbm_is_marked_discovery_and_lgg_registered_replication():
    """The pre-specification says so explicitly; presenting them as equals overstates GBM."""
    t = MS.read_text()
    assert "discovery" in t and "registered replication" in t, (
        "the manuscript presents GBM and LGG as equivalent cohorts. The pre-specification "
        "says 'This is LGG, not GBM. The anomaly was measured in GBM.' — GBM is discovery, "
        "LGG is the registered replication, and that ordering is the design's strength.")


@pytest.mark.skipif(not PRESPEC.exists(), reason="pre-specification absent")
def test_the_pre_registered_prediction_and_its_falsifier_still_exist():
    """Everything above rests on this file saying what the manuscript says it says."""
    t = PRESPEC.read_text()
    assert "DNA methylation will show T cells > B cells" in t, (
        "the pre-registered prediction is no longer in prespecified/; the manuscript's "
        "claim to have registered it in advance cannot stand")
    assert re.search(r"B\s*[>=≥]\s*T", t), (
        "the named falsifier ('methylation showing B >= T') is gone. A prediction without a "
        "stated falsifier is not a pre-registration.")


@pytest.mark.skipif(not RQ.exists(), reason="RESEARCH_QUESTION.md absent")
def test_the_claims_table_keeps_the_inconclusive_result_inconclusive():
    t = RQ.read_text()
    assert "INCONCLUSIVE" in t, (
        "docs/RESEARCH_QUESTION.md no longer marks the ACS-versus-accuracy result as "
        "inconclusive. At n = 12 it cannot be reported as a refutation.")
