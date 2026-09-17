"""
Tests for the real ground-truth yardstick loaders.

The point of these is the distinction the module exists to preserve: "we could not
obtain ground truth" and "ground truth exists and covers 2 of 10 methods" must not
produce the same output, because only the second one tells you what to do next.
"""

from __future__ import annotations

import pytest

from ivygap import config
from ivygap.bench import real_yardsticks as ry
from ivygap.anatomic.agreement import MIN_METHODS_FOR_CORRELATION


def test_real_yardstick_is_wired_to_actual_data():
    scores, prov = ry.load_synthetic_mixtures()
    if not prov.get("available"):
        pytest.skip("TCGA benchmark not vendored")
    assert scores, "the file is present but produced no scores"
    assert prov["is_real_ground_truth"] is True
    assert prov["n_mixtures"] == 500


def test_units_are_converted_not_left_mixed():
    """
    The source reports MAE in percent; this project works in fractions. Leaving two unit
    systems in one table is harmless for a rank correlation and wrong for anything else.
    """
    scores, prov = ry.load_synthetic_mixtures()
    if not scores:
        pytest.skip("TCGA benchmark not vendored")
    assert all(0.0 < v < 1.0 for v in scores.values()), scores
    assert "percent" in prov["units"]


def test_coverage_limitation_is_stated_with_a_number_and_a_remedy():
    """
    A limitation without a number is a shrug. This one has to say how many methods it
    covers, how many are needed, and what would extend it.
    """
    scores, prov = ry.load_synthetic_mixtures()
    if not scores:
        pytest.skip("TCGA benchmark not vendored")
    if len(scores) < MIN_METHODS_FOR_CORRELATION:
        lim = prov["limitation"]
        assert str(len(scores)) in lim
        assert str(MIN_METHODS_FOR_CORRELATION) in lim
        assert "single-cell" in lim


def test_method_names_are_mapped_explicitly_not_by_lowercasing():
    """
    A silent name-matching rule is how one method's score ends up attached to another.
    The mapping must be a declared table.
    """
    assert ry.METHOD_NAME_MAP == {"NNLS": "nnls", "SVR": "svr"}
    scores, prov = ry.load_synthetic_mixtures()
    if scores:
        assert not prov["methods_not_mapped"], (
            f"unmapped source labels would be silently dropped: "
            f"{prov['methods_not_mapped']}")


def test_absolute_purity_is_never_silent_about_its_state():
    """
    The invariant is NOT "ABSOLUTE purity is blocked" — that was merely true when this test
    was written, and `scripts/absolute_purity_yardstick.py` can now compute it. The invariant
    is that the yardstick never returns silently: unavailable means a specific reason naming
    the blocker, available means real scores plus the provenance a reader needs to judge them.

    Written state-agnostically on purpose. A test that pins the current state forces the next
    person to weaken it, and a weakened test is worse than a state-agnostic one.
    """
    scores, prov = ry.load_absolute_purity()
    assert prov["is_real_ground_truth"] is True, "this yardstick is measured, not simulated"

    if prov["available"]:
        assert scores, "available with no scores is the silence this test exists to prevent"
        assert all(isinstance(v, float) for v in scores.values())
        assert all(-1.0 <= v <= 1.0 for v in scores.values()), "these are correlations"
        # The independence claim must travel WITH the numbers. This yardstick is
        # independent on the ground-truth side and not on the reference side, and a reader
        # who gets the scores without that sentence will overclaim.
        assert prov.get("not_independent"), (
            "an independent-yardstick claim must carry its own limit; the reference is "
            "GBmap-derived and the provenance has to say so")
        assert prov.get("n_samples")
    else:
        assert scores == {}
        assert "BLOCKED" in prov["reason"]


def test_all_three_protocol_yardsticks_are_represented():
    scores, prov = ry.load_all()
    assert set(scores) == {"synthetic_mixtures", "absolute_purity", "sc_pseudobulk"}
    for name, p in prov.items():
        assert "available" in p, name
        if not p["available"]:
            assert p.get("reason"), f"{name} is unavailable with no reason given"


def test_vendored_evidence_matches_its_recorded_hashes():
    import json
    path = ry.TCGA_BENCHMARK_DIR / "PROVENANCE.json"
    if not path.exists():
        pytest.skip("TCGA benchmark not vendored")
    prov = json.loads(path.read_text())
    for name, info in prov["files"].items():
        f = ry.TCGA_BENCHMARK_DIR / name
        if f.exists():
            assert config.sha256_file(f) == info["sha256"], f"{name} changed"
