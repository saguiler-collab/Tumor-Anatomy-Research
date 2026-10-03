"""The power analysis of the registered agreement test (scripts/agreement_power.py)."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("ap", ROOT / "scripts" / "agreement_power.py")
ap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ap)


def test_vectorised_bootstrap_reproduces_the_registered_one():
    assert ap.check_against_registered()["max_abs_difference"] <= 0.01


def test_pearson_to_spearman_mapping():
    rng = np.random.default_rng(0)
    r = ap.pearson_for_spearman(0.6)
    x = rng.multivariate_normal([0, 0], [[1, r], [r, 1]], size=200_000)
    assert stats.spearmanr(x[:, 0], x[:, 1])[0] == pytest.approx(0.6, abs=0.01)


def test_no_true_correlation_rarely_meets_the_criterion():
    """Negative control: with rho = 0 the criterion must almost never be met."""
    assert ap.simulate(0.0, 12, 200, 300, seed=3)["p_criterion_met"] < 0.05


ART = ROOT / "results" / "agreement_power.json"


@pytest.mark.skipif(not ART.exists(), reason="run scripts/agreement_power.py first")
def test_power_rises_with_effect_and_panel_size():
    rows = json.loads(ART.read_text())["rows"]
    get = {(r["n_methods"], r["true_spearman"]): r["p_criterion_met"] for r in rows}
    assert get[(12, 0.9)] > get[(12, 0.6)] > get[(12, 0.0)]
    assert get[(50, 0.7)] > get[(12, 0.7)]


@pytest.mark.skipif(not ART.exists(), reason="run scripts/agreement_power.py first")
def test_a_true_correlation_at_the_bar_is_met_about_half_the_time_at_any_size():
    """The criterion thresholds the ESTIMATE at the bar, so a true value exactly at the bar
    cannot reach high power by adding methods -- the manuscript must not say otherwise."""
    rows = json.loads(ART.read_text())["rows"]
    at_bar = [r["p_criterion_met"] for r in rows if r["true_spearman"] == 0.6]
    assert at_bar and max(at_bar) < 0.65
