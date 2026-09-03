"""
Survival tests. The subject here is not the Cox model — lifelines is tested upstream —
but the guards around it: the power bound, the refusal to rank on an underpowered
cohort, and the firewall that keeps outcomes out of method selection.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from ivygap.survival.run_survival import (MIN_EVENTS_FOR_RANKING,
                                          approximate_detectable_delta, assess_power,
                                          assert_selection_frozen, rank_methods)


def _clinical(n_patients: int, n_events: int) -> pd.DataFrame:
    return pd.DataFrame({
        "time": np.linspace(100, 900, n_patients),
        "event": [1.0] * n_events + [0.0] * (n_patients - n_events),
    }, index=[f"P{i:02d}" for i in range(n_patients)])


def test_detectable_delta_shrinks_with_events():
    assert approximate_detectable_delta(400) < approximate_detectable_delta(40)


def test_ivygap_sized_cohort_is_flagged_underpowered():
    """The finding that justified moving the headline off survival."""
    c = assess_power(_clinical(40, 34), n_covariates=8)
    assert not c.adequately_powered
    assert c.verdict.startswith("UNDERPOWERED")
    assert c.detectable_delta_c_index > 0.4


def test_tcga_sized_cohort_is_adequate():
    c = assess_power(_clinical(115, 93), n_covariates=8)
    assert c.adequately_powered
    assert c.verdict.startswith("ADEQUATE")


def test_zero_events_is_blocked_not_reported():
    c = assess_power(_clinical(40, 0), n_covariates=8)
    assert c.verdict.startswith("BLOCKED")
    assert not c.adequately_powered


def test_underpowered_cohort_refuses_to_rank():
    """
    The guard that matters. Two methods differing by 0.04 C-index on 34 events must not
    be reported as one beating the other.
    """
    contract = assess_power(_clinical(40, 34), n_covariates=8)
    results = {
        "a": {"status": "OK", "c_index_with_composition": 0.68, "n_events": 34},
        "b": {"status": "OK", "c_index_with_composition": 0.64, "n_events": 34},
    }
    table = rank_methods(results, contract)
    assert table["interpretation"].str.startswith("INCONCLUSIVE").all()


def test_powered_cohort_separates_clearly_different_methods():
    contract = assess_power(_clinical(400, 350), n_covariates=8)
    results = {
        "a": {"status": "OK", "c_index_with_composition": 0.75},
        "b": {"status": "OK", "c_index_with_composition": 0.50},
    }
    table = rank_methods(results, contract)
    assert table.loc["a", "interpretation"].startswith("tied with best")
    assert table.loc["b", "interpretation"] == "worse than best"


def test_selection_firewall_rejects_outcome_based_criterion(tmp_path):
    """The rule that outcomes may not select methods, enforced rather than remembered."""
    path = tmp_path / "decision.json"
    path.write_text(json.dumps({
        "selected_method": "nnls",
        "criterion": "highest out-of-fold survival c_index",
    }))
    with pytest.raises(AssertionError, match="outcome-blind"):
        assert_selection_frozen(path)


def test_selection_firewall_accepts_outcome_blind_criterion(tmp_path):
    path = tmp_path / "decision.json"
    path.write_text(json.dumps({
        "selected_method": "svr",
        "criterion": "mean absolute error on primary cell types, donor-held-out "
                     "pseudobulk, aggregated donor-equally",
    }))
    assert_selection_frozen(path)          # must not raise


def test_missing_selection_file_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="before survival"):
        assert_selection_frozen(tmp_path / "nope.json")
