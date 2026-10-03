"""The BayesPrism verdict rule (prespecified/bayesprism_authors_prediction.md), as code."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("bv", ROOT / "scripts" / "bayesprism_verdict.py")
bv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bv)


def test_lower_b_share_under_authors_is_supported():
    assert bv.verdict_from(0.30, 0.45) == "SUPPORTED"


def test_higher_or_equal_is_not_supported():
    """The negative control: the rule must be able to say no."""
    assert bv.verdict_from(0.50, 0.45) == "NOT SUPPORTED"
    assert bv.verdict_from(0.45, 0.45) == "NOT SUPPORTED"


def test_missing_run_is_blocked_not_guessed():
    assert bv.verdict_from(None, 0.4) == "BLOCKED"
    assert bv.verdict_from(float("nan"), 0.4) == "BLOCKED"


def _frame(b, t, n):
    idx = [f"TCGA-00-{i:04d}-01A" for i in range(n)]
    cols = {c: 0.0 for c in ("Tumor", "T_cell", "NK_cell", "B_cell")}
    f = pd.DataFrame([cols] * n, index=idx)
    f["B_cell"], f["T_cell"] = b, t
    return f


def test_per_sample_comparison_is_inconclusive_below_ten_paired_samples():
    a, p = _frame(0.01, 0.02, 6), _frame(0.02, 0.01, 6)
    truth_idx = [bv.k4s(i) for i in a.index]
    assert bv.paired(a, p, truth_idx)["per_sample"] == "INCONCLUSIVE"


def test_per_sample_comparison_counts_direction():
    n = 30
    a, p = _frame(np.linspace(.01, .02, n), .03, n), _frame(.03, .03, n)
    truth_idx = [bv.k4s(i) for i in a.index]
    out = bv.paired(a, p, truth_idx)
    assert out["n_samples_with_signal_in_both"] == n
    assert out["n_lower_under_authors"] == n and out["median_difference_authors_minus_pipeline"] < 0
