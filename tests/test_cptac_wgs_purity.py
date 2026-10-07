"""
The CPTAC whole-genome secondary (scripts/cptac_wgs_purity.py; prespecified/cptac_wgs_purity_secondary.md)
on planted data: the per-sample correlation must be able to fail (an inverted truth scores negative, a
constant one is not scored), the coverage screen picks exactly the complete samples, and the S2.2
reading credits coverage only when the complete-sample run passes and the all-sample run does not.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import cptac_wgs_purity as w  # noqa: E402


def test_inverted_truth_is_rejected_and_planted_truth_recovered():
    rng = np.random.default_rng(0)
    purity = np.linspace(0.4, 0.95, 15)
    good = w._rho(purity + rng.normal(0, 0.02, 15), purity, rng)
    bad = w._rho(1 - purity, purity, rng)
    assert good["rho"] > 0.9 and good["perm_p"] < 0.01
    assert bad["rho"] < -0.9                       # an inverted truth must never read as agreement
    assert w._rho(np.full(15, 0.8), purity, rng)["rho"] is None        # constant: not scored
    assert w._rho(purity[:3], purity[:3], rng)["rho"] is None          # too few cases: not scored


def test_coverage_screen_selects_the_complete_samples_at_any_threshold_in_the_gap():
    rng = np.random.default_rng(1)
    m = pd.DataFrame(rng.random((1000, 6)), columns=[f"s{i}" for i in range(6)])
    for col, miss in (("s0", 0.45), ("s1", 0.35), ("s2", 0.03), ("s3", 0.01)):
        m.loc[rng.choice(1000, int(miss * 1000), replace=False), col] = np.nan
    expect = ["s2", "s3", "s4", "s5"]
    assert all(w.covered(m, t) == expect for t in (0.68, 0.90, 0.95))
    assert w.covered(m, 0.99) == ["s3", "s4", "s5"]                   # the screen is a real threshold


def test_s22_reading_credits_coverage_only_when_the_rule_says_so():
    yes = "failed primary control attributed to methylation coverage"
    assert w.reading_s22(0.15, 0.92) == yes
    assert w.reading_s22(None, 0.92) == yes
    assert w.reading_s22(0.55, 0.92) != yes        # the all-sample run already passed: nothing to explain
    assert w.reading_s22(0.15, 0.30) != yes        # the complete-sample run fails too
    assert w.reading_s22(0.15, None) != yes
