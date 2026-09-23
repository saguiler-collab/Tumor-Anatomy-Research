"""A yardstick's score direction must match what its loader actually returns.

WHY THIS EXISTS. `run_anatomic` passed `higher_is_better=False` for every yardstick under
the comment "every yardstick here is an error metric". That held for `synthetic_mixtures`
(mean absolute error) and not for `absolute_purity`, which returns a Spearman correlation of
estimated tumour content against DNA purity — higher is better. `test_agreement` negates the
truth vector when told lower is better, so the orthogonal yardstick's rho was published with
its sign inverted (-0.0810 where the data give +0.0810). |rho|, the p-value and every
pre-registered threshold were unaffected; the sign was not.

A wrong sign here inverts the study's central claim about whether ACS tracks accuracy, so it
is pinned by a test that FAILS if a direction is changed without changing the loader.
"""
from __future__ import annotations

import numpy as np
import pytest

from ivygap.anatomic import agreement
from ivygap.bench import real_yardsticks as ry


def test_every_yardstick_declares_a_direction():
    assert set(ry.HIGHER_IS_BETTER) >= set(agreement.YARDSTICKS), (
        "a yardstick in agreement.YARDSTICKS has no declared direction: "
        f"{sorted(set(agreement.YARDSTICKS) - set(ry.HIGHER_IS_BETTER))}")


def test_absolute_purity_is_a_correlation_so_higher_is_better():
    """Pinned against the loader's own source, not against a remembered convention."""
    assert ry.HIGHER_IS_BETTER["absolute_purity"] is True
    src = ry.load_absolute_purity.__code__.co_consts
    assert any(isinstance(c, str) and "spearman_vs_purity" in c for c in src), (
        "load_absolute_purity no longer reads `spearman_vs_purity`. If it now returns an "
        "ERROR metric, flip HIGHER_IS_BETTER['absolute_purity'] to False; if it returns a "
        "different correlation, update this test to name it.")


def test_synthetic_mixtures_is_an_error_metric_so_lower_is_better():
    assert ry.HIGHER_IS_BETTER["synthetic_mixtures"] is False
    src = ry.load_synthetic_mixtures.__code__.co_consts
    assert any(isinstance(c, str) and c == "mae" for c in src), (
        "load_synthetic_mixtures no longer reads `mae`. Re-derive the direction before "
        "changing this test.")


@pytest.mark.parametrize("higher_better,expect_sign", [(True, +1), (False, -1)])
def test_the_flag_actually_controls_the_sign(higher_better, expect_sign):
    """The negative control for the two tests above.

    If the flag did not move the sign, pinning it would prove nothing.
    """
    methods = [f"m{i}" for i in range(8)]
    acs = {m: float(i) for i, m in enumerate(methods)}
    truth = {m: float(i) for i, m in enumerate(methods)}      # perfectly concordant
    res = agreement.test_agreement(acs, truth, yardstick="synthetic_mixtures",
                                   higher_truth_is_better=higher_better, n_boot=200)
    assert np.sign(res.rho) == expect_sign, (
        f"with higher_truth_is_better={higher_better} a perfectly concordant pair gave "
        f"rho={res.rho}; the flag is not controlling the sign.")


def test_the_reported_absolute_purity_rho_matches_a_direct_recomputation():
    """End-to-end: the artefact's sign must match Spearman on the underlying per-method data.

    This is the check that would have caught the defect. It reads the two artefacts and
    correlates them itself rather than trusting either.
    """
    import json
    from pathlib import Path

    import pandas as pd
    from scipy.stats import spearmanr

    # THE REAL results tree, derived from this file's location. `config.RESULTS_DIR` is
    # redirected to a temp root by conftest so tests cannot write over a real run -- which
    # is right, but it also meant this check silently SKIPPED instead of ever reading the
    # shipped artefacts. A test that always skips is not a test.
    results = Path(__file__).resolve().parent.parent / "results"
    ypath = results / "absolute_purity_yardstick.json"
    lpath = results / "anatomic" / "acs_leaderboard.csv"
    apath = results / "anatomic" / "agreement_report.json"
    if not (ypath.exists() and lpath.exists() and apath.exists()):
        pytest.skip("a real run's artefacts are not present")

    per = json.loads(ypath.read_text())["methods"]
    lb = pd.read_csv(lpath).set_index("method")
    pairs = [(float(lb.loc[m, "acs"]), float(d["spearman_vs_purity"]))
             for m, d in per.items()
             if m in lb.index and not bool(lb.loc[m, "is_control"])
             and d.get("spearman_vs_purity") is not None]
    if len(pairs) < 3:
        pytest.skip("too few methods with both scores")
    rho = spearmanr([a for a, _ in pairs], [t for _, t in pairs])[0]

    rep = json.loads(apath.read_text())
    row = next((y for y in rep.get("by_yardstick", [])
                if y.get("yardstick") == "absolute_purity"), None)
    if row is None or row.get("rho") is None or not np.isfinite(row.get("rho", np.nan)):
        pytest.skip("absolute_purity was not computed in this run")

    # The report restricts to COMPARABLE methods, so |rho| may differ slightly. The SIGN
    # is what this test pins, and it must not be opposite.
    assert np.sign(row["rho"]) == np.sign(rho), (
        f"agreement_report.json reports rho={row['rho']:+.4f} for absolute_purity, but "
        f"Spearman recomputed directly from acs_leaderboard.csv and "
        f"absolute_purity_yardstick.json over {len(pairs)} methods gives {rho:+.4f}. "
        f"Opposite signs mean the direction flag is wrong again.")
