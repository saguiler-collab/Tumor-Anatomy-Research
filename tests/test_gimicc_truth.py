"""
The GIMiCC truth-confirmation analysis (prespecified/gimicc_truth_confirmation.md, Addenda 1-2), on
planted data, before it is trusted with real estimates.

The arithmetic must not manufacture an ordering: a planted T > B is read as T > B, an inverted one
as B > T (the negative control), and a sample with an undefined (NaN) or empty lymphoid estimate is
dropped and counted, never zero-filled. The reading rule follows the document's text.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import gimicc_truth as gt  # noqa: E402


def _h4(t_cd4, t_cd8, nk, b, index=None) -> pd.DataFrame:
    n = len(b)
    idx = index or [f"TCGA-00-{i:04d}-01" for i in range(n)]
    return pd.DataFrame({"Tumor": 0.7, "CD4Tcell": t_cd4, "CD8Tcell": t_cd8, "Bcell": b, "NK": nk,
                         "Mono": 0.01, "Neu": 0.005, "Microglia": 0.05}, index=idx)


def test_planted_t_over_b_is_read_as_t_over_b():
    df = _h4([0.010] * 20, [0.005] * 20, [0.004] * 20, [0.002] * 20)
    r = gt.q1(df)
    assert r["T_exceeds_B"] is True and r["ordering"] == "T>NK>B" and r["n_samples_with_lymphoid"] == 20
    assert r["share_samples_T_over_B"] == 1.0


def test_inverted_signal_is_read_as_b_over_t():
    """The negative control: the same numbers with T and B swapped must flip the reading."""
    df = _h4([0.001] * 20, [0.001] * 20, [0.004] * 20, [0.015] * 20)
    r = gt.q1(df)
    assert r["T_exceeds_B"] is False and r["ordering"].startswith("B")
    assert r["share_samples_T_over_B"] == 0.0


def test_level3_and_level4_agree_when_cd4_cd8_is_defined():
    rng = np.random.default_rng(0)
    cd4, cd8, nk, b = (rng.uniform(0.0005, 0.01, 30) for _ in range(4))
    h4 = _h4(cd4, cd8, nk, b)
    h3 = h4.drop(columns=["CD4Tcell", "CD8Tcell"]).assign(Tcell=cd4 + cd8)
    pd.testing.assert_frame_equal(gt.lymph_shares(h4), gt.lymph_shares(h3))


def test_nan_and_empty_lymphoid_rows_are_dropped_not_zero_filled():
    """A NaN CD4/CD8 (GIMiCC's 'noisy' row) and an all-zero lymphoid row are both excluded. Had the
    NaN T been treated as 0, sample 0 would enter as pure B and pull the mean toward B."""
    df = _h4([np.nan, 0.0, 0.010], [np.nan, 0.0, 0.002], [0.003, 0.0, 0.003], [0.009, 0.0, 0.001])
    L = gt.lymph_shares(df)
    assert list(L.index) == [df.index[2]]
    assert np.allclose(L.sum(axis=1), 1.0)


def test_immune_total_excludes_nan_rows():
    df = _h4([np.nan, 0.01, 0.02], [np.nan, 0.01, 0.01], [0.01] * 3, [0.01] * 3)
    tot, n_ex = gt.immune_total(df, gt.IMMUNE)
    assert n_ex == 1 and list(tot.index) == list(df.index[1:])
    assert tot.iloc[0] == pytest.approx(0.01 + 0.01 + 0.01 + 0.01 + 0.01 + 0.005 + 0.05)


@pytest.mark.parametrize("ok, both, robust, reading, symmetric", [
    (True, [True, True], True, "CONFIRMED", "CONFIRMED"),
    (True, [True, True], False, "INCONCLUSIVE", "INCONCLUSIVE"),     # fragile cannot be CONFIRMED
    (True, [True, False], True, "NOT CONFIRMED", "NOT CONFIRMED"),
    (True, [False, False], False, "NOT CONFIRMED", "INCONCLUSIVE"),  # the text vs the symmetric reading
    (False, [True, True], True, "INCONCLUSIVE", "INCONCLUSIVE"),     # a failing control
    (False, [False, True], True, "INCONCLUSIVE", "INCONCLUSIVE"),
])
def test_reading_follows_the_document(ok, both, robust, reading, symmetric):
    assert gt.reading_of(ok, both, robust) == (reading, symmetric)


def _estimates_csv(tmp_path, t, b, nk=0.01, n=60) -> tuple[Path, pd.DataFrame]:
    """A one-method estimates_full-style CSV and GIMiCC-style shares on the same 60 samples."""
    ids = [f"TCGA-00-{i:04d}-01A" for i in range(n)]
    est = pd.DataFrame({"sample_id": ids, "method": "m", "Tumor": 0.9, "Macrophage_Microglia": 0.05,
                        "T_cell": t, "NK_cell": nk, "B_cell": b, "Endothelial": 0.0,
                        "Oligodendrocyte": 0.0, "Astrocyte": 0.0})
    p = tmp_path / "est.csv"
    est.to_csv(p, index=False)
    shares = pd.DataFrame({"T": 0.6, "NK": 0.3, "B": 0.1}, index=[i[:15] for i in ids])
    return p, shares


def test_per_sample_discordance_planted_agreement_and_inversion(tmp_path):
    import gimicc_secondary as gs
    p, shares = _estimates_csv(tmp_path, t=0.03, b=0.01)            # method T > B, GIMiCC T > B
    r = gs.per_sample_discordance(p, shares)["m"]
    assert r["n_scored"] == 60 and r["frac_same_direction"] == 1.0
    assert r["frac_method_B_over_T_while_gimicc_T_over_B"] == 0.0
    p, shares = _estimates_csv(tmp_path, t=0.01, b=0.03)            # inverted: method B > T
    r = gs.per_sample_discordance(p, shares)["m"]
    assert r["frac_same_direction"] == 0.0 and r["frac_method_B_over_T_while_gimicc_T_over_B"] == 1.0


def test_per_sample_discordance_counts_empty_lymphoid_rows_separately(tmp_path):
    import gimicc_secondary as gs
    p, shares = _estimates_csv(tmp_path, t=0.0, b=0.0, nk=0.0)     # no lymphoid signal at all
    r = gs.per_sample_discordance(p, shares)["m"]
    assert r["n_scored"] == 0 and r["n_no_lymphoid_signal"] == 60 and r["frac_same_direction"] is None


# ------------------------------------------------------------------ the real artefacts
REAL = Path(__file__).resolve().parent.parent / "results"      # never the conftest-redirected RESULTS_DIR


@pytest.fixture(scope="module")
def real():
    p = REAL / "gimicc_truth_confirmation.json"
    if not p.exists():
        pytest.skip("results/gimicc_truth_confirmation.json not produced yet")
    import json
    return json.loads(p.read_text())


def test_every_control_matched_real_samples(real):
    """The 2026-10-03 defect: a sample-key mismatch gave n = 0 and was scored as a FAIL. A control with
    no data is a defect, never a result."""
    for c, r in real["cohorts"].items():
        assert r["control_tumor_vs_absolute"]["n"] >= 100, c
        assert r["control_immune_vs_leukocyte_fraction"]["n"] >= 100, c
        assert r["negative_control_shuffled_cpgs"]["n"] >= 100, c
        assert r["structural_invariance"]["n_samples_compared"] >= 100, c


def test_recorded_reading_recomputes_from_recorded_controls(real):
    ok = all(r[k]["pass"] for r in real["cohorts"].values()
             for k in ("control_tumor_vs_absolute", "control_immune_vs_leukocyte_fraction",
                       "negative_control_shuffled_cpgs", "structural_invariance"))
    both = [real["cohorts"][c]["Q1_all_methylation_samples"]["T_exceeds_B"] for c in ("gbm", "lgg")]
    robust = all(r["probe_drop_robustness"]["robust"] for r in real["cohorts"].values())
    assert gt.reading_of(ok, both, robust) == (real["Q1_reading"],
                                               real["Q1_reading_if_robustness_required_in_both_directions"])
    for c, r in real["cohorts"].items():          # the pass flags follow the registered thresholds
        assert r["control_tumor_vs_absolute"]["pass"] == (r["control_tumor_vs_absolute"]["spearman"] >= 0.40)
        assert r["negative_control_shuffled_cpgs"]["pass"] == (abs(r["negative_control_shuffled_cpgs"]["tumor_vs_absolute_spearman"]) < 0.20)


def test_q1_shares_are_compositions(real):
    for c, r in real["cohorts"].items():
        for key in ("Q1_all_methylation_samples", "Q1_level3_sensitivity"):
            m = r[key]["mean_within_lymphoid_share"]
            assert abs(sum(m.values()) - 1) < 1e-3 and min(m.values()) >= 0, (c, key)
