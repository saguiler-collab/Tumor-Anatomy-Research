"""
Extension E2 (truth-free identifiability diagnostics) and its Addendum 1 robustness checks.

The exact test and the stability statistic are checked on cases with known answers. On the real
artefacts, the registered values are RECOMPUTED from the saved per-sample estimates rather than
compared with typed literals. The negative control matters most: with each cohort's truths shuffled
across samples, stability must almost never "predict" truth agreement. If it did, the H1 test would
be measuring the stability spread, not accuracy.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import identifiability_diagnostics as e2  # noqa: E402
import identifiability_robustness as rb  # noqa: E402


# ---------- the exact one-sided permutation test (no data) ----------

def test_exact_test_perfect_concordance_is_one_in_720():
    rho, p = e2.exact_one_sided([1, 2, 3, 4, 5, 6], [10, 20, 30, 40, 50, 60])
    assert rho == pytest.approx(1.0)
    assert p == pytest.approx(1 / 720)


def test_exact_test_two_adjacent_swaps_is_twelve_in_720():
    """E2's own D1 configuration: sum of squared rank differences 4 -> rho 0.886, 12 of 720."""
    rho, p = e2.exact_one_sided([1, 2, 3, 4, 5, 6], [2, 1, 4, 3, 5, 6])
    assert rho == pytest.approx(1 - 6 * 4 / 210)
    assert p == pytest.approx(12 / 720)


def test_exact_test_rejects_an_inverted_relation():
    rho, p = e2.exact_one_sided([1, 2, 3, 4, 5, 6], [60, 50, 40, 30, 20, 10])
    assert rho == pytest.approx(-1.0)
    assert p == pytest.approx(1.0)


# ---------- the stability statistic ----------

def test_mean_pairwise_drops_constant_columns_and_needs_three():
    rng = np.random.default_rng(0)
    base = rng.normal(size=50)
    X = np.column_stack([base, base + rng.normal(scale=.1, size=50), base * 2, np.zeros(50)])
    with_const = rb.mean_pairwise(X)
    without = rb.mean_pairwise(X[:, :3])
    assert with_const == pytest.approx(without)            # the constant column is excluded
    assert np.isnan(rb.mean_pairwise(X[:, [0, 3]]))         # fewer than 3 usable -> not estimable


def test_identical_fits_are_collapsed(tmp_path, monkeypatch):
    """Addendum 1 check 5: byte-identical setting files count once."""
    monkeypatch.setattr(e2, "OUT_DIR", tmp_path)
    for i, (label, _, _) in enumerate(e2.SETTINGS):
        body = "same" if label in ("predeclared", "s1") else f"fit {i}"
        (tmp_path / f"unmix_gbm_{label}.csv").write_text(body)
    kept = rb.distinct_settings("gbm")
    assert kept[0] == "predeclared" and "s1" not in kept
    assert len(kept) == len(e2.SETTINGS) - 1


# ---------- the real artefacts (skipped where results/ is absent) ----------

def _have_real() -> bool:
    need = [e2.OUT, *(e2.OUT_DIR / f"unmix_{h}_{lab}.csv" for h in ("gbm", "lgg") for lab, _, _ in e2.SETTINGS)]
    return all(p.exists() for p in need)


real = pytest.mark.skipif(not _have_real(), reason="E2 artefacts absent (results/ is not tracked)")


#: The real results directory, captured at import -- before conftest's per-test isolation
#: re-points config.RESULTS_DIR at an empty temporary tree (the same pattern as
#: test_anatomic_auc.ARTEFACT). The loaders read config.RESULTS_DIR at call time.
REAL_RESULTS = e2.config.RESULTS_DIR


@pytest.fixture(scope="module")
def data():
    """Read-only load of the real E2 inputs; config is restored whatever happens."""
    saved = e2.config.RESULTS_DIR
    e2.config.RESULTS_DIR = REAL_RESULTS
    try:
        return {h: rb.load(h) for h in ("gbm", "lgg")}
    finally:
        e2.config.RESULTS_DIR = saved


@real
def test_registered_values_recompute_from_saved_estimates(data):
    reg = json.loads(e2.OUT.read_text())
    assert all(c["pass"] for c in reg["controls"].values())
    for (c, h) in [(c, h) for c in e2.COMPARTMENTS for h in ("gbm", "lgg")]:
        s = rb.unit_stats(data[h]["comp"][c], None, list(range(len(e2.SETTINGS))))
        r = reg["units"][f"{c}|{h}"]
        assert s["d1"] == pytest.approx(r["d1_stability"], abs=1e-9), (c, h)
        assert s["d2"] == pytest.approx(r["d2_agreement"], abs=1e-9), (c, h)
        if np.isfinite(r["truth_rho_unmix"]):
            assert s["t1"] == pytest.approx(r["truth_rho_unmix"], abs=1e-9), (c, h)
            assert s["t2"] == pytest.approx(r["truth_rho_median_methods"], abs=1e-9), (c, h)


@real
def test_shuffled_truths_almost_never_support_h1(data):
    """The negative control. Shuffled truths keep every unit's stability; only accuracy is destroyed."""
    rng = np.random.default_rng(20261002)
    supported, tumour_rho = 0, []
    n_shuffles = 20
    for _ in range(n_shuffles):
        st = {}
        for h in ("gbm", "lgg"):
            for c in e2.PRIMARY:
                u = dict(data[h]["comp"][c])
                u["truth"] = rng.permutation(u["truth"])
                st[(c, h)] = rb.unit_stats(u, None, list(range(len(e2.SETTINGS))))
        rho, p = rb.h1(st, "d1", "t1")
        supported += rho > 0 and p < 0.05
        tumour_rho.append(st[("Tumor", "lgg")]["t1"])
    assert supported <= n_shuffles // 4, f"{supported} of {n_shuffles} shuffles came out SUPPORTED"
    assert abs(float(np.median(tumour_rho))) < 0.1           # shuffling did destroy accuracy


@real
def test_robustness_artefact_is_internally_consistent():
    p = rb.OUT
    if not p.exists():
        pytest.skip("robustness artefact absent")
    r = json.loads(p.read_text())
    assert r["control_pass"] is True
    assert r["check1_contrast_probability"] == pytest.approx(1 / 15)
    assert r["check5_distinct_fits"]["fits_used"]["lgg"] == [lab for lab, _, _ in e2.SETTINGS]
    assert "s1" not in r["check5_distinct_fits"]["fits_used"]["gbm"]     # GBM's pre-declared shift is 1
    for h in ("gbm", "lgg"):                                              # degraded modes are labelled
        lk = r["check4_leukocyte_truth_by_method"][h]
        assert set(lk["degraded_mode"]) <= set(lk["methods"])
        assert set(lk["degraded_mode"]).isdisjoint(r["check6_comparable_panel"]["methods"][h])


# ---------- Addendum 1, check 7: matched denominators for the EpiDISH units ----------

import identifiability_denominators as dn  # noqa: E402
import pandas as pd  # noqa: E402


def _perfect_estimates(n=40, seed=3):
    """Estimates whose lymphoid share of leukocytes IS the truth, with leukocyte levels varying 20-fold."""
    rng = np.random.default_rng(seed)
    ids = [f"TCGA-AA-{i:04d}-01A" for i in range(n)]
    lym_share = rng.uniform(0.05, 0.6, n)                     # what EpiDISH's blood reference measures
    leuk = rng.uniform(0.01, 0.20, n)                         # varies independently across samples
    est = pd.DataFrame({"Tumor": 1 - leuk, "Macrophage_Microglia": leuk * (1 - lym_share),
                        "T_cell": leuk * lym_share / 2, "NK_cell": leuk * lym_share / 4,
                        "B_cell": leuk * lym_share / 4}, index=ids)
    truth = pd.Series(lym_share, index=[dn.e2.k4s(i) for i in ids])
    return est, truth


def test_share_is_undefined_without_leukocytes():
    est = pd.DataFrame({"Macrophage_Microglia": [0.10, 0.0], "T_cell": [0.02, 0.0], "NK_cell": [0.01, 0.0],
                        "B_cell": [0.02, 0.0], "Tumor": [0.85, 1.0]}, index=["a", "b"])
    s = dn.share(est, "Lymphoid")
    assert s["a"] == pytest.approx(0.05 / 0.15)
    assert np.isnan(s["b"])                                   # undefined, never 0


def test_matched_denominator_recovers_a_perfect_estimate_that_the_registered_one_attenuates():
    est, truth = _perfect_estimates()
    matched, _ = dn.truth_rho(dn.share(est, "Lymphoid"), truth)
    registered, _ = dn.truth_rho(dn.e2.comp(est, "Lymphoid"), truth)   # tissue fraction vs immune share
    assert matched == pytest.approx(1.0)
    assert registered < 0.9                                   # the attenuation check 7 corrects


def test_matched_denominator_rejects_a_shuffled_estimate():
    """Negative control: the matched comparison must not reward an estimate unrelated to the truth."""
    est, truth = _perfect_estimates()
    shuffled = est.copy()
    shuffled.index = np.random.default_rng(9).permutation(est.index)
    rho, _ = dn.truth_rho(dn.share(shuffled, "Lymphoid"), truth)
    assert abs(rho) < 0.4


# ---------- prespecified/agreement_selection_test.md: DECEPTICON's rule as implemented ----------

import agreement_selection_test as ag  # noqa: E402


def test_decepticon_rule_picks_the_two_most_correlated_pairs_and_weights_repeats():
    rng = np.random.default_rng(0)
    base = rng.normal(size=60)
    X = pd.DataFrame({"a": base, "b": base + rng.normal(scale=.01, size=60),
                      "c": base + rng.normal(scale=.2, size=60), "d": rng.normal(size=60),
                      "z": np.zeros(60)})                       # constant: excluded (DECEPTICON rule 5)
    pairs = ag.select_pairs(X)
    assert pairs[0] == ("a", "b") and "z" not in {m for p in pairs for m in p}
    w = ag.weights([("a", "b"), ("a", "c")], ["a", "b", "c", "d"])
    assert list(w) == [0.5, 0.25, 0.25, 0.0] and w.sum() == pytest.approx(1.0)


def test_agreement_selects_a_shared_bias_over_the_one_accurate_method():
    """Why the premise needs testing: two methods sharing a bias agree with each other, the accurate
    one agrees with nobody, and the rule picks the bias."""
    rng = np.random.default_rng(1)
    truth = rng.uniform(.2, .9, 80)
    bias = rng.uniform(.2, .9, 80)                             # a shared, wrong signal
    X = pd.DataFrame({"accurate": truth + rng.normal(scale=.02, size=80),
                      "biased1": bias + rng.normal(scale=.02, size=80),
                      "biased2": bias + rng.normal(scale=.02, size=80),
                      "biased3": .5 * bias + .5 * rng.uniform(.2, .9, 80)})
    pairs = ag.select_pairs(X)
    ens = X.to_numpy() @ ag.weights(pairs, list(X.columns))
    assert ("biased1", "biased2") in pairs
    assert ag.spearman(ens, truth) < ag.spearman(X["accurate"].to_numpy(), truth) - 0.5


def test_lymphoid_units_are_shares_of_the_leukocyte_total():
    c = {"Tumor": np.array([.8, 1.0]), "Macrophage_Microglia": np.array([.1, 0.]),
         "T_cell": np.array([.05, 0.]), "NK_cell": np.array([.0, 0.]), "B_cell": np.array([.05, 0.])}
    v = ag.unit_values(c)
    assert v["Lymphoid"][0] == pytest.approx(0.5) and np.isnan(v["Lymphoid"][1])
    assert v["Leukocytes"][0] == pytest.approx(0.2)


# ---------- Addendum 2 (S2): the anatomy arm's guard ----------

def test_s2_mode_refuses_any_method_but_unmix_before_loading_anything():
    """--e2-s2 is defined for deseq2_unmix only; for any other method it must stop before the
    multi-GB atlas load, never fall through to a re-measurement under the ablation's name."""
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "remeasure_method.py"), "--method", "dwls",
                        "--e2-s2"], capture_output=True, text=True, timeout=120, cwd=ROOT)
    assert r.returncode == 2
    assert "BLOCKED: --e2-s2 is defined for deseq2_unmix only" in r.stdout
