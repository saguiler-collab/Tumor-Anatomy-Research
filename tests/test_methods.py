"""
Every method must satisfy the same output contract, and must actually recover
composition better than chance. A method that returns a valid-looking simplex of noise
would pass a contract test alone, so accuracy is checked too.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.bench import metrics
from ivygap.deconv.base import project_to_simplex, to_cell_fractions
from ivygap.deconv.bayesian import BayesianDeconvolution
from ivygap.deconv.classical import NNLSDeconvolution
from ivygap.deconv.elastic_net import ElasticNetDeconvolution
from ivygap.deconv.reference_based import (BisqueDeconvolution, DWLSDeconvolution,
                                           MuSiCDeconvolution, SCDCDeconvolution,
                                           SCDCEnsembleDeconvolution)

FAST_METHODS = [NNLSDeconvolution, MuSiCDeconvolution, BisqueDeconvolution,
                SCDCDeconvolution, ElasticNetDeconvolution]


@pytest.mark.parametrize("cls", FAST_METHODS)
def test_output_contract(cls, problem):
    data, _ = problem
    est = cls().fit_predict(data)

    assert list(est.columns) == list(config.CELL_TYPES)
    assert list(est.index) == data.samples
    assert (est.to_numpy() >= -1e-9).all(), "fractions must be non-negative"
    sums = est.sum(axis=1)
    assert np.allclose(sums[np.isfinite(sums)], 1.0, atol=1e-6), "rows must sum to 1"


@pytest.mark.parametrize("cls", FAST_METHODS)
def test_beats_chance(cls, problem):
    """
    A uniform guess is the null. Any working method must beat it on MAE; one that does
    not is not deconvolving, whatever else it is doing.
    """
    data, test_set = problem
    est = cls().fit_predict(data)

    uniform = pd.DataFrame(1.0 / len(config.CELL_TYPES),
                           index=est.index, columns=est.columns)
    m_est = metrics.overall_metrics(est, test_set.truth, groups=test_set.donors)
    m_null = metrics.overall_metrics(uniform, test_set.truth, groups=test_set.donors)
    assert m_est["mae_primary"] < m_null["mae_primary"], (
        f"{cls.__name__} MAE {m_est['mae_primary']:.4f} did not beat the uniform "
        f"guess {m_null['mae_primary']:.4f}"
    )


def test_methods_cannot_mutate_shared_input(problem):
    """
    Every method receives the same object. One that mutated it would corrupt every
    method scheduled after it, and the corruption would look like a method difference.
    """
    data, _ = problem

    class Vandal(NNLSDeconvolution):
        name = "vandal"

        def _solve_all(self, d):
            d.bulk.iloc[0, 0] = -999.0
            return super()._solve_all(d)

    with pytest.raises(RuntimeError, match="mutated"):
        Vandal().fit_predict(data)


def test_bayesian_reports_credible_intervals(problem):
    data, _ = problem
    m = BayesianDeconvolution(n_iter=60, n_burnin=20)
    est = m.fit_predict(data)

    for key in ("mean", "sd", "ci_lower", "ci_upper"):
        assert key in m.posterior_
    lo, hi = m.posterior_["ci_lower"], m.posterior_["ci_upper"]
    assert (hi.to_numpy() >= lo.to_numpy() - 1e-9).all()
    assert (hi.to_numpy() > lo.to_numpy()).any(), "posteriors must have non-zero width"


def test_bayesian_prior_is_actually_used(problem):
    """
    A 'Bayesian' method whose answer does not move with the prior is not using one.
    A very large alpha forces the posterior toward uniform.
    """
    data, _ = problem
    weak = BayesianDeconvolution(alpha=0.5, n_iter=60, n_burnin=20).fit_predict(data)
    strong = BayesianDeconvolution(alpha=5000.0, n_iter=60, n_burnin=20).fit_predict(data)

    uniform = 1.0 / len(config.CELL_TYPES)
    assert (strong - uniform).abs().mean().mean() < (weak - uniform).abs().mean().mean()


def test_elastic_net_records_its_hyperparameters(problem):
    """Auditability: the chosen penalty must be inspectable, not buried in the fit."""
    data, _ = problem
    m = ElasticNetDeconvolution()
    m.fit_predict(data)
    report = m.selection_report(data.samples)
    assert len(report) == len(data.samples)
    assert (report["alpha"].dropna() > 0).all()


def test_scdc_ensemble_reports_degeneracy(problem):
    """With one reference there is no ensemble, and it must say so rather than imply one."""
    data, _ = problem
    m = SCDCEnsembleDeconvolution()
    m.fit_predict(data)
    assert m.degenerate_ is True
    assert m.ensemble_weights_ == {data.primary.name: 1.0}


def test_simplex_projection_refuses_to_invent_a_composition():
    """An all-zero solution has no composition; returning uniform would be a fabrication."""
    assert np.isnan(project_to_simplex(np.zeros(5))).all()
    assert project_to_simplex(np.array([1.0, 3.0])).tolist() == [0.25, 0.75]
    assert project_to_simplex(np.array([-1.0, 1.0])).tolist() == [0.0, 1.0]


def test_cell_size_correction_moves_the_estimate():
    """
    A type with twice the mRNA per cell must end up with a smaller *cell* fraction than
    its RNA fraction. If this is a no-op the estimand is wrong for every method at once.
    """
    rna = np.array([0.5, 0.5])
    out = to_cell_fractions(rna, np.array([2.0, 1.0]))
    assert out[0] < rna[0] and out[1] > rna[1]
    assert out.sum() == pytest.approx(1.0)


def test_calibration_is_measured_and_overconfidence_is_flagged(problem):
    """
    Realised coverage must be MEASURED and, when poor, flagged in words.

    Note what this does not assert. Credible intervals here describe uncertainty about
    the composition *given the reference*. The dominant error against known truth is
    systematic per-type bias from reference mismatch — on this fixture roughly -0.13
    for Tumor and +0.15 for Macrophage_Microglia — and no posterior over w conditional
    on the reference can cover a bias. Demanding nominal coverage would be demanding
    something the model does not claim; the honest requirement is that the shortfall is
    measured and stated, which is what calibration.py exists to do.
    """
    from ivygap.bench.calibration import measure_coverage

    data, test_set = problem
    m = BayesianDeconvolution(n_iter=300, n_burnin=100)
    m.fit_predict(data)
    result = measure_coverage(m, test_set.truth)

    assert result["uncertainty_reported"] is True
    assert 0.0 <= result["coverage_95"] <= 1.0
    assert result["mean_ci_width"] > 0
    assert np.isfinite(result["interval_score"])
    # Whatever the coverage, the verdict must describe it rather than stay silent.
    if result["coverage_95"] < 0.70:
        assert "OVERCONFIDENT" in result["verdict"]


def test_reference_uncertainty_propagation_widens_intervals(problem):
    """
    Drawing the signature matrix each sweep must widen the posterior, not leave it
    unchanged. The effect is modest — averaging over hundreds of genes cancels most
    per-gene reference noise — but a propagation step that changes nothing is not
    propagating anything.
    """
    data, _ = problem
    kw = dict(n_iter=300, n_burnin=100)

    off = BayesianDeconvolution(propagate_reference_uncertainty=False, **kw)
    off.fit_predict(data)
    w_off = (off.posterior_["ci_upper"] - off.posterior_["ci_lower"]).mean().mean()

    on = BayesianDeconvolution(propagate_reference_uncertainty=True, **kw)
    on.fit_predict(data)
    w_on = (on.posterior_["ci_upper"] - on.posterior_["ci_lower"]).mean().mean()

    assert w_on > w_off


def test_point_methods_are_listed_not_dropped(problem):
    """'This method does not quantify uncertainty' is a finding, not a gap in the table."""
    from ivygap.bench.calibration import measure_coverage

    data, test_set = problem
    m = NNLSDeconvolution()
    m.fit_predict(data)
    assert measure_coverage(m, test_set.truth) == {
        "method": "nnls", "uncertainty_reported": False}


# --- is the benchmark's winner actually winning? -----------------------------

def test_tie_detection_refuses_to_resolve_with_too_few_groups():
    """
    With two held-out donors a paired bootstrap has nothing to resample. Reporting a
    winner anyway would dress a coin flip as a decision — the same error the survival
    stage already refuses to make, applied to the criterion that actually selects.
    """
    from ivygap.bench.metrics import paired_bootstrap_ties

    rng = np.random.default_rng(0)
    idx = [f"s{i}" for i in range(20)]
    truth = pd.DataFrame(rng.dirichlet(np.ones(8), 20), index=idx,
                         columns=config.CELL_TYPES)
    ests = {f"m{k}": truth + rng.normal(0, 0.01 * (k + 1), truth.shape)
            for k in range(3)}
    donors = pd.Series(["D0"] * 10 + ["D1"] * 10, index=idx)

    out = paired_bootstrap_ties(ests, truth, groups=donors, n_boot=50)
    assert out["tied_with_best"].all()
    assert "cannot resolve" in out["note"].iloc[0]


def test_tie_detection_separates_clearly_different_methods():
    """The flip side: with enough groups and a real gap, it must NOT call everything tied."""
    from ivygap.bench.metrics import paired_bootstrap_ties

    rng = np.random.default_rng(1)
    idx = [f"s{i}" for i in range(60)]
    truth = pd.DataFrame(rng.dirichlet(np.ones(8), 60), index=idx,
                         columns=config.CELL_TYPES)
    ests = {
        "good": truth + rng.normal(0, 0.005, truth.shape),
        "bad": truth + rng.normal(0, 0.10, truth.shape),
    }
    donors = pd.Series([f"D{i % 8}" for i in range(60)], index=idx)

    out = paired_bootstrap_ties(ests, truth, groups=donors, n_boot=300)
    assert out.index[0] == "good"
    assert bool(out.loc["good", "tied_with_best"])
    assert not bool(out.loc["bad", "tied_with_best"])


def test_tie_detection_is_paired():
    """
    Two methods differing by a small constant on every sample must be distinguishable
    even when the samples themselves vary hugely. An unpaired bootstrap would drown that
    difference in between-sample variance that cancels in the comparison.
    """
    from ivygap.bench.metrics import paired_bootstrap_ties

    rng = np.random.default_rng(2)
    idx = [f"s{i}" for i in range(60)]
    truth = pd.DataFrame(rng.dirichlet(np.ones(8) * 0.4, 60), index=idx,
                         columns=config.CELL_TYPES)
    # Same noise realisation for both, scaled up for b. b's per-sample error is then
    # exactly 1.5x a's on EVERY sample — a perfectly consistent paired difference — while
    # the between-sample spread is large. An unpaired bootstrap drowns the former in the
    # latter; a paired one sees it immediately.
    shared_noise = rng.normal(0, 0.08, truth.shape)
    ests = {
        "a": truth + shared_noise,
        "b": truth + shared_noise * 1.5,
    }
    donors = pd.Series([f"D{i % 8}" for i in range(60)], index=idx)

    out = paired_bootstrap_ties(ests, truth, groups=donors, n_boot=300)
    assert out.index[0] == "a"
    assert not bool(out.loc["b", "tied_with_best"]), (
        "a consistent per-sample difference was lost — the bootstrap is not paired"
    )
