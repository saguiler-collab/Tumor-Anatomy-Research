"""
Tests for the simulated ground-truth yardstick.

The load-bearing test is `test_zero_perturbation_recovers_truth`. It caught a real
defect: the generative model originally divided the signature columns by their totals,
which put `n_k * L_k / colsum_k` into the mixture while the solver estimates `n_k * L_k`.
Since colsum varies by up to 1.33x across cell types, every method inherited a systematic
error of that size — and because it hit every method roughly equally, the resulting
leaderboard still looked entirely plausible.

A yardstick that is wrong in a way that looks reasonable is the worst kind, so the
zero-perturbation control is the thing to pin.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.bench import metrics
from ivygap.bench.pseudobulk import NICHE_TEMPLATES
from ivygap.bench.simulated_ground_truth import generate
from ivygap.deconv.base import DeconvolutionInput
from ivygap.deconv.classical import NNLSDeconvolution


@pytest.fixture(scope="module")
def reference():
    from ivygap.data.reference import load_frozen_reference, select_signature_genes
    if not config.FROZEN_SIGNATURE_PATH.exists():
        pytest.skip("frozen signature not present")
    ref = load_frozen_reference()
    return ref.subset_genes(select_signature_genes(ref, n_per_type=60))


def _fit(sim, reference):
    mf = pd.DataFrame({"patient_id": sim.donors.astype(str)},
                      index=sim.expression.columns)
    data = DeconvolutionInput(bulk=sim.expression, references=(reference,), manifest=mf)
    est = NNLSDeconvolution().fit_predict(data)
    return metrics.overall_metrics(est, sim.truth, sim.donors)


def test_output_contract(reference):
    sim = generate(reference, n_mixtures=30, n_donors=6)
    assert len(sim) > 0
    assert list(sim.truth.columns) == list(config.CELL_TYPES)
    assert np.allclose(sim.truth.sum(axis=1), 1.0)
    assert (sim.expression.to_numpy() >= 0).all()


def test_zero_perturbation_recovers_truth(reference):
    """
    With no donor perturbation the mixture is built from exactly the signature the
    solver uses, so recovery must be near-perfect. Anything above ~0.01 MAE means the
    generative model and the estimand disagree, and every downstream ranking is
    measuring that disagreement rather than the methods.
    """
    sim = generate(reference, n_mixtures=40, n_donors=4, donor_log_sd=0.0)
    m = _fit(sim, reference)
    assert m["mae_primary"] < 0.01, (
        f"zero-perturbation MAE is {m['mae_primary']:.4f}; the generative model does "
        f"not match the estimand the methods report"
    )


def test_perturbation_makes_the_problem_harder(reference):
    """
    The whole point: error must rise with donor mismatch. A yardstick where it does not
    is measuring sampling noise, and its ranking would say nothing about robustness.
    """
    easy = _fit(generate(reference, n_mixtures=40, n_donors=6, donor_log_sd=0.0),
                reference)["mae_primary"]
    mid = _fit(generate(reference, n_mixtures=40, n_donors=6, donor_log_sd=0.15),
               reference)["mae_primary"]
    hard = _fit(generate(reference, n_mixtures=40, n_donors=6, donor_log_sd=0.35),
                reference)["mae_primary"]
    assert easy < mid < hard


def test_niche_templates_match_the_structure_roster():
    """
    Regression: the templates kept the archive's CTmvp/CTpan names after the roster was
    renamed to MVP/PAN, so the per-niche breakdown was labelled with structures that no
    longer existed anywhere else in the project.
    """
    assert set(NICHE_TEMPLATES) == set(config.PRIMARY_STRUCTURES)


def test_provenance_declares_it_is_a_simulation(reference):
    """
    This is NOT the protocol's yardstick 1, and every artefact must say so. A simulated
    ranking presented as the real one would misrepresent the study's primary result.
    """
    prov = generate(reference, n_mixtures=20, n_donors=4).provenance()
    assert prov["is_simulation"] is True
    assert prov["yardstick"] == "simulated_donor_mismatch"
    assert "NOT the protocol's yardstick 1" in prov["warning"]


def test_donor_perturbation_preserves_column_totals(reference):
    """
    Renormalising after perturbation keeps reference mismatch from being confounded with
    cell size — otherwise a donor whose noise ran high would also have systematically
    larger cells and the yardstick would measure two things at once.
    """
    from ivygap.bench.simulated_ground_truth import _perturb
    rng = np.random.default_rng(0)
    profile = reference.profile.to_numpy(dtype="float64")
    perturbed = _perturb(profile, rng, log_sd=0.35)
    np.testing.assert_allclose(perturbed.sum(axis=0), profile.sum(axis=0), rtol=1e-9)
    assert not np.allclose(perturbed, profile), "perturbation did nothing"
