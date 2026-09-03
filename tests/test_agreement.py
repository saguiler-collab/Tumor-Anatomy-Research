"""
Tests for the study's primary result and for the two negative controls.

The most important property here is that the NULL outcome is reported as a finding
rather than as a failure. A codebase that can only express "it worked" will quietly
pressure whoever runs it toward the answer that renders cleanly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.anatomic import agreement as ag
from ivygap.deconv.controls import (CONTROL_NAMES, RandomFractionsControl,
                                    ShuffledSignatureControl)


# --- the agreement test -----------------------------------------------------

def _scores(order):
    """ACS scores descending, so `order` is the ACS ranking best-to-worst."""
    return {m: 1.0 - i * 0.05 for i, m in enumerate(order)}


METHODS = [f"m{i}" for i in range(8)]


def test_perfect_agreement_is_the_headline():
    acs = _scores(METHODS)
    mae = {m: i * 0.01 for i, m in enumerate(METHODS)}      # same order, lower is better
    r = ag.test_agreement(acs, mae, "synthetic_mixtures")
    assert r.rho == pytest.approx(1.0)
    assert r.meets_threshold and r.ci_excludes_zero
    assert r.verdict.startswith("HEADLINE")


def test_no_agreement_is_reported_as_a_null_result_not_a_failure():
    """The branch the protocol says is 'arguably more useful'."""
    acs = _scores(METHODS)
    mae = {m: i * 0.01 for i, m in enumerate(reversed(METHODS))}   # exactly inverted
    r = ag.test_agreement(acs, mae, "synthetic_mixtures")
    assert r.rho == pytest.approx(-1.0)
    assert not r.meets_threshold
    assert r.verdict.startswith("NULL RESULT")
    assert "NOT evidence" in r.verdict


def test_sign_convention_is_not_inverted():
    """
    MAE is an error metric: lower is better. Getting `higher_truth_is_better` wrong would
    turn a perfect result into a perfectly negative one, and the mistake would look like
    a finding.
    """
    acs = _scores(METHODS)
    mae = {m: i * 0.01 for i, m in enumerate(METHODS)}
    assert ag.test_agreement(acs, mae, "y", higher_truth_is_better=False).rho > 0
    assert ag.test_agreement(acs, mae, "y", higher_truth_is_better=True).rho < 0


def test_threshold_met_but_ci_spanning_zero_does_not_pass():
    """The pre-registered bar required BOTH conditions; meeting one is not meeting it."""
    v = ag.interpret(0.7, ci_low=-0.2, ci_high=0.95, n_methods=8)
    assert v.startswith("BELOW BAR")


def test_too_few_methods_is_flagged_not_silently_reported():
    acs = _scores(METHODS[:4])
    mae = {m: i * 0.01 for i, m in enumerate(METHODS[:4])}
    r = ag.test_agreement(acs, mae, "synthetic_mixtures")
    assert r.verdict.startswith("NOT INTERPRETABLE")


def test_controls_are_excluded_from_the_correlation():
    """
    Leaving controls in would add two points every yardstick agrees are terrible, which
    inflates rho without any evidence that ACS discriminates among real methods.
    """
    acs = {**_scores(METHODS), "control_random": 0.05}
    mae = {**{m: i * 0.01 for i, m in enumerate(METHODS)}, "control_random": 0.9}
    with_ctrl = ag.test_agreement(acs, mae, "y", exclude_controls=False)
    without = ag.test_agreement(acs, mae, "y", exclude_controls=True)
    assert with_ctrl.n_methods == without.n_methods + 1
    assert "control_random" not in without.method_ranks


def test_missing_yardstick_is_unavailable_not_disagreeing():
    """'We could not obtain it' and 'it disagreed' must not look the same."""
    table = ag.run_all_yardsticks(_scores(METHODS), {"synthetic_mixtures": {}})
    assert table.loc["absolute_purity", "verdict"].startswith("UNAVAILABLE")
    assert table.loc["synthetic_mixtures", "verdict"].startswith("UNAVAILABLE")


def test_preregistered_threshold_is_a_constant_not_an_argument():
    """If the bar could be passed in per call, it would not be pre-registered."""
    import inspect
    assert ag.PREREGISTERED_RHO_THRESHOLD == 0.60
    assert "threshold" not in inspect.signature(ag.interpret).parameters


# --- the negative controls --------------------------------------------------

def test_random_control_ignores_the_expression_data(problem):
    """If it consulted the data it would not be a floor."""
    data, _ = problem
    est_a = RandomFractionsControl().fit_predict(data)

    perturbed = data.bulk * 3.7
    from ivygap.deconv.base import DeconvolutionInput
    other = DeconvolutionInput(bulk=perturbed, references=data.references,
                               manifest=data.manifest)
    est_b = RandomFractionsControl().fit_predict(other)
    pd.testing.assert_frame_equal(est_a, est_b)


def test_shuffled_signature_control_still_produces_valid_output(problem):
    """
    The control must be a well-formed composition table — otherwise it fails for the
    wrong reason and tests nothing about the constraint set.
    """
    data, _ = problem
    est = ShuffledSignatureControl().fit_predict(data)
    assert list(est.columns) == list(config.CELL_TYPES)
    sums = est.sum(axis=1)
    assert np.allclose(sums[np.isfinite(sums)], 1.0, atol=1e-6)


def test_shuffled_signature_differs_from_the_real_solver(problem):
    """A 'shuffle' that changes nothing is not a control."""
    from ivygap.deconv.classical import NNLSDeconvolution
    data, _ = problem
    real = NNLSDeconvolution().fit_predict(data)
    shuffled = ShuffledSignatureControl().fit_predict(data)
    assert (real - shuffled).abs().mean().mean() > 0.01


def test_control_names_are_registered():
    assert {"control_random", "control_shuffled_signature"} == set(CONTROL_NAMES)
