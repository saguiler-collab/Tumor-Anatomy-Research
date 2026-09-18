"""Negative controls for the methylation-vs-estimate lymphoid comparison.

`docs/IMMUNE_ARM.md` §6 claims that comparing a methylation leukocyte SUB-composition against a
deconvolution ALL-CELL fraction is legitimate because "the denominator difference is removed by
construction". That is a claim about code, so it is tested here rather than asserted in prose.

The inverted control matters more than the planted one: a comparison that reports agreement on
both truth and its reverse is measuring nothing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from lymphoid_ordering import COMPARE, k4s, ordering_of, renormalise  # noqa: E402

TRUTH = pd.DataFrame({"T_cell": [0.50, 0.60, 0.40],
                      "B_cell": [0.20, 0.10, 0.25],
                      "NK_cell": [0.30, 0.30, 0.35]},
                     index=["s1", "s2", "s3"])


def test_planted_ordering_is_recovered():
    assert ordering_of(renormalise(TRUTH).mean()) == "T>NK>B"


def test_inverted_ordering_is_rejected():
    """The negative control. Reversing the columns must reverse the verdict."""
    inv = TRUTH.copy()
    inv[["T_cell", "B_cell"]] = inv[["B_cell", "T_cell"]].to_numpy()
    assert ordering_of(renormalise(inv).mean()) == "B>NK>T"
    assert ordering_of(renormalise(TRUTH).mean()) != ordering_of(renormalise(inv).mean())


def test_denominator_is_removed_by_construction():
    """THE CLAIM UNDER TEST.

    A deconvolution estimate carries Tumor and Macrophage columns that the methylation
    reference cannot see, and the methylation estimate carries Mono/Neutro/Eosino that the
    deconvolution roster does not model. If `renormalise` truly cancels the denominator, adding
    an arbitrarily large foreign column must not move the lymphoid ordering AT ALL.
    """
    base = renormalise(TRUTH)
    for foreign in (0.001, 1.0, 50.0, 1e4):
        wide = TRUTH.copy()
        wide["Tumor"] = foreign
        wide["Macrophage_Microglia"] = foreign * 3
        got = renormalise(wide)
        assert np.allclose(got.to_numpy(), base.to_numpy()), f"moved at Tumor={foreign}"
        assert ordering_of(got.mean()) == "T>NK>B"


def test_per_method_scaling_does_not_change_the_verdict():
    """A method that is uniformly under- or over-confident must not be scored differently.

    The comparison is about WHICH lymphocyte dominates, not about absolute scale, so scaling a
    method's whole output by a constant is required to be a no-op.
    """
    base = renormalise(TRUTH)
    for k in (0.01, 0.5, 2.0, 1000.0):
        assert np.allclose(renormalise(TRUTH * k).to_numpy(), base.to_numpy())


def test_all_zero_lymphoid_row_becomes_nan_not_zero():
    """A method that assigns no lymphoid signal must not be silently scored as B>NK>T.

    0/0 must propagate as NaN so the row is excluded, rather than dividing by a replaced 0 and
    inventing an ordering out of nothing.
    """
    z = pd.DataFrame({c: [0.0] for c in COMPARE}, index=["s1"])
    assert renormalise(z).isna().all(axis=None)


def test_k4s_strips_the_aliquot_letter_only():
    assert k4s("TCGA-FG-6692-01A") == "TCGA-FG-6692-01"
    assert k4s("TCGA-CS-4938-01") == "TCGA-CS-4938-01"
    # a barcode with no sample field is returned untouched rather than truncated
    assert k4s("TCGA-CS-4938") == "TCGA-CS-4938"


def test_renormalise_rejects_a_missing_type():
    """Never impute a missing input: a roster without NK must raise, not fill a zero."""
    with pytest.raises(KeyError):
        renormalise(TRUTH.drop(columns=["NK_cell"]))


def test_recovery_is_bias_invariant():
    """The ranking must not be the calibration shift in disguise.

    `docs/EQUAL_FOOTING.md` reports that the rebuilt reference reshuffles the ranking
    (tau +0.214) AND worsens the median bias (-0.028 -> -0.472). Those are only separable
    claims if the recovery statistic ignores a constant offset. Recovery is `1 + slope` of
    (estimate - truth) on truth, so it must be exactly invariant to adding a constant.
    """
    rng = np.random.default_rng(0)
    truth = rng.uniform(0.2, 1.0, 400)
    est = 0.4 * truth + 0.1 + rng.normal(0, 0.02, 400)

    def recovery(e):
        return 1.0 + np.polyfit(truth, e - truth, 1)[0]

    base = recovery(est)
    for shift in (-0.5, -0.3, 0.0, 0.25, 0.9):
        assert np.isclose(recovery(est + shift), base, atol=1e-9), shift
    # and it is NOT invariant to a change of slope, or it would measure nothing
    assert not np.isclose(recovery(est * 2.0), base, atol=1e-6)
