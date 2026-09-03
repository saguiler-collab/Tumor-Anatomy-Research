"""
The equal-footing contract is the project's central claim. These tests exist to make
sure the guard actually trips — a guard that silently passes everything is worse than
no guard, because it converts an unverified assumption into a printed certificate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.bench.equal_footing import (EqualFootingViolation, assign_patient_folds,
                                        build_certificate, check_estimates_aligned,
                                        patient_equal_mean, require_comparable)


def test_certificate_is_reproducible(problem):
    data, _ = problem
    assert build_certificate(data).to_json() == build_certificate(data).to_json()


def test_certificate_detects_changed_input(problem):
    data, _ = problem
    a = build_certificate(data)

    perturbed = data.bulk.copy()
    perturbed.iloc[0, 0] += 1.0
    from ivygap.deconv.base import DeconvolutionInput
    b = build_certificate(DeconvolutionInput(
        bulk=perturbed, references=data.references, manifest=data.manifest))

    assert a.comparable_with(b), "a changed bulk matrix must break comparability"
    with pytest.raises(EqualFootingViolation):
        require_comparable({"a": a, "b": b})


def test_outcome_column_in_manifest_is_rejected(problem):
    """The one violation no downstream check could catch: a method seeing the outcome."""
    data, _ = problem
    leaky = data.manifest.copy()
    leaky["survival_days"] = 100.0
    from ivygap.deconv.base import DeconvolutionInput
    with pytest.raises(EqualFootingViolation, match="outcome"):
        build_certificate(DeconvolutionInput(
            bulk=data.bulk, references=data.references, manifest=leaky))


def test_misaligned_estimates_are_rejected(problem):
    data, _ = problem
    a = pd.DataFrame(0.1, index=data.samples, columns=list(data.cell_types))
    b = a.iloc[:-3]                                       # one method dropped samples
    with pytest.raises(EqualFootingViolation, match="different samples"):
        check_estimates_aligned({"a": a, "b": b})


def test_reordered_cell_types_are_rejected(problem):
    data, _ = problem
    a = pd.DataFrame(0.1, index=data.samples, columns=list(data.cell_types))
    b = a[list(reversed(a.columns))]
    with pytest.raises(EqualFootingViolation):
        check_estimates_aligned({"a": a, "b": b})


def test_folds_never_split_a_patient():
    """The single most important property of the fold assignment."""
    manifest = pd.DataFrame({
        "patient_id": ["P1"] * 12 + ["P2"] * 2 + ["P3"] * 7 + ["P4"] * 3 + ["P5"] * 5,
    })
    folds = assign_patient_folds(manifest, n_folds=3)
    per_patient = pd.DataFrame({"p": manifest["patient_id"], "f": folds})
    assert (per_patient.groupby("p")["f"].nunique() == 1).all()


def test_folds_are_deterministic():
    manifest = pd.DataFrame({"patient_id": [f"P{i%7}" for i in range(40)]})
    assert assign_patient_folds(manifest).equals(assign_patient_folds(manifest))


def test_patient_equal_mean_ignores_block_count():
    """
    The nested-design hazard, stated as a test.

    P1 contributes ten samples at 1.0; P2 contributes one at 0.0. A plain mean gives
    0.909 and is essentially P1's answer. Patient-equal weighting gives 0.5.
    """
    values = pd.Series([1.0] * 10 + [0.0])
    patients = pd.Series(["P1"] * 10 + ["P2"])
    assert values.mean() == pytest.approx(0.909, abs=1e-3)
    assert patient_equal_mean(values, patients) == pytest.approx(0.5)


def test_patient_equal_mean_drops_all_nan_patients():
    values = pd.Series([1.0, np.nan, np.nan])
    patients = pd.Series(["P1", "P2", "P2"])
    assert patient_equal_mean(values, patients) == pytest.approx(1.0)
