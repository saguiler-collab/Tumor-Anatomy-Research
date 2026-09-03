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


# =============================================================================
# donor leakage through the R cell-level export
# -----------------------------------------------------------------------------
# The R-backed methods do not read the ReferenceBundle. They read the cell-level export,
# which r_bridge writes from whatever cell source is registered. `run_benchmark` builds
# its reference from TRAINING donors only, but a caller that registered the full atlas
# would leave every donor visible — so MuSiC, Bisque and SCDC would see the held-out
# donors whose cells were pooled into the test mixtures, while NNLS and SVR saw only the
# training signature. The R methods would win on their own reference.
# =============================================================================

def test_benchmark_registers_training_donors_only(monkeypatch):
    """
    The guard. After the donor split, the registered cell source must contain no
    held-out donor — otherwise the mixtures are scored against their own cells.
    """
    from ivygap.data.reference import build_synthetic
    from ivygap.bench import run_benchmark
    from ivygap.bench import pseudobulk as pb
    from ivygap.deconv import r_bridge
    from ivygap import config as cfg

    # The gene floor guards real runs; a fixture only needs enough genes to solve.
    monkeypatch.setattr(cfg, "MIN_GENES_SHARED", 20)
    _, expression, meta = build_synthetic(n_genes=200, n_donors=6,
                                          cells_per_type_per_donor=8)
    # Register the FULL source, exactly as run_all does before the benchmark.
    r_bridge.set_cell_source(cfg.PRIMARY_REFERENCE, expression, meta)

    seen = {}
    real_split = pb.build_train_test

    def spy_split(expr, m, **kw):
        out = real_split(expr, m, **kw)
        seen["train"], seen["test"] = set(out[1]), set(out[2])
        return out

    monkeypatch.setattr(pb, "build_train_test", spy_split)

    # Record EVERY registration. The last one is the restore of the caller's full
    # source, which legitimately contains held-out donors; the one that matters is the
    # registration in force while the methods are fitted.
    calls: list[set] = []
    real_set = r_bridge.set_cell_source

    def spy_set(name, e, mm):
        calls.append(set(mm["donor"].astype(str)))
        return real_set(name, e, mm)

    monkeypatch.setattr(r_bridge, "set_cell_source", spy_set)

    run_benchmark.run(expression, meta, prefer_r=False, n_test=12,
                      n_signature_genes=10, verbose=False)

    assert seen["test"], "the fixture produced no held-out donors"
    assert len(calls) >= 2, "the benchmark did not register a cell source of its own"

    during_fit = calls[-2]                  # set by run_benchmark, before the restore
    leaked = during_fit & seen["test"]
    assert not leaked, (
        f"held-out donors {sorted(leaked)} were visible to the R cell-level export; "
        f"the mixtures are pooled from exactly those cells")
    assert during_fit <= seen["train"]
    # and the restore really does put the caller's source back
    assert calls[-1] == set(meta["donor"].astype(str))


def test_benchmark_restores_the_callers_cell_source(monkeypatch):
    """A benchmark must not leave its training-only subset registered for later stages,
    which score real tissue and legitimately want the whole reference."""
    from ivygap.data.reference import build_synthetic
    from ivygap.bench import run_benchmark
    from ivygap.deconv import r_bridge
    from ivygap import config as cfg

    monkeypatch.setattr(cfg, "MIN_GENES_SHARED", 20)
    _, expression, meta = build_synthetic(n_genes=200, n_donors=6,
                                          cells_per_type_per_donor=8)
    r_bridge.set_cell_source(cfg.PRIMARY_REFERENCE, expression, meta)
    before = r_bridge.get_cell_source(cfg.PRIMARY_REFERENCE)

    run_benchmark.run(expression, meta, prefer_r=False, n_test=12,
                      n_signature_genes=10, verbose=False)

    after = r_bridge.get_cell_source(cfg.PRIMARY_REFERENCE)
    assert after is not None
    assert after[1].shape == before[1].shape, "the caller's cell source was not restored"
    r_bridge.clear_cell_source(cfg.PRIMARY_REFERENCE)
