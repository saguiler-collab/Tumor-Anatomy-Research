"""
CIBERSORTx S-mode, against Supplementary Note 1 (p39).

S-mode adjusts the SIGNATURE rather than the mixtures. The properties tested here are
the ones that distinguish it from B-mode and from a wrong implementation of itself:
it must consume real single cells, it must actually change the signature, and it must
refuse to run rather than quietly degrade into something else.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap.deconv import r_bridge
from ivygap.deconv.base import DeconvolutionInput
from ivygap.deconv.classical import (CIBERSORTxDeconvolution,
                                     CIBERSORTxSModeDeconvolution)
from ivygap.data.reference import build_reference


@pytest.fixture
def cell_problem():
    """A small cell-level reference plus mixtures built from it."""
    from ivygap import config
    rng = np.random.default_rng(0)
    # build_reference always emits the full project roster, so the fixture must supply
    # cells for all of it — a partial roster silently yields zero columns.
    types = list(config.CELL_TYPES)
    genes = [f"G{i}" for i in range(len(types) * 10)]

    cells, meta = [], []
    for ti, t in enumerate(types):
        base = rng.random(len(genes)) * 0.2
        base[ti * 10:(ti + 1) * 10] += 3.0                   # a marker block per type
        for d in range(4):
            for j in range(25):
                cells.append(base * rng.uniform(0.8, 1.2, len(genes)))
                meta.append({"donor": f"D{d}", "cell_type": t})
    expr = pd.DataFrame(np.array(cells).T, index=genes,
                        columns=[f"c{i}" for i in range(len(cells))])
    expr = expr / expr.sum(axis=0) * 1e6
    cmeta = pd.DataFrame(meta, index=expr.columns)

    ref = build_reference(expr, cmeta, name="toy_sc")
    S = ref.profile.to_numpy()
    F = rng.dirichlet(np.ones(len(types)), size=12).T
    bulk = pd.DataFrame(S @ F, index=ref.profile.index,
                        columns=[f"S{i}" for i in range(12)])
    bulk = bulk / bulk.sum(axis=0) * 1e6
    manifest = pd.DataFrame({"patient_id": ["P0"] * 6 + ["P1"] * 6,
                             "structure": ["CT", "PAN"] * 6}, index=bulk.columns)
    data = DeconvolutionInput(bulk=bulk, references=(ref,), manifest=manifest)
    return data, expr, cmeta


def test_smode_refuses_to_run_without_the_cells(cell_problem):
    """
    The property that keeps S-mode honest. B-mode reconstructs mixtures from the
    signature alone; S-mode resamples real cells. Without them it must raise, not
    silently produce a B-mode result under S-mode's name.
    """
    data, _, _ = cell_problem
    r_bridge.clear_cell_source(data.primary.name)
    with pytest.raises(RuntimeError, match="no cell source is registered"):
        CIBERSORTxSModeDeconvolution().fit_predict(data)


def test_smode_runs_and_actually_adjusts_the_signature(cell_problem):
    data, expr, cmeta = cell_problem
    r_bridge.set_cell_source(data.primary.name, expr, cmeta)
    try:
        m = CIBERSORTxSModeDeconvolution()
        est = m.fit_predict(data)
    finally:
        r_bridge.clear_cell_source(data.primary.name)

    assert m.signature_adjusted_ is True
    assert m.n_artificial_ == m.N_ARTIFICIAL
    assert "S-mode" in (m.adjustment_note_ or "")
    assert est.shape == (len(data.samples), len(data.cell_types))
    assert np.isfinite(est.to_numpy()).all()
    assert np.allclose(est.sum(axis=1), 1.0)


def test_smode_and_bmode_are_different_methods(cell_problem):
    """
    If the two agree exactly, one of them is not doing what it says. They adjust
    opposite sides of the same equation.
    """
    data, expr, cmeta = cell_problem
    r_bridge.set_cell_source(data.primary.name, expr, cmeta)
    try:
        s = CIBERSORTxSModeDeconvolution().fit_predict(data)
    finally:
        r_bridge.clear_cell_source(data.primary.name)
    b = CIBERSORTxDeconvolution().fit_predict(data)
    assert not np.allclose(s.to_numpy(), b.to_numpy()), \
        "S-mode reproduced B-mode exactly; the signature was not adjusted"


def test_smode_needs_more_artificial_mixtures_than_cell_types(cell_problem):
    """
    Supp. Note 1: the per-gene NNLS in step 6 has a (k x c) design, so k must exceed c.
    The paper manufactures extra pseudo-mixtures when it does not; here k is fixed at
    100 against 8 roster types, so the guard is that the default cannot drift under it.
    """
    m = CIBERSORTxSModeDeconvolution()
    from ivygap import config
    assert m.N_ARTIFICIAL > len(config.CELL_TYPES)


def test_a_type_absent_from_the_cells_is_refused(cell_problem):
    """
    mu is a cell type's fractional abundance in R. A type with no cells has mu = 0, so
    it would never enter an artificial mixture and its column of the adjusted signature
    would be unconstrained. That must fail loudly, not produce a zero column.
    """
    data, expr, cmeta = cell_problem
    drop = cmeta["cell_type"] != "T_cell"
    r_bridge.set_cell_source(data.primary.name, expr.loc[:, drop.to_numpy()],
                             cmeta[drop])
    try:
        with pytest.raises(RuntimeError, match="no cells for"):
            CIBERSORTxSModeDeconvolution().fit_predict(data)
    finally:
        r_bridge.clear_cell_source(data.primary.name)


def test_smode_is_registered_under_its_own_name():
    """Both modes appear; neither stands in for the other."""
    from ivygap.deconv.registry import build_methods
    names = [m.name for m in build_methods(prefer_r=False)]
    assert "cibersortx" in names and "cibersortx_smode" in names


def test_smode_is_deterministic(cell_problem):
    """Seeded from config.RANDOM_SEED, so a rerun reproduces the same artefact."""
    data, expr, cmeta = cell_problem
    r_bridge.set_cell_source(data.primary.name, expr, cmeta)
    try:
        a = CIBERSORTxSModeDeconvolution().fit_predict(data)
        b = CIBERSORTxSModeDeconvolution().fit_predict(data)
    finally:
        r_bridge.clear_cell_source(data.primary.name)
    pd.testing.assert_frame_equal(a, b)
