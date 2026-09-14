"""
The central cell-size conversion must be skipped for a package that already applied it.

D1, resolved 2026-09-13 after three wrong turns. Measured on a probe mixing two cell types
differing 3x in mRNA at equal cell counts, so 0.750 is an mRNA share and 0.500 a cell share:
MuSiC, SCDC, EPIC and BayesPrism all return ~0.750; **Bisque returns 0.5000**. Applying this
project's conversion on top of Bisque's own applied it twice.

Two things these tests pin, and the second matters as much as the first:

  * the skip happens for the measured package;
  * it does NOT happen on the Python fallback path, which has never been measured on that
    probe. Assuming the reimplementation shares the R package's convention would be the same
    unmeasured inference D1 was corrected for twice.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap.deconv import r_bridge
from ivygap.deconv.base import DeconvolutionMethod, to_cell_fractions
from ivygap.deconv.registry import build_methods


class _Fixed(DeconvolutionMethod):
    """Returns a known mRNA-share row so the conversion's effect is visible."""
    name = "fixed"
    family = "test"

    def __init__(self, row, **kw):
        super().__init__(**kw)
        self.row = np.asarray(row, dtype="float64")

    def _solve_all(self, data):
        return np.tile(self.row, (len(data.samples), 1))


def test_default_is_to_apply_the_conversion(problem):
    data, _ = problem
    n = len(data.cell_types)
    row = np.full(n, 1.0 / n)
    m = _Fixed(row)
    assert m.returns_cell_fractions is False, "the default must be to apply ours"
    out = m.fit_predict(data)

    cs = data.primary.cell_size.reindex(list(data.cell_types)).to_numpy()
    expected = to_cell_fractions(row, cs)
    assert np.allclose(out.iloc[0].to_numpy(), expected), (
        "an unflagged method must receive the central conversion")


def test_flagged_method_is_left_alone(problem):
    data, _ = problem
    n = len(data.cell_types)
    row = np.full(n, 1.0 / n)

    class _AlreadyCells(_Fixed):
        name = "already_cells"
        returns_cell_fractions = True

    out = _AlreadyCells(row).fit_predict(data)
    assert np.allclose(out.iloc[0].to_numpy(), row), (
        "a method that already returns cell fractions must NOT be converted again")

    # And the two paths must actually differ on this reference, or the test proves nothing.
    cs = data.primary.cell_size.reindex(list(data.cell_types)).to_numpy()
    assert not np.allclose(to_cell_fractions(row, cs), row), (
        "the reference's cell sizes are too uniform for this test to discriminate")


def test_bisque_is_the_only_flagged_r_package():
    assert r_bridge.R_RETURNS_CELL_FRACTIONS == frozenset({"bisque"}), (
        "the flagged set changed. It may only be changed on evidence from "
        "scripts/verify_cell_size_semantics.py, never to move a score.")
    for m in ("music", "scdc", "scdc_ensemble", "epic", "bayesprism", "dwls", "quantiseq"):
        assert m not in r_bridge.R_RETURNS_CELL_FRACTIONS


def test_the_flag_follows_the_implementation_that_ran():
    """R:BisqueRNA is measured; the Python reimplementation is not, so it must not inherit."""
    wrapper = next(m for m in build_methods(prefer_r=True) if m.name == "bisque")

    wrapper.implementation_ = "R:BisqueRNA"
    assert wrapper.returns_cell_fractions is True, (
        "the genuine package is measured to return cell fractions and must be skipped")

    wrapper.implementation_ = "python-reimplementation"
    assert wrapper.returns_cell_fractions is False, (
        "the Python fallback has NOT been measured on the cell-size probe. Flagging it "
        "would assert a convention nobody measured, which is the error D1 records twice.")

    wrapper.implementation_ = "unknown"
    assert wrapper.returns_cell_fractions is False


def test_a_non_bisque_wrapper_is_never_flagged():
    for name in ("music", "scdc", "epic"):
        w = next((m for m in build_methods(prefer_r=True) if m.name == name), None)
        if w is None or not hasattr(w, "implementation_"):
            continue
        w.implementation_ = f"R:{name}"
        assert w.returns_cell_fractions is False, (
            f"{name} is measured to return an mRNA share; ours is its FIRST conversion")
