"""
`fit_predict` and `scripts/remeasure_method.py` must post-process identically.

They used not to. `remeasure_method.py` reimplemented the three steps -- simplex projection,
cell-size conversion, frame assembly -- and applied the conversion unconditionally. When the
conversion was made conditional for Bisque, the only package measured to convert cell size
itself, the fix reached `fit_predict` and not the copy. The re-measurement kept
double-correcting and wrote an estimate table BIT-IDENTICAL to the uncorrected one, so the
fix appeared to have been applied and had not been. Both now call `finalize_estimates`.
"""
from __future__ import annotations

import numpy as np
import pytest

from ivygap.deconv import r_bridge
from ivygap.deconv.base import DeconvolutionMethod, finalize_estimates


class _Echo(DeconvolutionMethod):
    name = "echo"
    family = "test"

    def __init__(self, raw, returns_cells=False, **kw):
        super().__init__(**kw)
        self.raw = np.asarray(raw, dtype="float64")
        self._rc = returns_cells

    @property
    def returns_cell_fractions(self) -> bool:
        return self._rc

    def _solve_all(self, data):
        return self.raw


@pytest.fixture
def raw(problem):
    data, _ = problem
    rng = np.random.default_rng(3)
    v = rng.random((len(data.samples), len(data.cell_types)))
    return data, v / v.sum(axis=1, keepdims=True)


@pytest.mark.parametrize("returns_cells", [False, True])
def test_the_two_paths_agree(raw, returns_cells):
    data, v = raw
    via_method = _Echo(v, returns_cells=returns_cells).fit_predict(data)
    via_helper = finalize_estimates(v, data, covered=None,
                                    returns_cell_fractions=returns_cells)
    assert np.allclose(via_method.to_numpy(), via_helper.to_numpy()), (
        "fit_predict and the shared helper disagree; a copy of the post-processing has "
        "reappeared")
    assert list(via_method.columns) == list(via_helper.columns)
    assert list(via_method.index) == list(via_helper.index)


def test_the_flag_actually_changes_the_output(raw):
    """If it did not, the Bisque fix would be undetectable -- which is how it was missed."""
    data, v = raw
    off = finalize_estimates(v, data, returns_cell_fractions=False)
    on = finalize_estimates(v, data, returns_cell_fractions=True)
    assert not np.allclose(off.to_numpy(), on.to_numpy()), (
        "skipping the conversion changed nothing, so this fixture cannot detect whether the "
        "flag is honoured. Either the reference's cell sizes are uniform or the flag is dead.")
    assert np.allclose(on.to_numpy(), v / v.sum(axis=1, keepdims=True)), (
        "with the flag set the output must be the method's own values, projected only")


def test_remeasure_reads_the_same_registry():
    """The re-measurement's decision must come from the measured registry, not a literal."""
    src = (__import__("pathlib").Path(__file__).resolve().parent.parent
           / "scripts" / "remeasure_method.py").read_text()
    assert "r_bridge.R_RETURNS_CELL_FRACTIONS" in src, (
        "remeasure_method.py must consult the measured registry")
    assert "finalize_estimates" in src, "it must use the shared post-processing"
    assert "bisque" not in src.lower().split("r_returns_cell_fractions")[0][-400:], (
        "no method should be hard-coded at the decision site")
    assert "bisque" in r_bridge.R_RETURNS_CELL_FRACTIONS
