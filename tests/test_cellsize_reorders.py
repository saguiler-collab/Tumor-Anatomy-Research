"""
The cell-size conversion CAN reorder a cell type across samples.

This exists because the opposite was asserted, repeatedly, in project documentation and
in conversation: that applying the correction twice "cannot change ACS, because a per-type
constant cannot reorder a column". The argument was taken from the shape of the formula
and never tested. It is wrong, and these tests pin why, so it cannot be re-asserted.

`to_cell_fractions` does two things:

    out = rna_fractions / cell_size     # a per-column constant — order-preserving
    return project_to_simplex(out)      # renormalise to sum 1 — NOT order-preserving

The renormalisation divisor is sum_j(w_j / c_j), which depends on the sample's own
composition. The effective scaling of type k is therefore (1/c_k) / sum_j(w_j/c_j), and
that denominator differs between samples.
"""
from __future__ import annotations

import numpy as np

from ivygap.deconv.base import to_cell_fractions


def test_the_division_alone_preserves_order():
    """The half of the argument that was right."""
    cs = np.array([4.0, 1.0, 2.0])
    a = np.array([0.5, 0.3, 0.2])
    b = np.array([0.2, 0.5, 0.3])
    for k in range(3):
        before = a[k] > b[k]
        after = (a / cs)[k] > (b / cs)[k]
        assert before == after, "dividing by a per-type constant must preserve order"


def test_the_renormalisation_can_reverse_an_ordering():
    """
    The half that was not. Two samples where cell type 0 is higher in A than in B before
    conversion, and lower after — because A and B renormalise by different totals.
    """
    # Found by searching the simplex rather than by hand — a first hand-picked pair did
    # NOT flip, which is itself worth knowing: flips are real but not ubiquitous.
    cs = np.array([10.0, 1.0, 1.0])
    a = np.array([0.4488, 0.0710, 0.4802])
    b = np.array([0.2014, 0.2377, 0.5609])

    assert a[2] < b[2], "before conversion, type 2 is LOWER in A"
    ca, cb = to_cell_fractions(a, cs), to_cell_fractions(b, cs)
    assert ca[2] > cb[2], (
        f"after conversion it must be HIGHER in A; got {ca[2]:.4f} vs {cb[2]:.4f}. "
        f"If this stops flipping, the renormalisation has changed and D1's severity "
        f"needs revisiting.")


def test_a_uniform_cell_size_is_a_no_op_on_ordering():
    """Sanity: with equal sizes there is nothing to reorder."""
    cs = np.array([2.0, 2.0, 2.0])
    a = np.array([0.5, 0.3, 0.2]); b = np.array([0.2, 0.5, 0.3])
    ca, cb = to_cell_fractions(a, cs), to_cell_fractions(b, cs)
    for k in range(3):
        assert (a[k] > b[k]) == (ca[k] > cb[k])


def test_double_correction_is_not_idempotent():
    """
    Applying the conversion twice is not the same as applying it once — which is the
    whole of defect D1. If this ever passes as equal, the conversion has been made
    idempotent and D1's severity needs revisiting.
    """
    cs = np.array([10.0, 1.0, 1.0])
    w = np.array([0.40, 0.55, 0.05])
    once = to_cell_fractions(w, cs)
    twice = to_cell_fractions(once, cs)
    assert not np.allclose(once, twice), "double correction must differ from single"
