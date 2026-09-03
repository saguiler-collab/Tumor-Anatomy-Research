"""
controls.py — the two negative controls the protocol requires.

WHY A BENCHMARK WITHOUT NEGATIVE CONTROLS PROVES NOTHING
--------------------------------------------------------
Suppose every real method scores ACS around 0.8. Is the constraint set detecting
competence, or is 0.8 simply what any non-degenerate composition table scores against
these constraints? Without something that is *known to contain no information*, there is
no way to tell, and a permissive constraint set would look like a validated one.

So two controls run through the identical pipeline, on the identical inputs, scored the
identical way:

  RandomFractions   Dirichlet draws. No expression data is consulted at all. This is the
                    floor: whatever it scores is what the constraint set gives away for
                    free.
  ShuffledSignature NNLS against a signature matrix whose gene labels have been permuted.
                    The solver, the bulk data and the optimisation are all real; only the
                    correspondence between genes and cell types is destroyed. This is the
                    sharper control, because it produces structured, non-uniform,
                    bulk-dependent output that merely means nothing.

If a control scores well, the protocol's instruction is explicit and is followed:
**report it, publish the constraint file anyway, and state what a harder set would need.
Do not retune and re-report.** Loosening a control or tightening a constraint after
seeing this result would convert the whole design back into the circular reasoning it
exists to replace.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import nnls

from ivygap import config
from .base import DeconvolutionInput, DeconvolutionMethod


class RandomFractionsControl(DeconvolutionMethod):
    """
    Composition drawn from a Dirichlet, ignoring the expression data entirely.

    The concentration is deliberately not 1.0. A flat Dirichlet produces compositions
    that are on average uniform across cell types, which is an unrealistically *easy*
    thing for a constraint set to reject. Drawing at 0.7 gives the lumpy, sometimes
    near-degenerate compositions that a real but wrong method might produce, which makes
    this a fairer floor.
    """

    name = "control_random"
    family = "negative-control"

    def __init__(self, concentration: float = 0.7, **params):
        super().__init__(concentration=concentration, **params)
        self.concentration = concentration

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        rng = np.random.default_rng(config.RANDOM_SEED)
        n_s, n_t = len(data.samples), len(data.cell_types)
        return rng.dirichlet(np.full(n_t, self.concentration), size=n_s)


class ShuffledSignatureControl(DeconvolutionMethod):
    """
    Real NNLS against a signature whose gene rows have been permuted.

    Every ingredient except meaning is preserved: the same bulk matrix, the same solver,
    the same non-negativity constraint, the same marginal distribution of expression
    values in the signature. What is destroyed is which gene belongs to which cell type.

    Rows are permuted rather than entries because permuting entries independently would
    also destroy each gene's overall expression level, making the signature obviously
    pathological and the control too easy to beat. Row permutation keeps the signature
    statistically plausible and only breaks the gene-to-type correspondence — which is
    precisely the information a deconvolution method is supposed to be using.
    """

    name = "control_shuffled_signature"
    family = "negative-control"

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        rng = np.random.default_rng(config.RANDOM_SEED)
        S, B = self._as_arrays(data)
        S_shuffled = S[rng.permutation(S.shape[0]), :]
        return np.vstack([nnls(S_shuffled, B[:, j])[0] for j in range(B.shape[1])])


def build_controls() -> list[DeconvolutionMethod]:
    return [RandomFractionsControl(), ShuffledSignatureControl()]


CONTROL_NAMES = frozenset({RandomFractionsControl.name, ShuffledSignatureControl.name})
