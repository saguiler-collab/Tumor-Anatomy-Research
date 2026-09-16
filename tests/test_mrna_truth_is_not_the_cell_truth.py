"""
The mRNA truth and the cell truth must not be the same numbers.

docs/OPEN_DEFECTS.md D16. `truth_mrna` was added so Problem 1 (can the RNA contributions be
inferred?) could be scored separately from Problem 2 (can they be converted to cell
abundance?). It was computed by summing the mixture's own expression matrix over each type's
cells — and every column of that matrix is normalised to 1e6, so the sum over a type's cells
is exactly `1e6 * n_cells`. The mRNA share therefore equalled the cell share by construction,
bit for bit, and the two problems were the same problem wearing different labels.

The test that catches this is the one asserting the two truths DIFFER. A test that only
checked `truth_mrna` had the right shape and summed to one passed happily throughout.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.bench import pseudobulk as pb


def _cells(n_donors=6, per=40, seed=0):
    """A normalised cell matrix — every column sums to 1e6, as the real one does."""
    rng = np.random.default_rng(seed)
    types = list(config.CELL_TYPES)
    rows, ids = [], []
    for d in range(n_donors):
        for ct in types:
            for j in range(per):
                ids.append(f"D{d}|{ct}|{j}")
                rows.append((f"D{d}", ct))
    meta = pd.DataFrame(rows, columns=["donor", "cell_type"], index=ids)
    genes = [f"G{i}" for i in range(120)]
    X = rng.gamma(2.0, 50.0, size=(len(genes), len(ids)))
    # each cell type gets a real marker block, so the mixtures are not noise
    for k, ct in enumerate(types):
        X[k * 10:(k + 1) * 10, (meta["cell_type"] == ct).to_numpy()] *= 25
    expr = pd.DataFrame(X, index=genes, columns=ids)
    expr = expr.div(expr.sum(axis=0), axis=1) * 1e6           # the pipeline's normalisation
    # library sizes that genuinely differ BY TYPE, which is the whole point of mRNA share
    size = pd.Series({i: 3000.0 * (1 + 3.0 * (meta.loc[i, "cell_type"] == "Tumor"))
                      for i in ids})
    return expr, meta, size


def test_without_library_sizes_there_is_no_mrna_truth():
    """A missing input is reported as missing, never imputed from the normalised matrix."""
    expr, meta, _ = _cells()
    s = pb.generate(expr, meta, sorted(meta["donor"].unique()), 25, niche_fraction=0.5)
    assert s.truth_mrna is None, \
        "an mRNA truth was fabricated without the pre-normalisation totals it needs"


def test_the_two_truths_differ_when_types_carry_different_mrna():
    expr, meta, size = _cells()
    s = pb.generate(expr, meta, sorted(meta["donor"].unique()), 25,
                    niche_fraction=0.5, cell_mrna=size)
    assert s.truth_mrna is not None
    assert s.truth_mrna.shape == s.truth.shape
    np.testing.assert_allclose(s.truth_mrna.sum(axis=1), 1.0, atol=1e-12)

    gap = float(np.abs(s.truth_mrna.to_numpy() - s.truth.to_numpy()).max())
    assert gap > 0.01, (
        f"the mRNA truth is numerically the cell truth (max gap {gap:.2e}); it is being "
        f"derived from the normalised matrix again")

    # and it moves in the right DIRECTION: Tumor carries 4x the mRNA per cell here, so its
    # mRNA share must exceed its cell share wherever it is present at all.
    have = s.truth["Tumor"] > 0.02
    assert have.sum() > 5, "too few mixtures contain Tumor to test the direction"
    assert (s.truth_mrna.loc[have, "Tumor"] > s.truth.loc[have, "Tumor"]).all(), \
        "the type with more mRNA per cell did not gain share in the mRNA truth"


def test_equal_library_sizes_collapse_the_two_truths():
    """
    NEGATIVE CONTROL for the test above. With no per-type difference in mRNA content the
    two truths SHOULD coincide — so a passing `test_the_two_truths_differ` is detecting the
    library sizes, not just any change of code path.
    """
    expr, meta, _ = _cells()
    flat = pd.Series(3000.0, index=expr.columns)
    s = pb.generate(expr, meta, sorted(meta["donor"].unique()), 25,
                    niche_fraction=0.5, cell_mrna=flat)
    np.testing.assert_allclose(s.truth_mrna.to_numpy(), s.truth.to_numpy(), atol=1e-12)
