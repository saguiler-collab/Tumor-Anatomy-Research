"""
The cell-level export cache must key on the VALUES, not only on gene and cell identity.

The key used to be sha256(gene list || cell ids). Re-registering the same cells and genes
at a different scaling -- raw library sizes instead of 1e6 CPM -- produced an identical key
and handed R the stale file. A cache that ignores the content it caches is a silent input
substitution, the same shape as D10, where a re-measurement ran on a gene space nobody had
recorded.

That scaling is not a contrived variation: it decides whether MuSiC, SCDC and Bisque can
estimate cell size at all, which is the whole subject of OPEN_DEFECTS D1.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap.deconv import r_bridge


@pytest.fixture
def cells():
    rng = np.random.default_rng(0)
    genes = [f"G{i:03d}" for i in range(12)]
    ids = [f"C{i:03d}" for i in range(20)]
    expr = pd.DataFrame(rng.lognormal(0, 1, (len(genes), len(ids))),
                        index=genes, columns=ids)
    meta = pd.DataFrame({"donor": ["D0"] * 10 + ["D1"] * 10,
                         "cell_type": ["Tumor"] * 20}, index=ids)
    return genes, expr, meta


def _export(monkeypatch, tmp_path, genes, expr, meta):
    monkeypatch.setattr(r_bridge.config, "REFERENCE_DIR", tmp_path)
    r_bridge.set_cell_source("probe", expr, meta)
    counts, _ = r_bridge.export_for_genes("probe", genes)
    return counts


def test_rescaling_the_same_cells_gets_a_different_export(monkeypatch, tmp_path, cells):
    genes, expr, meta = cells

    raw_path = _export(monkeypatch, tmp_path, genes, expr, meta)

    # Same genes, same cell ids, different values: normalise each cell to a common total.
    norm = expr.div(expr.sum(axis=0), axis=1) * 1e6
    norm_path = _export(monkeypatch, tmp_path, genes, norm, meta)

    assert raw_path != norm_path, (
        "the raw and CPM-normalised exports share a cache path, so whichever was written "
        "first is handed to R for both. This is the bug that made the cell-size probe "
        "unable to vary the one thing it needed to vary.")
    assert raw_path.exists() and norm_path.exists()

    a = pd.read_csv(raw_path, index_col=0, compression="gzip")
    b = pd.read_csv(norm_path, index_col=0, compression="gzip")
    assert not np.allclose(a.to_numpy(), b.to_numpy()), "the two exports hold equal values"
    # The normalised one must actually be normalised -- proves the right file came back.
    assert np.allclose(b.sum(axis=0).to_numpy(), 1e6, rtol=1e-6)


def test_identical_input_still_reuses_one_file(monkeypatch, tmp_path, cells):
    """The cache must still cache; a key that never hits would rewrite on every method."""
    genes, expr, meta = cells
    first = _export(monkeypatch, tmp_path, genes, expr, meta)
    mtime = first.stat().st_mtime_ns
    second = _export(monkeypatch, tmp_path, genes, expr, meta)
    assert first == second
    assert second.stat().st_mtime_ns == mtime, "an unchanged input was rewritten"


def test_one_changed_value_changes_the_key(monkeypatch, tmp_path, cells):
    genes, expr, meta = cells
    a = _export(monkeypatch, tmp_path, genes, expr, meta)
    bumped = expr.copy()
    bumped.iloc[3, 7] = float(bumped.iloc[3, 7]) + 1.0
    b = _export(monkeypatch, tmp_path, genes, bumped, meta)
    assert a != b, "a single altered count did not change the cache key"


def test_different_genes_still_separate(monkeypatch, tmp_path, cells):
    genes, expr, meta = cells
    a = _export(monkeypatch, tmp_path, genes, expr, meta)
    b = _export(monkeypatch, tmp_path, genes[:6], expr, meta)
    assert a != b
