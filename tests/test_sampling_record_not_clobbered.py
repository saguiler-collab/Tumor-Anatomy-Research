"""
A variant reference build must not overwrite the primary reference's sampling record.

docs/OPEN_DEFECTS.md D17. `build_from_h5ad` wrote to `reference_sampling_{name}.json` and
`name` defaults to "gbmap", so every sensitivity build — the assay subsets, the linear
`raw/X` build — landed on the primary reference's record and replaced it. That is not a
tidiness problem: `scripts/figure_data.py` reads exactly that path, so the last sensitivity
build's numbers would have been published as the leaderboard's, and `/data/` is gitignored so
there was no committed copy to fall back on.

The negative control is the one that matters here. A guard that refuses every second write
would also pass a "did not clobber" test while breaking honest re-runs, so
`test_an_identical_rerun_still_overwrites` is the test that shows the guard discriminates
rather than just blocks.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ivygap import config


def _tiny_atlas(tmp_path, n_cells=400, n_genes=60, n_donors=6):
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(0)
    types = ["AC-like", "MES-like", "TAM-BDM", "CD4/CD8", "Oligodendrocyte", "Endothelial"]
    obs = pd.DataFrame(
        {"annotation_level_3": rng.choice(types, n_cells),
         "donor_id": rng.choice([f"D{i}" for i in range(n_donors)], n_cells)},
        index=[f"c{i}" for i in range(n_cells)])
    counts = sparse.random(n_cells, n_genes, density=0.3, random_state=0, format="csr")
    counts.data = np.ceil(counts.data * 50)
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])
    a = anndata.AnnData(X=counts.copy(), obs=obs, var=var)
    # a second matrix, so a variant build has something different to read
    a.raw = anndata.AnnData(X=counts.multiply(3).tocsr(), obs=obs, var=var)
    path = tmp_path / "atlas.h5ad"
    a.write_h5ad(path)
    return path


def _isolate(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "REFERENCE_DIR", tmp_path)
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS", [tmp_path])


def _build(path, **kw):
    from ivygap.data.reference import build_from_h5ad
    return build_from_h5ad(path, name="gbmap", max_cells_per_donor_type=5,
                           max_total_cells=10_000, export=False, **kw)


def test_a_variant_build_does_not_touch_the_canonical_record(tmp_path, monkeypatch):
    _isolate(tmp_path, monkeypatch)
    path = _tiny_atlas(tmp_path)
    canonical = tmp_path / "reference_sampling_gbmap.json"

    _build(path)                                          # the primary build
    first = canonical.read_text()
    assert json.loads(first)["matrix"] == "X"

    _build(path, matrix="raw/X")                          # the variant

    assert canonical.read_text() == first, \
        "the variant build overwrote the primary reference's sampling record"

    variants = [p for p in tmp_path.glob("reference_sampling_gbmap__*.json")]
    assert len(variants) == 1, f"the variant wrote no file of its own: {variants}"
    rec = json.loads(variants[0].read_text())
    assert rec["matrix"] == "raw/X"
    assert "not_the_canonical_record" in rec, \
        "a variant record must say it is not the primary reference's provenance"


def test_an_identical_rerun_still_overwrites(tmp_path, monkeypatch):
    """NEGATIVE CONTROL: the guard must key on the BUILD, not on the file existing."""
    _isolate(tmp_path, monkeypatch)
    path = _tiny_atlas(tmp_path)

    _build(path)
    _build(path)

    assert not list(tmp_path.glob("reference_sampling_gbmap__*.json")), \
        "an identical re-run was misread as a variant; the guard blocks instead of keying"


def test_the_variant_key_ignores_measured_outputs(tmp_path, monkeypatch):
    """
    The fingerprint must cover inputs, not results. Keying on a measured count would make
    every honest re-run look like a new variant and defeat the guard by noise.
    """
    from ivygap.data.reference import _variant_key
    base = {"matrix": "X", "seed": 0, "sha256": "abc", "n_genes_kept": 80,
            "cell_filter_applied": False, "n_cells_kept": 100}
    same = dict(base, n_cells_kept=999, n_cells_after_roster_mapping=7)
    assert _variant_key(base) == _variant_key(same)
    assert _variant_key(base) != _variant_key(dict(base, matrix="raw/X"))
    assert _variant_key(base) != _variant_key(dict(base, seed=1))
