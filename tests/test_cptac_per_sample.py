"""
The CPTAC single-nucleus naming (scripts/cptac_per_sample.py; prespecified/cptac_per_sample_truth.md) on
planted nuclei: clusters with clear markers are named correctly, an ambiguous cluster is left
UNMAPPED rather than forced into a type, and a mixed T/NK cluster is split per nucleus.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import cptac_per_sample as cp  # noqa: E402
from build_abdelfattah_reference import PANELS  # noqa: E402


def _nuclei(groups: list[tuple[str, list[str], int]], seed: int = 0):
    """Each group: (cluster id, marker genes expressed strongly, number of nuclei)."""
    rng = np.random.default_rng(seed)
    genes = sorted({g for p in PANELS.values() for g in p} | {f"BG{i}" for i in range(50)})
    gi = {g: i for i, g in enumerate(genes)}
    cols, bcs, cl = [], [], []
    for cid, markers, n in groups:
        for k in range(n):
            v = rng.poisson(1.0, len(genes)).astype(float)          # background
            for g in markers:
                v[gi[g]] += rng.poisson(40)
            cols.append(v); bcs.append(f"{cid}_{k}"); cl.append(cid)
    mtx = sp.csc_matrix(np.array(cols).T)
    return mtx, genes, bcs, pd.Series(cl, index=bcs)


def test_clear_clusters_are_named_and_ambiguous_ones_are_not():
    mtx, genes, bcs, clusters = _nuclei([
        ("c0", PANELS["Tumor"], 60), ("c1", PANELS["Macrophage_Microglia"], 40),
        ("c2", PANELS["Oligodendrocyte"], 30),
        ("c3", PANELS["Tumor"] + PANELS["Macrophage_Microglia"], 20)])     # doublet-like: ambiguous
    lab, per = cp.label_sample(mtx, genes, bcs, clusters)
    assert per["c0"]["label"] == "Tumor" and per["c1"]["label"] == "Macrophage_Microglia"
    assert per["c2"]["label"] == "Oligodendrocyte"
    assert per["c3"]["label"] is None and lab[clusters == "c3"].isna().all()


def test_mixed_t_nk_cluster_is_split_per_nucleus():
    mtx, genes, bcs, clusters = _nuclei([("tn", [], 0)])
    # one cluster holding 30 T-like and 30 NK-like nuclei
    mtx_t, genes, bcs_t, _ = _nuclei([("tn", PANELS["T_cell"], 30)], seed=1)
    mtx_n, _, bcs_n, _ = _nuclei([("tn", PANELS["NK_cell"], 30)], seed=2)
    mtx = sp.hstack([mtx_t, mtx_n]).tocsc()
    bcs = [f"t{i}" for i in range(30)] + [f"n{i}" for i in range(30)]
    clusters = pd.Series("tn", index=bcs)
    lab, per = cp.label_sample(mtx, genes, bcs, clusters)
    assert per["tn"]["label"] == "T/NK split"
    assert (lab.iloc[:30] == "T_cell").mean() > 0.9 and (lab.iloc[30:] == "NK_cell").mean() > 0.9
