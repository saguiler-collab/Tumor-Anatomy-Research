"""
Planted-truth gate for Linseed (prespecified/linseed_config.md), before any real data.

Linseed is reference-free: it must find the cell-type profiles as well as the proportions, and
it names nothing. Mixtures of four roster profiles from the frozen signature are deconvolved with
k = 4; each true type is matched to an estimated component by maximum total correlation
(Hungarian), and the weakest matched correlation must reach 0.80. The negative control permutes
each gene's values independently across samples -- the marginal distributions survive, the
mixing structure does not -- and must FAIL that bar. A run that errors, or returns no estimate,
counts as recovering nothing.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import linear_sum_assignment

from ivygap import config

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "R" / "run_linseed.R"
TYPES = ["Tumor", "Macrophage_Microglia", "Oligodendrocyte", "Endothelial"]
BAR = 0.80


def _has_linseed() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "cat(requireNamespace('linseed', quietly = TRUE))"],
                       capture_output=True, text=True, timeout=120)
    return r.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _has_linseed(), reason="R package linseed not installed (vendor/linseed)")


def _run(bulk: pd.DataFrame, tmp: Path, k: int = 4) -> tuple[pd.DataFrame | None, str]:
    bulk.to_csv(tmp / "bulk.csv")
    args = {"bulk": str(tmp / "bulk.csv"), "k": k, "top_genes": 10000, "sig_iters": 100,
            "pval": 0.01, "seed": 0, "out_prop": str(tmp / "prop.csv"), "out_sig": str(tmp / "sig.csv")}
    (tmp / "args.json").write_text(json.dumps(args))
    r = subprocess.run(["Rscript", str(SCRIPT), str(tmp / "args.json")],
                       capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not (tmp / "prop.csv").exists():
        return None, (r.stdout[-800:] + r.stderr[-800:])
    return pd.read_csv(tmp / "prop.csv", index_col=0), r.stdout[-300:]


def _matched_min_r(est: pd.DataFrame | None, truth: pd.DataFrame) -> float:
    """Weakest correlation after matching estimated components to true types (Hungarian)."""
    if est is None or est.shape[1] != truth.shape[1]:
        return 0.0
    k = truth.shape[1]
    C = np.corrcoef(est.to_numpy(dtype=float).T, truth.to_numpy(dtype=float).T)[:k, k:]
    C = np.nan_to_num(C)                                  # an undefined correlation recovers nothing
    rows, cols = linear_sum_assignment(-C)
    return float(C[rows, cols].min())


@pytest.fixture(scope="module")
def planted():
    sig = pd.read_csv(config.FROZEN_SIGNATURE_PATH, sep="\t", index_col=0)[TYPES]
    sig = sig.loc[sig.var(axis=1).sort_values(ascending=False).index[:3000]]
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(len(TYPES)), size=40)
    noise = rng.lognormal(0.0, 0.05, size=(len(sig), 40))
    bulk = pd.DataFrame((sig.to_numpy() @ p.T) * noise, index=sig.index,
                        columns=[f"s{i}" for i in range(40)])
    return bulk, pd.DataFrame(p, index=bulk.columns, columns=TYPES)


def test_planted_mixtures_are_recovered(planted, tmp_path):
    bulk, truth = planted
    est, log = _run(bulk, tmp_path)
    assert est is not None, log
    assert [str(s) for s in est.index] == list(bulk.columns)
    r = _matched_min_r(est, truth)
    assert r >= BAR, (r, log)


def test_destroying_the_mixing_structure_fails_the_bar(planted, tmp_path):
    """The negative control: per-gene independent permutation across samples."""
    bulk, truth = planted
    rng = np.random.default_rng(1)
    shuffled = pd.DataFrame(np.vstack([rng.permutation(row) for row in bulk.to_numpy()]),
                            index=bulk.index, columns=bulk.columns)
    est, log = _run(shuffled, tmp_path)
    r = _matched_min_r(est, truth)
    assert r < BAR, (r, log)
