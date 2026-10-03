"""
Planted-truth gate for MIXTURE (prespecified/mixture_config.md), before any real data.

Mixtures are built from the frozen reference profiles with known proportions; the genuine
package must recover them. The negative control -- the same mixtures against a reference whose
gene labels are shuffled -- must FAIL the same bar, or the gate is measuring nothing. The
Astrocyte sidecar is excluded: it is unidentifiable against Tumor by construction (CLAUDE.md),
which would test the signature, not the method.
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ivygap import config

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "R" / "run_mixture.R"
def _has_mixture() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "cat(requireNamespace('MIXTURE', quietly = TRUE))"],
                       capture_output=True, text=True, timeout=120)
    return r.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _has_mixture(), reason="R package MIXTURE not installed (vendor/MIXTURE)")


def _run(sig: pd.DataFrame, bulk: pd.DataFrame, tmp: Path) -> pd.DataFrame:
    sig.to_csv(tmp / "sig.csv")
    bulk.to_csv(tmp / "bulk.csv")
    args = {"bulk": str(tmp / "bulk.csv"), "signature": str(tmp / "sig.csv"),
            "out": str(tmp / "out.csv"), "seed": 0, "cell_types": list(sig.columns)}
    (tmp / "args.json").write_text(json.dumps(args))
    r = subprocess.run(["Rscript", str(SCRIPT), str(tmp / "args.json")],
                       capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    return pd.read_csv(tmp / "out.csv", index_col=0)


def _bar(est: pd.DataFrame, truth: pd.DataFrame) -> tuple[float, float]:
    """(min per-type r, max abs error). A missing estimate recovers nothing: NaN correlations
    count as 0 and NaN errors as infinite. Comparisons with NaN are all False, so without this a
    run that returned NO estimate would satisfy `min(r) < 0.90 or err > 0.10` as False and look
    like it passed the bar -- which is how this test first failed on MIXTURE (2026-10-02)."""
    est = est[truth.columns]
    r = [np.corrcoef(est[c], truth[c])[0, 1] for c in truth.columns]
    r = [x if np.isfinite(x) else 0.0 for x in r]
    diff = np.abs(est.to_numpy() - truth.to_numpy())
    return float(min(r)), (float(diff.max()) if np.isfinite(diff).all() else float("inf"))


@pytest.fixture(scope="module")
def planted():
    sig = pd.read_csv(config.FROZEN_SIGNATURE_PATH, sep="\t", index_col=0)
    types = [t for t in config.PRIMARY_CELL_TYPES if t in sig.columns]
    sig = sig[types]
    sig = sig.loc[sig.var(axis=1).sort_values(ascending=False).index[:800]]   # informative genes
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(len(types)), size=30)                           # samples x types
    noise = rng.lognormal(0.0, 0.1, size=(len(sig), 30))
    bulk = pd.DataFrame((sig.to_numpy() @ p.T) * noise, index=sig.index,
                        columns=[f"s{i}" for i in range(30)])
    return sig, bulk, pd.DataFrame(p, index=bulk.columns, columns=types)


def test_planted_proportions_are_recovered(planted, tmp_path):
    sig, bulk, truth = planted
    r, err = _bar(_run(sig, bulk, tmp_path), truth)
    assert r >= 0.90 and err <= 0.10, (r, err)


def test_shuffled_reference_fails_the_same_bar(planted, tmp_path):
    """The negative control: a reference with its gene labels shuffled must not pass."""
    sig, bulk, truth = planted
    rng = np.random.default_rng(1)
    shuffled = sig.copy()
    shuffled.index = rng.permutation(sig.index)
    r, err = _bar(_run(shuffled.loc[sig.index], bulk, tmp_path), truth)
    assert r < 0.90 or err > 0.10, (r, err)
