"""
FARDEEP's driver on its DEFAULT path (one worker), which crashed until 2026-10-02.

`R/run_fardeep.R` splits samples into blocks with `cut()`, and `cut()` rejects a single interval.
Every run before 2026-10-02 exported IVYGAP_FARDEEP_CORES=2, so the one-worker default was never
executed; the regeneration queue left it unset and both cohorts failed. The only earlier test
read the call's arguments with a regex and could not see this.

Three checks: the default path runs and recovers planted proportions; one and two workers give
identical estimates (the driver's comment claims blocking cannot change an estimate); and a
reference with shuffled gene labels fails the same bar (the negative control).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ivygap import config

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "R" / "run_fardeep.R"


def _has_fardeep() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "cat(requireNamespace('FARDEEP', quietly = TRUE))"],
                       capture_output=True, text=True, timeout=120)
    return r.stdout.strip().endswith("TRUE")


pytestmark = pytest.mark.skipif(not _has_fardeep(), reason="R package FARDEEP not available")


def _run(sig: pd.DataFrame, bulk: pd.DataFrame, tmp: Path, cores: str | None) -> pd.DataFrame:
    sig.to_csv(tmp / "sig.csv")
    bulk.to_csv(tmp / "bulk.csv")
    args = {"bulk": str(tmp / "bulk.csv"), "signature": str(tmp / "sig.csv"),
            "out": str(tmp / "out.csv"), "seed": 0, "cell_types": list(sig.columns)}
    (tmp / "args.json").write_text(json.dumps(args))
    env = {k: v for k, v in os.environ.items() if k != "IVYGAP_FARDEEP_CORES"}
    if cores is not None:
        env["IVYGAP_FARDEEP_CORES"] = cores
    r = subprocess.run(["Rscript", str(SCRIPT), str(tmp / "args.json")],
                       capture_output=True, text=True, timeout=1800, env=env)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    return pd.read_csv(tmp / "out.csv", index_col=0)


@pytest.fixture(scope="module")
def planted():
    sig = pd.read_csv(config.FROZEN_SIGNATURE_PATH, sep="\t", index_col=0)
    types = [t for t in config.PRIMARY_CELL_TYPES if t in sig.columns]
    sig = sig[types]
    sig = sig.loc[sig.var(axis=1).sort_values(ascending=False).index[:400]]   # informative genes
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(len(types)), size=20)                           # samples x types
    noise = rng.lognormal(0.0, 0.1, size=(len(sig), 20))
    bulk = pd.DataFrame((sig.to_numpy() @ p.T) * noise, index=sig.index,
                        columns=[f"s{i}" for i in range(20)])
    return sig, bulk, pd.DataFrame(p, index=bulk.columns, columns=types)


def _recovery(est: pd.DataFrame, truth: pd.DataFrame) -> tuple[float, float]:
    est = est[truth.columns]
    # a constant estimate has no correlation; count it as no recovery rather than let NaN
    # slip through min() and both comparisons
    r = min(float(np.nan_to_num(np.corrcoef(est[c], truth[c])[0, 1])) for c in truth.columns)
    return float(r), float(np.abs(est.to_numpy() - truth.to_numpy()).max())


def test_default_one_worker_path_runs_and_recovers(planted, tmp_path):
    sig, bulk, truth = planted
    r, err = _recovery(_run(sig, bulk, tmp_path, None), truth)
    assert r >= 0.90 and err <= 0.10, (r, err)


def test_one_and_two_workers_give_identical_estimates(planted, tmp_path):
    sig, bulk, _ = planted
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    one = _run(sig, bulk, tmp_path / "a", None)
    two = _run(sig, bulk, tmp_path / "b", "2")
    assert list(one.index) == list(two.index)
    assert float(np.abs(one.to_numpy() - two.loc[one.index, one.columns].to_numpy()).max()) == 0.0


def test_shuffled_reference_fails_the_same_bar(planted, tmp_path):
    """The negative control: a reference with its gene labels shuffled must not pass."""
    sig, bulk, truth = planted
    shuffled = sig.copy()
    shuffled.index = np.random.default_rng(1).permutation(sig.index)
    r, err = _recovery(_run(shuffled.loc[sig.index], bulk, tmp_path, None), truth)
    assert r < 0.90 or err > 0.10, (r, err)
