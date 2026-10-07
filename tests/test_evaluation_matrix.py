"""
The evaluation matrix (scripts/evaluation_matrix.py) must not let a method go quietly missing or a join drift:
every recomputed tumour correlation equals the registered yardstick's, every method of every registered arm
appears, and every empty cell carries a reason.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402

# conftest redirects RESULTS_DIR to a temporary tree so that no test can WRITE a real artefact. This test only
# READS them (build() writes nothing), so it points the modules at the real tree for its duration.
REAL = config.PROJECT_ROOT / "results"
NEEDED = [REAL / f"absolute_purity_yardstick{t}.json" for t in ("", "_h5ad", "_lgg", "_lgg_h5ad")] + \
         [REAL / "gimicc" / "runs" / "gbm_GBM_h4.csv", REAL / "anatomic" / "acs_leaderboard.csv"]


@pytest.mark.skipif(not all(p.exists() for p in NEEDED), reason="registered artefacts absent")
def test_no_method_missing_no_join_drift_every_gap_explained(monkeypatch):
    # No module is reloaded: a reload re-binds module-level paths for every later test (it broke
    # test_identifiability once). Only the paths this test reads are patched, and monkeypatch restores them.
    monkeypatch.setattr(config, "RESULTS_DIR", REAL)
    import gimicc_truth
    import evaluation_matrix
    monkeypatch.setattr(evaluation_matrix, "R", REAL)
    monkeypatch.setattr(gimicc_truth, "RUNS", REAL / "gimicc" / "runs")
    o = evaluation_matrix.build()
    em = evaluation_matrix
    assert o["n_join_checks"] >= 40 and o["n_join_disagreements"] == 0
    for cohort, ref, tag in em.ARMS:
        reg = json.loads((REAL / f"absolute_purity_yardstick{tag}.json").read_text())["methods"]
        assert set(reg) <= set(o["arms"][f"{cohort}|{ref}"])          # no registered method dropped
    assert all(g["why"] for g in o["gaps"])
