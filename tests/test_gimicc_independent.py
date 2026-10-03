"""
Independent recomputation of the GIMiCC headline numbers (results/gimicc_truth_confirmation.json),
sharing no code with scripts/gimicc_truth.py: plain pandas on GIMiCC's raw output CSVs, and ABSOLUTE
matched through its own `array` column (TCGA-xx-xxxx-01) rather than through any project key function
-- the key path that failed on 2026-10-03 is not reused here.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RUNS = RES / "gimicc" / "runs"
ABS = ROOT / "data" / "raw" / "tcga" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
DEFAULT = {"gbm": "GBM", "lgg": "AST"}


@pytest.fixture(scope="module")
def art():
    p = RES / "gimicc_truth_confirmation.json"
    if not p.exists() or not RUNS.exists():
        pytest.skip("GIMiCC artefacts not produced yet")
    return json.loads(p.read_text())


def _raw(cohort: str, tag: str) -> pd.DataFrame:
    return pd.read_csv(RUNS / f"{cohort}_{tag}.csv", index_col=0) / 100.0


@pytest.mark.parametrize("cohort", ["gbm", "lgg"])
def test_within_lymphoid_means_recompute(art, cohort):
    df = _raw(cohort, f"{DEFAULT[cohort]}_h4")
    lym = pd.DataFrame({"T": df["CD4Tcell"] + df["CD8Tcell"], "NK": df["NK"], "B": df["Bcell"]})
    lym = lym[lym.notna().all(axis=1)]
    lym = lym[lym.sum(axis=1) > 0]
    mu = lym.div(lym.sum(axis=1), axis=0).mean()
    rec = art["cohorts"][cohort]["Q1_all_methylation_samples"]
    assert rec["n_samples_with_lymphoid"] == len(lym)
    for k in ("T", "NK", "B"):
        assert abs(rec["mean_within_lymphoid_share"][k] - mu[k]) < 1e-4, (k, mu[k])
    assert rec["T_exceeds_B"] == bool(mu["T"] > mu["B"])


@pytest.mark.parametrize("cohort", ["gbm", "lgg"])
def test_tumour_controls_recompute_through_absolutes_own_array_column(art, cohort):
    if not ABS.exists():
        pytest.skip("ABSOLUTE table not on disk")
    a = pd.read_csv(ABS, sep="\t").dropna(subset=["Cancer DNA fraction"])
    pur = a.drop_duplicates("array").set_index("array")["Cancer DNA fraction"].astype(float)
    for tag, key in ((f"{DEFAULT[cohort]}_h4", "control_tumor_vs_absolute"),
                     (f"{DEFAULT[cohort]}_h4_shuffled", "negative_control_shuffled_cpgs")):
        t = _raw(cohort, tag)["Tumor"]
        t = t[~t.index.duplicated()]
        idx = t.index.intersection(pur.index)
        rho = float(stats.spearmanr(t[idx], pur[idx]).statistic)
        rec = art["cohorts"][cohort][key]
        recorded = rec["spearman"] if "spearman" in rec else rec["tumor_vs_absolute_spearman"]
        assert rec["n"] == len(idx) and len(idx) >= 100
        assert abs(recorded - rho) < 1e-4, (tag, rho, recorded)
