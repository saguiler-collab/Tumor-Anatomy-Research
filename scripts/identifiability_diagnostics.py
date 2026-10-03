"""
identifiability_diagnostics.py -- Extension E2, implementing prespecified/identifiability_diagnostics.md.

Two truth-free diagnostics per compartment and cohort, then the test of whether they predict
agreement with DNA truths:

  D1  loss-scale stability: genuine DESeq2 unmix at seven loss settings on the TCGA frozen arm;
      mean pairwise Spearman of per-sample estimates across settings
  D2  method agreement: the registered panel's saved per-sample estimates; mean pairwise Spearman
      across methods (degenerate duplicates dropped)
  truth  ABSOLUTE purity (Tumor), Thorsson leukocyte fraction (Leukocytes), EpiDISH (Lymphoid;
         T, B, NK secondary)
  S1  data scale vs loss scale: NNLS on log1p data vs linear NNLS (this project's code)

    python3 scripts/identifiability_diagnostics.py [--recompute]
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
import extension_tcga  # noqa: E402
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402
from lymphoid_ordering import k4s  # noqa: E402

OUT_DIR = config.RESULTS_DIR / "identifiability"
OUT = config.RESULTS_DIR / "identifiability_diagnostics.json"
LF_PATH = config.IMMUNE_FRACTION_DIR / "TCGA_all_leuk_estimate.masked.20170107.tsv"
SETTINGS = [("predeclared", None, None), ("s1", 1, None), ("s10", 10, None), ("s100", 100, None),
            ("s1000", 1000, None), ("s10000", 10000, None), ("predeclared_p2", None, 2)]
COMPARTMENTS = {"Tumor": ["Tumor"],
                "Leukocytes": ["Macrophage_Microglia", "T_cell", "NK_cell", "B_cell"],
                "Myeloid": ["Macrophage_Microglia"], "Lymphoid": ["T_cell", "NK_cell", "B_cell"],
                "T_cell": ["T_cell"], "B_cell": ["B_cell"], "NK_cell": ["NK_cell"],
                "Endothelial": ["Endothelial"], "Oligodendrocyte": ["Oligodendrocyte"]}
PRIMARY = ["Tumor", "Leukocytes", "Lymphoid"]
SECONDARY = ["T_cell", "B_cell", "NK_cell"]
DROP = {"music", "scdc_ensemble"}          # degenerate duplicates on the frozen signature


def comp(est: pd.DataFrame, c: str) -> pd.Series:
    return est[COMPARTMENTS[c]].sum(axis=1)


def truths(cohort: str) -> dict[str, pd.Series]:
    base = "" if cohort == "gbm" else "_lgg"
    a = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    a["k"] = a["sample"].map(key4)
    purity = a.drop_duplicates("k").set_index("k")[PURITY_COL].astype(float)
    lf = pd.read_csv(LF_PATH, sep="\t", header=None, names=["study", "barcode", "lf"])
    lf["k"] = lf["barcode"].map(k4s)
    lf["lf"] = pd.to_numeric(lf["lf"], errors="coerce")
    lf = lf.dropna(subset=["lf"]).drop_duplicates("k").set_index("k")["lf"]
    m = pd.read_csv(config.RESULTS_DIR / f"methylation_celltypes{base}.csv", index_col=0)
    m.index = [k4s(i) for i in m.index]
    m = m[~m.index.duplicated()]
    t = m["CD4T"] + m["CD8T"]
    return {"purity": purity, "lf": lf, "epi_lymphoid": t + m["NK"] + m["B"],
            "epi_T": t, "epi_B": m["B"], "epi_NK": m["NK"]}


TRUTH_OF = {"Tumor": ("purity", key4), "Leukocytes": ("lf", k4s), "Lymphoid": ("epi_lymphoid", k4s),
            "T_cell": ("epi_T", k4s), "B_cell": ("epi_B", k4s), "NK_cell": ("epi_NK", k4s)}


def truth_rho(est_c: pd.Series, c: str, tr: dict) -> tuple[float, int]:
    if c not in TRUTH_OF:
        return float("nan"), 0
    name, key = TRUTH_OF[c]
    e = est_c.copy()
    e.index = [key(i) for i in e.index]
    e = e[~e.index.duplicated()]
    shared = e.index.intersection(tr[name].index)
    if len(shared) < 10 or e[shared].nunique() < 2:
        return float("nan"), len(shared)
    return float(stats.spearmanr(e[shared], tr[name][shared]).statistic), len(shared)


def mean_pairwise(frames: dict[str, pd.Series]) -> tuple[float, int, int]:
    usable = {k: v for k, v in frames.items() if v.nunique() > 1}
    if len(usable) < 3:
        return float("nan"), len(usable), len(frames) - len(usable)
    keys = list(usable)
    idx = usable[keys[0]].index
    for k in keys[1:]:
        idx = idx.intersection(usable[k].index)
    r = [stats.spearmanr(usable[a][idx], usable[b][idx]).statistic for a, b in itertools.combinations(keys, 2)]
    return float(np.nanmean(r)), len(usable), len(frames) - len(usable)


def d1_estimates(cohort: str, recompute: bool) -> tuple[dict, dict]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ests, inputs = {}, None
    for label, shift, power in SETTINGS:
        path = OUT_DIR / f"unmix_{cohort}_{label}.csv"
        if path.exists() and not recompute and label != "predeclared":
            ests[label] = pd.read_csv(path, index_col=0)
            continue
        for k in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):
            os.environ.pop(k, None)
        if shift is not None:
            os.environ["IVYGAP_UNMIX_SHIFT"] = str(shift)
        if power is not None:
            os.environ["IVYGAP_UNMIX_POWER"] = str(power)
        res = extension_tcga.run(cohort, "frozen", ["deseq2_unmix"], return_estimates=True)
        e = res["_estimates"]["deseq2_unmix"]
        e.to_csv(path)
        ests[label] = e
        if label == "predeclared":
            inputs = res["_inputs"]
    for k in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):
        os.environ.pop(k, None)
    return ests, inputs


def s1_scale_arm(inputs: dict) -> dict:
    """NNLS on log1p data (Avila Cobos's case) against linear NNLS, same inputs, no conversion."""
    S, B, pur = inputs["signature"], inputs["bulk"], inputs["purity"]
    S = S.loc[B.index]
    ti = list(S.columns).index("Tumor")
    out = {}
    for name, f in (("linear_nnls", lambda x: x), ("log1p_data_nnls", np.log1p)):
        X = f(S.to_numpy(dtype=float))
        t = []
        for j in range(B.shape[1]):
            coef = nnls(X, f(B.iloc[:, j].to_numpy(dtype=float)))[0]
            t.append(coef[ti] / coef.sum() if coef.sum() > 0 else np.nan)
        t = np.array(t)
        ok = np.isfinite(t)
        out[name] = round(float(stats.spearmanr(t[ok], pur.to_numpy()[ok]).statistic), 4)
    out["prediction_log_worse"] = out["log1p_data_nnls"] < out["linear_nnls"]
    return out


def exact_one_sided(x: list, y: list) -> tuple[float, float]:
    rho = float(stats.spearmanr(x, y).statistic)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    perms = list(itertools.permutations(range(len(x)))) if len(x) <= 8 else None
    if perms is not None:
        null = [stats.spearmanr(rx, ry[list(p)]).statistic for p in perms]
    else:
        rng = np.random.default_rng(config.RANDOM_SEED)
        null = [stats.spearmanr(rx, rng.permutation(ry)).statistic for _ in range(100_000)]
    return rho, float(np.mean(np.array(null) >= rho - 1e-12))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--recompute", action="store_true")
    a = ap.parse_args()
    units, controls, s1 = {}, {}, {}
    for cohort in ("gbm", "lgg"):
        base = "" if cohort == "gbm" else "_lgg"
        tr = truths(cohort)
        ests, inputs = d1_estimates(cohort, a.recompute)
        # control 1: unmix Tumor rho at the pre-declared setting = the reported row
        rep = json.loads((config.RESULTS_DIR / "extension" / f"tcga_{cohort}_frozen.json").read_text())[
            "methods"]["deseq2_unmix"]["spearman_vs_purity"]
        mine = round(truth_rho(comp(ests["predeclared"], "Tumor"), "Tumor", tr)[0], 4)
        # control 2: registered NNLS purity rho from the saved per-sample estimates
        full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
        sid = "sample_id" if "sample_id" in full.columns else "sample"
        meth = {m: g.set_index(sid)[list(config.CELL_TYPES)] for m, g in full.groupby("method")
                if m not in DROP}
        reg = json.loads((config.RESULTS_DIR / f"absolute_purity_yardstick{base}.json").read_text())[
            "methods"]["nnls"]["spearman_vs_purity"]
        mine_nnls = round(truth_rho(comp(meth["nnls"], "Tumor"), "Tumor", tr)[0], 4)
        controls[cohort] = {"unmix_tumor_rho": [mine, rep], "nnls_tumor_rho": [mine_nnls, reg],
                            "pass": abs(mine - rep) < 1e-4 and abs(mine_nnls - reg) < 1e-4}
        print(f"{cohort}: control unmix {mine} vs {rep}; nnls {mine_nnls} vs {reg}", flush=True)
        if not controls[cohort]["pass"]:
            print("ABORT -- a control did not reproduce"); OUT.write_text(json.dumps({"controls": controls}))
            return 2
        for c in COMPARTMENTS:
            d1 = mean_pairwise({lab: comp(e, c) for lab, e in ests.items()})
            d2 = mean_pairwise({m: comp(e, c) for m, e in meth.items()})
            t1 = truth_rho(comp(ests["predeclared"], c), c, tr)
            t2s = [truth_rho(comp(e, c), c, tr)[0] for e in meth.values()]
            t2 = float(np.nanmedian(t2s)) if np.isfinite(t2s).any() else float("nan")
            zero = float(np.mean([(comp(e, c) == 0).mean() for e in ests.values()]))
            units[f"{c}|{cohort}"] = {"compartment": c, "cohort": cohort,
                                      "d1_stability": d1[0], "d1_settings_used": d1[1], "d1_constant": d1[2],
                                      "d2_agreement": d2[0], "d2_methods_used": d2[1], "d2_constant": d2[2],
                                      "truth_rho_unmix": t1[0], "n_truth": t1[1], "truth_rho_median_methods": t2,
                                      "unmix_zero_share": round(zero, 3),
                                      "truth_tier": ("primary" if c in PRIMARY else "secondary" if c in SECONDARY
                                                     else "none")}
        if inputs is not None:
            s1[cohort] = s1_scale_arm(inputs)
    def test(names, x_key, y_key):
        pts = [(units[f"{c}|{h}"][x_key], units[f"{c}|{h}"][y_key]) for c in names for h in ("gbm", "lgg")]
        pts = [p for p in pts if all(np.isfinite(p))]
        if len(pts) < 4:
            return {"n_units": len(pts), "note": "too few estimable units"}
        rho, p = exact_one_sided([q[0] for q in pts], [q[1] for q in pts])
        return {"n_units": len(pts), "rho": round(rho, 4), "p_one_sided": round(p, 4)}
    h1 = {"D1": test(PRIMARY, "d1_stability", "truth_rho_unmix"),
          "D2": test(PRIMARY, "d2_agreement", "truth_rho_median_methods")}
    def read(t):
        return (t.get("rho", 0) > 0 and t.get("p_one_sided", 1) < 0.05)
    h1["reading"] = ("SUPPORTED" if read(h1["D1"]) or read(h1["D2"]) else
                     "NOT SUPPORTED" if all(h1[k].get("rho", 0) <= 0 for k in ("D1", "D2")) else "INCONCLUSIVE")
    h1s = {"D1": test(PRIMARY + SECONDARY, "d1_stability", "truth_rho_unmix"),
           "D2": test(PRIMARY + SECONDARY, "d2_agreement", "truth_rho_median_methods")}
    h2 = {h: {"d1_stability": units[f"Leukocytes|{h}"]["d1_stability"],
              "truth_rho_unmix": units[f"Leukocytes|{h}"]["truth_rho_unmix"],
              "pass": bool(units[f"Leukocytes|{h}"]["d1_stability"] >= 0.8
                           and units[f"Leukocytes|{h}"]["truth_rho_unmix"] >= 0.4)} for h in ("gbm", "lgg")}
    out = {"rule": "prespecified/identifiability_diagnostics.md", "controls": controls, "units": units,
           "H1_primary": h1, "H1s_secondary": h1s, "H2_leukocytes": h2, "S1_scale_arm": s1}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print("\nunit                    D1 stab  D2 agree | truth unmix  truth med-methods | zero share")
    for k, u in units.items():
        print(f"{k:22s}  {u['d1_stability']:+.3f}   {u['d2_agreement']:+.3f}  |  {u['truth_rho_unmix']:+.3f}"
              f"        {u['truth_rho_median_methods']:+.3f}          | {u['unmix_zero_share']:.2f}  [{u['truth_tier']}]")
    print("H1 (primary):", h1); print("H1s (secondary):", h1s); print("H2:", h2); print("S1:", s1)
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
