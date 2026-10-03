"""
identifiability_denominators.py -- Addendum 1, check 7, of prespecified/identifiability_diagnostics.md.

EpiDISH with the blood reference returns shares of the immune compartment (its seven columns sum
to 1); E2 compared them with tissue fractions. This recomputes every EpiDISH unit with matched
denominators -- each estimate's share of its own leukocyte total (Macrophage_Microglia + T + NK + B)
-- for D1 stability (the registered 7 settings), D2 agreement (the registered panel) and truth
agreement, then re-tests H1 and H1s exactly as registered. Tumor and Leukocytes are tissue fractions
on both sides already and keep their registered values.

    python3 scripts/identifiability_denominators.py
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
import identifiability_diagnostics as e2  # noqa: E402

OUT = config.RESULTS_DIR / "identifiability_denominators.json"
LEUK = e2.COMPARTMENTS["Leukocytes"]
MATCHED = {"Lymphoid": ("Lymphoid", "epi_lymphoid"), "T_cell": ("T_cell", "epi_T"),
           "B_cell": ("B_cell", "epi_B"), "NK_cell": ("NK_cell", "epi_NK")}
COHORTS = ("gbm", "lgg")


def share(est: pd.DataFrame, c: str) -> pd.Series:
    """The compartment's share of the estimate's own leukocyte total; undefined where that is 0."""
    leuk = est[LEUK].sum(axis=1)
    return (e2.comp(est, c) / leuk.where(leuk > 0)).astype(float)


def mean_pairwise_finite(frames: dict[str, pd.Series]) -> tuple[float, int, int, int]:
    """E2's mean pairwise Spearman, on the samples where every frame is defined."""
    df = pd.DataFrame(frames)
    ok = df.notna().all(axis=1)
    df = df[ok]
    usable = [k for k in df.columns if df[k].nunique() > 1]
    if len(usable) < 3:
        return float("nan"), len(usable), int((~ok).sum()), len(df)
    r = [stats.spearmanr(df[a], df[b]).statistic for a, b in itertools.combinations(usable, 2)]
    return float(np.mean(r)), len(usable), int((~ok).sum()), len(df)


def truth_rho(est_c: pd.Series, truth: pd.Series) -> tuple[float, int]:
    e = est_c.copy()
    e.index = [e2.k4s(i) for i in e.index]
    e = e[~e.index.duplicated()].dropna()
    shared = e.index.intersection(truth.index)
    if len(shared) < 10 or e[shared].nunique() < 2:
        return float("nan"), len(shared)
    return float(stats.spearmanr(e[shared], truth[shared]).statistic), len(shared)


def main() -> int:
    reg = json.loads(e2.OUT.read_text())
    units: dict[str, dict] = {}
    for h in COHORTS:
        base = "" if h == "gbm" else "_lgg"
        tr = e2.truths(h)
        d1 = {lab: pd.read_csv(e2.OUT_DIR / f"unmix_{h}_{lab}.csv", index_col=0) for lab, _, _ in e2.SETTINGS}
        full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
        sid = "sample_id" if "sample_id" in full.columns else "sample"
        meth = {m: g.set_index(sid)[list(config.CELL_TYPES)] for m, g in full.groupby("method") if m not in e2.DROP}
        for c, (comp, tname) in MATCHED.items():
            s1 = mean_pairwise_finite({lab: share(e, comp) for lab, e in d1.items()})
            s2 = mean_pairwise_finite({m: share(e, comp) for m, e in meth.items()})
            t1 = truth_rho(share(d1["predeclared"], comp), tr[tname])
            t2s = [truth_rho(share(e, comp), tr[tname])[0] for e in meth.values()]
            units[f"{c}|{h}"] = {
                "d1_stability": s1[0], "d1_settings_used": s1[1], "d1_samples_dropped_undefined": s1[2],
                "d1_samples": s1[3], "d2_agreement": s2[0], "d2_methods_used": s2[1],
                "d2_samples_dropped_undefined": s2[2], "d2_samples": s2[3],
                "truth_rho_unmix": t1[0], "n_truth": t1[1],
                "truth_rho_median_methods": float(np.nanmedian(t2s)) if np.isfinite(t2s).any() else float("nan"),
                "registered_truth_rho_unmix": reg["units"][f"{c}|{h}"]["truth_rho_unmix"],
                "registered_d1_stability": reg["units"][f"{c}|{h}"]["d1_stability"]}
        for c in ("Tumor", "Leukocytes"):                 # tissue fractions on both sides: unchanged
            r = reg["units"][f"{c}|{h}"]
            units[f"{c}|{h}"] = {k: r[k] for k in ("d1_stability", "d2_agreement", "truth_rho_unmix",
                                                    "truth_rho_median_methods", "n_truth")}

    def test(names, x, y):
        pts = [(units[f"{c}|{h}"][x], units[f"{c}|{h}"][y]) for c in names for h in COHORTS]
        pts = [p for p in pts if all(np.isfinite(p))]
        rho, p = e2.exact_one_sided([q[0] for q in pts], [q[1] for q in pts])
        return {"n_units": len(pts), "rho": round(rho, 4), "p_one_sided": round(p, 4)}

    prim = ["Tumor", "Leukocytes", "Lymphoid"]
    sec = prim + ["T_cell", "B_cell", "NK_cell"]
    out = {"rule": "prespecified/identifiability_diagnostics.md, Addendum 1, check 7",
           "units": units,
           "H1_matched": {"D1": test(prim, "d1_stability", "truth_rho_unmix"),
                          "D2": test(prim, "d2_agreement", "truth_rho_median_methods")},
           "H1s_matched": {"D1": test(sec, "d1_stability", "truth_rho_unmix"),
                           "D2": test(sec, "d2_agreement", "truth_rho_median_methods")}}
    lym = [units[f"Lymphoid|{h}"]["truth_rho_unmix"] for h in COHORTS]
    out["reading"] = ("WITHDRAWN: the lymphoid non-agreement was a denominator artefact" if max(lym) >= 0.30
                      else "STANDS on matched denominators: the lymphoid estimate does not track methylation")
    OUT.write_text(json.dumps(out, indent=2, default=float))
    for k, u in units.items():
        if "registered_truth_rho_unmix" in u:
            print(f"{k:16s} D1 {u['d1_stability']:.3f} (reg {u['registered_d1_stability']:.3f}; "
                  f"{u['d1_samples_dropped_undefined']} undefined)  D2 {u['d2_agreement']:.3f} "
                  f"({u['d2_samples_dropped_undefined']} undefined)  truth unmix {u['truth_rho_unmix']:+.3f} "
                  f"(reg {u['registered_truth_rho_unmix']:+.3f}, n {u['n_truth']})  median {u['truth_rho_median_methods']:+.3f}")
    print("H1 matched:", out["H1_matched"]); print("H1s matched:", out["H1s_matched"])
    print("reading:", out["reading"]); print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
