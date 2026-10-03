"""
agreement_selection_test.py -- prespecified/agreement_selection_test.md.

Does agreement between deconvolution methods pick the accurate estimate? DECEPTICON's selection
criterion (Deng et al., Brief Bioinform 2025;26:bbaf234) applied to this study's registered panel
and scored against DNA truths in TCGA-GBM and TCGA-LGG. THIS PROJECT'S CODE implementing the
published selection rule -- not the DECEPTICON package, which runs its own ten strategies.

  P1  the agreement-selected ensemble vs ensembles from randomly chosen pairs (2000 draws)
  P2  across methods, is centrality (mean correlation with the other methods) associated with
      truth agreement? Mean per-unit Spearman over 6 primary units; within-unit permutation null

Truths are aligned to the estimates' samples once, so each null draw is array arithmetic. (The
first version re-mapped sample keys inside every draw and was stopped at its 30-minute limit.)

    python3 scripts/agreement_selection_test.py [--draws 2000]
"""
from __future__ import annotations

import argparse
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
from absolute_purity_yardstick import key4  # noqa: E402

OUT = config.RESULTS_DIR / "agreement_selection_test.json"
COHORTS = ("gbm", "lgg")
LEUK = ["Macrophage_Microglia", "T_cell", "NK_cell", "B_cell"]
NEEDED = ["Tumor"] + LEUK
TRUTH_OF = {"Tumor": ("purity", key4), "Leukocytes": ("lf", e2.k4s), "Lymphoid": ("epi_lymphoid", e2.k4s),
            "T_cell": ("epi_T", e2.k4s), "B_cell": ("epi_B", e2.k4s), "NK_cell": ("epi_NK", e2.k4s)}
PRIMARY, SECONDARY = ["Tumor", "Leukocytes", "Lymphoid"], ["T_cell", "B_cell", "NK_cell"]
NAMES = PRIMARY + SECONDARY


def unit_values(c: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Compartments from per-cell-type estimates; lymphoid units are shares of the leukocyte total
    (matched to EpiDISH's blood-reference denominator, E2 Addendum 1 check 7)."""
    leuk = sum(c[k] for k in LEUK)
    d = np.where(leuk > 0, leuk, np.nan)
    return {"Tumor": c["Tumor"], "Leukocytes": leuk, "Lymphoid": (c["T_cell"] + c["NK_cell"] + c["B_cell"]) / d,
            "T_cell": c["T_cell"] / d, "B_cell": c["B_cell"] / d, "NK_cell": c["NK_cell"] / d}


def spearman(x: np.ndarray, t: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(t)
    if ok.sum() < 10 or np.ptp(x[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(stats.rankdata(x[ok]), stats.rankdata(t[ok]))[0, 1])


def select_pairs(X: pd.DataFrame) -> list[tuple[str, str]]:
    """DECEPTICON step 2: the two method pairs with the highest Pearson correlation (constant methods out)."""
    usable = [m for m in X.columns if X[m].nunique() > 1]
    C = X[usable].corr(method="pearson")
    return sorted(itertools.combinations(usable, 2), key=lambda p: -C.loc[p[0], p[1]])[:2]


def weights(pairs, methods: list[str]) -> np.ndarray:
    """DECEPTICON step 3: 0.25 per selection; a method selected twice carries 0.5."""
    w = np.zeros(len(methods))
    for a, b in pairs:
        w[methods.index(a)] += .25
        w[methods.index(b)] += .25
    return w


def prepare(h: str, panel: list[str]) -> dict:
    base = "" if h == "gbm" else "_lgg"
    tr = e2.truths(h)
    full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
    sid = "sample_id" if "sample_id" in full.columns else "sample"
    meth = {m: g.set_index(sid)[list(config.CELL_TYPES)] for m, g in full.groupby("method") if m in panel}
    methods = sorted(meth)
    idx = meth[methods[0]].index
    per_type = {ct: pd.DataFrame({m: meth[m].loc[idx, ct].astype(float) for m in methods}) for ct in NEEDED}
    truth = {}
    for u, (name, key) in TRUTH_OF.items():
        keys = [key(i) for i in idx]
        truth[u] = tr[name].reindex(keys).to_numpy(dtype=float)
    return {"methods": methods, "per_type": per_type, "X": {ct: per_type[ct].to_numpy() for ct in NEEDED},
            "truth": truth, "n": len(idx)}


def run_panel(panel_of: dict[str, list[str]], draws: int, rng: np.random.Generator) -> dict:
    res: dict = {"units": {}, "methods": {}, "chosen_pairs": {}}
    p2_rows: dict[str, dict] = {}
    for h in COHORTS:
        P = prepare(h, panel_of[h])
        ms, X, T = P["methods"], P["X"], P["truth"]
        res["methods"][h] = ms
        chosen = {ct: select_pairs(P["per_type"][ct]) for ct in NEEDED}
        res["chosen_pairs"][h] = {ct: [list(p) for p in chosen[ct]] for ct in NEEDED}
        ens = unit_values({ct: X[ct] @ weights(chosen[ct], ms) for ct in NEEDED})
        obs = {u: spearman(ens[u], T[u]) for u in NAMES}
        singles = {m: unit_values({ct: X[ct][:, j] for ct in NEEDED}) for j, m in enumerate(ms)}
        single_rho = {m: {u: spearman(singles[m][u], T[u]) for u in NAMES} for m in ms}
        usable = {ct: [m for m in ms if np.ptp(X[ct][:, ms.index(m)]) > 0] for ct in NEEDED}
        allp = {ct: list(itertools.combinations(usable[ct], 2)) for ct in NEEDED}
        null = {u: np.empty(draws) for u in NAMES}
        for d in range(draws):
            pick = {ct: [allp[ct][i] for i in rng.choice(len(allp[ct]), 2, replace=False)] for ct in NEEDED}
            v = unit_values({ct: X[ct] @ weights(pick[ct], ms) for ct in NEEDED})
            for u in NAMES:
                null[u][d] = spearman(v[u], T[u])
        for u in NAMES:
            nv = null[u][np.isfinite(null[u])]
            sv = {m: single_rho[m][u] for m in ms}
            fin = {m: v for m, v in sv.items() if np.isfinite(v)}
            best = max(fin, key=fin.get) if fin else None
            res["units"][f"{u}|{h}"] = {
                "tier": "primary" if u in PRIMARY else "secondary", "decepticon_rho": obs[u],
                "null_median": float(np.median(nv)) if len(nv) else float("nan"),
                "null_percentile": float((nv < obs[u]).mean()) if len(nv) and np.isfinite(obs[u]) else float("nan"),
                "beats_null_median": bool(len(nv) and np.isfinite(obs[u]) and obs[u] > np.median(nv)),
                "median_single_method": float(np.nanmedian(list(sv.values()))),
                "best_single_method": best, "best_single_rho": fin.get(best, float("nan")) if best else float("nan"),
                "single_method_rho": sv,
                "ensemble_samples_undefined": int(np.isnan(ens[u]).sum())}
            # P2 inputs: centrality (mean Pearson with the other methods) and truth agreement, per method
            df = pd.DataFrame({m: singles[m][u] for m in ms})
            df = df[[m for m in ms if df[m].nunique() > 1]]
            C = df.corr(method="pearson", min_periods=10)
            cen = {m: float(C.loc[m].drop(m).mean()) for m in C.columns}
            ok = [m for m in cen if np.isfinite(cen[m]) and np.isfinite(sv[m])]
            p2_rows[f"{u}|{h}"] = {"centrality": {m: cen[m] for m in ok}, "truth": {m: sv[m] for m in ok}}
    prim = [k for k, v in res["units"].items() if v["tier"] == "primary"]
    wins = sum(res["units"][k]["beats_null_median"] for k in prim)
    res["P1"] = {"n_units": len(prim), "beats_null_median": wins,
                 "sign_test_p_one_sided": float(stats.binomtest(wins, len(prim), .5, alternative="greater").pvalue),
                 "reading": "SUPPORTED" if wins == len(prim) else "NOT SUPPORTED" if wins <= 3 else "INCONCLUSIVE"}

    def per_unit(keys, shuffle=None):
        out = []
        for k in keys:
            ms_k = list(p2_rows[k]["centrality"])
            if len(ms_k) < 4:
                continue
            c = np.array([p2_rows[k]["centrality"][m] for m in ms_k])
            t = np.array([p2_rows[k]["truth"][m] for m in ms_k])
            if shuffle is not None:
                t = shuffle.permutation(t)
            out.append(float(np.corrcoef(stats.rankdata(c), stats.rankdata(t))[0, 1]))
        return out
    pk = [k for k in p2_rows if k.split("|")[0] in PRIMARY]
    sk = [k for k in p2_rows if k.split("|")[0] in SECONDARY]
    obs_r = per_unit(pk)
    stat = float(np.mean(obs_r))
    prng = np.random.default_rng(config.RANDOM_SEED + 1)
    nulls = np.array([np.mean(per_unit(pk, prng)) for _ in range(10_000)])
    p = float((nulls >= stat - 1e-12).mean())
    res["P2"] = {"per_unit_rho": dict(zip(pk, [round(x, 4) for x in obs_r])), "mean_rho": round(stat, 4),
                 "p_one_sided": round(p, 4), "n_methods_per_unit": {k: len(p2_rows[k]["centrality"]) for k in pk},
                 "reading": "SUPPORTED" if stat > 0 and p < 0.05 else "NOT SUPPORTED" if stat <= 0 else "INCONCLUSIVE"}
    res["P2_secondary_units"] = dict(zip(sk, [round(x, 4) for x in per_unit(sk)]))
    res["p2_inputs"] = p2_rows
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--draws", type=int, default=2000)
    a = ap.parse_args()
    comparable, full10 = {}, {}
    for h in COHORTS:
        base = "" if h == "gbm" else "_lgg"
        comparable[h] = list(json.loads((config.RESULTS_DIR / f"immune_arm{base}.json").read_text())["comparable_methods"])
        full10[h] = sorted(set(pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv", usecols=["method"])["method"])
                           - e2.DROP)
    rng = np.random.default_rng(config.RANDOM_SEED)
    out = {"rule": "prespecified/agreement_selection_test.md",
           "implementation": "this project's code implementing DECEPTICON's published selection rule; not the package",
           "draws": a.draws,
           "primary_panel": run_panel(comparable, a.draws, rng),
           "secondary_panel_10_methods": run_panel(full10, a.draws, rng)}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    P = out["primary_panel"]
    for k, u in P["units"].items():
        print(f"{k:18s} DECEPTICON {u['decepticon_rho']:+.3f} | null median {u['null_median']:+.3f} "
              f"(pct {u['null_percentile']:.2f}) | median single {u['median_single_method']:+.3f} | "
              f"best {u['best_single_method']} {u['best_single_rho']:+.3f}")
    print("chosen pairs:", json.dumps(P["chosen_pairs"]))
    print("P1:", P["P1"])
    print("P2:", {k: v for k, v in P["P2"].items() if k != "per_unit_rho"}, "| per unit:", P["P2"]["per_unit_rho"])
    S = out["secondary_panel_10_methods"]
    print("10-method panel  P1:", S["P1"]["reading"], S["P1"]["beats_null_median"], "| P2:", S["P2"]["mean_rho"],
          S["P2"]["p_one_sided"], S["P2"]["reading"])
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
