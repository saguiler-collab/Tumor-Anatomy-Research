"""
identifiability_robustness.py -- Addendum 1 of prespecified/identifiability_diagnostics.md.

Post hoc checks on Extension E2's H1, written down before they were computed. None changes a
registered reading except check 5, whose consequence the addendum fixed in advance.

  1  contrast vs fine ordering: P(the two lowest-truth units are the two least stable) = 1/C(6,2)
  2  near-tie swap: H1 with the GBM Tumor / GBM Leukocytes truth values exchanged
  3  bootstrap over samples (B = 2000): percentiles of each H1 rho; how often the lymphoid units are
     the two lowest
  4  each registered method's all-leukocyte truth agreement (for the H2 wording)
  5  duplicate setting: in GBM `predeclared` and `s1` are the same fit; D1 over distinct fits only
  6  degraded methods: D2 and its truth agreement on the registered immune arm's comparable panel
     (Bisque and EPIC, which run in a degraded mode here, left out as that arm leaves them out)

    python3 scripts/identifiability_robustness.py [--B 2000]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
import identifiability_diagnostics as e2  # noqa: E402

OUT = config.RESULTS_DIR / "identifiability_robustness.json"
COHORTS = ("gbm", "lgg")
UNITS = [(c, h) for c in e2.PRIMARY for h in COHORTS]          # the order e2's test() uses


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distinct_settings(cohort: str) -> list[str]:
    """Setting labels with byte-identical estimate files collapsed to their first occurrence."""
    seen, keep = {}, []
    for label, _, _ in e2.SETTINGS:
        h = sha(e2.OUT_DIR / f"unmix_{cohort}_{label}.csv")
        if h not in seen:
            seen[h] = label
            keep.append(label)
    return keep


def load(cohort: str) -> dict:
    """Per-compartment sample-by-setting and sample-by-method matrices, with truth vectors aligned."""
    base = "" if cohort == "gbm" else "_lgg"
    tr = e2.truths(cohort)
    d1 = {lab: pd.read_csv(e2.OUT_DIR / f"unmix_{cohort}_{lab}.csv", index_col=0) for lab, _, _ in e2.SETTINGS}
    full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
    sid = "sample_id" if "sample_id" in full.columns else "sample"
    meth = {m: g.set_index(sid)[list(config.CELL_TYPES)] for m, g in full.groupby("method") if m not in e2.DROP}
    idx = d1["predeclared"].index
    for frame in list(d1.values()) + list(meth.values()):
        if set(frame.index) != set(idx):
            raise SystemExit(f"{cohort}: the arms do not cover the same samples")
    out = {"labels": [lab for lab, _, _ in e2.SETTINGS], "methods": list(meth), "n": len(idx), "comp": {}}
    for c in e2.COMPARTMENTS:
        X1 = np.column_stack([e2.comp(d1[lab], c).loc[idx].to_numpy(float) for lab in out["labels"]])
        X2 = np.column_stack([e2.comp(meth[m], c).loc[idx].to_numpy(float) for m in out["methods"]])
        t = np.full(len(idx), np.nan)
        if c in e2.TRUTH_OF:
            name, key = e2.TRUTH_OF[c]
            keys = pd.Index([key(i) for i in idx])
            if keys.duplicated().any():
                raise SystemExit(f"{cohort}/{c}: duplicate sample keys; alignment would differ from E2")
            t = tr[name].reindex(keys).to_numpy(float)
        out["comp"][c] = {"X1": X1, "X2": X2, "truth": t}
    return out


def rank_cols(X: np.ndarray) -> np.ndarray:
    return stats.rankdata(X, axis=0)


def mean_pairwise(X: np.ndarray, cols: list[int] | None = None) -> float:
    """E2's D1/D2: mean pairwise Spearman over non-constant columns (nan when fewer than 3)."""
    X = X if cols is None else X[:, cols]
    X = X[:, np.ptp(X, axis=0) > 0]
    if X.shape[1] < 3:
        return float("nan")
    C = np.corrcoef(rank_cols(X), rowvar=False)
    return float(C[np.triu_indices_from(C, 1)].mean())


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10 or np.ptp(x[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(stats.rankdata(x[ok]), stats.rankdata(y[ok]))[0, 1])


def nanmedian(v: list[float]) -> float:
    return float(np.nanmedian(v)) if np.isfinite(v).any() else float("nan")


def unit_stats(u: dict, rows: np.ndarray | None, d1_cols: list[int], d2_cols: list[int] | None = None) -> dict:
    X1, X2, t = u["X1"], u["X2"], u["truth"]
    if rows is not None:
        X1, X2, t = X1[rows], X2[rows], t[rows]
    t2 = [spearman(X2[:, j], t) for j in range(X2.shape[1])]
    out = {"d1": mean_pairwise(X1), "d1_distinct": mean_pairwise(X1, d1_cols),
           "d2": mean_pairwise(X2), "t1": spearman(X1[:, 0], t), "t2": nanmedian(t2)}
    if d2_cols is not None:                                     # check 6, point estimates only
        out["d2_comparable"] = mean_pairwise(X2, d2_cols)
        out["t2_comparable"] = nanmedian([t2[j] for j in d2_cols])
    return out


def immune_arm_panel(cohort: str) -> tuple[list[str], dict[str, str]]:
    """The registered immune arm's comparable methods and its reasons for excluding the others."""
    base = "" if cohort == "gbm" else "_lgg"
    arm = json.loads((config.RESULTS_DIR / f"immune_arm{base}.json").read_text())
    return list(arm["comparable_methods"]), {m: r.split(";")[0] for m, r in arm["excluded_from_comparison"].items()}


def h1(st: dict, x: str, y: str, names=e2.PRIMARY) -> tuple[float, float]:
    pts = [(st[(c, h)][x], st[(c, h)][y]) for c in names for h in COHORTS]
    return e2.exact_one_sided([p[0] for p in pts], [p[1] for p in pts])


def lymphoid_lowest(st: dict, x: str) -> bool:
    order = sorted(UNITS, key=lambda k: st[k][x])
    return {order[0][0], order[1][0]} == {"Lymphoid"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args()
    reg = json.loads(e2.OUT.read_text())
    data = {h: load(h) for h in COHORTS}
    keep = {h: distinct_settings(h) for h in COHORTS}
    cols = {h: [data[h]["labels"].index(lab) for lab in keep[h]] for h in COHORTS}
    print("distinct fits:", {h: f"{len(keep[h])} of {len(e2.SETTINGS)}" for h in COHORTS}, flush=True)
    panel = {h: immune_arm_panel(h) for h in COHORTS}
    comparable = {h: [j for j, m in enumerate(data[h]["methods"]) if m in panel[h][0]] for h in COHORTS}
    degraded = {h: {m: panel[h][1][m] for m in data[h]["methods"] if m in panel[h][1]} for h in COHORTS}

    # point estimates; control: the registered E2 values reproduce from these matrices
    point = {(c, h): unit_stats(data[h]["comp"][c], None, cols[h], comparable[h])
             for c in e2.COMPARTMENTS for h in COHORTS}
    worst = 0.0
    for (c, h), s in point.items():
        r = reg["units"][f"{c}|{h}"]
        for mine, theirs in ((s["d1"], r["d1_stability"]), (s["d2"], r["d2_agreement"]),
                             (s["t1"], r["truth_rho_unmix"]), (s["t2"], r["truth_rho_median_methods"])):
            if np.isfinite(theirs):
                worst = max(worst, abs(mine - theirs))
    control_ok = worst < 1e-9
    print(f"control: registered E2 values reproduce, max |diff| = {worst:.2e}", flush=True)
    if not control_ok:
        OUT.write_text(json.dumps({"control_pass": False, "max_abs_diff": worst}, indent=2))
        return 2

    out: dict = {"rule": "prespecified/identifiability_diagnostics.md, Addendum 1", "control_pass": True,
                 "control_max_abs_diff": worst}
    # 1 -- contrast vs fine ordering
    out["check1_contrast_probability"] = 1 / comb(len(UNITS), 2)
    # 2 -- near-tie swap of the GBM Tumor / Leukocytes truth values
    sw = {k: dict(v) for k, v in point.items()}
    sw[("Tumor", "gbm")]["t1"], sw[("Leukocytes", "gbm")]["t1"] = point[("Leukocytes", "gbm")]["t1"], point[("Tumor", "gbm")]["t1"]
    r_sw, p_sw = h1(sw, "d1", "t1")
    out["check2_near_tie_swap"] = {"truth_gap": abs(point[("Tumor", "gbm")]["t1"] - point[("Leukocytes", "gbm")]["t1"]),
                                   "rho": round(r_sw, 4), "p_one_sided": round(p_sw, 4)}
    # 5 -- distinct fits only (computed before the bootstrap so it is reported regardless)
    r5, p5 = h1(point, "d1_distinct", "t1")
    r5s, p5s = h1(point, "d1_distinct", "t1", e2.PRIMARY + e2.SECONDARY)
    out["check5_distinct_fits"] = {
        "fits_used": keep,
        "duplicate_sha256": {h: sha(e2.OUT_DIR / f"unmix_{h}_predeclared.csv")[:16] for h in COHORTS},
        "d1_distinct": {f"{c}|{h}": round(point[(c, h)]["d1_distinct"], 4) for c in e2.PRIMARY + e2.SECONDARY
                        for h in COHORTS},
        "H1_D1": {"rho": round(r5, 4), "p_one_sided": round(p5, 4),
                  "reading": "SUPPORTED" if r5 > 0 and p5 < 0.05 else "not significant"},
        "H1s_D1": {"rho": round(r5s, 4), "p_one_sided": round(p5s, 4)},
        "lymphoid_two_least_stable": lymphoid_lowest(point, "d1_distinct")}
    print("check 5:", json.dumps(out["check5_distinct_fits"]["H1_D1"]), flush=True)
    # 6 -- D2 on the registered immune arm's comparable panel (degraded-mode methods left out)
    r6, p6 = h1(point, "d2_comparable", "t2_comparable")
    r6s, p6s = h1(point, "d2_comparable", "t2_comparable", e2.PRIMARY + e2.SECONDARY)
    out["check6_comparable_panel"] = {
        "methods": {h: [data[h]["methods"][j] for j in comparable[h]] for h in COHORTS},
        "left_out_degraded": degraded,
        "units": {f"{c}|{h}": {"d2": round(point[(c, h)]["d2_comparable"], 4),
                               "truth_median": round(point[(c, h)]["t2_comparable"], 4)}
                  for c in e2.PRIMARY + e2.SECONDARY for h in COHORTS},
        "H1_D2": {"rho": round(r6, 4), "p_one_sided": round(p6, 4),
                  "reading": "SUPPORTED" if r6 > 0 and p6 < 0.05 else "not significant"},
        "H1s_D2": {"rho": round(r6s, 4), "p_one_sided": round(p6s, 4)},
        "lymphoid_two_lowest": lymphoid_lowest(point, "d2_comparable")}
    print("check 6:", json.dumps(out["check6_comparable_panel"]["H1_D2"]), flush=True)
    # 4 -- each registered method's all-leukocyte truth agreement (degraded modes labelled)
    out["check4_leukocyte_truth_by_method"] = {}
    for h in COHORTS:
        u = data[h]["comp"]["Leukocytes"]
        per = {m: round(spearman(u["X2"][:, j], u["truth"]), 4) for j, m in enumerate(data[h]["methods"])}
        out["check4_leukocyte_truth_by_method"][h] = {
            "unmix": round(point[("Leukocytes", h)]["t1"], 4), "methods": per,
            "degraded_mode": degraded[h],
            "unmix_exceeds_all": bool(all(point[("Leukocytes", h)]["t1"] > v for v in per.values()))}
    # 3 -- bootstrap over samples within cohort
    rng = np.random.default_rng(config.RANDOM_SEED)
    draws = {k: [] for k in ("d1", "d1_distinct", "d2")}
    low = {k: 0 for k in draws}
    for _ in range(a.B):
        st = {}
        for h in COHORTS:
            rows = rng.integers(0, data[h]["n"], data[h]["n"])
            for c in e2.PRIMARY:
                st[(c, h)] = unit_stats(data[h]["comp"][c], rows, cols[h])
        for k, y in (("d1", "t1"), ("d1_distinct", "t1"), ("d2", "t2")):
            pts = [(st[u][k], st[u][y]) for u in UNITS]
            draws[k].append(float(stats.spearmanr([p[0] for p in pts], [p[1] for p in pts]).statistic))
            low[k] += lymphoid_lowest(st, k)
    out["check3_bootstrap"] = {"B": a.B, "seed": config.RANDOM_SEED}
    for k in draws:
        q = np.nanpercentile(draws[k], [2.5, 50, 97.5])
        out["check3_bootstrap"][k] = {"h1_rho_2.5_50_97.5": [round(float(v), 4) for v in q],
                                      "share_lymphoid_two_lowest": round(low[k] / a.B, 4),
                                      "contrast_robust": low[k] / a.B >= 0.95,
                                      "fine_ordering_robust": bool(q[0] > 0)}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({k: v for k, v in out.items() if k.startswith(("control", "check1", "check2", "check3"))},
                     indent=1, default=float))
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
