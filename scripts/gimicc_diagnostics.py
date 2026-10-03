"""
gimicc_diagnostics.py -- prespecified/gimicc_truth_confirmation.md, Addendum 5: what carries GIMiCC's
B > T in TCGA-LGG, and why its CpG-shuffle negative control fails there. EXPLORATORY; nothing here
changes the registered reading (INCONCLUSIVE).

  D-a  the lymphoid split of the shuffled-CpG runs (Q1's computation): does B > T survive shuffling?
  D-b  per sample, GIMiCC level-3 B/(T+B) against the sample's mean beta over the complete-case library
       CpGs (a CIMP-like global summary), per cohort
  D-c  that mean beta against ABSOLUTE purity, per cohort: can label shuffling break the tumour layer?
  D-d  Q4's reference-free RNA marker index (T minus B), compared between the samples GIMiCC calls
       B-dominant and T-dominant (Mann-Whitney), per cohort
  D-e  (Addendum 6) GIMiCC's level-3 B/(T+B) against ABSOLUTE purity and against its own Tumor: is B
       an overflow for tumour signal, as in the RNA arm?
  D-f  (Addendum 6) partial Spearman of unmix's (and the methods' median) lymphoid total against
       GIMiCC's lymphoid tissue fraction, controlling for the Thorsson leukocyte fraction

    python3 scripts/gimicc_diagnostics.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import gimicc_truth as gt  # noqa: E402
from lymphoid_ordering import k4s  # noqa: E402

OUT = config.RESULTS_DIR / "gimicc_diagnostics.json"


def mean_beta(cohort: str) -> pd.Series:
    """Per-sample mean beta over the complete-case library CpGs -- exactly the CpGs GIMiCC saw."""
    m = pd.read_csv(gt.G / f"tcga_{cohort}_gimicc_cpgs.tsv", sep="\t", index_col=0, na_values=["", "NA", "NaN"])
    m = m.dropna(axis=0, how="any")
    s = m.mean(axis=0)
    s.index = [k4s(i) for i in s.index]
    return s[~s.index.duplicated()]


def b_share(df: pd.DataFrame) -> pd.Series:
    L = gt.lymph_shares(df)
    s = (L["B"] / (L["T"] + L["B"]).replace(0, np.nan)).dropna()
    s.index = [k4s(i) for i in s.index]
    return s[~s.index.duplicated()]


def sp(a: pd.Series, b: pd.Series) -> dict:
    idx = a.dropna().index.intersection(b.dropna().index)
    r = stats.spearmanr(a[idx], b[idx])
    return {"rho": round(float(r.statistic), 4), "p": float(r.pvalue), "n": int(len(idx))}


def partial_spearman(x: pd.Series, y: pd.Series, z: pd.Series) -> dict:
    """Spearman of x and y controlling for z: Pearson correlation of rank residuals."""
    idx = x.dropna().index.intersection(y.dropna().index).intersection(z.dropna().index)
    if len(idx) < 20:
        return {"rho": None, "n": int(len(idx))}
    rx, ry, rz = (stats.rankdata(v[idx]) for v in (x, y, z))
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)  # noqa: E731
    r = stats.pearsonr(res(rx), res(ry))
    return {"rho": round(float(r.statistic), 4), "p": float(r.pvalue), "n": int(len(idx)),
            "plain_rho": round(float(stats.spearmanr(x[idx], y[idx]).statistic), 4)}


def main() -> int:
    from lymphoid_tracking import marker_index  # noqa: PLC0415
    out = {"rule": "prespecified/gimicc_truth_confirmation.md, Addendum 5", "cohorts": {}}
    for c in gt.COHORTS:
        d = gt.DEFAULT_TYPE[c]
        tr = gt.truths(c)
        real4, shuf4, real3 = gt.load(c, f"{d}_h4"), gt.load(c, f"{d}_h4_shuffled"), gt.load(c, f"{d}_h3")
        mb = mean_beta(c)
        bs = b_share(real3)
        rec = {
            "D-a_shuffled_lymphoid_split": {"real": gt.q1(real4), "shuffled": gt.q1(shuf4)},
            "D-b_B_share_vs_mean_beta": sp(bs, mb),
            "D-c_mean_beta_vs_ABSOLUTE": sp(mb, tr["purity_meth"]),
            "D-c_real_tumour_layer_vs_ABSOLUTE": sp(gt.rekey(real4["Tumor"], k4s), tr["purity_meth"]),
            "D-c_shuffled_tumour_layer_vs_mean_beta": sp(gt.rekey(shuf4["Tumor"], k4s), mb),
        }
        idx, _ = marker_index(c)
        common = bs.index.intersection(idx.index)
        b_dom = [k for k in common if bs[k] > 0.5]
        t_dom = [k for k in common if bs[k] < 0.5]
        if len(b_dom) >= 5 and len(t_dom) >= 5:
            mw = stats.mannwhitneyu(idx[t_dom], idx[b_dom], alternative="greater")
            rec["D-d_marker_index_by_gimicc_dominance"] = {
                "n_gimicc_B_dominant": len(b_dom), "n_gimicc_T_dominant": len(t_dom),
                "median_index_T_dominant": round(float(idx[t_dom].median()), 4),
                "median_index_B_dominant": round(float(idx[b_dom].median()), 4),
                "mannwhitney_p_T_dominant_higher": float(mw.pvalue)}
        else:
            rec["D-d_marker_index_by_gimicc_dominance"] = {"n_gimicc_B_dominant": len(b_dom),
                                                           "n_gimicc_T_dominant": len(t_dom), "note": "too few"}
        # D-e: overflow
        rec["D-e_B_share_vs_ABSOLUTE"] = sp(bs, tr["purity_meth"])
        rec["D-e_B_share_vs_GIMiCC_Tumor"] = sp(bs, gt.rekey(real3["Tumor"], k4s))
        e = rec["D-e_B_share_vs_ABSOLUTE"]
        rec["D-e_reading"] = "consistent with overflow" if (e["rho"] > 0 and e["p"] < 0.05) else "not consistent with overflow"
        # D-f: is Q3's tracking the leukocyte level?
        import identifiability_diagnostics as e2  # noqa: PLC0415
        lym_g = gt.rekey(real3[["Tcell", "NK", "Bcell"]].sum(axis=1, min_count=3), k4s)
        unmix = pd.read_csv(config.RESULTS_DIR / "identifiability" / f"unmix_{c}_predeclared.csv", index_col=0)
        u = gt.rekey(e2.comp(unmix, "Lymphoid"), k4s)
        base = "" if c == "gbm" else "_lgg"
        full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
        sid = "sample_id" if "sample_id" in full.columns else "sample"
        per_m = {}
        for m, g in full.groupby("method"):
            if m in e2.DROP:
                continue
            em = gt.rekey(g.set_index(sid)[list(config.CELL_TYPES)].pipe(lambda d: e2.comp(d, "Lymphoid")), k4s)
            if em.nunique() > 2:
                per_m[m] = partial_spearman(em, lym_g, tr["lf"])
        vals = [v["rho"] for v in per_m.values() if v.get("rho") is not None]
        pu = partial_spearman(u, lym_g, tr["lf"])
        rec["D-f_partial_on_leukocyte_fraction"] = {
            "unmix": pu, "methods_median_partial_rho": round(float(np.median(vals)), 4) if vals else None,
            "per_method": per_m,
            "reading": ("tracking is the leukocyte level (partial < 0.20)" if (pu.get("rho") is not None and pu["rho"] < 0.20)
                        else "tracking survives the leukocyte level (partial >= 0.20)")}
        out["cohorts"][c] = rec
        print(f"--- {c}")
        print("  D-a real    ", rec["D-a_shuffled_lymphoid_split"]["real"]["ordering"],
              rec["D-a_shuffled_lymphoid_split"]["real"]["mean_within_lymphoid_share"])
        print("  D-a shuffled", rec["D-a_shuffled_lymphoid_split"]["shuffled"]["ordering"],
              rec["D-a_shuffled_lymphoid_split"]["shuffled"]["mean_within_lymphoid_share"],
              "n", rec["D-a_shuffled_lymphoid_split"]["shuffled"]["n_samples_with_lymphoid"])
        for k in ("D-b_B_share_vs_mean_beta", "D-c_mean_beta_vs_ABSOLUTE", "D-c_real_tumour_layer_vs_ABSOLUTE",
                  "D-c_shuffled_tumour_layer_vs_mean_beta"):
            print(f"  {k}: {rec[k]}")
        print("  D-d", rec["D-d_marker_index_by_gimicc_dominance"])
        print("  D-e B share vs ABSOLUTE", rec["D-e_B_share_vs_ABSOLUTE"], "| vs GIMiCC Tumor",
              rec["D-e_B_share_vs_GIMiCC_Tumor"], "->", rec["D-e_reading"])
        f_ = rec["D-f_partial_on_leukocyte_fraction"]
        print("  D-f unmix", f_["unmix"], "| methods' median partial", f_["methods_median_partial_rho"], "->", f_["reading"])
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
