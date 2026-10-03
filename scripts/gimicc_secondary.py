"""
gimicc_secondary.py -- Q2 and Q3 of prespecified/gimicc_truth_confirmation.md (Addendum 3), and Q4
(Addendum 4).

POST-REGISTRATION, EXPLORATORY. Nothing here changes a registered statistic.

  Q2  The headline re-scored against GIMiCC. Each method's recorded cohort ordering (registered panel,
      both reference builds, both cohorts: results/lymphoid_ordering{,_lgg}{,_h5ad}.json, read as
      recorded, not recomputed) against GIMiCC's cohort direction (Q1's level-4 computation on the
      RNA-matched samples), beside EpiDISH's. Per-sample discordance against GIMiCC's level-3 shares
      (level 4 beside). The extension panel (results/extension/tcga_*.json) is reported alongside,
      outside the reading.
  Q3  E2's units against GIMiCC tissue fractions (matched denominators by construction). D1 and D2 are
      truth-free and taken from results/identifiability_diagnostics.json; only the truth changes.
      Reproduction control first (saved predeclared unmix -> E2's recorded Tumor rho). H1 recomputed
      with GIMiCC's Lymphoid truth, beside the registered H1.
  Q4  D23's positive control on GIMiCC's per-sample T/(T+B): the reference-free T-minus-B marker index
      (scripts/lymphoid_tracking.py, panels fixed 2026-09-30) against GIMiCC level 3 (level 4 beside),
      10,000 permutations. Then every method's per-sample tracking of GIMiCC's T/(T+B) -- informative
      only if the control passes in both cohorts.

Needs results/gimicc_truth_confirmation.json (scripts/gimicc_truth.py --stage analyse).

    python3 scripts/gimicc_secondary.py
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
import identifiability_diagnostics as e2  # noqa: E402
from absolute_purity_yardstick import key4  # noqa: E402
from lymphoid_ordering import COMPARE, k4s, renormalise  # noqa: E402

OUT = config.RESULTS_DIR / "gimicc_secondary.json"
BAR = 0.40                      # E2's own "tracks" bar (H2)
Q3_UNITS = ["Tumor", "Leukocytes", "Lymphoid", "T_cell", "B_cell", "NK_cell"]


def rho(est_c: pd.Series, truth: pd.Series, key=k4s) -> tuple[float, int]:
    """E2's truth_rho, with an arbitrary truth series (keys already normalised)."""
    e = gt.rekey(est_c, key)
    shared = e.index.intersection(truth.dropna().index)
    if len(shared) < 10 or e[shared].nunique() < 2:
        return float("nan"), len(shared)
    return float(stats.spearmanr(e[shared], truth[shared]).statistic), len(shared)


def gimicc_truths(cohort: str) -> dict:
    """Per-sample GIMiCC tissue fractions at level 3 (Layer 3A's T) and level 4 (CD4 + CD8).
    A sample with a NaN in any lymphoid column is excluded from every lymphoid truth (counted)."""
    d = gt.DEFAULT_TYPE[cohort]
    out = {}
    for lvl, tag, imm_cols in (("level3", f"{d}_h3", gt.IMMUNE_H3), ("level4", f"{d}_h4", gt.IMMUNE)):
        df = gt.load(cohort, tag)
        t = df["Tcell"] if "Tcell" in df.columns else df["CD4Tcell"] + df["CD8Tcell"]
        lym = pd.DataFrame({"T_cell": t, "NK_cell": df["NK"], "B_cell": df["Bcell"]})
        ok = lym.notna().all(axis=1)
        tr = {k: lym.loc[ok, k] for k in lym.columns}
        tr["Lymphoid"] = lym[ok].sum(axis=1)
        tr["Leukocytes"], n_ex_imm = gt.immune_total(df, imm_cols)
        out[lvl] = {"truth": {k: gt.rekey(v, k4s) for k, v in tr.items()}, "n_samples": int(len(df)),
                    "n_excluded_nan_lymphoid": int((~ok).sum()), "n_excluded_nan_immune": n_ex_imm,
                    "shares": gt.rekey(gt.lymph_shares(df), k4s)}
    return out


# ------------------------------------------------------------------------------------------- Q2
def per_sample_discordance(est_csv: Path, shares: pd.DataFrame) -> dict:
    full = pd.read_csv(est_csv)
    scol = next(c for c in ("sample", "sample_id") if c in full.columns)
    full["k"] = full[scol].map(k4s)
    res = {}
    for m, g in full.groupby("method"):
        g = g.drop_duplicates("k").set_index("k")
        idx = [k for k in g.index if k in shares.index]
        if len(idx) < 50:
            continue
        e = renormalise(g.loc[idx])
        live = e.notna().all(axis=1)
        e, s = e[live], shares.loc[e.index[live]]
        m_bt, m_tb = (e["B_cell"] > e["T_cell"]).to_numpy(), (e["T_cell"] > e["B_cell"]).to_numpy()
        g_tb, g_bt = (s["T"] > s["B"]).to_numpy(), (s["B"] > s["T"]).to_numpy()
        n = int(live.sum())
        res[m] = {"n_scored": n, "n_no_lymphoid_signal": int((~live).sum()),
                  "frac_method_B_over_T_while_gimicc_T_over_B": round(float((m_bt & g_tb).mean()), 4) if n else None,
                  "frac_method_T_over_B_while_gimicc_B_over_T": round(float((m_tb & g_bt).mean()), 4) if n else None,
                  "frac_same_direction": round(float(((m_tb & g_tb) | (m_bt & g_bt)).mean()), 4) if n else None}
    return res


def q2(gim: dict, gtr: dict) -> dict:
    out = {}
    for c in gt.COHORTS:
        g = gim["cohorts"][c]
        gdir = bool(g["Q1_rna_matched_samples"]["T_exceeds_B"])
        for ref in ("frozen", "h5ad"):
            base = ("" if c == "gbm" else "_lgg") + ("_h5ad" if ref == "h5ad" else "")
            rec = json.loads((config.RESULTS_DIR / f"lymphoid_ordering{base}.json").read_text())
            tm = rec["truth_mean"]
            epi_dir = bool(tm["T_cell"] > tm["B_cell"])
            meths = rec["methods"]
            comp = {m: v for m, v in meths.items() if v.get("comparable", True)}
            ext = json.loads((config.RESULTS_DIR / "extension" / f"tcga_{c}_{ref}.json").read_text())["methods"]
            ext_dir = {m: v.get("T_exceeds_B") for m, v in ext.items()}
            out[f"{c}|{ref}"] = {
                "gimicc_T_exceeds_B_rna_matched_level4": gdir,
                "gimicc_T_exceeds_B_rna_matched_level3": bool(g["Q1_level3_rna_matched"]["T_exceeds_B"]),
                "epidish_T_exceeds_B": epi_dir,
                "n_methods": len(meths), "n_comparable": len(comp),
                "n_methods_T_exceeds_B": int(sum(v["T_exceeds_B"] for v in meths.values())),
                "n_match_gimicc": int(sum(v["T_exceeds_B"] == gdir for v in meths.values())),
                "n_match_epidish": int(sum(v["T_exceeds_B"] == epi_dir for v in meths.values())),
                "n_comparable_match_gimicc": int(sum(v["T_exceeds_B"] == gdir for v in comp.values())),
                "n_comparable_match_epidish": int(sum(v["T_exceeds_B"] == epi_dir for v in comp.values())),
                "methods_T_exceeds_B": sorted(m for m, v in meths.items() if v["T_exceeds_B"]),
                "per_sample_vs_gimicc_level3": per_sample_discordance(
                    config.RESULTS_DIR / f"estimates_full{base}.csv", gtr[c]["level3"]["shares"]),
                "per_sample_vs_gimicc_level4": per_sample_discordance(
                    config.RESULTS_DIR / f"estimates_full{base}.csv", gtr[c]["level4"]["shares"]),
                "extension_panel_outside_reading": {
                    "n_methods_with_ordering": int(sum(v is not None for v in ext_dir.values())),
                    "n_match_gimicc": int(sum(v is not None and bool(v) == gdir for v in ext_dir.values())),
                    "n_match_epidish": int(sum(v is not None and bool(v) == epi_dir for v in ext_dir.values())),
                    "T_exceeds_B": ext_dir}}
    dirs = {c: bool(gim["cohorts"][c]["Q1_rna_matched_samples"]["T_exceeds_B"]) for c in gt.COHORTS}
    reading = ("UNCHANGED" if all(dirs.values()) else
               "REVERSED for " + ", ".join(f"TCGA-{c.upper()}" for c, d in dirs.items() if not d))
    return {"by_cohort_reference": out, "gimicc_direction_T_exceeds_B": dirs, "reading": reading}


# ------------------------------------------------------------------------------------------- Q3
def q3(gtr: dict) -> dict:
    rec = json.loads((config.RESULTS_DIR / "identifiability_diagnostics.json").read_text())
    units, controls = {}, {}
    for c in gt.COHORTS:
        base = "" if c == "gbm" else "_lgg"
        tr_e2 = e2.truths(c)
        unmix = pd.read_csv(config.RESULTS_DIR / "identifiability" / f"unmix_{c}_predeclared.csv", index_col=0)
        full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
        sid = "sample_id" if "sample_id" in full.columns else "sample"
        meth = {m: g.set_index(sid)[list(config.CELL_TYPES)] for m, g in full.groupby("method") if m not in e2.DROP}
        # reproduction control: the saved estimates must give E2's recorded Tumor rhos
        t_unmix = rho(e2.comp(unmix, "Tumor"), tr_e2["purity"], key4)[0]
        t_med = float(np.nanmedian([rho(e2.comp(e, "Tumor"), tr_e2["purity"], key4)[0] for e in meth.values()]))
        r_unmix = rec["units"][f"Tumor|{c}"]["truth_rho_unmix"]
        r_med = rec["units"][f"Tumor|{c}"]["truth_rho_median_methods"]
        controls[c] = {"unmix_tumor_rho": [round(t_unmix, 6), round(r_unmix, 6)],
                       "median_methods_tumor_rho": [round(t_med, 6), round(r_med, 6)],
                       "pass": abs(t_unmix - r_unmix) < 1e-4 and abs(t_med - r_med) < 1e-4}
        for u in Q3_UNITS:
            r = rec["units"][f"{u}|{c}"]
            row = {"d1_stability": r["d1_stability"], "d2_agreement": r["d2_agreement"],
                   "e2_truth": e2.TRUTH_OF[u][0], "e2_truth_rho_unmix": r["truth_rho_unmix"],
                   "e2_truth_rho_median_methods": r["truth_rho_median_methods"]}
            for lvl in ("level3", "level4"):
                if u == "Tumor":
                    continue
                truth = gtr[c][lvl]["truth"][u]
                tu, n = rho(e2.comp(unmix, u), truth)
                tm = [rho(e2.comp(e, u), truth)[0] for e in meth.values()]
                row[f"gimicc_{lvl}_truth_rho_unmix"] = tu
                row[f"gimicc_{lvl}_truth_rho_median_methods"] = float(np.nanmedian(tm)) if np.isfinite(tm).any() else float("nan")
                row[f"gimicc_{lvl}_n"] = n
            units[f"{u}|{c}"] = row
    if not all(v["pass"] for v in controls.values()):
        return {"controls": controls, "note": "reproduction control failed; Q3 not reported"}

    def h1(y_unmix: dict, y_med: dict) -> dict:
        """E2's exact one-sided test over the primary units, with the given truth columns."""
        res = {}
        for name, xk, ys in (("D1", "d1_stability", y_unmix), ("D2", "d2_agreement", y_med)):
            pts = [(units[k][xk], ys[k]) for k in ys if np.isfinite(units[k][xk]) and np.isfinite(ys[k])]
            r, p = e2.exact_one_sided([q[0] for q in pts], [q[1] for q in pts])
            res[name] = {"n_units": len(pts), "rho": round(r, 4), "p_one_sided": round(p, 4)}
        return res

    prim = [f"{u}|{c}" for u in ("Tumor", "Leukocytes", "Lymphoid") for c in gt.COHORTS]

    def ycol(k: str, which: str, lymph_from: str, leuk_from: str) -> float:
        u = units[k]
        comp_ = k.split("|")[0]
        src = {"Lymphoid": lymph_from, "Leukocytes": leuk_from}.get(comp_, "e2")
        return u[f"e2_truth_rho_{which}"] if src == "e2" else u[f"gimicc_level3_truth_rho_{which}"]

    h1_lymph = h1({k: ycol(k, "unmix", "gimicc", "e2") for k in prim},
                  {k: ycol(k, "median_methods", "gimicc", "e2") for k in prim})
    h1_both = h1({k: ycol(k, "unmix", "gimicc", "gimicc") for k in prim},
                 {k: ycol(k, "median_methods", "gimicc", "gimicc") for k in prim})
    lym = {c: (units[f"Lymphoid|{c}"]["gimicc_level3_truth_rho_unmix"],
               units[f"Lymphoid|{c}"]["gimicc_level3_truth_rho_median_methods"]) for c in gt.COHORTS}
    holds = all(np.isfinite(v).all() and max(v) < BAR for v in lym.values())
    return {"controls": controls, "units": units,
            "H1_gimicc_lymphoid_truth": h1_lymph, "H1_gimicc_lymphoid_and_leukocyte_truth": h1_both,
            "H1_registered_for_reference": rec["H1_primary"],
            "lymphoid_truth_rho_level3_unmix_median": lym,
            "reading": ("HOLDS" if holds else "FAILS") + " (the lymphoid estimate does not track methylation; "
                        f"bar {BAR}, GIMiCC level 3)"}


# ------------------------------------------------------------------------------------------- Q4
def _tb(shares: pd.DataFrame) -> pd.Series:
    return (shares["T"] / (shares["T"] + shares["B"]).replace(0, np.nan)).dropna()


def q4(gtr: dict) -> dict:
    from lymphoid_tracking import MIN_PAIRS, marker_index, perm_spearman  # noqa: PLC0415
    rng = np.random.default_rng(config.RANDOM_SEED)
    control = {}
    for c in gt.COHORTS:
        idx, prov = marker_index(c)
        rec = {"marker_panels": {k: v for k, v in prov.items() if k != "source"}}
        for lvl in ("level3", "level4"):
            tb = _tb(gtr[c][lvl]["shares"])
            common = [k for k in idx.index if k in tb.index]
            r, p = perm_spearman(idx.loc[common].to_numpy(), tb.loc[common].to_numpy(), rng)
            rec[lvl] = {"n": len(common), "rho": round(r, 4), "p": round(p, 6),
                        "n_gimicc_samples": gtr[c][lvl]["n_samples"],
                        "n_gimicc_samples_without_T_over_TplusB": gtr[c][lvl]["n_samples"] - len(tb)}
        # the original's exploratory PTPRC-tertile stratification, level 3
        tb = _tb(gtr[c]["level3"]["shares"])
        path = config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz"
        pt = pd.concat([ch.loc[ch.index.intersection(["PTPRC"])]
                        for ch in pd.read_csv(path, index_col=0, chunksize=2000)]).iloc[0]
        pt.index = [k4s(x) for x in pt.index]
        pt = np.log2(pt[~pt.index.duplicated()] + 1)
        cs = [k for k in idx.index if k in tb.index and k in pt.index]
        dd = pd.DataFrame({"idx": idx.loc[cs], "tb": tb.loc[cs], "ptprc": pt.loc[cs]})
        dd["tert"] = pd.qcut(dd["ptprc"], 3, labels=["low", "mid", "high"])
        rec["level3_ptprc_tertiles"] = {}
        for tt, g in dd.groupby("tert", observed=True):
            r_, p_ = perm_spearman(g["idx"].to_numpy(), g["tb"].to_numpy(), rng)
            rec["level3_ptprc_tertiles"][str(tt)] = {"n": int(len(g)), "rho": round(r_, 4), "p": round(p_, 6)}
        r_m, p_m = perm_spearman(dd["ptprc"].to_numpy(), dd["tb"].to_numpy(), rng)
        rec["level3_ptprc_vs_gimicc_TB"] = {"rho": round(r_m, 4), "p": round(p_m, 6)}
        control[c] = rec
    passes = all(control[c]["level3"]["rho"] > 0 and control[c]["level3"]["p"] < 0.05 for c in gt.COHORTS)

    tracking = {}
    for c in gt.COHORTS:
        tb = _tb(gtr[c]["level3"]["shares"])
        for ref in ("frozen", "h5ad"):
            base = ("" if c == "gbm" else "_lgg") + ("_h5ad" if ref == "h5ad" else "")
            for panel, f in (("registered", config.RESULTS_DIR / f"estimates_full{base}.csv"),
                             ("extension", config.RESULTS_DIR / "extension" / f"estimates_full{base}_extension.csv")):
                if not f.exists():
                    continue
                full = pd.read_csv(f)
                scol = next(x for x in ("sample", "sample_id") if x in full.columns)
                full["k"] = full[scol].map(k4s)
                for m, g in full.groupby("method"):
                    g = g.drop_duplicates("k").set_index("k")
                    ids = [k for k in g.index if k in tb.index]
                    e = renormalise(g.loc[ids, COMPARE])
                    etb = (e["T_cell"] / (e["T_cell"] + e["B_cell"]).replace(0, np.nan)).dropna()
                    ids = list(etb.index)
                    key = f"{c}|{ref}|{panel}|{m}"
                    if len(ids) < MIN_PAIRS or etb.nunique() < 3:
                        tracking[key] = {"n": len(ids), "rho": None, "p": None, "why_none": "too few pairs"}
                        continue
                    r, p = perm_spearman(etb.to_numpy(), tb.loc[ids].to_numpy(), rng)
                    tracking[key] = {"n": len(ids), "rho": round(r, 4), "p": round(p, 6)}
    return {"positive_control": control, "positive_control_passes": passes,
            "reading": "PASSES" if passes else "FAILS",
            "method_tracking_of_gimicc_TB_level3": tracking,
            "method_tracking_status": "informative (control passed)" if passes else
                                      "INCONCLUSIVE (control failed; as D23 for EpiDISH)"}


def main() -> int:
    src = config.RESULTS_DIR / "gimicc_truth_confirmation.json"
    if not src.exists():
        print("BLOCKED: run scripts/gimicc_truth.py --stage analyse first"); return 2
    gim = json.loads(src.read_text())
    gtr = {c: gimicc_truths(c) for c in gt.COHORTS}
    out = {"rule": "prespecified/gimicc_truth_confirmation.md, Addendum 3",
           "gimicc_Q1_reading": gim["Q1_reading"],
           "nan_exclusions": {c: {lvl: {k: v for k, v in gtr[c][lvl].items() if k.startswith("n_")}
                                  for lvl in ("level3", "level4")} for c in gt.COHORTS},
           "Q2": q2(gim, gtr), "Q3": q3(gtr), "Q4": q4(gtr)}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    for k, v in out["Q2"]["by_cohort_reference"].items():
        print(f"Q2 {k:12s} GIMiCC T>B {v['gimicc_T_exceeds_B_rna_matched_level4']!s:5s} | EpiDISH T>B "
              f"{v['epidish_T_exceeds_B']!s:5s} | methods T>B {v['n_methods_T_exceeds_B']}/{v['n_methods']} | "
              f"match GIMiCC {v['n_match_gimicc']} | match EpiDISH {v['n_match_epidish']}")
    print("Q2 reading:", out["Q2"]["reading"])
    q = out["Q3"]
    print("Q3 controls:", q["controls"])
    if "units" in q:
        for k, u in q["units"].items():
            print(f"Q3 {k:16s} D1 {u['d1_stability']:+.3f} | E2 truth ({u['e2_truth']}) unmix {u['e2_truth_rho_unmix']:+.3f}"
                  + (f" | GIMiCC L3 unmix {u['gimicc_level3_truth_rho_unmix']:+.3f} med {u['gimicc_level3_truth_rho_median_methods']:+.3f}"
                     f" (n {u['gimicc_level3_n']})" if "gimicc_level3_n" in u else ""))
        print("Q3 H1 (GIMiCC lymphoid truth):", q["H1_gimicc_lymphoid_truth"])
        print("Q3 reading:", q["reading"])
    q = out["Q4"]
    for c, r in q["positive_control"].items():
        print(f"Q4 {c}: marker index vs GIMiCC T/(T+B) level 3 rho {r['level3']['rho']:+.3f} p {r['level3']['p']:.4g} "
              f"(n {r['level3']['n']}) | level 4 rho {r['level4']['rho']:+.3f} | tertiles {r['level3_ptprc_tertiles']}")
    print("Q4 reading:", q["reading"], "|", q["method_tracking_status"])
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
