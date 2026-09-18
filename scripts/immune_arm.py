#!/usr/bin/env python3
"""
immune_arm.py — the second ground truth: does deconvolution over-call immune content on tissue?

Executes `prespecified/immune_failure_factors.md`, committed before any immune error existed.

THE POINT
---------
This project's most-quoted clinical claim is that every method over-calls T cells by 4.6-7.9x.
That comes from SYNTHETIC mixtures built from GBmap's own cells, so it could be a simulator
artefact. Leukocyte fraction from DNA METHYLATION (Thorsson et al. 2018, Immunity 48(4)) tests it
on tissue with a molecule independent of both the RNA the methods consume and the copy number the
tumour arm uses.

THE DENOMINATORS DO NOT MATCH, and that is declared rather than hidden. Leukocyte fraction is a
DNA-derived cell fraction of the whole sample; the estimate is an mRNA proportion of an eight-type
roster (D12). The roster drops 7.0% of the atlas including some leukocytes, and immune cells carry
less mRNA per cell than tumour cells (2.245x measured) -- BOTH push the estimate DOWNWARD relative
to the truth. So the over-call prediction is conservative: a positive error has to overcome them.

    python scripts/immune_arm.py --cohort gbm
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

LF_PATH = ROOT / "Immune_Fraction" / "TCGA_all_leuk_estimate.masked.20170107.tsv"
IMMUNE = ["Macrophage_Microglia", "T_cell", "NK_cell", "B_cell"]
STUDY = {"gbm": "GBM", "lgg": "LGG"}


def k4s(s):
    f = str(s).split("-")
    return "-".join(f[:3] + [f[3][:2]]) if len(f) > 3 else str(s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["gbm", "lgg"], required=True)
    args = ap.parse_args()
    tag = "" if args.cohort == "gbm" else f"_{args.cohort}"

    full_p = config.RESULTS_DIR / f"estimates_full{tag}.csv"
    tum_p = config.RESULTS_DIR / f"absolute_purity_per_sample{tag}.csv"
    if not full_p.exists():
        print(f"BLOCKED: {full_p.name} missing. Re-run absolute_purity_yardstick.py "
              f"--cohort {args.cohort}; it now persists all eight types.")
        return 2

    full = pd.read_csv(full_p)
    full["sample"] = full["sample"].astype(str)
    lf = pd.read_csv(LF_PATH, sep="\t", header=None,
                     names=["study", "barcode", "leukocyte_fraction"])
    lf = lf[lf["study"] == STUDY[args.cohort]].copy()
    lf["k"] = lf["barcode"].map(k4s)
    lf["leukocyte_fraction"] = pd.to_numeric(lf["leukocyte_fraction"], errors="coerce")
    lf = lf.dropna(subset=["leukocyte_fraction"]).drop_duplicates("k").set_index("k")

    full["k"] = full["sample"].map(k4s)
    shared = sorted(set(full["k"]) & set(lf.index))
    print(f"{args.cohort.upper()}: {full['method'].nunique()} methods, "
          f"{full['sample'].nunique()} samples; {len(shared)} join a leukocyte fraction")
    if len(shared) < 50:
        print("BLOCKED: too few joined samples."); return 2
    truth = lf.loc[shared, "leukocyte_fraction"]
    print(f"  leukocyte fraction: median {truth.median():.4f}  "
          f"range {truth.min():.4f}-{truth.max():.4f}")

    # --- the outcome -------------------------------------------------------------------
    rows, per_method = {}, {}
    for m, g in full.groupby("method"):
        g = g.drop_duplicates("k").set_index("k")
        if not set(IMMUNE) <= set(g.columns):
            continue
        idx = [k for k in shared if k in g.index]
        imm = g.loc[idx, IMMUNE].sum(axis=1)
        t = truth.loc[idx]
        ok = np.isfinite(imm) & np.isfinite(t)
        if ok.sum() < 50:
            print(f"  {m:26s} SKIPPED — {int(ok.sum())} finite"); continue
        e = (imm[ok] - t[ok]).to_numpy(dtype="float64")
        rho = float(stats.spearmanr(imm[ok], t[ok]).statistic)
        slope = float(np.polyfit(t[ok].to_numpy(dtype="float64"), e, 1)[0])
        per_method[m] = {
            "n": int(ok.sum()),
            "mean_immune_estimate": round(float(imm[ok].mean()), 4),
            "mean_leukocyte_fraction": round(float(t[ok].mean()), 4),
            "mean_error": round(float(np.mean(e)), 4),
            "median_error": round(float(np.median(e)), 4),
            "fold_over_call": (round(float(imm[ok].mean() / t[ok].mean()), 3)
                               if t[ok].mean() > 0 else None),
            "spearman_vs_truth": round(rho, 4),
            "recovered_fraction": round(1 + slope, 4)}
        rows[m] = per_method[m]
        print(f"  {m:26s} est {imm[ok].mean():.3f} vs LF {t[ok].mean():.3f}  "
              f"err {np.mean(e):+.4f}  {per_method[m]['fold_over_call']}x  "
              f"rho {rho:+.4f}  recovers {per_method[m]['recovered_fraction']:.1%}")

    errs = np.array([v["mean_error"] for v in per_method.values()])
    n_over = int((errs > 0).sum())
    p1 = float(stats.binomtest(n_over, len(errs), 0.5).pvalue)
    print(f"\nP1 — do methods OVER-call immune on tissue?")
    print(f"  {n_over} of {len(errs)} methods have a positive mean error; "
          f"sign-test p = {p1:.5f}")
    print(f"  median across methods: {np.median(errs):+.4f}   "
          f"-> {'SUPPORTED' if n_over >= 0.75 * len(errs) and p1 < 0.05 else 'NOT SUPPORTED'}")

    # --- P2, the asymmetry -------------------------------------------------------------
    p2 = None
    if tum_p.exists():
        tum = pd.read_csv(tum_p, index_col=0)
        pur = tum.pop("absolute_purity")
        tb = {m: float((tum[m] - pur).mean()) for m in tum.columns if m in per_method}
        both = [(m, tb[m], per_method[m]["mean_error"]) for m in tb]
        opp = sum(1 for _, a, b in both if a < 0 < b)
        print(f"\nP2 — are the signs OPPOSITE between arms?")
        print(f"  {opp} of {len(both)} methods: tumour UNDER-called and immune OVER-called")
        for m, a, b in sorted(both, key=lambda x: x[1]):
            print(f"    {m:24s} tumour {a:+.4f}   immune {b:+.4f}"
                  f"{'   <- opposite' if a < 0 < b else ''}")
        p2 = {"n_methods": len(both), "n_opposite_signs": opp,
              "per_method": {m: {"tumour_bias": round(a, 4), "immune_bias": round(b, 4)}
                             for m, a, b in both}}

    out = {"cohort": args.cohort,
           "prespecification": "prespecified/immune_failure_factors.md",
           "ground_truth": "methylation-derived leukocyte fraction, Thorsson et al. 2018 "
                           "Immunity 48(4); 2,000 differentially methylated loci, mixture model",
           "outcome": "sum(Macrophage_Microglia, T_cell, NK_cell, B_cell) - leukocyte_fraction",
           "denominator_caveat": "an mRNA proportion of an eight-type roster against a "
                                 "DNA-derived cell fraction of the whole sample. The roster "
                                 "drops 7.0% of the atlas including some leukocytes, and immune "
                                 "cells carry ~2.245x less mRNA than tumour cells; both bias the "
                                 "estimate DOWNWARD, so an over-call finding is conservative.",
           "n_samples": len(shared),
           "P1_over_call": {"n_methods": len(errs), "n_positive": n_over,
                            "sign_test_p": round(p1, 5),
                            "median_mean_error": round(float(np.median(errs)), 4)},
           "P2_asymmetry": p2, "methods": rows}
    (config.RESULTS_DIR / f"immune_arm{tag}.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote results/immune_arm{tag}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
