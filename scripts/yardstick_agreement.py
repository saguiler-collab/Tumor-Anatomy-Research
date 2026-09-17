#!/usr/bin/env python3
"""
yardstick_agreement.py — the registered primary outcome, computed against BOTH yardsticks.

THE QUESTION, AND WHY THIS IS THE STUDY'S CLOSING MEASUREMENT
--------------------------------------------------------------
The registration asks: does the ACS ranking of deconvolution methods agree with a ranking from
real ground truth? It reported Spearman **rho = 0.7501**, above the pre-registered bar of 0.60,
against synthetic pseudobulk.

Correction **C4** then observed that BOTH arms of that comparison use GBmap — the ACS arm
deconvolves Ivy GAP against it, and the pseudobulk arm builds its mixtures from its cells — so
the agreement may be two GBmap-based rankings agreeing about GBmap rather than about tissue.
That was an argument. This makes it a measurement.

**ABSOLUTE tumour purity is measured from DNA copy number.** It shares no failure mode with any
RNA deconvolution method under test. Running the same comparison against it asks the registered
question with the shared dependence removed.

WHAT DIFFERS BETWEEN THE TWO ARMS, STATED BECAUSE IT QUALIFIES THE ANSWER
-------------------------------------------------------------------------
This is not a like-for-like re-run and must not be presented as one. The ACS scores come from
Ivy GAP scored against the h5ad-built GBmap reference; the purity correlations come from
TCGA-GBM scored against the vendored frozen signature. **Cohort and reference build both
differ.** What is held constant is the thing being ranked: the same named methods, each scored
once per yardstick — which is the same form the registered outcome itself takes, since its two
arms are also different cohorts.

Degenerate methods are reported and also excluded in a second statistic, because a method that
has silently reduced to another one contributes a duplicate point to a rank correlation.

    python scripts/yardstick_agreement.py
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
from ivygap import config                                          # noqa: E402

YARD = config.RESULTS_DIR / "absolute_purity_yardstick.json"
LB = config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv"
REGISTERED_BAR = 0.60
REGISTERED_RHO = 0.7501


def boot_spearman(x, y, n=10000, seed=config.RANDOM_SEED):
    rng = np.random.default_rng(seed)
    k = len(x)
    out = []
    for _ in range(n):
        i = rng.integers(0, k, k)
        if len(set(i)) < 3:
            continue
        r = stats.spearmanr(x[i], y[i]).statistic
        if np.isfinite(r):
            out.append(r)
    a = np.array(out)
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def main() -> int:
    if not YARD.exists():
        print(f"BLOCKED: {YARD.name} missing. Run scripts/absolute_purity_yardstick.py.")
        return 2
    rep = json.loads(YARD.read_text())
    lb = pd.read_csv(LB).set_index("method")

    rows = []
    for m, v in rep.get("methods", {}).items():
        if "spearman_vs_purity" not in v or v.get("is_control"):
            continue
        if m not in lb.index:
            continue
        rows.append({"method": m, "acs": float(lb.loc[m, "acs"]),
                     "rho_purity": float(v["spearman_vs_purity"]),
                     "bias": float(v.get("bias_tumor_minus_purity", np.nan)),
                     "degenerate": bool(v.get("degenerate", False))})
    d = pd.DataFrame(rows).set_index("method").sort_values("acs", ascending=False)
    if len(d) < 6:
        print(f"BLOCKED: only {len(d)} comparable methods; the registration requires >= 6.")
        return 2
    print(d.round(4).to_string())

    res = {}
    for label, sub in (("all_methods", d), ("excluding_degenerate", d[~d["degenerate"]])):
        if len(sub) < 3:
            continue
        r = stats.spearmanr(sub["acs"], sub["rho_purity"])
        lo, hi = boot_spearman(sub["acs"].to_numpy(), sub["rho_purity"].to_numpy())
        passes = bool(r.statistic >= REGISTERED_BAR and lo > 0)
        res[label] = {"n_methods": int(len(sub)), "spearman": round(float(r.statistic), 4),
                      "p": round(float(r.pvalue), 4), "ci": [round(lo, 4), round(hi, 4)],
                      "meets_registered_bar": passes}
        print(f"\n{label}: n={len(sub)}  rho = {r.statistic:+.4f}  p = {r.pvalue:.4f}  "
              f"95% CI [{lo:+.4f}, {hi:+.4f}]  -> "
              f"{'MEETS' if passes else 'FAILS'} the registered bar (rho >= {REGISTERED_BAR})")

    print(f"\nFor comparison, the SAME test against the GBmap-based pseudobulk yardstick "
          f"reported rho = {REGISTERED_RHO}, which met the bar.")

    out = {
        "what_this_is": "The registered primary outcome computed against a yardstick that "
                        "shares no failure mode with the methods under test.",
        "registered_bar": REGISTERED_BAR,
        "arm_1_pseudobulk": {"rho": REGISTERED_RHO, "yardstick": "synthetic pseudobulk",
                             "shares_with_acs_arm": "GBmap — the mixtures are built from its "
                                                    "cells and the ACS arm deconvolves "
                                                    "against it (correction C4)",
                             "meets_registered_bar": True},
        "arm_2_absolute_purity": {"yardstick": "ABSOLUTE tumour purity from DNA copy number",
                                  "shares_with_acs_arm": "nothing",
                                  **res},
        "what_differs_between_the_arms": (
            "ACS is scored on Ivy GAP against the h5ad-built GBmap reference; the purity "
            "correlations are scored on TCGA-GBM against the vendored frozen signature. "
            "Cohort and reference build both differ. Constant across them is the set of named "
            "methods being ranked, which is the same form the registered outcome takes."),
        "per_method": {m: {k: (None if pd.isna(v) else v) for k, v in r.items()}
                       for m, r in d.iterrows()},
    }
    (config.RESULTS_DIR / "yardstick_agreement.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/yardstick_agreement.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
