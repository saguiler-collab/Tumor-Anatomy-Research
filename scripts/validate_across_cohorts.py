#!/usr/bin/env python3
"""
validate_across_cohorts.py — did the GBM findings reproduce in LGG?

This is the question the whole two-cohort design exists to answer, and it is deliberately NOT
"is the LGG result significant". A factor can be significant in both cohorts and still have
failed to replicate, if it points the other way. So the comparison is made on three things,
in this order:

  1. DIRECTION   — does the sign agree with GBM? A factor that reverses has NOT replicated,
                   whatever its p-value.
  2. CONSISTENCY — is it still consistently signed across the twelve methods?
  3. SIGNIFICANCE— does it survive Holm in LGG on its own?

A factor replicates only if all three hold. Reporting them separately makes a partial
replication legible instead of collapsing it into a yes/no.

The headline quantity — the fraction of true purity variation each method recovers — is
compared directly, because it is a measurement rather than a test and a reader wants the two
numbers side by side.

    python scripts/validate_across_cohorts.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402


def main() -> int:
    g_p = config.RESULTS_DIR / "failure_factors_gbm.json"
    l_p = config.RESULTS_DIR / "failure_factors_lgg.json"
    for p in (g_p, l_p):
        if not p.exists():
            print(f"BLOCKED: {p.name} missing. Run failure_factors.py for that cohort.")
            return 2
    g, l = json.loads(g_p.read_text()), json.loads(l_p.read_text())

    # --- the headline measurement, side by side ---------------------------------------
    gr = g["factors"]["purity"]["interpretable_statistic"]
    lr = l["factors"]["purity"]["interpretable_statistic"]
    print("RECOVERY OF TRUE TUMOUR-CONTENT VARIATION (1 = perfect, 0 = no information)\n")
    print(f"{'method':18s} {'GBM':>9s} {'LGG':>9s} {'change':>9s}")
    print("-" * 49)
    rows = {}
    for m in sorted(gr["per_method"], key=lambda k: -gr["per_method"][k]["recovered_fraction"]):
        if m not in lr["per_method"]:
            continue
        a = gr["per_method"][m]["recovered_fraction"]
        b = lr["per_method"][m]["recovered_fraction"]
        rows[m] = {"gbm": a, "lgg": b, "change": b - a}
        print(f"{m:18s} {a:8.1%} {b:8.1%} {b - a:+8.1%}")
    med_g, med_l = gr["median_recovered_fraction"], lr["median_recovered_fraction"]
    print(f"\n{'median':18s} {med_g:8.1%} {med_l:8.1%} {med_l - med_g:+8.1%}")

    # --- factor-by-factor replication ---------------------------------------------------
    print("\n\nFACTOR REPLICATION\n")
    print(f"{'factor':20s} {'GBM':>18s} {'LGG':>18s}  verdict")
    print("-" * 78)
    rep = {}
    for f, gf in g["factors"].items():
        lf = l["factors"].get(f)
        if lf is None:
            rep[f] = {"verdict": "NOT TESTABLE in LGG",
                      "why": ("Verhaak class does not exist for LGG; carried by the "
                              "secondary MES score" if f == "is_mesenchymal" else
                              "methylation phenotype not in MC3, and 0.868 correlated with "
                              "IDH1 in GBM so not separately informative"
                              if f == "gcimp" else "absent")}
            print(f"{f:20s} {'':>18s} {'':>18s}  NOT TESTABLE — {rep[f]['why'][:40]}")
            continue
        gs = np.sign(gf["median_beta"]); ls_ = np.sign(lf["median_beta"])
        same_dir = bool(gs == ls_ and gs != 0)
        cons = lf["n_in_predicted_direction"] >= 0.75 * lf["n_methods"]
        sig = bool(lf.get("holm_p") is not None and lf["holm_p"] < 0.05)
        if gf.get("supported"):
            verdict = ("REPLICATED" if (same_dir and cons and sig) else
                       "direction agrees, not significant" if same_dir and cons else
                       "direction agrees, inconsistent" if same_dir else
                       "REVERSED — did not replicate")
        else:
            verdict = ("null in GBM, SUPPORTED in LGG" if sig and cons else
                       "null in both")
        rep[f] = {"gbm_beta": gf["median_beta"], "lgg_beta": lf["median_beta"],
                  "gbm_supported": bool(gf.get("supported")),
                  "lgg_holm_p": lf.get("holm_p"),
                  "lgg_signs": f"{lf['n_in_predicted_direction']}/{lf['n_methods']}",
                  "same_direction": same_dir, "verdict": verdict}
        print(f"{f:20s} {gf['median_beta']:+10.5f}{'*' if gf.get('supported') else ' ':>2s}"
              f"{lf['n_in_predicted_direction']:3d}/{lf['n_methods']:<3d}"
              f"{lf['median_beta']:+10.5f} {str(lf.get('holm_p')):>7s}  {verdict}")

    out = {"what_this_is": "Did the GBM discovery reproduce in an independent glioma cohort?",
           "replication_rule": "direction first, then consistency across methods, then "
                               "significance. A factor that reverses has not replicated "
                               "whatever its p-value.",
           "n_gbm": g["n_samples_complete_cases"], "n_lgg": l["n_samples_complete_cases"],
           "recovery_fraction": {"gbm_median": med_g, "lgg_median": med_l,
                                 "per_method": rows},
           "factors": rep}
    (config.RESULTS_DIR / "cross_cohort_validation.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/cross_cohort_validation.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
