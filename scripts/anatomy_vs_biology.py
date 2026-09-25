"""Part 2 of the registered question, tested directly.

THE REGISTERED QUESTION, from `Anatomy_Test.md`, has two clauses:

    "Can a tumor's own anatomy stand in for ground truth when choosing a cell-type
     deconvolution method -- AND DOES A METHOD THAT GETS THE ANATOMY RIGHT ALSO GET THE
     BIOLOGY RIGHT?"

The second clause has been answered only indirectly, through the ACS-versus-purity rank
correlation (rho = +0.081, n = 12, p = 0.80), which is underpowered and therefore reported as
INCONCLUSIVE. That is a weak answer to a sharp question.

This script answers it directly, using the study's own pre-registered biological criterion --
"DNA methylation will show T cells > B cells", fixed in `prespecified/immune_failure_factors.md`
at 2026-09-17 23:37, ninety-five minutes before the LGG measurement existed.

It asks, per method: does a high anatomic concordance score go with getting that ordering right?

WHAT IS AND IS NOT PRE-REGISTERED HERE. The T > B criterion is. Relating it to the ACS ranking
is **not** -- it is a direct reading of the registered question against two artefacts that
already exist, and it is labelled as such in the output. It selects nothing: no method is
chosen, no parameter tuned, no constraint altered. It reports a relationship and its power.

Usage:  python3 scripts/anatomy_vs_biology.py [--write]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                    # noqa: E402

COHORTS = {"GBM": "", "LGG": "_lgg"}


def build() -> dict:
    R = config.RESULTS_DIR
    lb = pd.read_csv(R / "anatomic/acs_leaderboard.csv").set_index("method")
    out: dict = {
        "what_this_is":
            "Part 2 of the registered research question -- 'does a method that gets the "
            "anatomy right also get the biology right?' -- tested directly against the "
            "pre-registered T > B criterion rather than through the underpowered "
            "ACS-versus-purity rank correlation.",
        "registered_status":
            "The T > B prediction IS pre-registered (prespecified/immune_failure_factors.md, "
            "committed 2026-09-17T23:37, 95 minutes before the LGG methylation measurement "
            "existed). Relating it to the ACS ranking is NOT pre-registered; it is a direct "
            "reading of the registered question against two existing artefacts.",
        "selects_nothing":
            "No method is chosen, no parameter tuned, no constraint altered. ACS remains the "
            "object under test, never a selection criterion.",
        "cohorts": {},
    }

    for cohort, tag in COHORTS.items():
        p = R / f"lymphoid_ordering{tag}.json"
        if not p.exists():
            continue
        lo = json.loads(p.read_text())["methods"]
        rows = []
        for m, d in lo.items():
            if m not in lb.index:
                continue
            frac = d.get("frac_samples_est_B_over_T")
            rows.append({
                "method": m,
                "acs": float(lb.loc[m, "acs"]),
                "acs_rank": None,
                "gets_T_over_B_right": bool(d.get("T_exceeds_B")),
                "frac_samples_placing_B_over_T": None if frac is None else float(frac),
                "comparable": bool(d.get("comparable")),
            })
        df = pd.DataFrame(rows).sort_values("acs", ascending=False).reset_index(drop=True)
        df["acs_rank"] = df["acs"].rank(ascending=False, method="min").astype(int)

        top = df.iloc[0]
        sub = df.dropna(subset=["frac_samples_placing_B_over_T"])
        rho = p_val = None
        if len(sub) >= 3:
            rho, p_val = stats.spearmanr(sub["acs"], sub["frac_samples_placing_B_over_T"])

        block = {
            "n_methods": int(len(df)),
            "n_getting_T_over_B_right": int(df["gets_T_over_B_right"].sum()),
            "top_method_by_acs": str(top["method"]),
            "top_method_acs": round(float(top["acs"]), 4),
            "top_method_gets_T_over_B_right": bool(top["gets_T_over_B_right"]),
            "top_method_frac_placing_B_over_T":
                None if pd.isna(top["frac_samples_placing_B_over_T"])
                else round(float(top["frac_samples_placing_B_over_T"]), 4),
            # The correlation is reported WITH its power, because at n = 12 the point
            # estimate is the least interesting thing about it.
            "spearman_acs_vs_fraction_placing_B_over_T":
                None if rho is None else round(float(rho), 4),
            "p_value": None if p_val is None else round(float(p_val), 4),
            "n_in_correlation": int(len(sub)),
            "direction_vs_hypothesis":
                None if rho is None else
                ("OPPOSITE — a higher anatomy score goes with MORE samples placing B above T"
                 if rho > 0 else
                 "consistent with the hypothesis — a higher anatomy score goes with fewer"),
            "per_method": df.to_dict(orient="records"),
        }
        out["cohorts"][cohort] = block

    # The claim that does not depend on the correlation at all.
    cats = [b for b in out["cohorts"].values()]
    out["categorical_finding"] = (
        "In every cohort tested, ZERO methods reproduce the pre-registered T > B ordering, "
        "and the method ranked FIRST by anatomic concordance is among the failures. Anatomic "
        "concordance therefore cannot be used to choose a method that gets the lymphoid "
        "compartment right, because no such method is present at any concordance level. This "
        "claim is categorical and needs no correlation."
        if cats and all(b["n_getting_T_over_B_right"] == 0 for b in cats) else
        "At least one method reproduces T > B; re-read the per-method table before quoting.")
    out["power_caveat"] = (
        "The rank correlation between ACS and lymphoid failure is computed on 12 methods. At "
        "that size it cannot establish a trend, and the point estimate is reported with its "
        "p-value rather than interpreted. What carries weight is the categorical finding.")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="write the artefact")
    args = ap.parse_args()

    res = build()
    for cohort, b in res["cohorts"].items():
        print(f"--- {cohort}")
        print(f"    methods reproducing the registered T > B ordering : "
              f"{b['n_getting_T_over_B_right']} of {b['n_methods']}")
        print(f"    top method by anatomic concordance                : "
              f"{b['top_method_by_acs']} (ACS {b['top_method_acs']})")
        print(f"      ... does it get T > B right?                    : "
              f"{b['top_method_gets_T_over_B_right']}")
        print(f"      ... fraction of samples it places B above T     : "
              f"{b['top_method_frac_placing_B_over_T']}")
        print(f"    Spearman(ACS, fraction placing B over T)          : "
              f"{b['spearman_acs_vs_fraction_placing_B_over_T']}  "
              f"p={b['p_value']}  n={b['n_in_correlation']}")
        print(f"    direction                                          : "
              f"{b['direction_vs_hypothesis']}")
        print()
    print(res["categorical_finding"])
    print()
    print(res["power_caveat"])

    if args.write:
        out = config.RESULTS_DIR / "anatomy_vs_biology.json"
        out.write_text(json.dumps(res, indent=2) + "\n")
        print(f"\nwrote {out.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
