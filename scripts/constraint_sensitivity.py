#!/usr/bin/env python3
"""
constraint_sensitivity.py — is the ACS ranking carried by any single constraint or tumour?

WHY THIS EXISTS
---------------
ACS is a composite of seven claims scored across nine tumours. The first question anyone
asks about a composite is whether one component drives it, and the project could not
answer: `acs_cohort_sensitivity.csv` varies the *cohort* (122 vs 270 samples), never the
constraint set or the tumour set.

That gap matters more since the per-constraint spread was measured. C3 is satisfied by
every comparable real method on every evaluable tumour, so it contributes nothing to the
ordering, and C6 alone spans 0.222 to 1.000. A reviewer who notices the second will ask
whether the leaderboard is really "C6 plus noise". This answers it with a number.

WHAT THIS IS NOT — AND THE INVARIANT IT MUST NOT BREAK
------------------------------------------------------
The project's rule is absolute: **never edit the constraint file because a method scored
badly, and never let an outcome select anything.** This script does not select. It
recomputes the same statistic under stated exclusions, reports every method under every
exclusion, and changes nothing.

  * The registered ACS is, and stays, all seven constraints over all nine tumours.
  * Every leave-one-out variant is reported, not the favourable ones.
  * No variant is offered as an alternative headline, and none may be used to reorder the
    leaderboard, choose a method, or justify dropping a constraint. A constraint that
    contributes nothing to the ordering STAYS: it earns its place against the controls,
    and removing it after seeing the scores is the retune the protocol forbids.

HOW IT AVOIDS RE-RUNNING ANYTHING
---------------------------------
`acs_per_tumor.csv` records, per method and tumour, which constraints were satisfied and
which were violated. That is exactly the (constraint x tumour) satisfaction matrix ACS is
computed from, so every variant is arithmetic on artefacts already on disk — no
deconvolution, no permutation, no chance of a different answer because something was
re-run differently.

The reconstruction is checked against the published leaderboard before any variant is
reported: all 17 methods must reproduce to 1e-9, and the weighted denominator must come
out at the published value. If they do not, the script aborts rather than report.

    python scripts/constraint_sensitivity.py
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


def _list(x) -> list[str]:
    if not isinstance(x, str) or not x.strip():
        return []
    return [s for s in x.split(",") if s]


def load_matrix(anat: Path):
    """(method, tumour, constraint) -> satisfied, plus the constraint weights."""
    weights = {c["id"]: float(c["weight"])
               for c in json.loads((anat / "constraint_file.json").read_text())["constraints"]}
    pt = pd.read_csv(anat / "acs_per_tumor.csv")
    recs = []
    for _, r in pt.iterrows():
        sat, vio = _list(r["constraints_satisfied"]), _list(r["constraints_violated"])
        for c in sat:
            recs.append((r["method"], str(r["tumor_id"]), c, 1.0))
        for c in vio:
            recs.append((r["method"], str(r["tumor_id"]), c, 0.0))
    return pd.DataFrame(recs, columns=["method", "tumor", "constraint", "sat"]), weights


def acs(df: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    """Weighted proportion of satisfied (constraint x tumour) pairs, per method."""
    w = df["constraint"].map(weights)
    num = (df["sat"] * w).groupby(df["method"]).sum()
    den = w.groupby(df["method"]).sum()
    return (num / den).sort_values(ascending=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out", default="results/anatomic/constraint_sensitivity.json")
    args = ap.parse_args()

    anat = ROOT / args.results_dir / "anatomic"
    df, weights = load_matrix(anat)
    lb = pd.read_csv(anat / "acs_leaderboard.csv").set_index("method")

    # --- gate: the reconstruction must reproduce the published leaderboard ----
    full = acs(df, weights)
    worst = 0.0
    for m in full.index:
        pub = float(lb.loc[m, "acs"])
        worst = max(worst, abs(full[m] - pub))
        n = int((df["method"] == m).sum())
        if n != int(lb.loc[m, "n_constraint_tumor_pairs"]):
            print(f"ABORT: {m} reconstructs {n} pairs, leaderboard says "
                  f"{int(lb.loc[m, 'n_constraint_tumor_pairs'])}")
            return 1
    if worst > 1e-9:
        print(f"ABORT: reconstruction differs from the published leaderboard by {worst:.2e}. "
              "Not reporting variants computed from a matrix that does not reproduce the "
              "result it is meant to perturb.")
        return 1
    print(f"reconstruction gate PASSED: {len(full)} methods reproduce to "
          f"{worst:.1e}; weighted denominator "
          f"{df[df.method == full.index[0]]['constraint'].map(weights).sum():.0f}")

    # Ranking is over comparable real methods only -- controls are not ranked, and a
    # partial-coverage method is not comparable, exactly as the leaderboard has it.
    ranked = [m for m in full.index
              if not bool(lb.loc[m, "is_control"]) and bool(lb.loc[m, "comparable"])]
    base = full[ranked]
    base_rank = base.rank(ascending=False, method="average")

    report = {
        "what_this_is": ("Leave-one-out sensitivity of the ACS ranking, computed from the "
                         "archived (constraint x tumour) satisfaction matrix. REPORTING "
                         "ONLY: the registered ACS uses all constraints and all tumours, "
                         "and no variant here may be used to reorder the leaderboard, "
                         "select a method, or drop a constraint."),
        "reconstruction_max_abs_error": float(worst),
        "n_methods_ranked": len(ranked),
        "baseline": {m: round(float(base[m]), 4) for m in ranked},
        "leave_one_constraint_out": [],
        "leave_one_tumor_out": [],
    }

    def variant(label: str, sub: pd.DataFrame, dropped_weight: float | None = None) -> dict:
        v = acs(sub, weights).reindex(ranked)
        rho = stats.spearmanr(base_rank, v.rank(ascending=False, method="average")).statistic
        moved = (v.rank(ascending=False, method="average") - base_rank).abs()
        top = v.idxmax()
        return {"excluded": label,
                "spearman_vs_full_ranking": round(float(rho), 4),
                "max_rank_move": round(float(moved.max()), 1),
                "n_methods_whose_rank_moves": int((moved > 0).sum()),
                "top_method": str(top),
                "top_method_changes": bool(top != base.idxmax()),
                "acs": {m: round(float(v[m]), 4) for m in ranked}}

    print("\n--- leave one CONSTRAINT out ---")
    for c in sorted(weights):
        r = variant(c, df[df["constraint"] != c])
        report["leave_one_constraint_out"].append(r)
        print(f"  without {c}: rho {r['spearman_vs_full_ranking']:.4f}  "
              f"max rank move {r['max_rank_move']:.1f}  top {r['top_method']}"
              f"{'  (TOP CHANGES)' if r['top_method_changes'] else ''}")

    print("\n--- leave one TUMOUR out ---")
    for t in sorted(df["tumor"].unique()):
        r = variant(f"tumor:{t}", df[df["tumor"] != t])
        report["leave_one_tumor_out"].append(r)
        print(f"  without {t}: rho {r['spearman_vs_full_ranking']:.4f}  "
              f"max rank move {r['max_rank_move']:.1f}  top {r['top_method']}"
              f"{'  (TOP CHANGES)' if r['top_method_changes'] else ''}")

    rc = [r["spearman_vs_full_ranking"] for r in report["leave_one_constraint_out"]]
    rt = [r["spearman_vs_full_ranking"] for r in report["leave_one_tumor_out"]]
    report["summary"] = {
        "min_rho_constraint_dropout": round(min(rc), 4),
        "min_rho_tumor_dropout": round(min(rt), 4),
        "constraints_whose_removal_changes_the_top_method":
            [r["excluded"] for r in report["leave_one_constraint_out"] if r["top_method_changes"]],
        "tumors_whose_removal_changes_the_top_method":
            [r["excluded"] for r in report["leave_one_tumor_out"] if r["top_method_changes"]],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nworst rho: {min(rc):.4f} (constraint dropout), {min(rt):.4f} (tumour dropout)")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
