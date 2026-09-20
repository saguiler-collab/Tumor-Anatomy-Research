"""Did rebuilding the reference from raw counts change the ACS leaderboard?

OPEN_DEFECTS D16. Every archived result, including the ACS leaderboard and the registered
primary outcome, was computed against a reference built from GBmap's `X` -- which is
`log1p(counts * size factor)`, not linear expression. `run_all.py --matrix raw/X` rebuilds the
same atlas from the genuine counts.

This compares the two. It SELECTS NOTHING: both leaderboards are reported, the rank change per
method is shown, and no method, parameter or reference is chosen on the basis of the result.
The comparison exists because a ranking that moves when the same atlas is merely supplied in a
different space is measuring the harness rather than the methods -- which is the claim
`docs/ENDPOINT.md` §2b already makes about the atlas, tested here one level deeper.

    python scripts/compare_matrix_arms.py --before <dir> [--after results]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ivygap import config


def load(d: Path) -> pd.DataFrame | None:
    p = d / "anatomic" / "acs_leaderboard.csv"
    return pd.read_csv(p).set_index("method") if p.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True, help="results tree built from matrix X (log)")
    ap.add_argument("--after", default=str(config.RESULTS_DIR),
                    help="results tree built from raw/X (counts)")
    a = ap.parse_args()
    b, f = load(Path(a.before)), load(Path(a.after))
    if b is None or f is None:
        print(f"BLOCKED: need acs_leaderboard.csv in both trees "
              f"(before={b is not None}, after={f is not None})")
        return 2

    shared = [m for m in b.index if m in f.index]
    print(f"=== ACS: matrix X (log) vs raw/X (counts) — {len(shared)} methods in both ===\n")
    print(f"{'method':26s} {'X (log)':>9s} {'raw/X':>9s} {'delta':>8s}  "
          f"{'rank X':>7s} {'rank raw':>9s}  moved")
    print("-" * 82)
    rows = {}
    bb = b.loc[shared, "acs"]
    ff = f.loc[shared, "acs"]
    rb, rf = bb.rank(ascending=False), ff.rank(ascending=False)
    for m in sorted(shared, key=lambda x: -ff[x]):
        mv = int(rb[m]) - int(rf[m])
        ctrl = " [control]" if bool(b.loc[m].get("is_control", False)) else ""
        print(f"{m:26s} {bb[m]:9.4f} {ff[m]:9.4f} {ff[m] - bb[m]:+8.4f}  "
              f"{int(rb[m]):7d} {int(rf[m]):9d}  {'same' if mv == 0 else f'{mv:+d}'}{ctrl}")
        rows[m] = {"acs_log": round(float(bb[m]), 4), "acs_counts": round(float(ff[m]), 4),
                   "delta": round(float(ff[m] - bb[m]), 4),
                   "rank_log": int(rb[m]), "rank_counts": int(rf[m]),
                   "is_control": bool(b.loc[m].get("is_control", False))}

    tau = float(rb.corr(rf, method="kendall"))
    rho = float(stats.spearmanr(bb, ff).statistic)
    print(f"\nKendall tau between the two orderings : {tau:+.4f}")
    print(f"Spearman between the two ACS vectors  : {rho:+.4f}")
    print(f"median |delta|                        : {float(np.median(np.abs(ff - bb))):.4f}")
    print(f"largest move                          : "
          f"{max(rows, key=lambda m: abs(rows[m]['delta']))} "
          f"({max(abs(r['delta']) for r in rows.values()):+.4f})")

    # The control separation is what makes ACS a detector at all; it has to survive.
    ctl = [m for m in rows if rows[m]["is_control"]]
    real = [m for m in rows if not rows[m]["is_control"] and rows[m]["acs_counts"] > 0]
    if ctl and real:
        print(f"\nCONTROL SEPARATION (the property ACS must keep to be a detector)")
        print(f"  worst real method, counts : {min(rows[m]['acs_counts'] for m in real):.4f}")
        print(f"  best control, counts      : {max(rows[m]['acs_counts'] for m in ctl):.4f}")
        print(f"  still separated           : "
              f"{min(rows[m]['acs_counts'] for m in real) > max(rows[m]['acs_counts'] for m in ctl)}")

    # THE CONFOUND THAT CAUGHT THIS PROJECT ONCE ALREADY.
    #
    # In the equal-footing comparison, one method (`scdc`) switched from a Python
    # reimplementation to the R package between the two arms, and a rank change caused by
    # swapping implementations says nothing about the variable under test. Budgets are
    # unchanged here, but a method sitting near its timeout can fall back on one run and not
    # the other purely on machine load -- OPEN_DEFECTS D19. So the implementations are
    # compared explicitly rather than assumed identical.
    impl = {}
    if "implementation" in b.columns and "implementation" in f.columns:
        for m in shared:
            ib, if_ = str(b.loc[m, "implementation"]), str(f.loc[m, "implementation"])
            impl[m] = {"log": ib, "counts": if_, "same": ib == if_}
        switched = [m for m, v in impl.items() if not v["same"]]
        print(f"\nIMPLEMENTATION CHECK — did the software change under the method's name?")
        if not switched:
            print(f"  none of the {len(shared)} methods changed implementation. "
                  f"The comparison is clean.")
        else:
            for m in switched:
                print(f"  {m:24s} {impl[m]['log']:26s} -> {impl[m]['counts']}")
            keep = [m for m in shared if impl[m]["same"]]
            if len(keep) >= 3:
                tc = float(rb[keep].corr(rf[keep], method="kendall"))
                print(f"  {len(switched)} changed. Excluding them, tau on the remaining "
                      f"{len(keep)} is {tc:+.4f} (vs {tau:+.4f} overall).")
                print("  => read the EXCLUDED figure, not the overall one: a rank change from "
                      "swapping\n     a reimplementation for a package is not evidence about "
                      "the matrix.")

    out = config.RESULTS_DIR / "matrix_arm_comparison.json"
    out.write_text(json.dumps({
        "what_this_is": "ACS leaderboard built from GBmap X (log1p) vs raw/X (counts). "
                        "OPEN_DEFECTS D16. Selects nothing.",
        "kendall_tau": round(tau, 4), "spearman": round(rho, 4),
        "n_methods": len(shared), "methods": rows,
        "implementation_check": impl}, indent=2))
    print(f"\nwrote {out.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
