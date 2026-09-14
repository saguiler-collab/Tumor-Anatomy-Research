#!/usr/bin/env python3
"""
uncertainty_conformal.py — calibrated prediction intervals for EVERY method, no re-running.

THE CLINICAL GAP THIS CLOSES
----------------------------
`results/method_completeness.json` finds 0 of 15 methods complete, and one criterion fails
for all 15: **no per-sample uncertainty**. One method (`bayesian`) reports intervals and they
are severely overconfident — 9.25% coverage against a nominal 95%.

For a clinical purpose that is the single most disqualifying property in the panel. A number
without an honest error bar cannot support a decision, and a number with a dishonest one is
worse than no number at all: 9.25% coverage means an interval a clinician would read as
"almost certainly inside" contains the truth about one time in eleven.

WHY THIS CAN BE DONE WITHOUT RE-SOLVING ANYTHING
------------------------------------------------
The donor-held-out pseudobulk benchmark already produced 500 mixtures with **known truth**
and every method's predictions on them. Split-conformal prediction turns exactly that into
intervals with a coverage guarantee: take the absolute residuals on a calibration split, take
their (1-alpha) quantile `q`, and the interval `estimate +/- q` covers the truth with
probability at least 1-alpha on exchangeable data. No model, no distributional assumption, no
re-solving, and it works identically for an R package, a Python reimplementation and a
control.

It also inverts the usual relationship in a way that suits this study: the interval width IS
the method's measured reliability. A method needing +/-0.25 to reach 95% coverage on the
tumour compartment is telling you it cannot support a decision that turns on 10 percentage
points of purity, whatever its rank.

THE TWO LIMITATIONS, BOTH REAL, NEITHER HIDDEN
----------------------------------------------
1. **Exchangeability with real tissue is an assumption, not a result.** Calibration is on
   pseudobulk mixtures; the target is laser-captured Ivy GAP tissue. The intervals are
   honest error bars for the pseudobulk regime and transfer to real tissue only insofar as
   the error distributions match. They are labelled that way everywhere.
2. **Mixtures can share donors.** 500 mixtures are drawn from 22 held-out donors, and the
   benchmark does not persist which donor built which mixture, so a random calibration split
   is not donor-disjoint and measured coverage may be optimistic. The fix is to persist the
   per-mixture donor labels, which `run_benchmark` now does; until a run writes them this
   script reports the random-split coverage and says it is optimistic.

    python scripts/uncertainty_conformal.py
    python scripts/uncertainty_conformal.py --alpha 0.10    # 90% intervals
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PRIMARY = ["Tumor", "Macrophage_Microglia", "Endothelial", "Oligodendrocyte"]


def conformal_q(resid: np.ndarray, alpha: float) -> float:
    """
    Split-conformal quantile with the finite-sample correction.

    The level is ceil((n+1)(1-alpha))/n rather than the plain (1-alpha) empirical quantile.
    Without the correction the guarantee does not hold at small n, and n here is 250 per
    cell type after the split -- small enough for it to matter.
    """
    r = np.asarray(resid, dtype="float64")
    r = r[np.isfinite(r)]
    n = len(r)
    if n == 0:
        return float("nan")
    k = int(np.ceil((n + 1) * (1.0 - alpha)))
    if k > n:
        return float("inf")            # n too small to certify this level
    return float(np.sort(r)[k - 1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--alpha", type=float, default=0.05, help="1-alpha coverage; 0.05 -> 95%%")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/uncertainty_conformal.json")
    args = ap.parse_args()

    R = ROOT / args.results_dir
    est_dir = R / "estimates"
    truth = pd.read_csv(est_dir / "pseudobulk_truth.csv", index_col=0)

    donor_path = R / "benchmark" / "mixture_donors.csv"
    donors = None
    if donor_path.exists():
        donors = pd.read_csv(donor_path, index_col=0).iloc[:, 0].reindex(truth.index)
        print(f"donor-grouped split: {donors.nunique()} donors across {len(truth)} mixtures")
    else:
        print("per-mixture donor labels are NOT on disk -- using a random split. Coverage "
              "may be OPTIMISTIC because mixtures can share donors. run_benchmark now "
              "persists them; re-run this after the next full run.")

    rng = np.random.default_rng(args.seed)
    if donors is not None and donors.notna().all():
        uniq = np.array(sorted(donors.unique()))
        rng.shuffle(uniq)
        cal_donors = set(uniq[: len(uniq) // 2])
        cal_mask = donors.isin(cal_donors).to_numpy()
        split_kind = "donor-disjoint"
    else:
        cal_mask = np.zeros(len(truth), dtype=bool)
        cal_mask[rng.choice(len(truth), len(truth) // 2, replace=False)] = True
        split_kind = "random (mixtures may share donors)"

    rows = []
    for f in sorted(est_dir.glob("pseudobulk_*.csv")):
        m = f.stem.replace("pseudobulk_", "")
        if m == "truth":
            continue
        pred = pd.read_csv(f, index_col=0).reindex(index=truth.index, columns=truth.columns)
        resid = (pred - truth).abs()

        per_type, widths, covs = {}, [], []
        for ct in truth.columns:
            r = resid[ct].to_numpy()
            if not np.isfinite(r).any():
                per_type[ct] = {"q": None, "coverage": None,
                                "note": "method does not model this type"}
                continue
            q = conformal_q(r[cal_mask], args.alpha)
            held = r[~cal_mask]
            held = held[np.isfinite(held)]
            cov = float(np.mean(held <= q)) if len(held) and np.isfinite(q) else None
            q_full = conformal_q(r, args.alpha)
            per_type[ct] = {"q": None if not np.isfinite(q) else round(q, 4),
                            "q_full_data": None if not np.isfinite(q_full) else round(q_full, 4),
                            "coverage_on_heldout": None if cov is None else round(cov, 4),
                            "n_calibration": int(np.isfinite(r[cal_mask]).sum())}
            if np.isfinite(q):
                widths.append(q)
            if cov is not None:
                covs.append(cov)

        tum = per_type.get("Tumor", {})
        rows.append({
            "method": m,
            "mean_half_width": round(float(np.mean(widths)), 4) if widths else None,
            "tumour_half_width": tum.get("q"),
            "mean_coverage": round(float(np.mean(covs)), 4) if covs else None,
            "n_types_certified": len(widths),
            "per_type": per_type,
        })

    rows.sort(key=lambda r: (r["tumour_half_width"] is None, r["tumour_half_width"] or 9))
    print(f"\nsplit: {split_kind};  target coverage {1 - args.alpha:.0%}\n")
    print(f"{'method':28s} {'tumour +/-':>11s} {'mean +/-':>9s} {'coverage':>9s} {'types':>6s}")
    print("-" * 68)
    for r in rows:
        t = "n/a" if r["tumour_half_width"] is None else f"{r['tumour_half_width']:.4f}"
        w = "n/a" if r["mean_half_width"] is None else f"{r['mean_half_width']:.4f}"
        c = "n/a" if r["mean_coverage"] is None else f"{r['mean_coverage']:.3f}"
        print(f"{r['method']:28s} {t:>11s} {w:>9s} {c:>9s} {r['n_types_certified']:6d}")

    ok = [r for r in rows if r["mean_coverage"] is not None
          and r["mean_coverage"] >= (1 - args.alpha) - 0.03]
    report = {
        "what_this_is": ("Split-conformal prediction intervals for every method, calibrated "
                         "on the donor-held-out pseudobulk benchmark's known truth. Closes "
                         "the 'no per-sample uncertainty' gap that failed all 15 methods in "
                         "method_completeness.json."),
        "how_to_read": ("`tumour +/- q` is the half-width needed to cover the truth at the "
                        "target rate on the tumour compartment. It IS the method's measured "
                        "reliability: a method needing +/-0.25 cannot support a decision "
                        "turning on 10 points of purity, whatever its rank."),
        "limitations": [
            "Calibrated on pseudobulk mixtures; the target is laser-captured tissue. "
            "Transfer is an assumption, not a result.",
            f"Split is {split_kind}. A random split is not donor-disjoint, so coverage may "
            "be optimistic; run_benchmark now persists per-mixture donors so the next run "
            "enables the donor-grouped version.",
            "Intervals are marginal per cell type, not joint across the roster.",
        ],
        "target_coverage": 1 - args.alpha,
        "split": split_kind,
        "n_mixtures": int(len(truth)),
        "n_methods_achieving_target": len(ok),
        "methods": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\n{len(ok)} of {len(rows)} methods achieve the target coverage on held-out data")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
