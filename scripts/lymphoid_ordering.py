"""Do methods reproduce the lymphoid ordering that DNA methylation measures?

WHAT IS PRE-SPECIFIED AND WHAT IS NOT.
`prespecified/immune_failure_factors.md` (per-cell-type addendum) registered ONE directional
prediction before measurement: **T > B**, with the falsifier "methylation showing B >= T"
stated in advance. That prediction, and only that, is confirmatory.

The three-way ordering including NK reported here is **SECONDARY and EXPLORATORY**. It was not
registered, no direction was fixed for it, and it is labelled as such in every artefact this
writes. It is reported because the same measurement yields it at no extra cost, not because it
was predicted.

WHY THE COMPARISON IS LIKE-FOR-LIKE. EpiDISH returns a leukocyte SUB-composition; a
deconvolution method returns a fraction of all cells. Both sides are renormalised within
{T_cell, B_cell, NK_cell} before anything is compared, so the differing denominators cancel by
construction and only the relative lymphoid composition is compared.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.deconv.comparability import non_comparable

COMPARE = ["T_cell", "B_cell", "NK_cell"]
SHORT = {"T_cell": "T", "B_cell": "B", "NK_cell": "NK"}


def k4s(s: object) -> str:
    """TCGA barcode to patient+sample-type key, dropping the aliquot letter."""
    f = str(s).split("-")
    return "-".join(f[:3] + [f[3][:2]]) if len(f) > 3 else str(s)


def renormalise(frame: pd.DataFrame, types=COMPARE) -> pd.DataFrame:
    """Relative composition WITHIN `types`, row-wise.

    This is the step that makes a methylation leukocyte sub-composition and a
    deconvolution all-cell fraction comparable: after it, both are the same quantity --
    relative composition among lymphocytes. Anything outside `types` (Tumor, Macrophage,
    the methylation reference's Mono/Neutro/Eosino) is dropped BEFORE the division, so the
    differing denominators cancel rather than being assumed away.
    """
    sub = frame.loc[:, list(types)]
    return sub.div(sub.sum(axis=1).replace(0, np.nan), axis=0)


def ordering_of(mean_row) -> str:
    """Short ordering string, e.g. "T>NK>B", highest first."""
    return ">".join(SHORT[c] for c in mean_row.sort_values(ascending=False).index)


def _truth(meth_csv) -> pd.DataFrame:
    meth = pd.read_csv(meth_csv, index_col=0)
    ros = renormalise(pd.DataFrame({"T_cell": meth[["CD4T", "CD8T"]].sum(axis=1),
                                    "B_cell": meth["B"], "NK_cell": meth["NK"]}))
    ros.index = [k4s(i) for i in ros.index]
    return ros[~ros.index.duplicated()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="lgg", choices=["gbm", "lgg"])
    ap.add_argument("--reference", default="frozen", choices=["frozen", "h5ad"])
    a = ap.parse_args()
    base = "" if a.cohort == "gbm" else f"_{a.cohort}"
    tag = base + ("_h5ad" if a.reference == "h5ad" else "")

    meth_csv = config.RESULTS_DIR / f"methylation_celltypes{base}.csv"
    est_csv = config.RESULTS_DIR / f"estimates_full{tag}.csv"
    if not meth_csv.exists():
        print(f"BLOCKED: no methylation truth for {a.cohort} ({meth_csv.name})"); return 2
    if not est_csv.exists():
        print(f"BLOCKED: no estimates for {a.cohort}/{a.reference} ({est_csv.name})"); return 2

    ros = _truth(meth_csv)
    mu_t = ros.mean()
    truth_ord = ordering_of(mu_t)
    per_sample = ros.apply(ordering_of, axis=1)
    print(f"=== {a.cohort.upper()} / {a.reference} reference ===")
    print(f"\nMETHYLATION truth within {COMPARE}, n={len(ros)}:")
    print(f"  mean  T {mu_t.T_cell:.4f}  NK {mu_t.NK_cell:.4f}  B {mu_t.B_cell:.4f}"
          f"   -> {truth_ord}")
    print("  per-sample ordering frequency:")
    for o, c in per_sample.value_counts().items():
        print(f"    {o:10s} {c:4d}  {c / len(ros):6.1%}")
    t_first = float(per_sample.str.startswith("T>").mean())
    b_first = float(per_sample.str.startswith("B>").mean())
    print(f"  T ranked first in {t_first:.1%} of samples; B ranked first in {b_first:.1%}")

    full = pd.read_csv(est_csv)
    scol = next(c for c in ("sample", "sample_id") if c in full.columns)
    full["k"] = full[scol].map(k4s)
    nc = non_comparable(a.reference)

    print(f"\n{'method':18s} {'T':>7s} {'NK':>7s} {'B':>7s}  {'ordering':10s} "
          f"{'T>B (registered)':>17s}  {'full order':>10s}")
    print("-" * 82)
    res: dict[str, dict] = {}
    for m, g in full.groupby("method"):
        g = g.drop_duplicates("k").set_index("k")
        idx = [k for k in g.index if k in ros.index]
        if len(idx) < 50:
            continue
        mu = renormalise(g.loc[idx]).mean()
        o = ordering_of(mu)
        tb = bool(mu.T_cell > mu.B_cell)
        res[m] = {"n": len(idx), "T": round(float(mu.T_cell), 4),
                  "NK": round(float(mu.NK_cell), 4), "B": round(float(mu.B_cell), 4),
                  "ordering": o, "T_exceeds_B": tb, "full_ordering_correct": o == truth_ord,
                  "comparable": m not in nc}
        print(f"{m:18s} {mu.T_cell:7.4f} {mu.NK_cell:7.4f} {mu.B_cell:7.4f}  {o:10s} "
              f"{str(tb):>17s}  {str(o == truth_ord):>10s}"
              f"{'  [not comparable]' if m in nc else ''}")

    n_tb = sum(v["T_exceeds_B"] for v in res.values())
    n_full = sum(v["full_ordering_correct"] for v in res.values())
    comp = {m: v for m, v in res.items() if v["comparable"]}
    print(f"\nREGISTERED prediction (T > B): {n_tb} of {len(res)} methods agree with "
          f"methylation.")
    print(f"SECONDARY, not registered (full {truth_ord}): {n_full} of {len(res)} methods "
          f"reproduce it.")
    if comp:
        print(f"Restricted to the {len(comp)} methods comparable on this reference: "
              f"{sum(v['T_exceeds_B'] for v in comp.values())} agree on T > B.")

    out = config.RESULTS_DIR / f"lymphoid_ordering{tag}.json"
    out.write_text(json.dumps({
        "cohort": a.cohort, "reference": a.reference,
        "truth_mean": {k: round(float(v), 4) for k, v in mu_t.items()},
        "truth_ordering": truth_ord,
        "truth_T_first_fraction": round(t_first, 4),
        "truth_B_first_fraction": round(b_first, 4),
        "n_methylation_samples": int(len(ros)),
        "methods": res, "n_agree_T_over_B": int(n_tb),
        "n_reproduce_full_ordering": int(n_full), "n_methods": len(res),
        "registered": "Only T > B was pre-specified "
                      "(prespecified/immune_failure_factors.md, per-cell-type addendum).",
        "EXPLORATORY": "The three-way ordering including NK was NOT registered and is "
                       "secondary.",
    }, indent=2))
    print(f"\nwrote results/lymphoid_ordering{tag}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
