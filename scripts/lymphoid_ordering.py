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


def ordering_of(mean_row) -> str | None:
    """Short ordering string, e.g. "T>NK>B", highest first. None if not orderable.

    RETURNS None FOR ANY NaN, and that is the whole point. `sort_values(ascending=False)` puts
    NaN LAST, so a row where every lymphoid value is NaN -- which is what `renormalise` produces
    for a sample with zero T, zero B and zero NK -- sorted to the frame's own COLUMN ORDER and
    came back "T>B>NK". The sample was then counted as "T ranked first" purely because T_cell is
    the first column.

    Found 2026-09-18 by a full-project sweep: TCGA-14-0862-01 has EpiDISH returning exactly 0
    for all three lymphoid types, and it was inflating the GBM T-ranked-first figure from
    0.6323 to 0.6387. One sample in 155, and in the one direction that flatters the paper's
    own prediction, which is the kind of error that has to be impossible rather than small.
    """
    if mean_row.isna().any():
        return None
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
    n_unorderable = int(per_sample.isna().sum())
    if n_unorderable:
        print(f"  NOTE: {n_unorderable} sample(s) have zero signal in all of {COMPARE} and are "
              f"EXCLUDED from\n        the per-sample ordering counts rather than credited to "
              f"whichever column sorts first.")
    per_sample = per_sample.dropna()
    print(f"=== {a.cohort.upper()} / {a.reference} reference ===")
    print(f"\nMETHYLATION truth within {COMPARE}, n={len(ros)}:")
    print(f"  mean  T {mu_t.T_cell:.4f}  NK {mu_t.NK_cell:.4f}  B {mu_t.B_cell:.4f}"
          f"   -> {truth_ord}")
    print("  per-sample ordering frequency:")
    for o, c in per_sample.value_counts().items():
        print(f"    {o:10s} {c:4d}  {c / len(ros):6.1%}")
    t_first = float(per_sample.str.startswith("T>").mean())
    b_first = float(per_sample.str.startswith("B>").mean())
    # denominator is the ORDERABLE samples; n_unorderable is reported separately above
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

    # PER-SAMPLE, not per-method-mean. Comparing two averages over a cohort can hide the
    # case where a method is right on most samples and dragged over by a few. The honest
    # statistic pairs each sample's methylation call against the same sample's estimate.
    print(f"\nPER-SAMPLE AGREEMENT (paired within sample, not a comparison of two means)")
    print(f"{'method':18s} {'n':>5s} {'method says B>T':>16s} {'meth says T>B':>14s} "
          f"{'discordant':>11s} {'no lymph':>9s}")
    print("-" * 80)
    for m, g in full.groupby("method"):
        g = g.drop_duplicates("k").set_index("k")
        idx = [k for k in g.index if k in ros.index]
        if len(idx) < 50:
            continue
        e = renormalise(g.loc[idx])
        # DENOMINATOR. A row where the method assigned NOTHING to T, B and NK renormalises to
        # NaN. Counting those as "did not say B>T" mixes two different statements -- "put B
        # below T" and "found no lymphocytes at all" -- and silently disagrees with the mean
        # above, which skips NaN. They are separated here and counted separately.
        live = e.notna().all(axis=1).to_numpy()
        n_empty = int((~live).sum())
        el = e[live]
        meth_l = ros.loc[idx, :][live]
        est_bt = (el["B_cell"] > el["T_cell"]).to_numpy()
        meth_tb = (meth_l["T_cell"] > meth_l["B_cell"]).to_numpy()
        disc = float((est_bt & meth_tb).mean()) if live.sum() else float("nan")
        res[m]["n_scored"] = int(live.sum())
        res[m]["n_no_lymphoid_signal"] = n_empty
        res[m]["frac_samples_est_B_over_T"] = round(float(est_bt.mean()), 4) if live.sum() else None
        res[m]["frac_samples_discordant"] = round(disc, 4) if live.sum() else None
        print(f"{m:18s} {int(live.sum()):5d} {est_bt.mean():15.1%} {meth_tb.mean():14.1%} "
              f"{disc:10.1%} {n_empty:9d}")
    scored = [v["frac_samples_discordant"] for v in res.values()
              if v.get("frac_samples_discordant") is not None]
    if scored:
        print(f"\nEven the BEST method disagrees with methylation on {min(scored):.1%} of the "
              f"samples it scored.")
    empties = {m: v["n_no_lymphoid_signal"] for m, v in res.items()
               if v.get("n_no_lymphoid_signal")}
    if empties:
        print("Methods returning NO lymphoid signal at all on some samples (excluded above, "
              "not counted as agreement):")
        for m, n0 in sorted(empties.items(), key=lambda kv: -kv[1]):
            print(f"  {m:18s} {n0} sample(s)")

    n_tb = sum(v["T_exceeds_B"] for v in res.values())
    n_full = sum(v["full_ordering_correct"] for v in res.values())
    comp = {m: v for m, v in res.items() if v["comparable"]}
    # DISCLOSE THE DENOMINATOR. "N of 12 methods put B above T" is a statement about MEANS,
    # and for a method that returns exactly zero lymphoid content on most samples that mean
    # rests on a small minority of the cohort. Reporting it without the sample count invites
    # exactly the misreading this project keeps catching elsewhere.
    n_tot = len(res)
    substantive = {m: v for m, v in res.items()
                   if v.get("n_scored", v["n"]) >= 0.5 * v["n"]}
    empty_major = {m: v for m, v in res.items() if m not in substantive}
    print(f"\nTWO DISTINCT FAILURE MODES, counted separately:")
    print(f"\n  ABSENCE — {len(empty_major)} of {n_tot} methods return EXACTLY ZERO T, B and NK "
          f"on\n  the majority of samples. For these the B-vs-T question is vacuous: they are "
          f"not\n  placing B above T, they are reporting no lymphocytes at all.")
    for m, v in sorted(empty_major.items(), key=lambda kv: -kv[1]["n_no_lymphoid_signal"]):
        print(f"    {m:18s} zero lymphoid in {v['n_no_lymphoid_signal']} of {v['n']} "
              f"({v['n_no_lymphoid_signal'] / v['n']:.1%})")
    sub_bt = sum(1 for v in substantive.values()
                 if (v.get("frac_samples_est_B_over_T") or 0) > 0.5)
    print(f"\n  MISASSIGNMENT — of the {len(substantive)} methods that DO return lymphoid signal "
          f"on most\n  samples, {sub_bt} place B above T on a majority of the samples they score.")
    for m, v in sorted(substantive.items(),
                       key=lambda kv: -(kv[1].get("frac_samples_est_B_over_T") or 0)):
        print(f"    {m:18s} B>T on {v['frac_samples_est_B_over_T']:.1%} of {v['n_scored']} "
              f"scored;  discordant with methylation {v['frac_samples_discordant']:.1%}")
    print(f"\nREGISTERED prediction (T > B): {n_tb} of {n_tot} methods agree with methylation "
          f"ON THE MEAN\n  (that mean skips samples with no lymphoid signal -- see above).")
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
        "n_orderable_samples": int(len(per_sample)),
        "n_unorderable_zero_lymphoid": n_unorderable,
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
