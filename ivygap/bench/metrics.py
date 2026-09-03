"""
metrics.py — how a composition estimate is scored against known truth.

EVERY METRIC HERE IS PATIENT-EQUAL WHERE A PATIENT EXISTS
---------------------------------------------------------
On pseudobulk the grouping variable is the reference donor the mixture was pooled
from; on real Ivy GAP it is the tumour. Either way, error is averaged within group
first and across groups second, so a donor or tumour that contributed more samples
does not count for more. `ivygap.bench.equal_footing` explains why this is condition 7
and not a stylistic preference.

WHY THESE METRICS AND NOT JUST CORRELATION
------------------------------------------
Correlation is the most-quoted deconvolution metric and the least informative one. A
method that reports every tumour's macrophage content as exactly twice the truth
correlates at r = 1.0 and is wrong about every sample. So correlation is reported
*alongside* absolute error and bias, never instead of them, and the composite score
used for method selection is built on error rather than on correlation.

Two failure modes get their own metrics because they matter specifically here:

  structural-zero rate   When a cell type is genuinely absent, does the method say
                         zero, or does it sprinkle a few percent there? Ivy GAP's
                         leading-edge samples really do lack lymphocytes, and a method
                         that cannot represent absence will manufacture immune
                         infiltrate at the tumour margin.
  limit of detection     The smallest true fraction at which a type is still recovered
                         above noise. This is the number that decides whether a 1-2%
                         T-cell population is measurable at all.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ivygap import config
from .equal_footing import patient_equal_mean

#: Below this true fraction a cell type counts as structurally absent.
ABSENT_THRESHOLD = 1e-6


def _aligned(estimate: pd.DataFrame, truth: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    shared = [s for s in truth.index if s in estimate.index]
    if not shared:
        raise ValueError("estimate and truth share no samples")
    types = [c for c in truth.columns if c in estimate.columns]
    return estimate.loc[shared, types], truth.loc[shared, types]


def per_type_metrics(estimate: pd.DataFrame, truth: pd.DataFrame,
                     groups: pd.Series | None = None) -> pd.DataFrame:
    """
    One row per cell type: error, bias, correlation, and the two absence metrics.

    `groups` is the donor or patient each sample belongs to. Supplying it switches
    error aggregation to group-equal; omitting it falls back to a plain sample mean and
    the returned frame records which was used, so a table can never be misread as
    patient-equal when it is not.
    """
    est, tru = _aligned(estimate, truth)
    grp = groups.reindex(est.index).astype(str) if groups is not None else None

    rows = []
    for ctype in tru.columns:
        e, t = est[ctype], tru[ctype]
        ok = np.isfinite(e) & np.isfinite(t)
        if ok.sum() == 0:
            rows.append({"cell_type": ctype, "n": 0})
            continue
        e, t = e[ok], t[ok]
        g = grp[ok] if grp is not None else None

        abs_err = (e - t).abs()
        sq_err = (e - t) ** 2
        signed = e - t

        mae = patient_equal_mean(abs_err, g) if g is not None else float(abs_err.mean())
        mse = patient_equal_mean(sq_err, g) if g is not None else float(sq_err.mean())
        bias = patient_equal_mean(signed, g) if g is not None else float(signed.mean())

        present = t > ABSENT_THRESHOLD
        absent = ~present

        # Correlation is undefined without spread, and reporting 0.0 there would read
        # as "no relationship" when the truth is "not enough variation to tell".
        if present.sum() >= 3 and t[present].std() > 0 and e[present].std() > 0:
            pearson = float(stats.pearsonr(t[present], e[present])[0])
            spearman = float(stats.spearmanr(t[present], e[present])[0])
        else:
            pearson = spearman = np.nan

        rows.append({
            "cell_type": ctype,
            "n": int(ok.sum()),
            "n_groups": int(g.nunique()) if g is not None else np.nan,
            "mae": mae,
            "rmse": float(np.sqrt(mse)),
            "bias": bias,
            "pearson_r": pearson,
            "spearman_rho": spearman,
            "n_absent": int(absent.sum()),
            # Of the samples where this type is truly absent, the mean fraction the
            # method invented.
            "false_positive_fraction": (float(e[absent].mean()) if absent.any() else np.nan),
            "limit_of_detection": _limit_of_detection(e, t),
        })
    out = pd.DataFrame(rows)
    out.attrs["aggregation"] = "group-equal" if grp is not None else "sample-mean"
    return out


def _limit_of_detection(estimate: pd.Series, truth: pd.Series,
                        n_bins: int = 8) -> float:
    """
    Smallest true fraction at which the estimate still tracks the truth.

    True fractions are binned; a bin counts as detected if the mean estimate there
    exceeds the mean estimate in the truly-absent samples by more than that group's
    spread. Returns the lower edge of the smallest detected bin, or NaN when there is
    no absent group to compare against.
    """
    present = truth > ABSENT_THRESHOLD
    absent = ~present
    if absent.sum() < 3 or present.sum() < n_bins:
        return np.nan

    noise_mu = float(estimate[absent].mean())
    noise_sd = float(estimate[absent].std()) or 1e-9

    t_present = truth[present]
    try:
        bins = pd.qcut(t_present, q=n_bins, duplicates="drop")
    except ValueError:
        return np.nan

    for interval, idx in t_present.groupby(bins, observed=True).groups.items():
        idx = list(idx)
        if len(idx) < 2:
            continue
        if float(estimate[idx].mean()) > noise_mu + noise_sd:
            return float(interval.left)
    return np.nan


def overall_metrics(estimate: pd.DataFrame, truth: pd.DataFrame,
                    groups: pd.Series | None = None) -> dict:
    """
    Whole-table summary.

    `mae_primary` deliberately excludes the Astrocyte sidecar. Astrocyte is known to be
    unidentifiable against Tumor, so including it would let a method's handling of an
    acknowledged-unresolvable column move the headline number.
    """
    est, tru = _aligned(estimate, truth)
    grp = groups.reindex(est.index).astype(str) if groups is not None else None

    primary = [c for c in tru.columns if c in config.PRIMARY_CELL_TYPES]

    def _agg(per_sample: pd.Series) -> float:
        """Group-equal when we know the grouping, plain mean when we do not."""
        return (patient_equal_mean(per_sample, grp) if grp is not None
                else float(per_sample.mean()))

    def _mae(cols):
        return _agg((est[cols] - tru[cols]).abs().mean(axis=1))

    # Per-sample correlation across cell types: does the method get the *shape* of each
    # individual sample right? Averaged patient-equally like everything else.
    per_sample_r = []
    for s in est.index:
        e, t = est.loc[s], tru.loc[s]
        if np.isfinite(e).all() and e.std() > 0 and t.std() > 0:
            per_sample_r.append(stats.pearsonr(t, e)[0])
        else:
            per_sample_r.append(np.nan)
    r_series = pd.Series(per_sample_r, index=est.index)

    n_failed = int((~np.isfinite(est.to_numpy())).any(axis=1).sum())

    return {
        "n_samples": int(len(est)),
        "n_groups": int(grp.nunique()) if grp is not None else None,
        "aggregation": "group-equal" if grp is not None else "sample-mean",
        "mae_all_types": _mae(list(tru.columns)),
        "mae_primary": _mae(primary),
        # Aggregate the squared error first, take the root last: rooting per sample
        # and then averaging would compute a mean of RMSEs, which is a different and
        # smaller quantity than the RMSE.
        "rmse_primary": float(np.sqrt(
            _agg(((est[primary] - tru[primary]) ** 2).mean(axis=1)))),
        "mean_per_sample_r": _agg(r_series),
        "n_failed_samples": n_failed,
    }


def summarise_methods(estimates: dict[str, pd.DataFrame], truth: pd.DataFrame,
                      groups: pd.Series | None = None) -> pd.DataFrame:
    """One row per method, sorted by primary-type MAE (lower is better)."""
    rows = []
    for name, est in estimates.items():
        m = overall_metrics(est, truth, groups)
        m["method"] = name
        rows.append(m)
    out = pd.DataFrame(rows).set_index("method")
    cols = ["mae_primary", "rmse_primary", "mae_all_types", "mean_per_sample_r",
            "n_failed_samples", "n_samples", "n_groups", "aggregation"]
    return out[[c for c in cols if c in out.columns]].sort_values("mae_primary")


# =============================================================================
# IS THE WINNER ACTUALLY WINNING?
# =============================================================================

def paired_bootstrap_ties(estimates: dict[str, pd.DataFrame], truth: pd.DataFrame,
                          groups: pd.Series | None = None, n_boot: int = 2000,
                          seed: int = 0) -> pd.DataFrame:
    """
    Which methods are statistically indistinguishable from the best one?

    WHY THIS IS NOT OPTIONAL
    -----------------------
    A ranking table always has a first row. On a small benchmark the top several methods
    can be separated by less than the noise, and reporting the top row as "the selected
    method" then dresses up a coin flip as a decision. The survival stage already
    refuses to rank when its cohort cannot support it; applying a weaker standard to the
    criterion that actually *selects* the method would be backwards.

    The bootstrap is PAIRED and resamples GROUPS (donors, or tumours), not samples. Two
    methods scored on the same mixtures share all the difficulty of those mixtures, so
    resampling them independently would swamp a real difference with variance that
    cancels in the comparison. Groups rather than samples because samples within a donor
    are not independent — the same reason every other aggregation here is group-equal.

    Returns one row per method with the mean MAE difference from the best method, a 95%
    CI on that difference, and `tied_with_best` when the interval contains zero.
    """
    names = list(estimates)
    if not names:
        return pd.DataFrame()

    primary = [c for c in truth.columns if c in config.PRIMARY_CELL_TYPES]
    ref_index = [s for s in truth.index if s in estimates[names[0]].index]

    # Per-sample MAE for every method, on identical samples.
    per_sample = pd.DataFrame({
        name: (estimates[name].loc[ref_index, primary]
               - truth.loc[ref_index, primary]).abs().mean(axis=1)
        for name in names
    })

    if groups is not None:
        g = groups.reindex(ref_index).astype(str)
    else:
        g = pd.Series(ref_index, index=ref_index).astype(str)   # each sample its own group

    unique_groups = list(pd.unique(g))
    by_group = {k: [s for s in ref_index if g[s] == k] for k in unique_groups}

    def _score(sample_ids) -> pd.Series:
        return per_sample.loc[sample_ids].groupby(g[sample_ids]).mean().mean()

    observed = _score(ref_index)
    best = str(observed.idxmin())

    # A single group cannot be resampled: every draw is the same data, the interval
    # collapses to zero width, and every method looks significantly different. Say so
    # instead of emitting a confident answer.
    if len(unique_groups) < 3:
        out = pd.DataFrame({
            "mae": observed,
            "delta_vs_best": observed - observed[best],
            "ci_low": np.nan, "ci_high": np.nan,
            "tied_with_best": True,
        })
        out["note"] = (
            f"only {len(unique_groups)} group(s); a paired bootstrap over groups cannot "
            f"resolve differences. Every method is reported as tied and the selection "
            f"is a documented tie-break, not a measured win."
        )
        return out.sort_values("mae")

    rng = np.random.default_rng(seed)
    deltas = {n: np.empty(n_boot) for n in names}
    for i in range(n_boot):
        pick = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        ids = [s for k in pick for s in by_group[k]]
        sub = per_sample.loc[ids]
        sc = sub.groupby(g[ids]).mean().mean()
        for n in names:
            deltas[n][i] = sc[n] - sc[best]

    rows = []
    for n in names:
        lo, hi = np.quantile(deltas[n], [0.025, 0.975])
        rows.append({
            "method": n,
            "mae": float(observed[n]),
            "delta_vs_best": float(observed[n] - observed[best]),
            "ci_low": float(lo),
            "ci_high": float(hi),
            # Zero inside the interval means "not distinguishable from the best".
            "tied_with_best": bool(lo <= 0.0 <= hi),
            "note": "",
        })
    return pd.DataFrame(rows).set_index("method").sort_values("mae")
