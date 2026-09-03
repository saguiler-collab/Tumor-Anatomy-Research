"""
calibration.py — does a method's stated uncertainty match its actual error?

WHY THIS EXISTS
---------------
A credible interval is a claim: "the truth falls inside this range 95% of the time."
Nothing about producing intervals makes that claim true, and an interval that covers
the truth 23% of the time while advertising 95% is worse than reporting no interval at
all — it converts a known unknown into false confidence.

This module measures the claim against known-truth pseudobulk. It found a real defect
during development: the Bayesian sampler's default treated FPKM as 100,000 independent
multinomial reads, which made intervals roughly five times too narrow (coverage 0.23
against a nominal 0.95). The default was corrected on the strength of this measurement,
and the measurement now runs on every benchmark so the same drift cannot recur
unnoticed.

WHAT IS REPORTED
----------------
    coverage_95      fraction of (sample, cell type) pairs whose truth lies inside the
                     95% credible interval. Nominal is 0.95.
    mean_ci_width    average interval width. Coverage bought by making every interval
                     span [0, 1] is not calibration, so width is reported beside it.
    calibration_gap  coverage_95 - 0.95. Negative means overconfident, the dangerous
                     direction.
    interval_score   proper scoring rule (Gneiting & Raftery 2007) that rewards narrow
                     intervals and penalises misses. Lower is better; unlike coverage
                     alone it cannot be gamed by widening everything.

Only methods that expose a `posterior_` with `ci_lower`/`ci_upper` are scored. Point
methods are listed as `uncertainty_reported: false` rather than being given a default
or silently omitted — "this method does not quantify its uncertainty" is itself a
finding worth putting in the results table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

NOMINAL = 0.95


def interval_score(truth: np.ndarray, lo: np.ndarray, hi: np.ndarray,
                   alpha: float = 1 - NOMINAL) -> float:
    """
    Negatively-oriented interval score for a central (1-alpha) interval.

        (hi - lo) + (2/alpha)(lo - y) if y < lo + (2/alpha)(y - hi) if y > hi

    A proper scoring rule: it is minimised by the true predictive distribution, so a
    method cannot improve it by widening its intervals to buy coverage.
    """
    width = hi - lo
    below = (lo - truth) * (truth < lo)
    above = (truth - hi) * (truth > hi)
    return float(np.mean(width + (2.0 / alpha) * (below + above)))


def measure_coverage(method, truth: pd.DataFrame) -> dict:
    """
    Score one already-fitted method's intervals against known composition.

    `method` must have been run via `fit_predict` already; this reads its stored
    `posterior_` rather than re-fitting, so calibration never costs a second run.
    """
    posterior = getattr(method, "posterior_", None) or {}
    if "ci_lower" not in posterior or "ci_upper" not in posterior:
        return {"method": method.name, "uncertainty_reported": False}

    lo, hi = posterior["ci_lower"], posterior["ci_upper"]
    shared = [s for s in truth.index if s in lo.index]
    types = [c for c in truth.columns if c in lo.columns]
    if not shared or not types:
        return {"method": method.name, "uncertainty_reported": False,
                "reason": "posterior and truth share no samples or cell types"}

    y = truth.loc[shared, types].to_numpy(dtype="float64")
    l = lo.loc[shared, types].to_numpy(dtype="float64")
    h = hi.loc[shared, types].to_numpy(dtype="float64")

    ok = np.isfinite(y) & np.isfinite(l) & np.isfinite(h)
    if not ok.any():
        return {"method": method.name, "uncertainty_reported": False,
                "reason": "no finite interval/truth pairs"}
    y, l, h = y[ok], l[ok], h[ok]

    coverage = float(np.mean((y >= l) & (y <= h)))
    return {
        "method": method.name,
        "uncertainty_reported": True,
        "n_intervals": int(ok.sum()),
        "coverage_95": coverage,
        "nominal": NOMINAL,
        "calibration_gap": coverage - NOMINAL,
        "mean_ci_width": float(np.mean(h - l)),
        "interval_score": interval_score(y, l, h),
        "verdict": _verdict(coverage),
    }


def _verdict(coverage: float) -> str:
    """
    Plain-language reading. The asymmetry is deliberate: an interval that is too wide
    is uninformative, but one that is too narrow is actively misleading, so the
    threshold for calling out overconfidence is tighter.
    """
    gap = coverage - NOMINAL
    if gap < -0.25:
        return ("SEVERELY OVERCONFIDENT — intervals are far too narrow. Do not quote "
                "them as credible intervals without this caveat.")
    if gap < -0.10:
        return "OVERCONFIDENT — intervals are narrower than the errors justify."
    if gap > 0.04:
        return "CONSERVATIVE — intervals are wider than needed; check mean_ci_width."
    return "WELL CALIBRATED — realised coverage is close to nominal."


def calibration_table(methods: list, truth: pd.DataFrame) -> pd.DataFrame:
    """One row per method. Methods without uncertainty are listed, not dropped."""
    rows = [measure_coverage(m, truth) for m in methods]
    table = pd.DataFrame(rows).set_index("method")
    if "coverage_95" in table.columns:
        table = table.sort_values("coverage_95", ascending=False, na_position="last")
    return table
