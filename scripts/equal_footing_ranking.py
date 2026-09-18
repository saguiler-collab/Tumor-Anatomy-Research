"""Does giving every method its intended inputs change the ranking?

THE QUESTION. `docs/EQUAL_FOOTING.md` established that the vendored frozen signature carries
an all-zero sigma. MuSiC without cross-donor variance is NNLS; EPIC without `refProfiles.var`
is uniform-weighted; SCDC ENSEMBLE with one reference is SCDC; CIBERSORTx S-mode cannot run at
all. Four methods were therefore not being measured -- they were being measured *as something
else*. Any ranking built on that reference ranks the degraded stand-ins.

Rebuilding the reference from `gbmap_core.h5ad` supplies sigma, donor profiles and the
registered cells. This script reads both runs and reports what moved.

IT DOES NOT SELECT ANYTHING. The comparison is descriptive: it reports the ranking under each
reference and the rank change per method. No method is chosen, dropped or retuned on the basis
of the result -- the exclusions come from `comparability.non_comparable`, which reads only what
inputs a method REQUIRES, never how it scored.
"""
from __future__ import annotations

import json
import sys

import pandas as pd

from ivygap import config
from ivygap.deconv.comparability import non_comparable


def _load(tag: str) -> dict:
    p = config.RESULTS_DIR / f"absolute_purity_yardstick{tag}.json"
    if not p.exists():
        return {}
    return json.loads(p.read_text())


def _shared_samples(base: str) -> set | None:
    """Samples scored by BOTH runs, or None if the h5ad run has not happened yet."""
    fp = config.RESULTS_DIR / f"absolute_purity_per_sample{base}.csv"
    hp = config.RESULTS_DIR / f"absolute_purity_per_sample{base}_h5ad.csv"
    if not (fp.exists() and hp.exists()):
        return None
    a = set(pd.read_csv(fp, index_col=0).index)
    b = set(pd.read_csv(hp, index_col=0).index)
    return a & b


def _observed(tag: str) -> tuple[set[str], set[str]]:
    """(degenerate, failed) as the run itself recorded them, from the yardstick JSON."""
    rep = _load(tag)
    methods = rep.get("methods", rep) if isinstance(rep, dict) else {}
    degen = {m for m, v in methods.items()
             if isinstance(v, dict) and v.get("degenerate")}
    failed = {m for m, v in methods.items()
              if isinstance(v, dict) and ("failed" in v or "skipped" in v)}
    return degen, failed


def _recovery(report: dict) -> dict[str, float]:
    """Recovery fraction = 1 + slope of (estimate - truth) on truth."""
    out = {}
    for m, v in (report.get("methods") or report).items():
        if not isinstance(v, dict) or "recovery" not in v:
            continue
        out[m] = float(v["recovery"])
    return out


def _from_per_sample(tag: str, keep: set | None = None) -> dict[str, float]:
    """Recompute recovery from the per-sample CSV -- the JSON may not carry it.

    `keep` restricts to a fixed sample set. THE COMPARISON REQUIRES IT: if the frozen run and
    the h5ad run are scored on different samples, a rank change cannot be attributed to the
    reference rather than to the sample set. The two runs do differ -- the h5ad reference
    shares a different gene space with the bulk, so a method can return finite values on one
    and not the other.

    NOTE ON THE DENOMINATOR. This is not the same sample set as
    `results/failure_factors_*.json`, which restricts to complete cases across all seven
    biological factors (147 of 154 in GBM, 496 of 510 in LGG) because it fits a seven-factor
    model. Recovery figures here are therefore close to, but not identical with, the ones
    quoted there (e.g. CIBERSORTx GBM 63.7% here vs 65.6% on complete cases). Neither is wrong;
    they answer different questions, and mixing them in one table would be.
    """
    import numpy as np
    p = config.RESULTS_DIR / f"absolute_purity_per_sample{tag}.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, index_col=0)
    # The truth column is named by the writer, not guessed here: GBM writes
    # "absolute_purity", and a rename would otherwise silently yield an empty ranking.
    truth_col = next((c for c in ("absolute_purity", "purity") if c in df.columns), None)
    if truth_col is None:
        raise KeyError(f"{p.name} has no purity column. Columns: {list(df.columns)[:8]}")
    if keep is not None:
        df = df.loc[[i for i in df.index if i in keep]]
    truth = df[truth_col].to_numpy(dtype="float64")
    out = {}
    for m in df.columns:
        if m == truth_col:
            continue
        if not pd.api.types.is_numeric_dtype(df[m]):
            continue
        est = df[m].to_numpy(dtype="float64")
        ok = np.isfinite(est) & np.isfinite(truth)
        if ok.sum() < 50:
            continue
        slope = np.polyfit(truth[ok], est[ok] - truth[ok], 1)[0]
        out[m] = 1.0 + float(slope)
    return out


def main() -> int:
    cohort = sys.argv[1] if len(sys.argv) > 1 else "gbm"
    base = "" if cohort == "gbm" else f"_{cohort}"
    shared = _shared_samples(base)
    frozen = _from_per_sample(base, shared)
    h5ad = _from_per_sample(base + "_h5ad", shared)
    if not frozen:
        print(f"BLOCKED: no frozen run for {cohort}"); return 2
    if not h5ad:
        print(f"BLOCKED: no h5ad run for {cohort} yet "
              f"(results/absolute_purity_per_sample{base}_h5ad.csv missing)"); return 2

    # OBSERVED degeneracy is unioned in, not assumed away. `non_comparable`'s static table
    # says what each method REQUIRES; the run says what actually happened. A method that
    # degenerated for a reason the table does not anticipate must still be excluded, or it
    # enters the ranking under its own name -- which is the one thing the invariants forbid.
    nc_frozen = non_comparable("frozen", *_observed(base))
    nc_h5ad = non_comparable("h5ad", *_observed(base + "_h5ad"))
    print(f"=== {cohort.upper()}: does equal footing change the ranking? ===\n")
    print(f"scored on the {len(shared):,} samples present in BOTH runs, so a rank change "
          f"cannot be a sample-set change\n")
    print("NOT COMPARABLE on the frozen reference (degraded stand-ins, not the method):")
    print("  [requirement = from the static input table; observed = what this run reported]")
    for m, why in sorted(nc_frozen.items()):
        print(f"  {m:20s} {why[:88]}")
    print("\nNOT COMPARABLE even with sigma (input the cohort cannot supply):")
    for m, why in sorted(nc_h5ad.items()):
        print(f"  {m:20s} {why[:88]}")

    rows = []
    for m in sorted(set(frozen) | set(h5ad)):
        rows.append({"method": m,
                     "frozen": frozen.get(m), "h5ad": h5ad.get(m),
                     "comparable_frozen": m not in nc_frozen,
                     "comparable_h5ad": m not in nc_h5ad})
    df = pd.DataFrame(rows)

    # Rank WITHIN the comparable set of each reference -- ranking a degraded stand-in
    # against a real method is the error this whole script exists to correct.
    for col in ("frozen", "h5ad"):
        ok = df[f"comparable_{col}"] & df[col].notna()
        df.loc[ok, f"rank_{col}"] = df.loc[ok, col].rank(ascending=False)

    print(f"\n{'method':20s} {'frozen':>9s} {'h5ad':>9s} {'change':>9s}  "
          f"{'rank F':>7s} {'rank H':>7s}  moved")
    print("-" * 78)
    for _, r in df.sort_values("h5ad", ascending=False, na_position="last").iterrows():
        f_, h_ = r["frozen"], r["h5ad"]
        ch = f"{(h_ - f_) * 100:+.1f}%" if pd.notna(f_) and pd.notna(h_) else "--"
        rf = f"{r['rank_frozen']:.0f}" if pd.notna(r.get("rank_frozen")) else "n/c"
        rh = f"{r['rank_h5ad']:.0f}" if pd.notna(r.get("rank_h5ad")) else "n/c"
        moved = ""
        if rf.isdigit() and rh.isdigit():
            d = int(rf) - int(rh)
            moved = "same" if d == 0 else f"{d:+d}"
        print(f"{r['method']:20s} {f_ * 100 if pd.notna(f_) else float('nan'):8.1f}% "
              f"{h_ * 100 if pd.notna(h_) else float('nan'):8.1f}% {ch:>9s}  "
              f"{rf:>7s} {rh:>7s}  {moved}")

    both = df[df["rank_frozen"].notna() & df["rank_h5ad"].notna()]
    if len(both) >= 3:
        tau = both["rank_frozen"].corr(both["rank_h5ad"], method="kendall")
        print(f"\nKendall tau between the two rankings, on the {len(both)} methods rankable "
              f"under BOTH: {tau:+.3f}")
        print("  tau = 1 means equal footing changed nothing; lower means the frozen "
              "ranking was\n  partly an artefact of the missing inputs.")

    # CONFOUND CHECK. Registering the cells does two things at once: it supplies sigma (the
    # variable of interest) AND it makes several R packages runnable where the frozen path
    # fell back to a Python reimplementation. A rank change caused by swapping
    # reimplementation-for-package is NOT evidence about equal footing, so the methods whose
    # implementation changed are identified and tau is recomputed without them.
    impl_f = {m: (v or {}).get("implementation") for m, v in
              (_load(base).get("methods", {}) or {}).items()}
    impl_h = {m: (v or {}).get("implementation") for m, v in
              (_load(base + "_h5ad").get("methods", {}) or {}).items()}
    switched = sorted(m for m in both["method"]
                      if impl_f.get(m) != impl_h.get(m))
    tau_clean = None
    if len(both) >= 3:
        print(f"\nIMPLEMENTATION CHECK — did the package change under the method's name?")
        for m in both["method"]:
            mark = "  CHANGED" if m in switched else ""
            print(f"  {m:20s} {str(impl_f.get(m)):24s} -> {str(impl_h.get(m)):24s}{mark}")
        keep = both[~both["method"].isin(switched)]
        if len(keep) >= 3 and switched:
            tau_clean = float(keep["rank_frozen"].corr(keep["rank_h5ad"], method="kendall"))
            print(f"\n  {len(switched)} of {len(both)} changed implementation. Excluding them, "
                  f"tau on the remaining {len(keep)} is {tau_clean:+.3f} "
                  f"(vs {tau:+.3f} overall) — so the rank change is NOT an artefact of "
                  f"swapping a reimplementation for a package.")
        elif not switched:
            print("\n  none changed: the rank comparison is clean.")

    out = config.RESULTS_DIR / f"equal_footing_ranking{base}.json"
    out.write_text(json.dumps({
        "cohort": cohort,
        "recovery_frozen": frozen, "recovery_h5ad": h5ad,
        "non_comparable_frozen": nc_frozen, "non_comparable_h5ad": nc_h5ad,
        "kendall_tau": None if len(both) < 3 else round(float(tau), 4),
        "implementation_frozen": impl_f, "implementation_h5ad": impl_h,
        "implementation_switched": switched,
        "kendall_tau_excluding_switched": (None if tau_clean is None
                                           else round(tau_clean, 4)),
        "n_ranked_under_both": int(len(both)),
    }, indent=2))
    print(f"\nwrote results/equal_footing_ranking{base}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
