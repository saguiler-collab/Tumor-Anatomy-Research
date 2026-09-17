#!/usr/bin/env python3
"""
epic_close_out.py — settle ROAD_TO_PAPER 0.4, the last item marked BLOCKING.

THE DECISION RULE, STATED BEFORE THE NUMBERS ARE COMPUTED
---------------------------------------------------------
This project's invariant is *never report a degenerate method under its own name* — MuSiC
without cross-donor variance is NNLS, SCDC ENSEMBLE with one reference is SCDC. EPIC's
published contribution is weighting genes by their variability across the reference, supplied
as `refProfiles.var`. The harness never passed it, and the package said so:

    'refProfiles.var' not defined; using identical weights for all genes

**So EPIC-with-`refProfiles.var` is EPIC, and EPIC-without is a constrained least squares with
uniform weights. That follows from the invariant and from the package's own warning, and it is
true whichever version scores higher.** The ACS numbers below cannot change it; they are
computed to report the size of the correction, not to choose it.

The other two findings are recording gaps, also score-independent:
  * `otherCells` — EPIC's ninth column, dropped by the harness and never archived;
  * `fit.gof` — convergence, never checked.

    python scripts/epic_close_out.py
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
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

D = config.RESULTS_DIR / "diagnostics"


def main() -> int:
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    man = meta.loc[anat]

    out: dict = {"decision_rule": (
        "EPIC's published contribution is variance weighting via refProfiles.var. The harness "
        "never supplied it and the package warned. By this project's invariant -- never report "
        "a degenerate method under its own name -- the variance-weighted run IS EPIC and the "
        "uniform-weight run is not, whichever scores higher. Stated before the scores were "
        "computed.")}

    scored = {}
    for tag, f in (("EPIC_variance_weighted", "epic_variance_weighted_estimates.csv"),
                   ("EPIC_uniform_weights", "epic_uniform_weights_estimates.csv")):
        p = D / f
        if not p.exists():
            print(f"BLOCKED: {p} missing"); return 2
        est = pd.read_csv(p, index_col=0)
        est.index = est.index.astype(str)          # numeric Ivy GAP ids round-trip as int64
        keep = [s for s in est.index if s in man.index]
        if len(keep) != len(anat):
            print(f"  {tag}: {len(keep)} of {len(anat)} samples match the manifest")
        r = acs_score(est.loc[keep], man, method=tag, n_permutations=10000, n_boot=2000)
        s = r.summary()
        scored[tag] = {"acs": round(float(s["acs"]), 4),
                       "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
                       "null_p": float(s["null_p"]),
                       "n_pairs": int(s["n_constraint_tumor_pairs"]),
                       "n_samples": len(keep)}
        print(f"{tag:26s} ACS {s['acs']:.4f}  CI [{s['ci_low']:.4f}, {s['ci_high']:.4f}]  "
              f"p={s['null_p']:.4f}  n={len(keep)}")

    a, b = scored["EPIC_variance_weighted"]["acs"], scored["EPIC_uniform_weights"]["acs"]
    print(f"\ndelta (EPIC as published - EPIC as run) = {a - b:+.4f}")
    out["acs"] = scored
    out["delta_variance_minus_uniform"] = round(a - b, 4)

    # --- how different are the estimates themselves, not just their rank? ------------
    v = pd.read_csv(D / "epic_variance_weighted_estimates.csv", index_col=0)
    u = pd.read_csv(D / "epic_uniform_weights_estimates.csv", index_col=0)
    shared = v.index.intersection(u.index)
    diff = (v.loc[shared] - u.loc[shared]).abs()
    print(f"\nper-cell-type mean |difference| between the two EPIC runs:")
    print(diff.mean().round(4).to_string())
    out["mean_abs_difference_per_type"] = {k: round(float(x), 4)
                                           for k, x in diff.mean().items()}
    out["max_abs_difference"] = round(float(diff.to_numpy().max()), 4)
    print(f"max |difference| anywhere: {out['max_abs_difference']:.4f}")

    # --- the two recording gaps -------------------------------------------------------
    c = pd.read_csv(D / "epic_convergence.csv", index_col=0)
    o = pd.read_csv(D / "epic_othercells.csv", index_col=0)
    codes = c["convergeCode"].value_counts().to_dict()
    out["convergence"] = {
        "n_samples_probed": int(len(c)),
        "n_samples_in_cohort": len(anat),
        "converge_code_counts": {str(k): int(x) for k, x in codes.items()},
        "n_failed": int(sum(x for k, x in codes.items() if k != 0)),
        "spearman_fit_min": round(float(c["spearmanR"].min()), 4),
        "caveat": "PARTIAL. The diagnostic covers %d of %d anatomic samples, so the "
                  "non-convergence rate is an estimate on a subset, not a cohort figure."
                  % (len(c), len(anat))}
    out["other_cells"] = {
        "n_samples_probed": int(len(o)),
        "max": float(o["otherCells"].max()), "median": float(o["otherCells"].median()),
        "reading": "EPIC's ninth column is essentially empty on this cohort, so the harness "
                   "dropping it changes almost nothing -- but it was never recorded until "
                   "now, and 'small' is a measurement rather than an assumption."}
    print(f"\nconvergence: {out['convergence']['n_failed']} of {len(c)} probed samples "
          f"failed (codes {codes}); fit spearman min {c['spearmanR'].min():.4f}")
    print(f"otherCells: max {o['otherCells'].max():.3e}, median "
          f"{o['otherCells'].median():.3e}, on {len(o)} probed samples")

    (config.RESULTS_DIR / "epic_close_out.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote results/epic_close_out.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
