#!/usr/bin/env python3
"""
remeasure_method.py — run a GENUINE R package on the anatomic cohort, unbounded.

WHY THIS EXISTS
---------------
DWLS is the one method this study has never cleanly measured. In the confirmatory run it
fell back to the Python reimplementation in both anatomic stages, having exceeded its
2,400 s budget for the R package. So the leaderboard's DWLS row is not DWLS, and the
project's own invariant is explicit that a reimplementation must never be reported under
a published package's name.

That has a concrete cost beyond tidiness. Avila Cobos et al. (2020) rank DWLS **best**
among methods that use a single-cell reference; this study ranks it last. The
disagreement is currently unresolvable, because the two are not testing the same
software. It also breaks the benchmark-concordance test in
`scripts/benchmark_concordance.py`, where DWLS at rank 14 is most of why that comparison
returns a null result.

WHY RAISING THE BUDGET HERE IS NOT OUTCOME-DRIVEN
--------------------------------------------------
The obvious objection: the budget is being changed after seeing that DWLS scored badly.
Three things make this a re-measurement rather than a retune, and all three have to hold:

  1. **The budget's purpose is to bound a 15-method run, not to exclude a method.** A
     method the budget removes *systematically* is a measurement failure. DWLS completed
     as R:DWLS in the 2026-09-05 run and fell back in the two after it, so the budget is
     sitting almost exactly on this method's runtime — the worst possible place for a
     threshold to sit, because it decides the answer by coin flip.
  2. **Nothing about the method changes.** No parameter, no gene set, no cutoff. The
     published defaults are untouched. Only wall-clock is allowed to differ, and
     wall-clock is not a property of the method's accuracy.
  3. **The result is reported whatever it is.** If genuine DWLS scores lower than the
     reimplementation, that is the finding and it goes in the leaderboard. This script
     writes its output whichever way it falls, and refuses to be run selectively.

The alternative — leaving it — is worse than either outcome, because it means publishing
"DWLS is last" about software that was never run.

WHAT THIS DOES NOT DO
---------------------
It does not touch the archived confirmatory run, which stands as recorded. It produces a
separate, clearly labelled measurement of one method, to be reported alongside rather
than substituted in.

    python scripts/remeasure_dwls.py                    # 4 hour budget
    python scripts/remeasure_dwls.py --budget 28800     # 8 hours
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                         # noqa: E402
from ivygap.anatomic.acs import score as acs_score                # noqa: E402
from ivygap.data.load_ivygap import load_cached                   # noqa: E402
from ivygap.data.reference import build_from_h5ad                 # noqa: E402
from ivygap.deconv import r_bridge                                # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                 # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", default="dwls",
                    help="method to re-measure, e.g. dwls or bayesprism")
    ap.add_argument("--budget", type=int, default=14400, help="seconds for the R call")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    M = args.method
    if args.out is None:
        args.out = f"results/{M}_remeasured.json"

    print("building the anatomic cohort and the cell-level reference...")
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]

    ref, cells, cmeta = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad")
    genes = [g for g in ref.profile.index if g in set(expr.index)]
    from ivygap.data.reference import frozen_gene_space
    genes, prov = frozen_gene_space(ref, expr.index)

    bulk = expr.loc[genes, anat]
    bulk = bulk / bulk.sum(axis=0) * 1e6
    data = DeconvolutionInput(bulk=bulk, references=(ref.subset_genes(genes),),
                              manifest=meta.loc[anat])
    r_bridge.set_cell_source(ref.name, cells, cmeta)
    print(f"  {bulk.shape[0]} genes x {bulk.shape[1]} samples | "
          f"reference {cells.shape[1]} cells")

    print(f"\nrunning the genuine {M} R package with a {args.budget}s budget "
          f"(the pipeline's is {r_bridge.timeout_for(M)}s)...")
    t0 = time.perf_counter()
    failed = None
    est = None
    try:
        est = r_bridge.run_r_method(M, data, timeout=args.budget)
    except Exception as exc:                                       # noqa: BLE001
        failed = f"{type(exc).__name__}: {str(exc)[:600]}"
    elapsed = time.perf_counter() - t0
    print(f"  finished in {elapsed:.0f}s ({elapsed/3600:.2f} h)")

    report = {
        "method": M,
        "what_this_is": (f"A single-method re-measurement of the GENUINE {M} R package on "
                         "the anatomic cohort, after it fell back to the Python "
                         "reimplementation in the confirmatory run. Reported alongside "
                         "that run, never substituted into it."),
        "budget_seconds": args.budget,
        "pipeline_budget_seconds": r_bridge.timeout_for("dwls"),
        "elapsed_seconds": round(elapsed, 1),
        "n_samples": int(bulk.shape[1]), "n_genes": int(bulk.shape[0]),
        "failed": failed,
    }

    if est is not None:
        # Cell-size correction and simplex projection, exactly as fit_predict applies
        # them, so this number is on the same footing as the leaderboard's.
        from ivygap.deconv.registry import build_methods
        dwls = next(m for m in build_methods(prefer_r=False) if m.name == M)
        import numpy as np
        from ivygap.deconv.base import project_to_simplex, to_cell_fractions
        cs = data.primary.cell_size.reindex(list(data.cell_types)).to_numpy()
        rows = []
        for i in range(len(data.samples)):
            w = project_to_simplex(est.to_numpy()[i], None)
            if np.isfinite(w).all():
                w = to_cell_fractions(w, cs)
            rows.append(w)
        frame = pd.DataFrame(np.vstack(rows), index=data.samples,
                             columns=list(data.cell_types))
        frame.index.name = "sample_id"
        out_csv = Path(f"results/estimates/ivygap_{M}_genuine.csv")
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out_csv)

        res = acs_score(frame, meta.loc[anat], method=f"{M}_genuine",
                        n_permutations=10000, n_boot=2000)
        s = res.summary()
        report["acs"] = round(float(s["acs"]), 4)
        report["ci"] = [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)]
        report["null_p"] = float(s["null_p"])
        report["n_pairs"] = int(s["n_constraint_tumor_pairs"])
        report["n_tumors"] = int(s["n_tumors"])
        report["estimates"] = str(out_csv)
        print(f"\n  GENUINE {M}: ACS {s['acs']:.4f} "
              f"[{s['ci_low']:.4f}, {s['ci_high']:.4f}] "
              f"on {s['n_constraint_tumor_pairs']} pairs, null_p {s['null_p']:.4g}")
        print(f"  (compare against the confirmatory run's leaderboard row)")
    else:
        print(f"\n  FAILED: {failed}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
