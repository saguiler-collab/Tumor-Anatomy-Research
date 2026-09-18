#!/usr/bin/env python3
"""
methylation_celltypes.py — per-sample, per-cell-type immune composition from DNA methylation.

Executes the per-cell-type addendum of `prespecified/immune_failure_factors.md`, committed
before this ran.

WHY IT EXISTS
-------------
Leukocyte fraction is a single aggregate number, so the immune arm could corroborate the
DIRECTION of an error and never a per-type magnitude. EpiDISH's reference resolves individual
immune types from EpiDISH's 333 reference HM450 CpGs -- of which only the COMPLETE-CASE
subset is actually used, because `epidish()` is given `b[complete.cases(b), ]`. On TCGA-LGG
that is 255 of 333 (64 CpGs are all-NA across every sample). The count that belongs in a
results sentence is the one EpiDISH consumed, not the one the package ships.

It exists specifically to test one unpredicted observation: methods place MORE B cells than
T cells, while the reference atlas holds 5x more T than B.

WHAT IT DOES NOT DO
-------------------
`Macrophage_Microglia` is not mapped. Microglia are brain-resident and appear in no blood
reference; `Mono` is a different population and substituting it would invent a correspondence.
And the reference is a BLOOD reference applied to brain tumour tissue, so what it estimates is
the leukocyte SUB-COMPOSITION -- relative within the immune compartment, not an absolute
fraction of the sample.

    python scripts/methylation_celltypes.py --cohort lgg
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

#: EpiDISH type -> roster type. Deliberately partial; see the module docstring.
MAP = {"CD4T": "T_cell", "CD8T": "T_cell", "B": "B_cell", "NK": "NK_cell"}
COMPARE = ["T_cell", "B_cell", "NK_cell"]


def k4s(s):
    f = str(s).split("-")
    return "-".join(f[:3] + [f[3][:2]]) if len(f) > 3 else str(s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["gbm", "lgg"], default="lgg")
    args = ap.parse_args()
    tag = "" if args.cohort == "gbm" else f"_{args.cohort}"
    beta_p = config.PROCESSED_DIR / f"{args.cohort}_methylation_epidish_probes.csv.gz"
    if not beta_p.exists():
        print(f"BLOCKED: {beta_p.name} missing."); return 2

    out_p = config.RESULTS_DIR / f"methylation_celltypes{tag}.csv"
    r = subprocess.run(["Rscript", "-e", f'''
        suppressMessages(library(EpiDISH)); data(centDHSbloodDMC.m)
        b <- as.matrix(read.csv(gzfile("{beta_p}"), row.names=1, check.names=FALSE))
        b <- b[complete.cases(b), , drop=FALSE]
        cat("beta:", nrow(b), "probes x", ncol(b), "samples\\n")
        res <- epidish(beta.m = b, ref.m = centDHSbloodDMC.m, method = "RPC")
        write.csv(res$estF, "{out_p}")
        cat("estF:", paste(dim(res$estF), collapse=" x "), "\\n")
    '''], capture_output=True, text=True)
    print(r.stdout.strip())
    # NEVER VALIDATE AGAINST A HARD-CODED LITERAL: read what EpiDISH actually consumed out of
    # its own stdout, rather than asserting 333 and being wrong in the results sentence.
    m_used = re.search(r"beta:\s*(\d+)\s*probes", r.stdout or "")
    n_used = int(m_used.group(1)) if m_used else -1
    n_supplied = int(pd.read_csv(beta_p, index_col=0).shape[0])
    if n_used > 0 and n_used < n_supplied:
        print(f"  NOTE: {n_supplied - n_used} of {n_supplied} reference CpGs dropped by "
              f"complete-case filtering; EpiDISH ran on {n_used}.")
    if r.returncode != 0:
        print("EpiDISH FAILED:\n" + r.stderr[-1200:]); return 2

    est = pd.read_csv(out_p, index_col=0)
    print(f"\nmethylation-derived leukocyte sub-composition, {args.cohort.upper()} "
          f"({len(est)} samples):")
    print(est.mean().round(4).sort_values(ascending=False).to_string())

    # roster view: the three types this reference can speak to
    ros = pd.DataFrame({
        "T_cell": est[["CD4T", "CD8T"]].sum(axis=1),
        "B_cell": est["B"], "NK_cell": est["NK"]})
    ros = ros.div(ros.sum(axis=1), axis=0)          # within the three compared types
    print(f"\nrenormalised within {COMPARE} (what the comparison uses):")
    print(ros.mean().round(4).to_string())
    tb = float((ros["T_cell"] > ros["B_cell"]).mean())
    print(f"\nPREDICTION: T > B.  Observed in {tb:.1%} of samples; "
          f"mean T {ros['T_cell'].mean():.4f} vs B {ros['B_cell'].mean():.4f}  "
          f"-> {'HELD' if ros['T_cell'].mean() > ros['B_cell'].mean() else 'FAILED'}")

    # --- does the anomaly even exist in this cohort's estimates? ----------------------
    full_p = config.RESULTS_DIR / f"estimates_full{tag}.csv"
    comp = {}
    if full_p.exists():
        full = pd.read_csv(full_p)
        scol = next(c for c in ("sample", "sample_id") if c in full.columns)
        full["k"] = full[scol].map(k4s)
        ros.index = [k4s(i) for i in ros.index]
        ros = ros[~ros.index.duplicated()]
        print(f"\n{'method':18s} {'est T':>8s} {'est B':>8s} {'est T>B?':>9s} "
              f"{'meth T':>8s} {'meth B':>8s}")
        print("-" * 64)
        for m, g in full.groupby("method"):
            g = g.drop_duplicates("k").set_index("k")
            idx = [k for k in g.index if k in ros.index]
            if len(idx) < 50:
                continue
            e = g.loc[idx, COMPARE]
            e = e.div(e.sum(axis=1).replace(0, np.nan), axis=0)
            mt, mb = ros.loc[idx, "T_cell"].mean(), ros.loc[idx, "B_cell"].mean()
            et, eb = e["T_cell"].mean(), e["B_cell"].mean()
            comp[m] = {"n": len(idx), "est_T": round(float(et), 4),
                       "est_B": round(float(eb), 4),
                       "est_T_exceeds_B": bool(et > eb),
                       "meth_T": round(float(mt), 4), "meth_B": round(float(mb), 4)}
            print(f"{m:18s} {et:8.4f} {eb:8.4f} {str(et > eb):>9s} {mt:8.4f} {mb:8.4f}")
        n_wrong = sum(1 for v in comp.values() if not v["est_T_exceeds_B"])
        print(f"\n{n_wrong} of {len(comp)} methods put B above T; "
              f"methylation puts T above B in {tb:.1%} of samples")

    json.dump({"cohort": args.cohort,
               "prespecification": "prespecified/immune_failure_factors.md (per-cell-type "
                                   "addendum)",
               "reference": (f"EpiDISH centDHSbloodDMC.m, RPC; {n_used} of {n_supplied} "
                             f"reference CpGs used after complete-case filtering"),
               "n_reference_cpgs_supplied": int(n_supplied),
               "n_reference_cpgs_used": int(n_used),
               "not_mapped": {"Macrophage_Microglia": "microglia are brain-resident and appear "
                                                      "in no blood reference; Mono is a "
                                                      "different population"},
               "is_relative_within_compartment": True,
               "n_samples": int(len(est)),
               "mean_subcomposition": {k: round(float(v), 4) for k, v in est.mean().items()},
               "prediction_T_exceeds_B": {"held": bool(ros["T_cell"].mean() >
                                                       ros["B_cell"].mean()),
                                          "fraction_of_samples": round(tb, 4)},
               "methods": comp},
              open(config.RESULTS_DIR / f"methylation_celltypes{tag}.json", "w"), indent=2)
    print(f"\nwrote results/methylation_celltypes{tag}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
