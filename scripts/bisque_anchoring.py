"""Is Bisque's cohort-mean composition the TISSUE's, or its REFERENCE's? A truth-free test.

Under the h5ad (raw/X, donor-level) reference Bisque is the one method that orders the lymphoid
compartment T > NK > B in both TCGA cohorts, as DNA methylation does (`lymphoid_ordering_*h5ad`).
This asks whether that agreement is a measurement.

MECHANISM, FROM THE PACKAGE SOURCE (BisqueRNA::ReferenceBasedDecomposition, use.overlap = FALSE --
the only mode TCGA permits, since no TCGA subject is in GBmap):

    Y.train <- sc.ref %*% sc.props                         # one pseudobulk per atlas DONOR
    SemisupervisedTransformBulk(gene, Y.train, X.pred):
        X.pred.scaled <- scale(X.pred[gene, ])             # z-score the bulk across samples
        Y.pred <- X.pred.scaled * shrink.scale + Y.center  # ...re-centred on the DONOR mean

Every gene of the transformed cohort therefore has mean `sc.ref %*% rowMeans(sc.props)`, and a
sum-to-one non-negative least-squares fit returns, on average, `rowMeans(sc.props)` -- the
reference's donor-mean composition -- whatever the tissue holds. The per-sample deviations come
from the bulk; the LEVEL comes from the reference.

THE TEST NEEDS NO GROUND TRUTH. Reproduce the exact cells the TCGA h5ad arm handed Bisque (the
project's own `_balanced_cell_sample`, same caps, same seed -- obs only, no expression), compute
the donor-mean composition Bisque anchors to, and measure how far every method's cohort mean lies
from it. A method measuring tissue has no reason to land on the reference's composition; an
anchored one lands on it in every cohort.

Selects nothing; uses no outcome.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from ivygap import config
from ivygap.data.reference import (MAX_CELLS_PER_DONOR_TYPE, MAX_TOTAL_CELLS,
                                   _balanced_cell_sample, _collapse_to_roster)

LYMPH = ["T_cell", "NK_cell", "B_cell"]


def donor_mean(obs: pd.DataFrame) -> pd.Series:
    """BisqueRNA CalculateSCCellProportions (per-donor table / n), averaged over donors."""
    p = obs.groupby("donor")["cell_type"].value_counts(normalize=True).unstack(fill_value=0)
    return p.reindex(columns=config.CELL_TYPES, fill_value=0).mean()


def lymph(s: pd.Series) -> dict:
    v = s[LYMPH]
    return {k: round(float(x), 4) for k, x in (v / v.sum()).items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(config.RESULTS_DIR / "bisque_anchoring.json"))
    a = ap.parse_args()

    with h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r") as f:
        obs_raw = read_elem(f["obs"])
    obs = pd.DataFrame({
        "donor": obs_raw["donor_id"].astype(str).to_numpy(),
        "cell_type": _collapse_to_roster(obs_raw["annotation_level_3"].astype(str),
                                         config.GBMAP_CELL_TYPE_MAP).to_numpy()},
        index=obs_raw.index.astype(str))
    obs["_row"] = np.arange(len(obs))
    mapped = obs.dropna(subset=["cell_type"])
    keep, _ = _balanced_cell_sample(mapped, max_cells_per_donor_type=MAX_CELLS_PER_DONOR_TYPE,
                                    max_total_cells=MAX_TOTAL_CELLS, seed=config.RANDOM_SEED)
    ref_cells = mapped.loc[keep]
    prior = donor_mean(ref_cells)

    report = {
        "what_this_is": ("Truth-free test of whether each method's TCGA cohort-mean composition "
                         "is the reference's donor-mean composition (Bisque no-overlap "
                         "anchoring). Selects nothing."),
        "reference_cells": {"n_atlas_cells_mapped": int(len(mapped)),
                            "n_reference_cells": int(len(ref_cells)),
                            "n_donors": int(ref_cells["donor"].nunique()),
                            "max_cells_per_donor_type": MAX_CELLS_PER_DONOR_TYPE,
                            "max_total_cells": MAX_TOTAL_CELLS, "seed": config.RANDOM_SEED},
        "lymphoid_composition": {
            "atlas_pooled_cells": lymph(mapped["cell_type"].value_counts(normalize=True)
                                        .reindex(config.CELL_TYPES, fill_value=0)),
            "atlas_donor_mean": lymph(donor_mean(mapped)),
            "capped_reference_donor_mean_THE_PRIOR": lymph(prior),
        },
        "reference_prior_all_types": {k: round(float(v), 4) for k, v in prior.items()},
        "cohorts": {},
    }
    for c, est_f, lo_f in (("gbm", "estimates_full_h5ad.csv", "lymphoid_ordering_h5ad.json"),
                           ("lgg", "estimates_full_lgg_h5ad.csv", "lymphoid_ordering_lgg_h5ad.json")):
        est = pd.read_csv(config.RESULTS_DIR / est_f)
        truth = json.load(open(config.RESULTS_DIR / lo_f))["truth_mean"]
        rows = {}
        for m, g in est.groupby("method"):
            mu = g[config.CELL_TYPES].mean()
            ml = lymph(mu) if mu[LYMPH].sum() > 0 else None
            rows[m] = {
                "L1_to_reference_prior_all_types": round(float(np.abs(mu - prior).sum()), 4),
                "lymphoid": ml,
                "L1_lymphoid_to_prior": (None if ml is None else round(sum(
                    abs(ml[k] - report["lymphoid_composition"]["capped_reference_donor_mean_THE_PRIOR"][k])
                    for k in LYMPH), 4)),
                "L1_lymphoid_to_methylation": (None if ml is None else round(sum(
                    abs(ml[k] - truth[k]) for k in LYMPH), 4)),
                "cohort_mean_all_types": {k: round(float(v), 4) for k, v in mu.items()},
            }
        report["cohorts"][c] = {"methylation_lymphoid_mean": truth,
                                "methods": dict(sorted(rows.items(),
                                                       key=lambda kv: kv[1]["L1_to_reference_prior_all_types"]))}
    Path(a.out).write_text(json.dumps(report, indent=2))

    lc = report["lymphoid_composition"]
    print("lymphoid T/NK/B -- atlas pooled", lc["atlas_pooled_cells"])
    print("                   capped reference donor-mean (Bisque's prior)",
          lc["capped_reference_donor_mean_THE_PRIOR"])
    for c, v in report["cohorts"].items():
        b = v["methods"]["bisque"]
        print(f"\n{c.upper()}: bisque lymphoid {b['lymphoid']}  | L1 to prior {b['L1_lymphoid_to_prior']}"
              f"  vs to methylation {b['L1_lymphoid_to_methylation']}")
        print("  L1 of cohort mean from the reference prior, all 8 types:")
        print("   " + "  ".join(f"{m} {r['L1_to_reference_prior_all_types']:.3f}"
                                for m, r in v["methods"].items()))
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
