#!/usr/bin/env python3
"""
mes_by_niche.py — is the mesenchymal program spatially concentrated in glioblastoma?

Executes addendum 3 of `prespecified/biological_failure_factors.md`, whose prediction — MES
highest in pseudopalisading cells around necrosis — was committed before this ran.

WHY THIS IS INTERPRETATION AND NOT VALIDATION
----------------------------------------------
Ivy GAP has no per-sample ground truth, so deconvolution error cannot be measured here at all.
What can be measured is whether the biological property that predicts error in GBM is
spatially concentrated. The link from niche to error is inferential, through the score.

THE DESIGN OBEYS THE PROJECT'S INVARIANTS
------------------------------------------
Samples are nested in tumours, so everything is collapsed to one value per (tumour, structure)
before aggregation. The null is a within-tumour permutation of structure labels, not a t-test
across samples, which would treat several blocks from one tumour as independent observations.

    python scripts/mes_by_niche.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

GBMD = (ROOT / "pipeline_packages " / " repos" / "GBMDeconvoluteR"
        / "GBMDeconvoluteR-main" / "data")
RDS = "Neftel_et_al_2019_four_state_neoplastic_markers.rds"
N_PERM = 10000


def main() -> int:
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    print(f"Ivy GAP anatomic samples: {len(anat)} "
          f"(H&E-labelled only; the 148 ISH samples are excluded as circular)")

    tsv = config.RESULTS_DIR / "_neftel_ivy.tsv"
    r = subprocess.run(["Rscript", "-e", f'''
        x <- readRDS("{GBMD / RDS}")
        write.table(data.frame(gene = unlist(x, use.names = FALSE),
                               state = rep(names(x), lengths(x))),
                    "{tsv}", sep = "\\t", row.names = FALSE, quote = FALSE)
    '''], capture_output=True, text=True)
    if r.returncode != 0:
        print("BLOCKED: could not read the marker object."); return 2
    mk = pd.read_csv(tsv, sep="\t"); tsv.unlink(missing_ok=True)

    # IDENTICAL score to GBM and LGG: within-sample percentile, MES minus the other states.
    sub = expr.loc[:, anat]
    pct = sub.rank(axis=0, pct=True)
    sc = {}
    for state, g in mk.groupby("state"):
        genes = [x for x in g["gene"] if x in pct.index]
        if len(genes) >= 5:
            sc[state] = pct.loc[genes].mean(axis=0)
        print(f"  {state:4s} {len(genes):3d} of {len(g):3d} markers present")
    s = pd.DataFrame(sc)
    score = s["MES"] - s[[c for c in s.columns if c != "MES"]].mean(axis=1)

    # COLLAPSE TO (tumour, structure) FIRST. Several blocks per tumour are not independent.
    d = pd.DataFrame({"score": score,
                      "tumor": meta.loc[anat, "patient_id"].astype(str),
                      "structure": meta.loc[anat, "structure"].astype(str)})
    cell = d.groupby(["tumor", "structure"])["score"].mean().reset_index()
    print(f"\ncollapsed {len(d)} samples -> {len(cell)} (tumour, structure) cells "
          f"across {cell['tumor'].nunique()} tumours")

    order = [x for x in config.PRIMARY_STRUCTURES if x in set(cell["structure"])]
    means = cell.groupby("structure")["score"].agg(["mean", "median", "size"]).reindex(order)
    print("\nMES score by structure, one value per (tumour, structure):")
    print(means.round(5).to_string())
    top = means["mean"].idxmax()
    print(f"\nHIGHEST: {top}  (pre-specified prediction: PAN)")

    # WITHIN-TUMOUR PERMUTATION NULL for "structure X is the maximum".
    rng = np.random.default_rng(config.RANDOM_SEED)
    obs_rank = {st: float(means.loc[st, "mean"]) for st in order}
    hits = {st: 0 for st in order}
    piv = cell.pivot(index="tumor", columns="structure", values="score").reindex(columns=order)
    arr = piv.to_numpy(dtype="float64")
    for _ in range(N_PERM):
        sh = np.array([rng.permutation(row) if np.isfinite(row).all()
                       else _shuffle_partial(row, rng) for row in arr])
        m = np.nanmean(sh, axis=0)
        hits[order[int(np.nanargmax(m))]] += 1
    print(f"\nwithin-tumour permutation null, {N_PERM:,} draws — how often each structure "
          f"comes out highest by chance:")
    for st in order:
        print(f"  {st:5s} {hits[st] / N_PERM:7.2%}"
              + ("   <- observed maximum" if st == top else ""))
    p = hits[top] / N_PERM
    print(f"\np(observed maximum is {top} by chance) = {p:.4f}")

    prediction_held = top == "PAN"
    print(f"\nPRE-SPECIFIED PREDICTION: {'HELD' if prediction_held else 'FAILED'} — "
          f"predicted PAN, observed {top}")

    # --- POST HOC, AND LABELLED AS SUCH ------------------------------------------------
    #
    # The pre-specification fixed the PREDICTION ("MES highest in PAN") but did NOT fix the
    # test statistic. I chose "which structure is the maximum", and it is a weak statistic:
    # with five structures the chance rate is ~20% before any signal, so PAN at 33% is barely
    # elevated. A monotone-trend or paired PAN-vs-LE test would be far more powerful.
    #
    # Switching to one of those NOW, having seen p = 0.33, would be choosing the test after
    # the result. So they are computed and reported as POST HOC, and no significance claim for
    # the spatial layer rests on them. This is the honest maximum, not the flattering one.
    post = {}
    if {"PAN", "LE"} <= set(order):
        both = piv[["LE", "PAN"]].dropna()
        diff = (both["PAN"] - both["LE"]).to_numpy()
        if len(diff) >= 4:
            cnt = sum(1 for _ in range(N_PERM)
                      if np.mean(rng.choice([-1, 1], len(diff)) * np.abs(diff)) >= np.mean(diff))
            post["paired_PAN_minus_LE"] = {
                "n_tumors": int(len(diff)), "mean_difference": round(float(np.mean(diff)), 5),
                "n_tumors_positive": int((diff > 0).sum()),
                "sign_flip_permutation_p": round(cnt / N_PERM, 4)}
            print(f"\nPOST HOC (statistic not pre-specified): paired PAN - LE across "
                  f"{len(diff)} tumours, mean {np.mean(diff):+.5f}, "
                  f"{int((diff > 0).sum())}/{len(diff)} positive, "
                  f"sign-flip p = {cnt / N_PERM:.4f}")
    ranks = np.arange(1, len(order) + 1)
    obs_t = float(np.corrcoef(ranks, [means.loc[st, "mean"] for st in order])[0, 1])
    post["monotone_trend_pearson_on_ranks"] = round(obs_t, 4)
    post["caveat"] = ("The statistic was NOT pre-specified. These are reported so a reader "
                      "sees the full picture, and no significance claim rests on them.")
    print(f"POST HOC: trend along LE->IT->CT->MVP->PAN, r = {obs_t:+.4f} on structure ranks")

    out = {"post_hoc_not_prespecified": post,
           "what_this_is": "Is the mesenchymal program spatially concentrated in GBM? "
                           "Interpretation only -- Ivy GAP has no per-sample ground truth, so "
                           "error itself cannot be measured here.",
           "prespecification": "prespecified/biological_failure_factors.md addendum 3",
           "predicted_maximum": "PAN", "observed_maximum": top,
           "prediction_held": bool(prediction_held),
           "n_samples": len(anat), "n_cells": int(len(cell)),
           "n_tumors": int(cell["tumor"].nunique()),
           "score": "MES - mean(AC, NPC, OPC), within-sample percentiles; identical to the "
                    "GBM and LGG score",
           "by_structure": {st: {"mean": round(float(means.loc[st, "mean"]), 5),
                                 "median": round(float(means.loc[st, "median"]), 5),
                                 "n_tumors": int(means.loc[st, "size"])} for st in order},
           "null": {"draws": N_PERM, "design": "within-tumour permutation of structure labels",
                    "p_observed_max_by_chance": round(p, 4),
                    "chance_rate_per_structure": {st: hits[st] / N_PERM for st in order}},
           "cannot_show": "that deconvolution error is higher in this niche; the link from "
                          "niche to error is inferential, through the score."}
    (config.RESULTS_DIR / "mes_by_niche.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/mes_by_niche.json")
    return 0


def _shuffle_partial(row, rng):
    """Permute only the observed entries of a tumour's row, leaving missing structures missing."""
    out = row.copy()
    idx = np.where(np.isfinite(row))[0]
    out[idx] = rng.permutation(row[idx])
    return out


if __name__ == "__main__":
    sys.exit(main())
