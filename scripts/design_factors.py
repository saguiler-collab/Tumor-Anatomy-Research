#!/usr/bin/env python3
"""
design_factors.py — which design choices actually buy accuracy, and which buy nothing.

THE QUESTION, AND WHY THE OBVIOUS VERSION OF IT IS UNANSWERABLE
---------------------------------------------------------------
"Why did MuSiC win?" cannot be answered from this study, and pretending otherwise would be
the post-hoc storytelling the protocol exists to prevent. MuSiC leads ACS by **one
constraint-tumour pair out of 57** (1.000 against 0.9846); six methods have confidence
intervals containing its score; ACS takes only 66 distinct values on a weighted denominator
of 65. Attributing a margin that thin to a design feature is a story, not a finding.

What IS answerable is better posed and better powered: **which design property, isolated,
changes accuracy?** Several pairs of methods in this panel differ in exactly one property
while sharing the reference, the gene space, the harness and the cell-size handling. Those
are near-controlled contrasts, and they are the only causal-ish evidence available here.

THE CONTRASTS, AND WHY EACH IS CLEAN
------------------------------------
  MuSiC - NNLS          cross-donor variance weighting. The project's own invariant states
                        that MuSiC without cross-subject variance IS NNLS, so this pair
                        differs in that weighting and nothing else.
  Elastic Net - NNLS    L1/L2 regularisation on the same least-squares objective.
  CIBERSORTx B - SVR    batch correction. CIBERSORTx's base solver IS this project's SVR,
                        so B-mode is the correction layer alone.
  CIBERSORTx S - B      S-mode's per-cell-type adjustment against B-mode's.
  Hierarchical - Bayes   a donor-level hierarchy over the same likelihood.
  SCDC ENSEMBLE - SCDC  ensembling across references. With one reference this MUST be
                        exactly zero, which makes it a null control on the contrast method
                        itself: any non-zero value here means the machinery is reading noise.

WHAT THIS IS NOT
----------------
Not a selection procedure. ACS is the object under test in this study and may not choose a
method; the frozen method was chosen by the outcome-blind pseudobulk benchmark. Nothing here
reorders the leaderboard, and a property that looks good here is a hypothesis for the second
tissue, not a conclusion.

Not powered for a regression across methods. Fourteen comparable rows contain four exact
ties and several near-duplicate pairs, so the effective number of independent designs is
closer to nine. Per-attribute correlations are reported for completeness and should be read
as descriptive only.

    python scripts/design_factors.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

#: Near-controlled contrasts: (label, method_with, method_without, the property isolated).
CONTRASTS = [
    ("cross-donor variance weighting", "music", "nnls",
     "MuSiC weights genes by consistency across the 88 reference donors; without that "
     "weighting it reduces to NNLS, which is the comparison row."),
    ("L1/L2 regularisation", "elastic_net", "nnls",
     "Same least-squares objective, penalised."),
    ("batch correction (B-mode)", "cibersortx", "svr",
     "CIBERSORTx's base solver is this project's nu-SVR, so the difference is the "
     "batch-correction layer alone."),
    ("S-mode vs B-mode", "cibersortx_smode", "cibersortx",
     "Per-cell-type expression adjustment against whole-matrix batch correction."),
    ("donor-level hierarchy", "bayesian_hierarchical", "bayesian",
     "The same likelihood with a donor-level prior added."),
    ("ensembling across references", "scdc_ensemble", "scdc",
     "NULL CONTROL. With one reference, SCDC ENSEMBLE reduces exactly to SCDC, so this "
     "contrast must be 0.0000 on every outcome. A non-zero value means the machinery "
     "below is reading noise rather than a design property."),
]

#: Design attributes, coded from each method's PUBLISHED description before any score was
#: read. Coding notes live beside each key so a reader can disagree with a specific cell.
ATTRIBUTES = {
    "uses_individual_cells":  {"music", "scdc", "scdc_ensemble", "bisque", "bayesprism", "dwls"},
    "uses_cross_donor_variance": {"music", "scdc", "scdc_ensemble"},
    "regularised":            {"elastic_net", "svr", "cibersortx", "cibersortx_smode", "dwls"},
    "models_other_compartment": {"epic", "quantiseq"},
    "probabilistic":          {"bayesprism", "bayesian", "bayesian_hierarchical"},
    "platform_correction":    {"cibersortx", "cibersortx_smode", "bisque"},
    "weights_genes_by_informativeness": {"music", "dwls", "epic", "scdc", "scdc_ensemble"},
    "brings_own_signature":   {"quantiseq"},
}

OUTCOMES = ["acs", "mae_primary", "tumour_mae", "tumour_abs_bias", "tumour_r", "worst_niche_mae"]
#: Direction in which a LOWER number is better.
LOWER_IS_BETTER = {"mae_primary", "tumour_mae", "tumour_abs_bias", "worst_niche_mae"}


def build_table(R: Path) -> pd.DataFrame:
    lb = pd.read_csv(R / "anatomic/acs_leaderboard.csv").set_index("method")
    bench = pd.read_csv(R / "benchmark/benchmark_summary.csv").set_index("method")
    pt = pd.read_csv(R / "benchmark/benchmark_per_type.csv")
    pn = pd.read_csv(R / "benchmark/benchmark_per_niche.csv")

    tum = pt[pt.cell_type == "Tumor"].set_index("method")
    worst = pn.pivot_table(index="method", columns="niche", values="mae_primary").max(axis=1)

    rows = []
    for m in lb.index:
        if bool(lb.loc[m, "is_control"]):
            continue
        rows.append({
            "method": m,
            "acs": float(lb.loc[m, "acs"]),
            "comparable": bool(lb.loc[m, "comparable"]),
            "mae_primary": float(bench.loc[m, "mae_primary"]) if m in bench.index else np.nan,
            "tumour_mae": float(tum.loc[m, "mae"]) if m in tum.index else np.nan,
            "tumour_abs_bias": abs(float(tum.loc[m, "bias"])) if m in tum.index else np.nan,
            "tumour_r": float(tum.loc[m, "pearson_r"]) if m in tum.index else np.nan,
            "worst_niche_mae": float(worst.get(m, np.nan)),
        })
    return pd.DataFrame(rows).set_index("method")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out", default="results/design_factors.json")
    args = ap.parse_args()

    R = ROOT / args.results_dir
    t = build_table(R)

    print("=== near-controlled contrasts: one design property, everything else shared ===\n")
    print(f"{'property':34s} {'ACS':>8s} {'MAE':>9s} {'TumMAE':>9s} {'|TumBias|':>10s} "
          f"{'TumR':>8s} {'WorstMAE':>9s}")
    print("-" * 92)
    contrast_rows = []
    for label, a, b, note in CONTRASTS:
        if a not in t.index or b not in t.index:
            continue
        d = {o: float(t.loc[a, o] - t.loc[b, o]) for o in OUTCOMES}
        contrast_rows.append({"property": label, "with": a, "without": b,
                              "note": note, "deltas": {k: round(v, 5) for k, v in d.items()},
                              "n_outcomes_improved": sum(
                                  1 for o in OUTCOMES
                                  if (d[o] < 0) == (o in LOWER_IS_BETTER) and abs(d[o]) > 1e-9)})
        print(f"{label:34s} {d['acs']:+8.4f} {d['mae_primary']:+9.4f} "
              f"{d['tumour_mae']:+9.4f} {d['tumour_abs_bias']:+10.4f} "
              f"{d['tumour_r']:+8.4f} {d['worst_niche_mae']:+9.4f}")
    print("\nnegative is better for MAE, |bias| and worst-niche; positive is better for "
          "ACS and TumR.")

    null = next((c for c in contrast_rows if c["with"] == "scdc_ensemble"), None)
    if null:
        mx = max(abs(v) for v in null["deltas"].values())
        print(f"\nNULL CONTROL (ensembling with one reference): max |delta| = {mx:.2e} "
              f"— {'PASSES, exactly zero as it must be' if mx < 1e-9 else 'FAILS'}")

    # Descriptive attribute association. Low power, stated as such.
    comp = t[t["comparable"]]
    print(f"\n=== attribute association across the {len(comp)} comparable methods "
          f"(DESCRIPTIVE ONLY) ===\n")
    print(f"{'attribute':36s} {'n':>3s} {'mean MAE with':>14s} {'without':>9s} {'delta':>8s}")
    print("-" * 74)
    attr_rows = []
    for attr, members in ATTRIBUTES.items():
        have = [m for m in comp.index if m in members]
        lack = [m for m in comp.index if m not in members]
        if not have or not lack:
            continue
        a_mae, b_mae = comp.loc[have, "mae_primary"].mean(), comp.loc[lack, "mae_primary"].mean()
        attr_rows.append({"attribute": attr, "n_with": len(have), "methods_with": have,
                          "mean_mae_with": round(float(a_mae), 4),
                          "mean_mae_without": round(float(b_mae), 4),
                          "delta_mae": round(float(a_mae - b_mae), 4),
                          "mean_acs_with": round(float(comp.loc[have, "acs"].mean()), 4),
                          "mean_acs_without": round(float(comp.loc[lack, "acs"].mean()), 4)})
        print(f"{attr:36s} {len(have):3d} {a_mae:14.4f} {b_mae:9.4f} {a_mae - b_mae:+8.4f}")

    report = {
        "what_this_is": ("Which design properties change accuracy, measured by "
                         "near-controlled contrasts between methods differing in one "
                         "property. NOT a selection procedure: ACS is the object under test "
                         "and may not choose a method, and nothing here reorders the "
                         "leaderboard."),
        "why_not_why_did_music_win": (
            "MuSiC leads ACS by one constraint-tumour pair out of 57 (1.000 vs 0.9846); six "
            "methods' CIs contain its score; ACS takes only 66 distinct values. A margin "
            "that thin cannot be attributed to a design feature, so the question is posed "
            "as 'which property, isolated, changes accuracy' instead."),
        "power": (f"{len(comp)} comparable methods, containing four exact ACS ties and "
                  "several near-duplicate pairs, so the effective number of independent "
                  "designs is closer to nine. Attribute associations are descriptive."),
        "contrasts": contrast_rows,
        "attribute_association_descriptive_only": attr_rows,
        "per_method": json.loads(t.round(5).reset_index().to_json(orient="records")),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
