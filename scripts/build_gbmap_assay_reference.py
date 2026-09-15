#!/usr/bin/env python3
"""
build_gbmap_assay_reference.py — separate ATLAS from PLATFORM (OPEN_DEFECTS D14).

THE QUESTION THIS SETTLES
-------------------------
D14 found that the ACS ordering under GBmap is reproduced by neither of two independent
references, while those two reproduce each other (rho 0.710 between them, 0.138 and 0.038
against GBmap). Two readings survive that, and the design could not separate them:

    "the ordering depends on the ATLAS"    vs    "it depends on the atlas's PLATFORM"

Both alternatives are Smart-seq2. GBmap is 87% 10x. So the two explanations are perfectly
confounded in every comparison made so far.

**GBmap contains its own Smart-seq2 subset: 9,275 cells across 24 donors.** Building a
reference from those, and another from its 10x cells alone, completes a 2x2 that breaks the
confound:

                        same atlas              different atlas
    same platform       --                      GBmap_SS2 vs Neftel   (both Smart-seq2)
    diff platform       GBmap_10x vs GBmap_SS2  GBmap_10x vs Neftel   (already: +0.038)

    If PLATFORM drives it:  GBmap_10x vs GBmap_SS2 is LOW,  GBmap_SS2 vs Neftel is HIGH.
    If ATLAS drives it:     GBmap_10x vs GBmap_SS2 is HIGH, GBmap_SS2 vs Neftel is LOW.

Two cells of that table are new and this script builds what they need.

WHY THE FILTER GOES INSIDE THE LOADER
-------------------------------------
`build_from_h5ad` gained a `cell_filter` applied BEFORE `_balanced_cell_sample`, and the order
is the point. The sampler caps cells per (donor, cell type) knowing nothing about assay, so
filtering afterwards would leave whatever 2.7% survived by chance — a handful per type and
often none per (donor, type). Filtering first lets the sampler balance within the chosen cells.

    python scripts/build_gbmap_assay_reference.py --assay smartseq2
    python scripts/build_gbmap_assay_reference.py --assay tenx
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.data.reference import build_from_h5ad                  # noqa: E402

ASSAYS = {
    "smartseq2": (lambda s: s.str.contains("Smart", case=False, na=False),
                  "GBmap's own Smart-seq2 cells"),
    "tenx": (lambda s: s.str.startswith("10x"), "GBmap's 10x cells only"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assay", choices=sorted(ASSAYS), required=True)
    args = ap.parse_args()
    pred, why = ASSAYS[args.assay]
    out = ROOT / f"data/reference/gbmap_{args.assay}"

    from ivygap.data.load_ivygap import load_cached
    expr, _ = load_cached()
    print(f"building a GBmap reference from {why}")
    print(f"  restricted to the {len(expr.index):,} genes the Ivy GAP bulk carries")

    def cell_filter(obs_raw: pd.DataFrame) -> np.ndarray:
        if "assay" not in obs_raw.columns:
            raise KeyError("the atlas has no `assay` column")
        return pred(obs_raw["assay"].astype(str)).to_numpy()

    ref, cells, meta = build_from_h5ad(
        config.REFERENCE_DIR / "gbmap_core.h5ad",
        restrict_to_genes=expr.index, export=False, cell_filter=cell_filter)

    n_don = meta.groupby("cell_type")["donor"].nunique().to_dict()
    n_cells = meta["cell_type"].value_counts().to_dict()
    print(f"\n=== {ref.profile.shape[0]:,} genes x {len(ref.cell_types)} types ===")
    print(f"{'type':24s} {'cells':>7s} {'donors':>7s} {'cell_size':>12s}")
    for t in ref.cell_types:
        print(f"{t:24s} {n_cells.get(t, 0):7d} {n_don.get(t, 0):7d} "
              f"{ref.cell_size.get(t, float('nan')):12,.0f}")
    spread = float(ref.cell_size.max() / ref.cell_size.min())
    print(f"\ncell_size spread max/min: {spread:.4f}")

    out.mkdir(parents=True, exist_ok=True)
    ref.profile.to_csv(out / "profile.csv")
    ref.sigma.to_csv(out / "sigma.csv")
    ref.cell_size.to_csv(out / "cell_size.csv", header=["cell_size"])
    (out / "provenance.json").write_text(json.dumps({
        "what_this_is": f"A GBmap reference built from {why}, to separate ATLAS from "
                        f"PLATFORM in docs/OPEN_DEFECTS.md D14.",
        "assay_subset": args.assay,
        "n_genes": int(ref.profile.shape[0]),
        "types": list(ref.cell_types),
        "cells_per_type": {k: int(v) for k, v in n_cells.items()},
        "donors_per_type": {k: int(v) for k, v in n_don.items()},
        "n_donors_total": int(meta["donor"].nunique()),
        "donors": sorted(meta["donor"].astype(str).unique().tolist()),
        "cell_size": {k: round(float(v), 1) for k, v in ref.cell_size.items()},
        "cell_size_spread_max_over_min": round(spread, 4),
        "filter_applied_before_subsampling": True,
        "why_that_order": ("_balanced_cell_sample caps per (donor, cell type) knowing nothing "
                           "about assay; filtering afterwards would leave whatever survived by "
                           "chance, a handful per type and often none per (donor, type)."),
    }, indent=2))
    print(f"wrote {out.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
