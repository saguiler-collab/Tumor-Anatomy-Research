#!/usr/bin/env python3
"""
build_gbmap_linear_reference.py — the same atlas, the same cells, LINEAR expression.

Why this exists (docs/OPEN_DEFECTS.md D16 and D14)
--------------------------------------------------
`gbmap_core.h5ad` carries two matrices. The one the pipeline reads, `X`, is not linear
expression: it is `log1p(counts * s_i)` with one size factor per cell. That was verified to
float32 precision — the implied `s_i` is constant within a cell to 3.7e-07, and
`corr(expm1(X), raw counts)` within a cell is 1.000000. `raw/X` holds the genuine integer
counts.

Averaging log values produces a signature in log space while the Ivy GAP bulk is linear
FPKM, which is not the mixing model any of the fifteen methods assumes.

It also puts D14's conclusion in question. D14 attributed the collapse of the ACS ordering to
the ATLAS, on a 2x2 over atlas and platform. But every agreeing arm of that 2x2 was log
against log (GBmap vs GBmap), and every disagreeing arm was log against linear — the
Darmanis reference was built from raw counts and the Neftel one from inverted TPM. So "atlas"
and "expression space" are confounded in exactly the pattern D14 read as atlas.

This reference breaks that confound. It changes ONE thing relative to the published GBmap
reference: `matrix="raw/X"` instead of `"X"`. Same atlas, same annotation column, same donor
column, same roster mapping, same seed, same balanced sampler, same gene restriction. So

    GBmap_linear vs Neftel / Darmanis   isolates ATLAS with expression space held linear
    GBmap_linear vs GBmap (published)   isolates EXPRESSION SPACE with the atlas held fixed

and the pair of them answers which one D14 was actually measuring.

    python scripts/build_gbmap_linear_reference.py

Reads nothing that later selection depends on, and writes only under
`data/reference/gbmap_linear/`.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.data.reference import build_from_h5ad                  # noqa: E402


def main() -> int:
    out = ROOT / "data/reference/gbmap_linear"
    from ivygap.data.load_ivygap import load_cached
    expr, _ = load_cached()
    print("building a GBmap reference from raw/X — the genuine integer counts")
    print(f"  restricted to the {len(expr.index):,} genes the Ivy GAP bulk carries")

    ref, cells, meta = build_from_h5ad(
        config.REFERENCE_DIR / "gbmap_core.h5ad",
        restrict_to_genes=expr.index, matrix="raw/X", export=False)

    n_don = meta.groupby("cell_type")["donor"].nunique().to_dict()
    n_cells = meta["cell_type"].value_counts().to_dict()
    print(f"\n=== {ref.profile.shape[0]:,} genes x {len(ref.cell_types)} types ===")
    print(f"{'type':24s} {'cells':>7s} {'donors':>7s} {'cell_size':>14s}")
    for t in ref.cell_types:
        print(f"{t:24s} {n_cells.get(t, 0):7d} {n_don.get(t, 0):7d} "
              f"{ref.cell_size.get(t, float('nan')):14,.0f}")
    spread = float(ref.cell_size.max() / ref.cell_size.min())
    print(f"\ncell_size spread max/min: {spread:.4f}")

    out.mkdir(parents=True, exist_ok=True)
    ref.profile.to_csv(out / "profile.csv")
    ref.sigma.to_csv(out / "sigma.csv")
    ref.cell_size.to_csv(out / "cell_size.csv", header=["cell_size"])
    (out / "provenance.json").write_text(json.dumps({
        "what_this_is": "The GBmap atlas read from raw/X (integer counts) rather than X "
                        "(log1p of normalised counts), to break the atlas/expression-space "
                        "confound in docs/OPEN_DEFECTS.md D14 and to measure D16.",
        "matrix": "raw/X",
        "the_one_thing_changed": "matrix='raw/X'; every other argument matches the "
                                 "published GBmap reference build.",
        "n_genes": int(ref.profile.shape[0]),
        "types": list(ref.cell_types),
        "cells_per_type": {k: int(v) for k, v in n_cells.items()},
        "donors_per_type": {k: int(v) for k, v in n_don.items()},
        "n_donors_total": int(meta["donor"].nunique()),
        "donors": sorted(meta["donor"].astype(str).unique().tolist()),
        "cell_size": {k: round(float(v), 1) for k, v in ref.cell_size.items()},
        "cell_size_spread_max_over_min": round(spread, 4),
        "caveat_cell_size_is_platform_confounded": (
            "GBmap pools 10x and Smart-seq2. Smart-seq2 cells carry ~100x the read count "
            "of 10x cells (measured: Tumor 785,065 vs 8,301), and the platform mix differs "
            "sharply by type, so a cell_size vector pooled over both measures platform, not "
            "mRNA content. See results/gbmap_cell_size_by_space.json; use the 10x-only "
            "vector for any mRNA-content claim."),
    }, indent=2))
    print(f"wrote {out.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
