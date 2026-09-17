#!/usr/bin/env python3
"""
build_tcga_bulk.py — assemble the TCGA-GBM bulk matrix for YARDSTICK 2 (ABSOLUTE purity).

WHY THIS EXISTS
---------------
`Anatomy_Test.md` names a second yardstick: ABSOLUTE tumour purity, measured from **DNA**, so
it shares no failure mode with any RNA method under test. `agreement.py` already evaluates it
and reports UNAVAILABLE. It needs data, not code — all methods run on TCGA-GBM, each method's
Tumor column correlated against purity.

THE TRAP, AND IT IS THE SAME ONE TWICE IN ONE DAY
--------------------------------------------------
The input is named `TCGA-GBM.star_counts.tsv`. It does not contain counts. Measured:

  * nonzero values run 0.584963 .. 19.1542, and 0.584963 = log2(1.5);
  * after `2^v - 1`, 67.9% of small values are EXACT integers and 1.7% are exact half-integers,
    i.e. the underlying quantity lies on a 0.5 grid.

So the file is **log2(x + 1)** with x count-like, and it is inverted here with `2**v - 1`
before anything else happens. Feeding it in as delivered would repeat `OPEN_DEFECTS.md` D16
exactly — a log-transformed matrix used as linear expression — in a study whose main finding
is that this kind of error is invisible to the evaluation.

The predecessor's own `bulk_expression_standardized.csv` is NOT used for the same reason: its
metadata records `"expression_scale_guess": "appears_log2_transformed"`.

The predecessor project is READ-ONLY from here. Nothing under it is written or modified.

    python scripts/build_tcga_bulk.py
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

PRED = Path("/Users/tatopro9130/Downloads/cancer_judging_machine-master 22"
            "/public/03_Reference_Free_TME/00_inputs")
STAR = PRED / "TCGA-GBM.star_counts.tsv"
GTF = PRED / "gencode.v36.annotation.gtf.gz"
OUT = config.PROCESSED_DIR / "tcga_gbm_bulk_cpm.csv.gz"
GRID = 0.5          # the underlying quantity's spacing, measured not assumed


def gene_map() -> dict[str, str]:
    """Ensembl gene id -> HUGO symbol, from the GENCODE annotation the counts were built on."""
    pat = re.compile(r'gene_id "([^"]+)".*?gene_name "([^"]+)"')
    m: dict[str, str] = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 2 and f[2] == "gene":
                g = pat.search(f[8])
                if g:
                    m[g.group(1)] = g.group(2)
    return m


def main() -> int:
    for p in (STAR, GTF):
        if not p.exists():
            print(f"BLOCKED: {p} not found."); return 2

    print(f"reading {STAR.name} ({STAR.stat().st_size / 1e6:.0f} MB) ...")
    df = pd.read_csv(STAR, sep="\t", index_col=0)
    print(f"  {df.shape[0]:,} gene rows x {df.shape[1]} samples")

    # --- GATE 1: is it really log2(x+1) with x on a 0.5 grid? ------------------------
    v = df.to_numpy(dtype="float64")
    nz = v[v > 0]
    lin = 2.0 ** nz - 1.0
    small = lin[lin < 200]
    off = np.abs(small / GRID - np.round(small / GRID))
    on_grid = float((off < 1e-6).mean())
    print(f"\nGATE 1 — scale: {on_grid:.4%} of small inverted values sit on a "
          f"{GRID} grid (min nonzero {nz.min():.6f} = log2({2**nz.min():.1f}))")
    if on_grid < 0.95:
        print("BLOCKED: the log2 hypothesis does not hold; do not guess the scale.")
        return 2

    expr = pd.DataFrame(2.0 ** v - 1.0, index=df.index, columns=df.columns).clip(lower=0.0)
    del df, v

    # --- GATE 2: gene identity, which is where this project has been bitten ----------
    gm = gene_map()
    base = [g.split(".")[0] for g in expr.index]
    sym = pd.Series([gm.get(g) or gm.get(b) for g, b in zip(expr.index, base)],
                    index=expr.index)
    mapped = sym.notna()
    print(f"GATE 2 — symbols: {mapped.sum():,} of {len(sym):,} gene rows mapped "
          f"({mapped.mean():.2%}) from {GTF.name}")
    if mapped.mean() < 0.90:
        print("BLOCKED: too few genes mapped; the annotation does not match the matrix.")
        return 2
    expr = expr[mapped.to_numpy()]
    expr.index = sym[mapped].to_numpy()
    n_before = len(expr)
    expr = expr.groupby(level=0).mean()          # predecessor's rule, stated
    print(f"  collapsed {n_before:,} -> {len(expr):,} symbols by arithmetic mean")

    # --- GATE 3: does it overlap the frozen signature at all? ------------------------
    sig = pd.read_csv(config.FROZEN_SIGNATURE_PATH, sep="\t", index_col=0)
    shared = expr.index.intersection(sig.index)
    print(f"GATE 3 — overlap with the frozen signature: {len(shared):,} of "
          f"{len(sig):,} signature genes")
    if len(shared) < config.MIN_GENES_SHARED:
        print(f"BLOCKED: only {len(shared)} shared genes."); return 2

    # CPM, on the LINEAR scale, which is the whole point of the inversion above.
    totals = expr.sum(axis=0).replace(0.0, np.nan)
    cpm = (expr.div(totals, axis=1) * 1e6).fillna(0.0)
    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        cpm.to_csv(fh, float_format="%.4f")

    prov = {
        "what_this_is": "TCGA-GBM bulk expression, linearised and CPM-normalised, for "
                        "yardstick 2 (ABSOLUTE purity).",
        "source": str(STAR), "annotation": str(GTF),
        "source_is_read_only": "the predecessor project is never written to",
        "scale_correction": {
            "as_delivered": "log2(x + 1) despite the filename saying counts",
            "evidence": f"min nonzero {float(nz.min()):.6f} = log2(1.5); after 2^v - 1, "
                        f"{on_grid:.2%} of values below 200 sit on a {GRID} grid",
            "applied": "2**v - 1, before any normalisation",
            "why_it_matters": "using it as delivered would repeat OPEN_DEFECTS D16 in a study "
                              "whose finding is that this error is invisible to the evaluation",
            "predecessor_standardized_file_not_used":
                "bulk_expression_standardized.csv records "
                "'expression_scale_guess: appears_log2_transformed'",
        },
        "n_gene_rows_in": int(len(sym)), "n_mapped": int(mapped.sum()),
        "n_symbols_after_collapse": int(len(expr)),
        "duplicate_symbol_rule": "arithmetic mean, matching the predecessor's documented rule",
        "n_samples": int(cpm.shape[1]),
        "n_shared_with_frozen_signature": int(len(shared)),
        "normalisation": "CPM on the linear scale; columns sum to 1e6",
        "output": str(OUT.relative_to(ROOT)),
    }
    (config.PROCESSED_DIR / "tcga_gbm_bulk_provenance.json").write_text(
        json.dumps(prov, indent=2))
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({cpm.shape[0]:,} genes x {cpm.shape[1]} samples)")
    print(f"column sums: min {cpm.sum(axis=0).min():,.0f} max {cpm.sum(axis=0).max():,.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
