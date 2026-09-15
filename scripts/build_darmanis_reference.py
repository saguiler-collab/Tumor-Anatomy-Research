#!/usr/bin/env python3
"""
build_darmanis_reference.py — a SECOND deconvolution reference, from Darmanis 2017.

WHY A SECOND REFERENCE
----------------------
Two things in this project need one and cannot have it from GBmap alone.

`SCDC ENSEMBLE` is **degenerate** (D2): with one reference there is nothing to weight across,
so it reduces exactly to SCDC and is reported as degraded rather than as an ENSEMBLE result.
And nothing currently distinguishes "this is how these methods behave" from "this is how these
methods behave *on GBmap*" — every number in the study rests on one atlas.

Darmanis 2017 (GSE84465) is the best candidate available, and the reason is slightly
counter-intuitive: it is **heavily FACS-panned**, which destroys it for composition (see
`scripts/darmanis_constraint_check.py`) and is **irrelevant to a reference**. A reference needs
per-cell-type expression *profiles*, not per-region abundances, and sorting helps by enriching
rare types.

WHAT IT COVERS, AND WHAT IT DOES NOT
------------------------------------
Five of this project's eight roster types, across **four donors** — so cross-donor variance is
estimable and MuSiC and SCDC can use it, which Albiach (one donor) cannot support.

    Tumor                 <- Neoplastic
    Macrophage_Microglia  <- Immune cell
    Oligodendrocyte       <- Oligodendrocyte
    Endothelial           <- Vascular
    Astrocyte             <- Astocyte      [sic: GEO's own spelling]

**T_cell, NK_cell and B_cell are NOT recoverable.** Darmanis pools every lymphoid and myeloid
cell into one `Immune cell` label. This reference therefore supports a **declared 5-type
sub-roster only**, and cannot be swapped in for GBmap on the full eight. Mapping `Immune cell`
onto Macrophage_Microglia is itself an approximation — that label contains the T cells too, and
in Darmanis's own data they are a minority of it, but not zero.

`OPC` and `Neuron` are EXCLUDED rather than forced. Folding OPC into Oligodendrocyte would
inflate exactly the type C2 is about, and the roster has no neuron column at all.

NORMALISATION — THE D12 LESSON, APPLIED
---------------------------------------
`cell_size` is computed from each cell's **raw** library size, captured before any
normalisation. D12 records what happens otherwise: `run_benchmark` rebuilt a reference from an
already-normalised matrix without passing `cell_totals`, `build_reference` fell back to
`expression.sum(axis=0)` = 1e6 per cell, and the cell-size conversion became the identity for
the entire study. This script takes two streaming passes specifically so the raw totals are
captured first, and reports the resulting spread so a reader can see it is not 1.0.

Streaming rather than loading: the matrix is 23,465 x 3,589, and accumulating per-group sums
gene-row by gene-row keeps peak memory in megabytes on a machine with under a gigabyte free.

    python scripts/build_darmanis_reference.py
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

COUNTS = ROOT / "data/raw/darmanis_2017/GSE84465_GBM_All_data.csv.gz"
META = ROOT / "data/raw/darmanis_2017/cell_metadata.csv"
OUT = ROOT / "data/reference/darmanis_2017"

#: Darmanis label -> this project's roster. Declared before any expression is read.
ROSTER = {
    "Neoplastic": "Tumor",
    "Immune cell": "Macrophage_Microglia",
    "Oligodendrocyte": "Oligodendrocyte",
    "Vascular": "Endothelial",
    "Astocyte": "Astrocyte",            # GEO's own spelling
}
EXCLUDED = {"OPC": "no roster column; folding it into Oligodendrocyte would inflate the type "
                   "C2 is about",
            "Neuron": "the roster has no neuron column"}


def cell_ids_from_meta() -> pd.DataFrame:
    m = pd.read_csv(META).rename(columns={"patient id": "donor", "cell type": "label",
                                          "plate id": "plate"})
    m["cell_id"] = m["plate"].astype(str) + "." + m["well"].astype(str)
    m["roster"] = m["label"].map(ROSTER)
    return m.set_index("cell_id")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--restrict-to-bulk", action="store_true", default=True,
                    help="keep only genes the Ivy GAP bulk carries, as run_all does")
    args = ap.parse_args()

    for p in (COUNTS, META):
        if not p.exists():
            print(f"BLOCKED: {p} not found. See docs/EXTERNAL_ACTIONS.md item 15.")
            return 2

    meta = cell_ids_from_meta()
    print(f"metadata: {len(meta)} cells, {meta['donor'].nunique()} donors")
    print(f"  mapped onto the roster: {meta['roster'].notna().sum()}")
    for lab, why in EXCLUDED.items():
        n = int((meta['label'] == lab).sum())
        print(f"  EXCLUDED {lab} ({n} cells): {why}")

    keep_genes = None
    if args.restrict_to_bulk:
        from ivygap.data.load_ivygap import load_cached
        expr, _ = load_cached()
        keep_genes = set(expr.index)
        print(f"  restricting to the {len(keep_genes):,} genes the Ivy GAP bulk carries")
    del_expr = None

    with gzip.open(COUNTS, "rt") as fh:
        header = fh.readline().split()
        cols = [c.strip('"') for c in header]
    print(f"matrix columns: {len(cols)}")
    matched = [c for c in cols if c in meta.index]
    print(f"  matched to metadata: {len(matched)}")
    if len(matched) < 0.9 * len(cols):
        print(f"BLOCKED: only {len(matched)} of {len(cols)} columns matched the metadata. "
              f"Nothing is imputed. Example column {cols[0]!r}; example metadata id "
              f"{meta.index[0]!r}.")
        return 2

    col_idx = np.array([i for i, c in enumerate(cols) if c in meta.index])
    col_ids = [cols[i] for i in col_idx]
    info = meta.loc[col_ids]
    usable = info["roster"].notna().to_numpy()

    # --- pass 1: raw library size per cell, BEFORE any normalisation ---------
    print("\npass 1/2: raw library sizes ...")
    totals = np.zeros(len(col_idx), dtype="float64")
    n_genes_seen = 0
    with gzip.open(COUNTS, "rt") as fh:
        fh.readline()
        for line in fh:
            parts = line.split()
            gene = parts[0].strip('"')
            if keep_genes is not None and gene not in keep_genes:
                continue
            v = np.asarray(parts[1:], dtype="float64")[col_idx]
            totals += v
            n_genes_seen += 1
    print(f"  {n_genes_seen:,} genes kept; median raw library size "
          f"{np.median(totals[usable]):,.0f}")

    # --- pass 2: per (donor, type) sums of CPM-normalised values -------------
    donors = sorted(info.loc[usable, "donor"].unique())
    types = [t for t in config.CELL_TYPES if t in set(info["roster"].dropna())]
    print(f"\npass 2/2: per-donor, per-type profiles for {len(types)} types x "
          f"{len(donors)} donors ...")
    key = {(d, t): i for i, (d, t) in enumerate((d, t) for d in donors for t in types)}
    grp = np.full(len(col_idx), -1, dtype="int64")
    for i in range(len(col_idx)):
        if usable[i]:
            grp[i] = key[(info["donor"].iat[i], info["roster"].iat[i])]
    counts_per_group = np.bincount(grp[grp >= 0], minlength=len(key)).astype("float64")

    safe_tot = np.where(totals > 0, totals, np.nan)
    genes, sums = [], []
    with gzip.open(COUNTS, "rt") as fh:
        fh.readline()
        for line in fh:
            parts = line.split()
            gene = parts[0].strip('"')
            if keep_genes is not None and gene not in keep_genes:
                continue
            v = np.asarray(parts[1:], dtype="float64")[col_idx]
            cpm = np.nan_to_num(v / safe_tot * 1e6)
            acc = np.zeros(len(key), dtype="float64")
            np.add.at(acc, grp[grp >= 0], cpm[grp >= 0])
            genes.append(gene); sums.append(acc)
    S = np.vstack(sums)                                 # genes x (donor,type)
    print(f"  accumulated {S.shape[0]:,} genes x {S.shape[1]} (donor,type) groups")

    with np.errstate(invalid="ignore", divide="ignore"):
        means = S / np.where(counts_per_group > 0, counts_per_group, np.nan)

    # profile = mean over DONORS of each donor's mean, so a donor contributing many cells
    # does not dominate. sigma = variance ACROSS donor means, which is the cross-donor
    # variance MuSiC and SCDC weight by.
    profile = np.zeros((S.shape[0], len(types)))
    sigma = np.zeros((S.shape[0], len(types)))
    n_don = {}
    for j, t in enumerate(types):
        idx = [key[(d, t)] for d in donors if counts_per_group[key[(d, t)]] > 0]
        n_don[t] = len(idx)
        block = means[:, idx]
        profile[:, j] = np.nanmean(block, axis=1)
        sigma[:, j] = np.nanvar(block, axis=1, ddof=1) if len(idx) > 1 else 0.0

    prof = pd.DataFrame(profile, index=genes, columns=types)
    sig = pd.DataFrame(sigma, index=genes, columns=types)

    cs = {}
    for t in types:
        sel = usable & (info["roster"].to_numpy() == t)
        cs[t] = float(totals[sel].mean())
    cell_size = pd.Series(cs, index=types)

    OUT.mkdir(parents=True, exist_ok=True)
    prof.to_csv(OUT / "profile.csv")
    sig.to_csv(OUT / "sigma.csv")
    cell_size.to_csv(OUT / "cell_size.csv", header=["cell_size"])

    prov = {
        "what_this_is": ("A SECOND deconvolution reference built from Darmanis 2017 "
                         "(GSE84465). Covers 5 of the 8 roster types across 4 donors, so "
                         "cross-donor variance is estimable. Supports a DECLARED 5-type "
                         "sub-roster only and cannot replace GBmap on the full eight."),
        "source": "GSE84465, Darmanis et al. 2017",
        "counts_sha256": config.sha256_file(COUNTS),
        "n_cells_in_matrix": int(len(cols)),
        "n_cells_matched_to_metadata": int(len(col_idx)),
        "n_cells_used": int(usable.sum()),
        "n_genes": int(prof.shape[0]),
        "restricted_to_ivygap_bulk_genes": bool(args.restrict_to_bulk),
        "roster_mapping": ROSTER,
        "excluded_labels": EXCLUDED,
        "types": types,
        "donors": donors,
        "n_donors_per_type": n_don,
        "cells_per_type": {t: int((usable & (info["roster"].to_numpy() == t)).sum())
                           for t in types},
        "cell_size": {k: round(v, 1) for k, v in cs.items()},
        "cell_size_spread_max_over_min": round(float(max(cs.values()) / min(cs.values())), 4),
        "cell_size_note": ("computed from RAW library sizes captured before normalisation. "
                           "D12 records what happens otherwise: a reference rebuilt from an "
                           "already-normalised matrix gets cell_size 1e6 for every type and "
                           "the conversion becomes the identity."),
        "known_limits": [
            "T_cell, NK_cell and B_cell are not recoverable — Darmanis pools all lymphoid and "
            "myeloid cells into one 'Immune cell' label, so mapping it onto "
            "Macrophage_Microglia is an approximation that includes T cells.",
            "Astrocyte (88 cells), Oligodendrocyte (85) and Vascular (51) are thin, and "
            "thinner per donor.",
            "Smart-seq2, where GBmap is 87% 10x. A platform confound, and also the regime "
            "CIBERSORTx's S-mode exists for.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2))

    print(f"\n=== reference built: {prof.shape[0]:,} genes x {len(types)} types ===")
    print(f"{'type':24s} {'cells':>7s} {'donors':>7s} {'cell_size':>12s}")
    for t in types:
        print(f"{t:24s} {prov['cells_per_type'][t]:7d} {n_don[t]:7d} {cs[t]:12,.0f}")
    print(f"\ncell_size spread max/min: {prov['cell_size_spread_max_over_min']}  "
          f"({'NOT the identity — good' if prov['cell_size_spread_max_over_min'] > 1.05 else 'SUSPICIOUS, see D12'})")
    print(f"\nwrote {OUT.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
