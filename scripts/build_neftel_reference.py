#!/usr/bin/env python3
"""
build_neftel_reference.py — a THIRD reference, from Neftel 2019's own annotations.

WHY A THIRD
-----------
`scripts/reference_sensitivity.py` found that swapping GBmap for Darmanis moves the ACS
ordering to a Spearman of **0.138** — 10 of 13 methods change rank, and the two Bayesian
methods go from near-bottom to top. That is the most consequential result in the study, and it
has an obvious alternative explanation: **Darmanis might simply be a bad reference.** Three of
its five types rest on 51 to 88 cells.

A third, independent reference decides between those readings. If GBmap versus Neftel also
disagrees, reference sensitivity is a general property of this cohort. If they agree, Darmanis
was the outlier and the finding is about thin cell types rather than about references.

WHAT NEFTEL PROVIDES
--------------------
The Broad Single Cell Portal files the user supplied — `GBM_Metadata` and
`IDHwtGBM.processed.SS2.logTPM.txt` — carry the **authors' own** `CellAssignment`, which is
what `docs/EXTERNAL_ACTIONS.md` item 9 was asking for. No derivation, no clustering, no
inferCNV, and no risk that our own labels biased a reference toward our own roster.

Adult cells only (5,742 of 7,930; the rest are paediatric and this study is adult GBM), four
types across **20 donors**:

    Tumor                 <- Malignant        4,916 cells, 20 donors
    Macrophage_Microglia  <- Macrophage         536 cells, 11 donors
    Oligodendrocyte       <- Oligodendrocyte    210 cells, 13 donors
    T_cell                <- T-cell              80 cells,  6 donors

**T_cell is the type Darmanis cannot give** — it pools all lymphoid and myeloid cells into one
label — so Neftel and Darmanis are complementary rather than redundant. Neither has
Endothelial or Astrocyte, so C3 and C4 are unscoreable on either.

TWO LIMITS THAT CANNOT BE ENGINEERED AWAY
-----------------------------------------
1. **`cell_size` is unrecoverable.** The published matrix is TPM, which is already normalised
   per cell, so per-cell mRNA content is gone at source. `cell_size` is therefore set uniform
   and **declared**, not estimated. This is the D12 situation arriving for an unavoidable
   reason rather than a fixable one, and it means this reference cannot drive a cell-size
   conversion. Since D12 shows the conversion is the identity in the main pipeline anyway, the
   arms stay comparable.
2. **`CrossSection` cannot test the constraints.** Only `MGH105` has one, its four sections are
   labelled A/B/C/D, and those are spatial slices rather than anatomic structures. There is no
   defensible mapping onto LE/IT/CT/MVP/PAN, so Neftel is a reference and not a second
   validation specimen.

The matrix is log2(TPM/10 + 1), so it is inverted to TPM as (2^E - 1) * 10 before averaging.
Averaging log values would give a geometric mean, which is not what a signature matrix is.

    python scripts/build_neftel_reference.py
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

NEF = ROOT / "pipeline_packages / repos/SCDC/neftel_data"
META = NEF / "GBM_Metadata"
MATRIX = NEF / "IDHwtGBM.processed.SS2.logTPM.txt"
OUT = ROOT / "data/reference/neftel_2019"

ROSTER = {
    "Malignant": "Tumor",
    "Macrophage": "Macrophage_Microglia",
    "Oligodendrocyte": "Oligodendrocyte",
    "T-cell": "T_cell",
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunk", type=int, default=2000, help="gene rows per read chunk")
    args = ap.parse_args()

    for p in (META, MATRIX):
        if not p.exists():
            print(f"BLOCKED: {p} not found.")
            return 2

    meta = pd.read_csv(META, sep="\t", skiprows=[1])
    meta = meta.set_index("NAME")
    adult = meta[(meta["GBMType"] == "Adult") & (meta["CellAssignment"].isin(ROSTER))].copy()
    adult["roster"] = adult["CellAssignment"].map(ROSTER)
    print(f"metadata: {len(meta)} cells; adult with a mapped type: {len(adult)}")
    print(f"  donors: {adult['Sample'].nunique()}")

    from ivygap.data.load_ivygap import load_cached
    expr, _ = load_cached()
    bulk_genes = set(expr.index)
    print(f"  restricting to the {len(bulk_genes):,} genes the Ivy GAP bulk carries")

    with open(MATRIX) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cells = header[1:]
    use = [c for c in cells if c in adult.index]
    print(f"matrix: {len(cells)} cells; usable: {len(use)}")
    if len(use) < 0.5 * len(adult):
        print(f"BLOCKED: only {len(use)} of {len(adult)} adult cells found in the matrix. "
              f"Nothing is imputed.")
        return 2

    info = adult.loc[use]
    donors = sorted(info["Sample"].unique())
    types = [t for t in config.CELL_TYPES if t in set(info["roster"])]
    key = {(d, t): i for i, (d, t) in enumerate((d, t) for d in donors for t in types)}
    grp = np.array([key[(info["Sample"].iat[i], info["roster"].iat[i])]
                    for i in range(len(use))])
    n_per_group = np.bincount(grp, minlength=len(key)).astype("float64")
    print(f"  {len(types)} types x {len(donors)} donors = {len(key)} groups; "
          f"{int((n_per_group > 0).sum())} non-empty")

    print(f"\nstreaming the matrix in {args.chunk}-row chunks "
          f"(log2(TPM/10+1) -> TPM before averaging) ...")
    genes, sums = [], []
    n_seen = 0
    reader = pd.read_csv(MATRIX, sep="\t", index_col=0, usecols=["GENE"] + use,
                         chunksize=args.chunk)
    for chunk in reader:
        n_seen += len(chunk)
        chunk = chunk[chunk.index.isin(bulk_genes)]
        if chunk.empty:
            continue
        chunk = chunk[use]
        tpm = (np.power(2.0, chunk.to_numpy(dtype="float64")) - 1.0) * 10.0
        np.clip(tpm, 0.0, None, out=tpm)
        acc = np.zeros((tpm.shape[0], len(key)), dtype="float64")
        for j in range(len(key)):
            sel = grp == j
            if sel.any():
                acc[:, j] = tpm[:, sel].sum(axis=1)
        genes.extend(chunk.index.tolist())
        sums.append(acc)
        if n_seen % 10000 < args.chunk:
            print(f"    {n_seen:,} rows read, {len(genes):,} kept")
    S = np.vstack(sums)
    print(f"  done: {n_seen:,} rows read, {S.shape[0]:,} genes kept")

    with np.errstate(invalid="ignore", divide="ignore"):
        means = S / np.where(n_per_group > 0, n_per_group, np.nan)

    profile = np.zeros((S.shape[0], len(types)))
    sigma = np.zeros((S.shape[0], len(types)))
    n_don = {}
    for j, t in enumerate(types):
        idx = [key[(d, t)] for d in donors if n_per_group[key[(d, t)]] > 0]
        n_don[t] = len(idx)
        block = means[:, idx]
        profile[:, j] = np.nanmean(block, axis=1)
        sigma[:, j] = np.nanvar(block, axis=1, ddof=1) if len(idx) > 1 else 0.0

    prof = pd.DataFrame(np.nan_to_num(profile), index=genes, columns=types)
    sig = pd.DataFrame(np.nan_to_num(sigma), index=genes, columns=types)
    # UNIFORM and DECLARED: the source matrix is TPM, so per-cell library size is gone.
    cell_size = pd.Series(1.0e6, index=types)

    OUT.mkdir(parents=True, exist_ok=True)
    prof.to_csv(OUT / "profile.csv")
    sig.to_csv(OUT / "sigma.csv")
    cell_size.to_csv(OUT / "cell_size.csv", header=["cell_size"])
    prov = {
        "what_this_is": ("A THIRD deconvolution reference, from Neftel 2019's OWN "
                         "CellAssignment labels (Broad Single Cell Portal). Adult cells only, "
                         "4 roster types, 20 donors. Built to decide whether the GBmap-vs-"
                         "Darmanis ordering disagreement is general or specific to Darmanis."),
        "source": "GSE131928 / Broad SCP: GBM_Metadata + IDHwtGBM.processed.SS2.logTPM.txt",
        "labels_are_the_authors_own": True,
        "n_cells_in_metadata": int(len(meta)),
        "n_cells_used": int(len(use)),
        "adult_only": True,
        "why_adult_only": "this study is adult glioblastoma; 2,188 paediatric cells excluded",
        "roster_mapping": ROSTER,
        "types": types, "donors": donors,
        "n_donors_per_type": n_don,
        "cells_per_type": {t: int((info["roster"] == t).sum()) for t in types},
        "n_genes": int(prof.shape[0]),
        "restricted_to_ivygap_bulk_genes": True,
        "log_inversion": "matrix is log2(TPM/10+1); inverted to TPM as (2^E - 1)*10 BEFORE "
                         "averaging, because averaging log values gives a geometric mean, "
                         "which is not what a signature matrix is",
        "cell_size": "UNIFORM 1e6, DECLARED not estimated",
        "cell_size_why": ("the published matrix is TPM, already normalised per cell, so "
                          "per-cell mRNA content is unrecoverable at source. This reference "
                          "cannot drive a cell-size conversion. D12 shows the conversion is "
                          "the identity in the main pipeline anyway, so the arms stay "
                          "comparable."),
        "known_limits": [
            "No Endothelial and no Astrocyte, so C3 and C4 are unscoreable on this reference.",
            "T_cell rests on 80 cells across 6 donors — thin, but it is the type Darmanis "
            "cannot give at all.",
            "CrossSection exists only for MGH105 and its A/B/C/D sections are spatial slices, "
            "not anatomic structures, so this is a reference and NOT a validation specimen.",
            "Smart-seq2, where GBmap is 87% 10x.",
        ],
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2))

    print(f"\n=== reference built: {prof.shape[0]:,} genes x {len(types)} types ===")
    print(f"{'type':24s} {'cells':>7s} {'donors':>7s}")
    for t in types:
        print(f"{t:24s} {prov['cells_per_type'][t]:7d} {n_don[t]:7d}")
    print(f"\nwrote {OUT.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
