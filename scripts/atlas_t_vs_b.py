"""
atlas_t_vs_b.py -- prespecified/atlas_t_vs_b.md: T against B counted cell by cell in an independent
single-cell atlas (Abdelfattah 2022, GSE182109), with no methylation and no deconvolution.

EXPLORATORY. C3 (T and NK together) is split per cell by the split already declared in
prespecified/abdelfattah_cluster_mapping.md: the sign of (mean T panel - mean NK panel), log1p-CPM.
Cells with no marker of either panel (difference exactly 0) are counted as unresolved; the
comparison uses strict T, and the declared rule's version (ties -> T) is reported beside it.
B = every C11 cell. Counts are recovered from the processed matrix with the reference builder's own
streaming code, and its integer-count identity must hold again.

One pass over the full matrix (about 20 minutes, CPU-bound): run it when no heavy job is running.

    python3 scripts/atlas_t_vs_b.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
from build_abdelfattah_reference import PANELS, SRC, matrix_path, per_cell_units, stream  # noqa: E402

OUT = config.RESULTS_DIR / "atlas_t_vs_b.json"
TNK, B, UNMAPPED = "C3", "C11", "C12"


def main() -> int:
    prov = json.loads((config.REFERENCE_DIR / "abdelfattah_2022" / "provenance.json").read_text())
    mp = prov["cluster_mapping"]
    if mp.get(TNK) != "T_cell" or mp.get(B) != "B_cell" or mp.get(UNMAPPED) is not None:
        print(f"BLOCKED: the recorded cluster mapping is not the one this rule names: {mp}"); return 2
    genes = [l.rstrip("\n").split("\t")[-1] for l in open(SRC / "genes.tsv")]
    cells = [l.strip() for l in open(SRC / "barcodes.tsv")]
    cl = pd.read_csv(SRC / "Cluster_GBM.txt", sep=",", skiprows=[1], header=None,
                     names=["NAME", "N2", "X", "Y", "Category"], low_memory=False)
    cl = cl[cl.NAME != "NAME"]
    assert list(cl.NAME) == cells, "barcodes.tsv and Cluster_GBM.txt disagree on cell order"
    clusters = cl.Category.to_numpy()
    gsm2s = dict(l.rstrip("\n").split("\t") for l in open(SRC / "gsm_to_sample.tsv"))
    patient = np.array(["-".join(gsm2s[n.split("_")[0]].split("-")[:2]) for n in cells])

    tgenes = [g for g in PANELS["T_cell"] if g in set(genes)]
    ngenes = [g for g in PANELS["NK_cell"] if g in set(genes)]
    panel = tgenes + ngenes
    gpos = {genes.index(g): k for k, g in enumerate(panel)}
    is_panel = np.zeros(len(genes), bool); is_panel[list(gpos)] = True
    n_cells = len(cells)
    lib = np.full(n_cells, np.nan); total = np.zeros(n_cells); maxdev = np.zeros(n_cells)
    counts = np.zeros((n_cells, len(panel)), np.float32)

    path = matrix_path(); print(f"streaming {path.name}", flush=True)
    it = stream(path); dims = next(it)
    assert dims[0] == len(genes) and dims[1] == n_cells, f"matrix dims {dims} vs files"
    carry, seen = None, 0

    def take(g, c, v):
        e = np.expm1(v)
        starts = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
        ids = c[starts]; units = np.minimum.reduceat(e, starts)
        cnt = e / np.repeat(units, np.diff(np.r_[starts, len(c)]))
        lib[ids] = 1e4 / units
        total[ids] = np.add.reduceat(cnt, starts)
        maxdev[ids] = np.maximum.reduceat(np.abs(cnt - np.round(cnt)) / np.maximum(cnt, 1), starts)
        m = is_panel[g]
        counts[c[m], [gpos[x] for x in g[m]]] = np.round(cnt[m])

    for g, c, v in it:
        seen += len(v)
        (g, c, v), carry = per_cell_units(g, c, v, carry)
        if len(c):
            take(g, c, v)
    if carry is not None and len(carry[1]):
        take(*carry)
    assert seen == dims[2], f"read {seen:,} entries, header says {dims[2]:,}"
    ok_int = bool(np.nanmax(maxdev) < 1e-6)
    lib_err = float(np.nanmax(np.abs(total - lib) / lib))
    print(f"identity: integral counts {ok_int} (max rel dev {np.nanmax(maxdev):.1e}); totals vs library "
          f"max rel error {lib_err:.1e}", flush=True)
    if not ok_int or lib_err > 1e-6:
        print("BLOCKED: the log-normalisation identity does not hold for every cell"); return 2

    lcpm = np.log1p(counts / lib[:, None] * 1e6)
    d = lcpm[:, :len(tgenes)].mean(1) - lcpm[:, len(tgenes):].mean(1)
    rows = []
    for p in sorted(set(patient)):
        inp = patient == p
        tnk = inp & (clusters == TNK)
        r = {"patient": p, "n_cells": int(inp.sum()),
             "T_strict": int((tnk & (d > 0)).sum()), "NK": int((tnk & (d < 0)).sum()),
             "unresolved_tie": int((tnk & (d == 0)).sum()), "B": int((inp & (clusters == B)).sum()),
             "C12_unmapped": int((inp & (clusters == UNMAPPED)).sum())}
        r["T_declared_rule"] = r["T_strict"] + r["unresolved_tie"]
        rows.append(r)
    df = pd.DataFrame(rows).set_index("patient")
    pooled = df.sum()
    lym = df[["T_strict", "NK", "B"]]
    ok = lym.sum(axis=1) > 0
    shares = lym[ok].div(lym[ok].sum(axis=1), axis=0)
    n_t_over_b = int((df["T_strict"] > df["B"]).sum())
    n_b_over_t = int((df["B"] >= df["T_strict"]).sum())
    n = len(df)
    if pooled["T_strict"] > pooled["B"] and n_t_over_b > n / 2:
        reading = "SUPPORTS T > B"
    elif pooled["B"] >= pooled["T_strict"] or n_b_over_t > n / 2:
        reading = "CONTRADICTS T > B"
    else:
        reading = "MIXED"
    epi = {}
    for c, base in (("gbm", ""), ("lgg", "_lgg")):
        epi[c] = json.loads((config.RESULTS_DIR / f"lymphoid_ordering{base}.json").read_text())["truth_mean"]
    out = {"rule": "prespecified/atlas_t_vs_b.md", "source": prov["source"], "matrix": path.name,
           "identity_check": {"integral_counts": ok_int, "library_max_rel_error": lib_err},
           "panels": {"T": tgenes, "NK": ngenes}, "n_patients": n,
           "per_patient": df.reset_index().to_dict(orient="records"),
           "pooled": {k: int(v) for k, v in pooled.items()},
           "pooled_within_lymphoid_share": (pooled[["T_strict", "NK", "B"]] / pooled[["T_strict", "NK", "B"]].sum()).round(4).to_dict(),
           "mean_patient_within_lymphoid_share": shares.mean().round(4).to_dict(),
           "n_patients_T_strict_over_B": n_t_over_b, "n_patients_B_at_or_over_T": n_b_over_t,
           "epidish_tcga_within_lymphoid_mean": epi, "reading": reading}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(df.to_string())
    print("pooled:", out["pooled"]); print("pooled shares:", out["pooled_within_lymphoid_share"])
    print(f"T (strict) > B in {n_t_over_b} of {n} patients | reading: {reading}")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
