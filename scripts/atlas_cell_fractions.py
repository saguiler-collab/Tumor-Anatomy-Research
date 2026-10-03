"""Per-cell ribosomal-protein and myelin UMI fractions, by cell type, in GBmap's raw counts (WHY_B_OVER_T §7h).

Two questions about the atlas the B profile comes from, measured from `raw/X` (genuine counts),
training donors only, up to 1,500 cells per type (seeded):
  * is the B profile's ribosomal excess a property of the atlas's B cells?  -> RP UMI fraction
  * is the B profile's tissue signal ambient myelin RNA?  -> fraction on MBP, PLP1, MOBP, MAG, MOG
    (gene set declared before the first run, 2026-10-01)
Also records median UMIs and mitochondrial fraction as quality controls. Selects nothing.
"""
from __future__ import annotations

import json
import re

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from ivygap import config
from ivygap.data.reference import _collapse_to_roster

MYELIN = ["MBP", "PLP1", "MOBP", "MAG", "MOG"]
RP = re.compile(r"^(RP[SL]\d|RPLP\d)")


def main() -> int:
    f = h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r")
    obs = read_elem(f["obs"])
    var = read_elem(f["raw/var"]) if "raw/var" in f else read_elem(f["var"])
    names = (var["feature_name"] if "feature_name" in var.columns else pd.Series(var.index)).astype(str).to_numpy()
    rp, mt, my = (np.array([bool(RP.match(n)) for n in names]), np.array([n.startswith("MT-") for n in names]),
                  np.isin(names, MYELIN))
    ct = _collapse_to_roster(obs["annotation_level_3"].astype(str), config.GBMAP_CELL_TYPE_MAP).to_numpy()
    sel = json.loads((config.BENCH_DIR / "method_selection_decision.json").read_text())
    train = obs["donor_id"].astype(str).isin(set(map(str, sel["train_donors"]))).to_numpy()
    X = f["raw/X"]; indptr = X["indptr"][:]; rng = np.random.default_rng(config.RANDOM_SEED)
    out = {"what_this_is": __doc__.split("\n")[0], "layer": "raw/X", "myelin_genes_found": list(names[my]),
           "n_rp_pattern_genes": int(rp.sum()), "by_type": {}}
    for t in config.CELL_TYPES:
        idx = np.where((ct == t) & train)[0]
        if len(idx) > 1500:
            idx = np.sort(rng.choice(idx, 1500, replace=False))
        tot, rpf, mtf, myf = [], [], [], []
        for i in idx:
            a, b = indptr[i], indptr[i + 1]
            d = X["data"][a:b]; j = X["indices"][a:b]; s = float(d.sum())
            tot.append(s)
            rpf.append(d[rp[j]].sum() / s if s else np.nan)
            mtf.append(d[mt[j]].sum() / s if s else np.nan)
            myf.append(d[my[j]].sum() / s if s else np.nan)
        out["by_type"][t] = {"n_cells": int(len(idx)), "median_umis": float(np.median(tot)),
                             "median_rp_fraction": round(float(np.nanmedian(rpf)), 4),
                             "median_mt_fraction": round(float(np.nanmedian(mtf)), 4),
                             "mean_myelin_fraction": round(float(np.nanmean(myf)), 6),
                             "frac_cells_with_myelin": round(float(np.mean(np.array(myf) > 0)), 4)}
        print(f"{t:<22} {out['by_type'][t]}", flush=True)
    (config.RESULTS_DIR / "atlas_cell_fractions.json").write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
