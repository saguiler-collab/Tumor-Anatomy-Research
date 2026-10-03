"""How tissue-like is each reference profile? (WHY_B_OVER_T §7h)

For each cell-type column of the frozen signature, the Pearson correlation (log1p scale) between
the column and the MEAN bulk profile of each TCGA cohort, over the TCGA arm's own gene space --
with and without ribosomal-protein genes. A minor column whose profile resembles generic bulk tissue
is the one best placed to absorb whatever the major columns leave unexplained. Also records the
donor concentration of GBmap's training B and T cells. Selects nothing.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from b_column_diagnostics import cohort_inputs   # noqa: E402

from ivygap import config                         # noqa: E402
from ivygap.data.reference import load_frozen_reference  # noqa: E402

RP = re.compile(r"^(RP[SL]\d|RPLP\d)")
IG = re.compile(r"^(IG[HKL][VDJC]|IGH[GAMDE]|IGKC|IGLC|IGLL|JCHAIN)")   # immunoglobulin genes


def main() -> int:
    ref = load_frozen_reference()
    out = {"what_this_is": __doc__.split("\n")[0], "cohorts": {}}
    for c in ("gbm", "lgg"):
        sub, genes, _ = cohort_inputs(c, ref)
        P = ref.profile.loc[genes]; mb = sub.mean(axis=1)
        nonrp = [g for g in genes if not RP.match(g)]
        ig = [g for g in genes if IG.match(g)]
        rp = [g for g in genes if RP.match(g)]
        out.setdefault("signature_gene_families", {})[c] = {
            "immunoglobulin_genes_in_space": ig,
            "immunoglobulin_mass_share": {t: round(float(P.loc[ig, t].sum() / P[t].sum()), 5) for t in P.columns},
            "ribosomal_genes_in_space": len(rp),
            "ribosomal_mass_share": {t: round(float(P.loc[rp, t].sum() / P[t].sum()), 4) for t in P.columns}}
        out["cohorts"][c] = {"n_genes": len(genes), "n_non_rp": len(nonrp), "r_with_mean_bulk": {
            t: {"all_genes": round(float(np.corrcoef(np.log1p(P[t]), np.log1p(mb))[0, 1]), 4),
                "without_rp": round(float(np.corrcoef(np.log1p(P.loc[nonrp, t]), np.log1p(mb[nonrp]))[0, 1]), 4)}
            for t in P.columns}}
    # donor concentration of training B and T cells (obs only)
    import h5py
    import pandas as pd
    from anndata.io import read_elem
    from ivygap.data.reference import _collapse_to_roster
    with h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r") as f:
        obs = read_elem(f["obs"])
    ct = _collapse_to_roster(obs["annotation_level_3"].astype(str), config.GBMAP_CELL_TYPE_MAP).to_numpy()
    donor = obs["donor_id"].astype(str).to_numpy()
    sel = json.loads((config.BENCH_DIR / "method_selection_decision.json").read_text())
    train = np.isin(donor, list(map(str, sel["train_donors"])))
    conc = {}
    for t in ("B_cell", "T_cell", "NK_cell"):
        vc = pd.Series(donor[(ct == t) & train]).value_counts()
        conc[t] = {"n_cells": int(vc.sum()), "n_donors": int(len(vc)),
                   "top3_donor_share": round(float(vc.head(3).sum() / vc.sum()), 4)}
    out["training_donor_concentration"] = conc
    (config.RESULTS_DIR / "b_profile_tissue_likeness.json").write_text(json.dumps(out, indent=2))
    for c, r in out["cohorts"].items():
        print(c.upper(), {t: v["all_genes"] for t, v in r["r_with_mean_bulk"].items() if t in ("T_cell", "NK_cell", "B_cell")})
    print("donor concentration:", conc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
