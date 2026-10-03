"""
gbmap_t_vs_b.py -- T against B cells counted directly in GBmap's core atlas (all glioblastoma), per donor.

DESCRIPTIVE, no choices beyond the atlas's own CELLxGENE labels: T = "mature T cell", NK = "natural
killer cell", B = "B cell"; plasma cells are counted separately and never added to B. Reads `obs` only
(no expression). Companion to scripts/atlas_t_vs_b.py (an independent atlas); quoted in the manuscript's
"is the truth itself right?" paragraph (§4.4) and WHY_B_OVER_T §7n.

    python3 scripts/gbmap_t_vs_b.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

OUT = config.RESULTS_DIR / "gbmap_t_vs_b.json"
LABELS = {"T": "mature T cell", "NK": "natural killer cell", "B": "B cell", "plasma": "plasma cell"}


def main() -> int:
    with h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r") as f:
        obs = read_elem(f["obs"])
    disease = obs["disease"].astype(str).value_counts().to_dict()
    ct = obs["cell_type"].astype(str)
    donor = obs["donor_id"].astype(str)
    per = pd.DataFrame({k: (ct == v).groupby(donor).sum() for k, v in LABELS.items()})
    lym = per[["T", "NK", "B"]].sum(axis=1)
    with_lym = per[lym > 0]
    pooled = per.sum()
    out = {"source": "GBmap core (gbmap_core.h5ad), obs only", "labels": LABELS, "disease": disease,
           "n_cells": int(len(obs)), "n_donors": int(per.shape[0]),
           "n_donors_with_any_lymphoid": int(len(with_lym)),
           "pooled": {k: int(v) for k, v in pooled.items()},
           "pooled_T_to_B": round(float(pooled["T"] / pooled["B"]), 2) if pooled["B"] else None,
           "pooled_within_lymphoid_share": (pooled[["T", "NK", "B"]] / pooled[["T", "NK", "B"]].sum()).round(4).to_dict(),
           "n_donors_T_over_B": int((with_lym["T"] > with_lym["B"]).sum()),
           "n_donors_B_at_or_over_T": int((with_lym["B"] >= with_lym["T"]).sum()),
           "n_donors_B_zero": int((with_lym["B"] == 0).sum())}
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k != "labels"}, indent=1))
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
