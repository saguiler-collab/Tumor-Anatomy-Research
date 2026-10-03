"""WHY_B_OVER_T §5.1 -- re-annotate NK and refit. The decisive test of suspect 4a.

Suspect 4a: the frozen reference's NK column carries the defining T-cell markers (CD3D/E/G, LCK),
so real T-cell signal could be attributed to NK. Summing fitted T and NK (scripts/
lymphoid_pooled_tnk.py) does not rescue the ordering, but a sum is not a refit: collinearity
between the T and NK columns could change how much the fit gives B. So the reference itself is
changed and the methods are refit:

  * as_registered   -- the eight-type frozen signature, unchanged (reproduces the headline);
  * nk_removed      -- the NK column deleted; nothing can take T-cell signal into NK;
  * t_nk_merged     -- T and NK replaced by one T/NK column, profile and cell size averaged with
                       weights from the atlas's TRAINING-donor cell counts (recorded).

THE GENE SPACE IS FIXED ACROSS VARIANTS -- chosen once, from the original eight-type reference,
exactly as absolute_purity_yardstick.py chooses it -- so only the roster changes. Methods: the
registered NNLS and nu-SVR (CIBERSORT core) implementations; fast enough to refit every variant.
Methylation is used at the cohort level only (results/lymphoid_tracking.json). Selects nothing.
"""
from __future__ import annotations

import gzip
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402
from lymphoid_ordering import _truth, k4s                      # noqa: E402

from ivygap import config                                      # noqa: E402
from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput              # noqa: E402
from ivygap.deconv.classical import NNLSDeconvolution, SVRDeconvolution  # noqa: E402


def training_weights() -> dict:
    import h5py                                                # noqa: PLC0415
    from anndata.io import read_elem                           # noqa: PLC0415
    from ivygap.data.reference import _collapse_to_roster      # noqa: PLC0415
    sel = json.loads((config.BENCH_DIR / "method_selection_decision.json").read_text())
    with h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r") as f:
        obs = read_elem(f["obs"])
    ct = _collapse_to_roster(obs["annotation_level_3"].astype(str), config.GBMAP_CELL_TYPE_MAP)
    keep = obs["donor_id"].astype(str).isin(set(map(str, sel["train_donors"]))).to_numpy()
    n = pd.Series(ct.to_numpy()[keep]).value_counts()
    t, nk = int(n.get("T_cell", 0)), int(n.get("NK_cell", 0))
    return {"n_T_cells": t, "n_NK_cells": nk, "w_T": t / (t + nk), "w_NK": nk / (t + nk)}


def variant(ref, kind: str, w: dict):
    if kind == "as_registered":
        return ref
    keep = [c for c in ref.cell_types if c != "NK_cell"]
    prof, sig, cs = ref.profile[keep].copy(), ref.sigma[keep].copy(), ref.cell_size[keep].copy()
    if kind == "t_nk_merged":
        prof["T_cell"] = w["w_T"] * ref.profile["T_cell"] + w["w_NK"] * ref.profile["NK_cell"]
        sig["T_cell"] = w["w_T"] * ref.sigma["T_cell"] + w["w_NK"] * ref.sigma["NK_cell"]
        cs["T_cell"] = w["w_T"] * ref.cell_size["T_cell"] + w["w_NK"] * ref.cell_size["NK_cell"]
    return replace(ref, name=f"{ref.name}_{kind}", profile=prof, sigma=sig, cell_size=cs,
                   donor_profiles=None)


def main() -> int:
    w = training_weights()
    base_ref = load_frozen_reference()
    out = {"what_this_is": __doc__.split("\n")[0], "merge_weights": w, "cohorts": {}}
    for c in ("gbm", "lgg"):
        base = "" if c == "gbm" else "_lgg"
        with gzip.open(config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz", "rt") as fh:
            bulk = pd.read_csv(fh, index_col=0)
        absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
        absol["k"] = absol["sample"].map(key4)
        absol = absol.drop_duplicates("k").set_index("k")
        bmap = {key4(x): x for x in bulk.columns}
        cols = [bmap[k] for k in sorted(set(absol.index) & set(bmap))]
        shared = [g for g in base_ref.profile.index if g in bulk.index]
        genes = [g for g in select_signature_genes(base_ref.subset_genes(shared),
                                                   n_per_type=config.SIGNATURE_GENES_PER_TYPE)
                 if g in bulk.index]
        sub = bulk.loc[genes, cols]
        sub = sub / sub.sum(axis=0) * 1e6
        man = pd.DataFrame({"patient_id": ["-".join(str(x).split("-")[:3]) for x in sub.columns]},
                           index=sub.columns)
        truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv").dropna()
        res = {"n_samples": len(cols), "n_genes": len(genes),
               "methylation": {"frac_T_over_B": round(float((truth.T_cell > truth.B_cell).mean()), 4),
                               "frac_TNK_over_B": round(float((truth.T_cell + truth.NK_cell > truth.B_cell).mean()), 4)},
               "variants": {}}
        for kind in ("as_registered", "nk_removed", "t_nk_merged"):
            ref = variant(base_ref, kind, w).subset_genes(genes)
            data = DeconvolutionInput(bulk=sub, references=(ref,), manifest=man,
                                      cell_types=tuple(ref.cell_types))
            vr = {}
            for m in (NNLSDeconvolution(), SVRDeconvolution()):
                print(f'  {c}/{kind}: {m.name} ...', flush=True)
                est = m.fit_predict(data)
                e = est.copy(); e.index = [k4s(i) for i in e.index]; e = e[~e.index.duplicated()]
                e = e.loc[[k for k in e.index if k in truth.index]]
                lymph = [x for x in ("T_cell", "NK_cell", "B_cell") if x in e.columns]
                tot = e[lymph].sum(axis=1)
                ok = tot > 0
                tt = e.loc[ok, "T_cell"] + (e.loc[ok, "NK_cell"] if kind == "as_registered" else 0)
                vr[m.name] = {
                    "n_with_lymphoid_signal": int(ok.sum()), "n_matched": int(len(e)),
                    "frac_B_over_T": round(float((e.loc[ok, "B_cell"] > e.loc[ok, "T_cell"]).mean()), 4),
                    "frac_B_over_T_or_TNK": round(float((e.loc[ok, "B_cell"] > tt).mean()), 4),
                    "mean_share_within_lymphoid": {x: round(float((e.loc[ok, x] / tot[ok]).mean()), 4)
                                                   for x in lymph}}
            res["variants"][kind] = vr
            out["cohorts"][c] = res
            (config.RESULTS_DIR / "nk_reannotation.json").write_text(json.dumps(out, indent=2))
            print(f"  saved {c}/{kind}", flush=True)
        out["cohorts"][c] = res
    (config.RESULTS_DIR / "nk_reannotation.json").write_text(json.dumps(out, indent=2))
    print(f"merge weights from training cells: T {w['n_T_cells']:,} / NK {w['n_NK_cells']:,} -> "
          f"w_T {w['w_T']:.3f}")
    for c, r in out["cohorts"].items():
        print(f"\n=== {c.upper()} (n={r['n_samples']}, {r['n_genes']} genes)  methylation: "
              f"T>B in {r['methylation']['frac_T_over_B']:.1%}, T+NK>B in {r['methylation']['frac_TNK_over_B']:.1%}")
        for kind, vr in r["variants"].items():
            for m, v in vr.items():
                print(f"   {kind:<14} {m:<5} lymphoid signal {v['n_with_lymphoid_signal']:>3}/{v['n_matched']:<3} "
                      f"B>T {v['frac_B_over_T']:>6.1%}  B>(T or T/NK) {v['frac_B_over_T_or_TNK']:>6.1%}  "
                      f"shares {v['mean_share_within_lymphoid']}")
    print("\nwrote results/nk_reannotation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
