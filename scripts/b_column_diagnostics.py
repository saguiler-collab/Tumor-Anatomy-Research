"""Why does the B profile win the lymphoid signal? Two exploratory measurements (WHY_B_OVER_T §7h).

The NK re-annotation (scripts/nk_reannotation.py) showed that however the roster treats NK, the
B column takes the lymphoid share and T gets almost nothing. Two questions follow, both declared
before running (2026-10-01), both exploratory, neither selecting anything:

1. WHERE DOES B's MASS GO WHEN B IS REMOVED? Refit with the B column deleted (gene space fixed,
   as-registered baseline reproduced). If B's share moves mainly to T, B was competing with T
   for genuine lymphoid signal. If it moves mainly to non-lymphoid columns, B was absorbing
   misfit from the populations that dominate the tissue.

2. WHICH GENES PULL THE FIT TOWARD B? Fit the five non-lymphoid columns alone (NNLS); the
   positive residual is bulk signal they cannot explain. For each gene, the contribution to
   preferring B over T is  sum_samples max(r_g, 0) * (B_g/||B|| - T_g/||T||). The top genes are
   reported with every column's value, flagging whether B is that gene's highest column (a B
   marker) or the gene is shared with another population.

Frozen reference, the TCGA arm's own gene space, registered NNLS (and SVR for question 1).
Writes results/b_column_diagnostics.json.
"""
from __future__ import annotations

import gzip
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402

from ivygap import config                                      # noqa: E402
from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput              # noqa: E402
from ivygap.deconv.classical import NNLSDeconvolution, SVRDeconvolution  # noqa: E402

LYMPH = ["T_cell", "NK_cell", "B_cell"]
OUT = config.RESULTS_DIR / "b_column_diagnostics.json"


def cohort_inputs(c, ref):
    with gzip.open(config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz", "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(x): x for x in bulk.columns}
    cols = [bmap[k] for k in sorted(set(absol.index) & set(bmap))]
    shared = [g for g in ref.profile.index if g in bulk.index]
    genes = [g for g in select_signature_genes(ref.subset_genes(shared),
                                               n_per_type=config.SIGNATURE_GENES_PER_TYPE)
             if g in bulk.index]
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6
    man = pd.DataFrame({"patient_id": ["-".join(str(x).split("-")[:3]) for x in sub.columns]},
                       index=sub.columns)
    return sub, genes, man


def main() -> int:
    base_ref = load_frozen_reference()
    out = json.loads(OUT.read_text()) if OUT.exists() else {"cohorts": {}}
    out["what_this_is"] = __doc__.split("\n")[0]
    for c in ("gbm", "lgg"):
        sub, genes, man = cohort_inputs(c, base_ref)
        ref = base_ref.subset_genes(genes)
        rec = out["cohorts"].setdefault(c, {"n_samples": int(sub.shape[1]), "n_genes": len(genes)})

        # ---------- 2. driver genes (NNLS on the five non-lymphoid columns) ----------
        if "drivers" not in rec:
            S = ref.profile
            non = [t for t in S.columns if t not in LYMPH]
            S5, Bm = S[non].to_numpy(float), sub.to_numpy(float)
            R = np.column_stack([Bm[:, j] - S5 @ nnls(S5, Bm[:, j])[0] for j in range(Bm.shape[1])])
            Rpos = np.clip(R, 0, None)
            bn = S["B_cell"] / np.linalg.norm(S["B_cell"]); tn = S["T_cell"] / np.linalg.norm(S["T_cell"])
            contrib = pd.Series(Rpos.sum(axis=1), index=S.index) * (bn - tn)
            top = contrib.sort_values(ascending=False).head(25)
            rows = []
            for g, v in top.items():
                col = S.loc[g]
                rows.append({"gene": g, "contribution": round(float(v), 2),
                             "highest_column": str(col.idxmax()),
                             "B_cell": round(float(col["B_cell"]), 1), "T_cell": round(float(col["T_cell"]), 1),
                             "NK_cell": round(float(col["NK_cell"]), 1),
                             "Macrophage_Microglia": round(float(col["Macrophage_Microglia"]), 1),
                             "Tumor": round(float(col["Tumor"]), 1),
                             "mean_positive_residual": round(float(Rpos[S.index.get_loc(g)].mean()), 1)})
            tot = float(contrib.clip(lower=0).sum())
            rec["drivers"] = {
                "top25": rows,
                "share_of_positive_B_preference_in_top25": round(float(top.clip(lower=0).sum()) / tot, 4) if tot else None,
                "n_top25_where_B_is_highest_column": sum(r["highest_column"] == "B_cell" for r in rows),
                "highest_column_counts_top25": pd.Series([r["highest_column"] for r in rows]).value_counts().to_dict(),
            }
            OUT.write_text(json.dumps(out, indent=2)); print(f"  saved {c}/drivers", flush=True)

        # ---------- 1. where does B's mass go when B is removed? ----------
        for mcls in (NNLSDeconvolution, SVRDeconvolution):
            key = f"b_removed_{mcls.name}"
            if key in rec:
                continue
            print(f"  {c}/{key} ...", flush=True)
            full = mcls().fit_predict(DeconvolutionInput(bulk=sub, references=(ref,), manifest=man,
                                                         cell_types=tuple(ref.cell_types)))
            keep = [t for t in ref.cell_types if t != "B_cell"]
            nob_ref = replace(ref, name=ref.name + "_b_removed", profile=ref.profile[keep],
                              sigma=ref.sigma[keep], cell_size=ref.cell_size[keep], donor_profiles=None)
            nob = mcls().fit_predict(DeconvolutionInput(bulk=sub, references=(nob_ref,), manifest=man,
                                                        cell_types=tuple(keep)))
            has_b = full["B_cell"] > 1e-6
            dB = float(full.loc[has_b, "B_cell"].mean()) if has_b.any() else 0.0
            delta = {t: float((nob.loc[has_b, t] - full.loc[has_b, t]).mean()) for t in keep}
            rec[key] = {
                "n_samples_with_B": int(has_b.sum()),
                "mean_B_fraction_removed": round(dB, 5),
                "mean_change_by_type": {t: round(v, 5) for t, v in delta.items()},
                "share_of_B_mass_to": {t: (round(v / dB, 4) if dB else None) for t, v in delta.items()},
                "share_to_T_or_NK": round((delta["T_cell"] + delta["NK_cell"]) / dB, 4) if dB else None,
            }
            OUT.write_text(json.dumps(out, indent=2)); print(f"  saved {c}/{key}", flush=True)

    for c, rec in out["cohorts"].items():
        d = rec["drivers"]
        print(f"\n=== {c.upper()} drivers: B is the highest column for {d['n_top25_where_B_is_highest_column']} "
              f"of the top 25; highest-column counts {d['highest_column_counts_top25']}; top 25 carry "
              f"{d['share_of_positive_B_preference_in_top25']:.1%} of the B-over-T preference")
        for r in d["top25"][:15]:
            print(f"   {r['gene']:<10} contrib {r['contribution']:>10.1f}  highest {r['highest_column']:<20} "
                  f"B {r['B_cell']:>8.1f}  T {r['T_cell']:>8.1f}  NK {r['NK_cell']:>8.1f}  Mac {r['Macrophage_Microglia']:>8.1f}  "
                  f"Tum {r['Tumor']:>7.1f}  resid+ {r['mean_positive_residual']:>7.1f}")
        for k in ("b_removed_nnls", "b_removed_svr"):
            if k in rec:
                v = rec[k]
                print(f"   {k}: {v['n_samples_with_B']} samples with B; mean B {v['mean_B_fraction_removed']:.4f}; "
                      f"to T+NK {v['share_to_T_or_NK']}; by type {v['share_of_B_mass_to']}")
    print(f"\nwrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
