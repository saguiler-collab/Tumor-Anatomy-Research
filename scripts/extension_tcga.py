"""The POST-REGISTRATION EXTENSION PANEL on TCGA: purity, lymphoid ordering and anchoring.

Runs extension methods (`ivygap.deconv.extension`) on exactly the inputs the registered TCGA arm
(`absolute_purity_yardstick.py`) gave the registered panel -- same samples (those with an
ABSOLUTE call), same reference-only marker selection, same CPM re-normalisation on the marker
space, same `fit_predict` post-processing -- and scores them on the same three questions:

  1. Spearman of the Tumor fraction against ABSOLUTE purity (DNA copy number);
  2. the lymphoid ordering against DNA methylation (`lymphoid_ordering.py`'s own functions);
  3. for the h5ad reference, the L1 distance of the cohort mean from the reference's donor-mean
     composition (`results/bisque_anchoring.json`).

Writes only under `results/extension/`. Never touches a registered artefact; never enters a
registered statistic.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4      # noqa: E402
from lymphoid_ordering import COMPARE, _truth, k4s, ordering_of, renormalise  # noqa: E402

from ivygap import config                                             # noqa: E402
from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                     # noqa: E402
from ivygap.deconv.extension import EXTENSION_SPECS, build_extension_methods  # noqa: E402

OUT = config.RESULTS_DIR / "extension"


def run(cohort: str, reference: str, methods: list[str], return_estimates: bool = False) -> dict:
    base = "" if cohort == "gbm" else f"_{cohort}"
    tag = base + ("_h5ad" if reference == "h5ad" else "")
    with gzip.open(config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz", "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(c): c for c in bulk.columns}
    shared = sorted(set(absol.index) & set(bmap))
    cols = [bmap[k] for k in shared]
    purity = absol.loc[shared, PURITY_COL].to_numpy(dtype="float64")

    if reference == "h5ad":
        from ivygap.data.reference import build_from_h5ad             # noqa: PLC0415
        from ivygap.deconv import r_bridge                            # noqa: PLC0415
        ref, sc_expr, sc_meta = build_from_h5ad(
            config.REFERENCE_DIR / "gbmap_core.h5ad", matrix="raw/X",
            restrict_to_genes=bulk.index, export=False)
        r_bridge.set_cell_source(config.PRIMARY_REFERENCE, sc_expr, sc_meta)
    else:
        ref = load_frozen_reference()
    shared_genes = [g for g in ref.profile.index if g in bulk.index]
    picked = select_signature_genes(ref.subset_genes(shared_genes),
                                    n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    genes = [g for g in picked if g in bulk.index]
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6
    manifest = pd.DataFrame({"patient_id": ["-".join(str(c).split("-")[:3]) for c in sub.columns]},
                            index=sub.columns)
    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),), manifest=manifest,
                              cell_types=tuple(ref.cell_types), bulk_full=bulk.loc[:, cols])

    truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
    truth_ord = ordering_of(truth.mean())
    prior = None
    anch = config.RESULTS_DIR / "bisque_anchoring.json"
    if reference == "h5ad" and anch.exists():
        prior = pd.Series(json.loads(anch.read_text())["reference_prior_all_types"])

    res = {"cohort": cohort, "reference": reference, "n_samples": len(cols), "n_genes": len(genes),
           "truth_lymphoid_ordering": truth_ord, "methods": {}}
    frames, ests = [], {}
    for m in build_extension_methods(methods):
        print(f"  {cohort}/{reference}: {m.name} ...", flush=True)
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                       # noqa: BLE001
            res["methods"][m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:300]}"}
            print(f"    FAILED {res['methods'][m.name]['failed'][:120]}")
            continue
        full = est.reindex(index=list(sub.columns), columns=list(config.CELL_TYPES))
        t = full["Tumor"].to_numpy(dtype="float64")
        ok = np.isfinite(t)
        rec = {"implementation": getattr(m, "implementation_", "?"),
               "spearman_vs_purity": round(float(stats.spearmanr(t[ok], purity[ok]).statistic), 4),
               "n_purity": int(ok.sum()),
               "bias_tumor_minus_purity": round(float(np.mean(t[ok] - purity[ok])), 4)}
        g = full.copy()
        g.index = [k4s(i) for i in g.index]
        g = g[~g.index.duplicated()]
        idx = [k for k in g.index if k in truth.index]
        ren = renormalise(g.loc[idx])
        mu = ren.mean()
        rec.update({"n_lymphoid": len(idx),
                    "lymphoid_mean": {c: round(float(mu[c]), 4) for c in COMPARE},
                    "lymphoid_ordering": ordering_of(mu),
                    "T_exceeds_B": bool(mu["T_cell"] > mu["B_cell"]),
                    "full_ordering_correct": ordering_of(mu) == truth_ord,
                    "n_zero_lymphoid": int(ren.isna().all(axis=1).sum()),
                    "frac_samples_B_over_T": round(float(
                        (ren.dropna()["B_cell"] > ren.dropna()["T_cell"]).mean()), 4)})
        if prior is not None:
            rec["L1_cohort_mean_to_reference_prior"] = round(float(
                (full[list(prior.index)].mean() - prior).abs().sum()), 4)
        res["methods"][m.name] = rec
        frames.append(full.assign(method=m.name).rename_axis("sample").reset_index())
        ests[m.name] = full
        print(f"    purity rho {rec['spearman_vs_purity']:+.4f} | lymphoid {rec['lymphoid_ordering']} "
              f"(T>B {rec['T_exceeds_B']}, B>T in {rec['frac_samples_B_over_T']:.1%}, "
              f"zero-lymphoid {rec['n_zero_lymphoid']}/{rec['n_lymphoid']})"
              + (f" | L1 to prior {rec['L1_cohort_mean_to_reference_prior']}" if prior is not None else ""))
    # The per-sample file is SHARED by every extension method. Until 2026-10-02 each call overwrote it
    # with only its own methods, and the unmix ablation stored override settings under the method's
    # name. Now: never written while an ablation override is set, and merged by method otherwise.
    ablation = any(os.environ.get(k) for k in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"))
    if frames and not ablation:
        path = OUT / f"estimates_full{tag}_extension.csv"
        new = pd.concat(frames)
        if path.exists():
            old = pd.read_csv(path)
            new = pd.concat([old[~old["method"].isin(new["method"].unique())], new])
        new.to_csv(path, index=False)
    elif frames:
        print("    (ablation override set: per-sample file NOT written)")
    if return_estimates:
        res["_estimates"] = ests
        res["_inputs"] = {"bulk": sub, "signature": ref.subset_genes(genes).profile,
                          "purity": pd.Series(purity, index=cols)}
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--methods", nargs="+", default=list(EXTENSION_SPECS))
    ap.add_argument("--references", nargs="+", default=["frozen", "h5ad"], choices=["frozen", "h5ad"])
    ap.add_argument("--cohorts", nargs="+", default=["gbm", "lgg"], choices=["gbm", "lgg"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for r in a.references:
        for c in a.cohorts:
            res = run(c, r, a.methods)
            res["panel"] = "extension (post-registration, 2026-09-30)"
            p = OUT / f"tcga_{c}_{r}.json"
            prior = json.loads(p.read_text()) if p.exists() else {}
            prior.update({k: v for k, v in res.items() if k != "methods"})
            prior.setdefault("methods", {}).update(res["methods"])
            p.write_text(json.dumps(prior, indent=2))
            print(f"  wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
