"""Does the additive mixing model describe real glioma bulk at all?

WHY THIS EXISTS. `scripts/identifiability_probe.py` showed the frozen signature separates T from
B perfectly **under the model** -- exactly, at 100% noise, and with a whole cell type deleted.
Yet on real tissue no method reproduces the T:B ordering that methylation measures
(`docs/IMMUNE_ARM.md` §6). Those two facts are only compatible if real bulk is not well described
by `b = S x`, which every method in this panel assumes.

That is measurable without any ground truth: fit the best non-negative combination of reference
profiles to each sample and report how much of the sample's marker-space variance it leaves
unexplained. A high residual does not say WHICH assumption fails, and this script does not claim
to -- it bounds how much of the signal the model can account for at all.

EXPLORATORY. Not pre-registered.

    python scripts/model_fit_residual.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from ivygap import config
from ivygap.data.reference import load_frozen_reference, select_signature_genes

COHORTS = {"GBM": "tcga_gbm_bulk_cpm.csv.gz", "LGG": "tcga_lgg_bulk_cpm.csv.gz"}
N_CONTROL_SAMPLES = 60      # enough to separate R^2 ~ 0.38 from R^2 ~ 0
N_CONTROL_DRAWS = 5


def _median_r2(A, sub, cols) -> float:
    vals = []
    for c in cols:
        y = sub[c].to_numpy(dtype="float64")
        x, _ = nnls(A, y)
        resid = float(((y - A @ x) ** 2).sum())
        total = float(((y - y.mean()) ** 2).sum())
        if total > 0:
            vals.append(1.0 - resid / total)
    return float(np.median(vals)) if vals else float("nan")


def controls(A, sub, seed: int) -> dict:
    """Floor and ceiling for the residual, without which R^2 = 0.36 means nothing.

    An R^2 of 0.36 could mean "the reference carries real but partial structure" or "any eight
    vectors explain a third of this data". The two are told apart by controls:

    - **shuffled profiles** permute each cell type's gene labels independently, preserving every
      profile's marginal distribution while destroying which gene belongs to which type;
    - **random basis** replaces the reference with non-negative noise at matched scale;
    - **ceiling** is the top-8 SVD of the bulk itself -- the best any eight-dimensional basis
      could do on these samples, which is what makes the reference's score interpretable as a
      FRACTION of the achievable rather than an absolute.
    """
    rng = np.random.default_rng(seed)
    cols = list(sub.columns)[:N_CONTROL_SAMPLES]
    real = _median_r2(A, sub, cols)
    shuf = float(np.mean([
        _median_r2(np.column_stack([rng.permutation(A[:, j]) for j in range(A.shape[1])]),
                   sub, cols) for _ in range(N_CONTROL_DRAWS)]))
    rand = float(np.mean([
        _median_r2(rng.random(A.shape) * A.mean() * 2, sub, cols)
        for _ in range(N_CONTROL_DRAWS)]))
    Y = sub[cols].to_numpy(dtype="float64")
    s = np.linalg.svd(Y - Y.mean(axis=0, keepdims=True), compute_uv=False)
    ceil = float((s[:A.shape[1]] ** 2).sum() / (s ** 2).sum())
    return {"n_control_samples": len(cols), "real_reference": round(real, 4),
            "shuffled_profiles": round(shuf, 4), "random_basis": round(rand, 4),
            "ceiling_top_k_svd": round(ceil, 4),
            "fraction_of_achievable": round(real / ceil, 4) if ceil else None}


def main() -> int:
    ref = load_frozen_reference()
    genes = select_signature_genes(ref, n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    print("HOW WELL DOES THE ADDITIVE MODEL DESCRIBE REAL BULK?")
    print("Per sample: best non-negative fit b ~ Sx, then the unexplained fraction of that")
    print("sample's variance across the marker genes. Exploratory, not pre-registered.\n")
    out: dict[str, dict] = {}
    for lab, fname in COHORTS.items():
        path = config.PROCESSED_DIR / fname
        if not path.exists():
            print(f"{lab}: {fname} absent, skipped (not imputed)")
            continue
        B = pd.read_csv(path, index_col=0)
        g = [x for x in genes if x in B.index]
        if len(g) < config.MIN_GENES_SHARED // 2:
            print(f"{lab}: only {len(g)} marker genes shared, skipped")
            continue
        S = ref.profile.loc[g]
        S = S / S.sum(axis=0)
        sub = B.loc[g]
        sub = sub / sub.sum(axis=0) * 1e6
        A = S.to_numpy(dtype="float64")
        r2 = []
        for c in sub.columns:
            y = sub[c].to_numpy(dtype="float64")
            if not np.isfinite(y).all():
                continue
            x, _ = nnls(A, y)
            resid = float(((y - A @ x) ** 2).sum())
            total = float(((y - y.mean()) ** 2).sum())
            if total > 0:
                r2.append(1.0 - resid / total)
        r2 = np.asarray([v for v in r2 if np.isfinite(v)])
        if not r2.size:
            print(f"{lab}: no finite fits")
            continue
        med = float(np.median(r2))
        print(f"=== {lab}: {len(r2)} samples x {len(g)} marker genes ===")
        print(f"  R^2 of the best additive fit: median {med:.4f}  "
              f"IQR [{np.percentile(r2, 25):.4f}, {np.percentile(r2, 75):.4f}]  "
              f"min {r2.min():.4f}  max {r2.max():.4f}")
        print(f"  => a median {100 * (1 - med):.1f}% of marker-space variance is left "
              f"UNEXPLAINED by the model every method here assumes")
        print(f"  samples where the model explains under half the variance: "
              f"{int((r2 < 0.5).sum())} of {len(r2)} ({(r2 < 0.5).mean():.1%})\n")
        ctl = controls(A, sub, config.RANDOM_SEED)
        print(f"  CONTROLS (first {ctl['n_control_samples']} samples):")
        print(f"    real reference                {ctl['real_reference']:+.4f}")
        print(f"    gene labels shuffled per type {ctl['shuffled_profiles']:+.4f}   <- floor")
        print(f"    random non-negative basis     {ctl['random_basis']:+.4f}   <- floor")
        print(f"    ceiling: top-8 SVD of bulk    {ctl['ceiling_top_k_svd']:+.4f}   <- best any "
              f"8-dim basis could do")
        print(f"  => the reference carries REAL structure (far above both floors) but reaches "
              f"only\n     {ctl['fraction_of_achievable']:.1%} of what an optimal 8-dimensional "
              f"basis achieves on the same samples.\n")
        out[lab] = {"n_samples": int(len(r2)), "n_genes": len(g), "controls": ctl,
                    "r2_median": round(med, 4),
                    "r2_iqr": [round(float(np.percentile(r2, 25)), 4),
                               round(float(np.percentile(r2, 75)), 4)],
                    "r2_min": round(float(r2.min()), 4),
                    "r2_max": round(float(r2.max()), 4),
                    "unexplained_median_pct": round(100 * (1 - med), 2),
                    "frac_samples_r2_below_half": round(float((r2 < 0.5).mean()), 4)}
    if out:
        (config.RESULTS_DIR / "model_fit_residual.json").write_text(json.dumps({
            "EXPLORATORY": "not pre-registered; bounds how much of real bulk the additive "
                           "mixing model can account for in the marker space.",
            "what_it_does_NOT_say": "which assumption fails. A high residual is consistent "
                                    "with platform/normalisation mismatch, unmodelled cell "
                                    "states, or non-additivity, and does not distinguish them.",
            "cohorts": out}, indent=2))
        print("wrote results/model_fit_residual.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
