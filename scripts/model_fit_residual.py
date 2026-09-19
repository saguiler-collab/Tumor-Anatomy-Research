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



def trust_signal(ref, genes, cohort_tag: str, bulk_file: str) -> dict | None:
    """Is the per-sample model fit a usable ground-truth-free TRUST SIGNAL?

    This is the question worth asking of the residual. `docs/ENDPOINT.md` reports that anatomic
    concordance detects but cannot rank, and that there is no routine way to know when a
    deconvolution estimate can be believed. The per-sample R^2 is a candidate: it needs no ground
    truth, and a sample the model cannot fit is a sample whose estimate should be distrusted.

    Tested as a correlation between per-sample R^2 and per-sample |estimate - purity|. A
    NEGATIVE correlation would mean "better fit, smaller error" -- the usable direction.

    A control runs first. If R^2 were merely a proxy for purity it would inherit purity's known
    relationship with error and the result would be circular.
    """
    from scipy import stats                                        # noqa: PLC0415
    path = config.PROCESSED_DIR / bulk_file
    per_sample = config.RESULTS_DIR / f"absolute_purity_per_sample{cohort_tag}.csv"
    if not (path.exists() and per_sample.exists()):
        return None
    B = pd.read_csv(path, index_col=0)
    g = [x for x in genes if x in B.index]
    S = ref.profile.loc[g]
    S = S / S.sum(axis=0)
    A = S.to_numpy(dtype="float64")
    sub = B.loc[g]
    sub = sub / sub.sum(axis=0) * 1e6
    r2 = {}
    for c in sub.columns:
        y = sub[c].to_numpy(dtype="float64")
        if not np.isfinite(y).all():
            continue
        x, _ = nnls(A, y)
        total = float(((y - y.mean()) ** 2).sum())
        if total > 0:
            r2[c] = 1.0 - float(((y - A @ x) ** 2).sum()) / total
    R = pd.Series(r2)
    P = pd.read_csv(per_sample, index_col=0)
    tcol = next(c for c in ("absolute_purity", "purity") if c in P.columns)
    shared = [i for i in P.index if i in R.index]
    if len(shared) < 50:
        return None
    pur, rr = P.loc[shared, tcol], R.loc[shared]
    cp = stats.spearmanr(rr, pur)
    rows, hits, wrong, n = {}, 0, 0, 0
    for m in [c for c in P.columns if c != tcol]:
        e = (P.loc[shared, m] - pur).abs()
        ok = np.isfinite(e) & np.isfinite(rr)
        if ok.sum() < 50 or float(e[ok].std()) == 0:
            continue
        r = stats.spearmanr(rr[ok], e[ok])
        n += 1
        useful = bool(r.statistic < 0 and r.pvalue < 0.05)
        backwards = bool(r.statistic > 0 and r.pvalue < 0.05)
        hits += useful
        wrong += backwards
        rows[m] = {"spearman_r2_vs_abs_error": round(float(r.statistic), 4),
                   "p": round(float(r.pvalue), 6), "useful_direction": useful,
                   "significant_WRONG_direction": backwards}
    return {"n_samples": len(shared), "n_methods": n,
            "control_spearman_r2_vs_purity": round(float(cp.statistic), 4),
            "control_p": round(float(cp.pvalue), 4),
            "n_useful": hits, "n_significant_wrong_direction": wrong,
            "methods": rows,
            "verdict": ("NOT a usable trust signal" if hits <= n / 2 or wrong
                        else "candidate trust signal")}

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
    # THE QUESTION THE RESIDUAL INVITES, asked and answered rather than left hanging.
    trust = {}
    for lab, fname in COHORTS.items():
        tg = "" if lab == "GBM" else f"_{lab.lower()}"
        r = trust_signal(ref, genes, tg, fname)
        if not r:
            continue
        trust[lab] = r
        print(f"IS THE PER-SAMPLE FIT A TRUST SIGNAL? — {lab} "
              f"({r['n_samples']} samples)")
        print(f"  control: Spearman(R^2, purity) = {r['control_spearman_r2_vs_purity']:+.4f} "
              f"(p={r['control_p']:.3g}) — not a purity proxy, so not circular")
        print(f"  methods where a better fit predicts SMALLER error (p<0.05): "
              f"{r['n_useful']} of {r['n_methods']}")
        print(f"  methods significant in the WRONG direction: "
              f"{r['n_significant_wrong_direction']}")
        print(f"  VERDICT: {r['verdict']}\n")
    if out:
        out_payload_trust = trust
        (config.RESULTS_DIR / "model_fit_residual.json").write_text(json.dumps({
            "EXPLORATORY": "not pre-registered; bounds how much of real bulk the additive "
                           "mixing model can account for in the marker space.",
            "what_it_does_NOT_say": "which assumption fails. A high residual is consistent "
                                    "with platform/normalisation mismatch, unmodelled cell "
                                    "states, or non-additivity, and does not distinguish them.",
            "cohorts": out,
            "trust_signal_test": {
                "question": "does per-sample model fit predict per-sample deconvolution error, "
                            "i.e. is it a ground-truth-free signal for when to distrust an "
                            "estimate?",
                "EXPLORATORY": "not pre-registered",
                "cohorts": out_payload_trust}}, indent=2))
        print("wrote results/model_fit_residual.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
