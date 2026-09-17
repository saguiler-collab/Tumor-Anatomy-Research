#!/usr/bin/env python3
"""
failure_factors.py — which biological properties of a GBM predict deconvolution error?

Executes `prespecified/biological_failure_factors.md` exactly. That document was committed
before any error-versus-factor relationship was examined, and nothing here deviates from it:
the seven factors, their directions, the joint OLS model, the sign test across methods, Holm
across seven, and the falsification criterion are all fixed there.

THE HYPOTHESIS
--------------
Biological properties of GBM that alter the relationship between malignant-cell abundance and
transcriptomic composition will predict systematic errors in bulk RNA deconvolution.

WHY A SIGN TEST ACROSS METHODS RATHER THAN ONE POOLED REGRESSION
-----------------------------------------------------------------
Twelve methods on the same 154 samples are twelve semi-independent replications of the same
measurement. A real biological effect should not depend on which solver produced the estimate.
Pooling would let one method with a large residual variance dominate; requiring a consistent
SIGN across methods asks the question the hypothesis actually poses.

    python scripts/failure_factors.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

PER_SAMPLE = config.RESULTS_DIR / "absolute_purity_per_sample.csv"
ABS_T = config.RAW_DIR / "tcga" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
CBIO = Path("/Users/tatopro9130/Downloads/cancer_judging_machine-master 22"
            "/public/04 Subtypes and Pathways/source")

#: factor -> (pre-specified direction on the stated outcome, which outcome)
#: "neg" = coefficient expected NEGATIVE on signed error; "pos" = POSITIVE on |error|.
FACTORS = {
    "purity":            ("neg", "signed"),
    "ploidy":            ("pos", "abs"),
    "genome_doublings":  ("pos", "abs"),
    "subclonal_fraction":("pos", "abs"),
    "is_mesenchymal":    ("neg", "signed"),
    "idh1_mutant":       ("pos", "abs"),
    "gcimp":             ("pos", "abs"),
}
ALPHA = 0.05


def k4(s): return "-".join(str(s).split("-")[:4])
def k3(s): return "-".join(str(s).split("-")[:3])


def build_covariates(samples: list[str]) -> pd.DataFrame:
    a = pd.read_csv(ABS_T, sep="\t").dropna(subset=["Cancer DNA fraction"])
    a["k"] = a["sample"].map(k4)
    a = a.drop_duplicates("k").set_index("k")
    s = pd.read_csv(CBIO / "gbm_tcga_pub2013_clinical_sample.txt", sep="\t", comment="#")
    s["k3"] = s["PATIENT_ID"].astype(str)
    s = s.drop_duplicates("k3").set_index("k3")

    rows = {}
    for smp in samples:
        a_k, c_k = k4(smp), k3(smp)
        if a_k not in a.index:
            continue
        r = {"purity": pd.to_numeric(a.loc[a_k, "purity"], errors="coerce"),
             "ploidy": pd.to_numeric(a.loc[a_k, "ploidy"], errors="coerce"),
             "genome_doublings": pd.to_numeric(a.loc[a_k, "Genome doublings"],
                                               errors="coerce"),
             "subclonal_fraction": pd.to_numeric(a.loc[a_k, "Subclonal genome fraction"],
                                                 errors="coerce")}
        if c_k in s.index:
            sub = str(s.loc[c_k, "EXPRESSION_SUBTYPE"]).strip().lower()
            idh = str(s.loc[c_k, "IDH1_MUTATION"]).strip().lower()
            gc = str(s.loc[c_k, "G_CIMP_METHYLATION"]).strip().lower()
            r["is_mesenchymal"] = 1.0 if sub.startswith("mesenchymal") else (
                np.nan if sub in ("nan", "", "na") else 0.0)
            # ANY R132 substitution counts. The table carries R132H (28), R132G (2) and
            # R132C (1); matching only R132H would silently call the other three wild-type.
            r["idh1_mutant"] = (1.0 if idh.startswith("r132")
                                else (0.0 if idh == "wt"
                                      else (np.nan if idh in ("nan", "", "na") else 0.0)))
            # "G-CIMP" vs "non-G-CIMP" — and note startswith would match BOTH, so the
            # negative case is tested first.
            r["gcimp"] = (0.0 if gc.startswith("non-g-cimp")
                          else (1.0 if gc.startswith("g-cimp")
                                else (np.nan if gc in ("nan", "", "na") else 0.0)))
        rows[smp] = r
    return pd.DataFrame(rows).T


def main() -> int:
    if not PER_SAMPLE.exists():
        print(f"BLOCKED: {PER_SAMPLE.name} missing. Re-run "
              f"scripts/absolute_purity_yardstick.py (it now persists per-sample estimates).")
        return 2
    est = pd.read_csv(PER_SAMPLE, index_col=0)
    purity = est.pop("absolute_purity")
    methods = list(est.columns)
    print(f"per-sample estimates: {est.shape[0]} samples x {len(methods)} methods")

    cov = build_covariates(list(est.index))
    cov = cov.reindex(est.index)
    print("\nfactor coverage (non-null of %d samples):" % len(cov))
    for f in FACTORS:
        print(f"  {f:20s} {int(cov[f].notna().sum()):4d}"
              + (f"   [{int(cov[f].sum())} positive]" if cov[f].dropna().isin([0, 1]).all()
                 else ""))

    keep = cov.dropna().index
    print(f"\ncomplete cases across all seven factors: {len(keep)} of {len(cov)}")
    if len(keep) < 80:
        print("BLOCKED: too few complete cases for a seven-factor joint model.")
        return 2
    X = cov.loc[keep].astype(float)
    Xz = (X - X.mean()) / X.std(ddof=0).replace(0.0, 1.0)     # standardised, comparable betas
    Xd = np.column_stack([np.ones(len(Xz)), Xz.to_numpy()])

    from ivygap.deconv.registry import build_methods            # noqa: PLC0415
    degen = set()
    try:
        yj = json.loads((config.RESULTS_DIR / "absolute_purity_yardstick.json").read_text())
        degen = {m for m, v in yj.get("methods", {}).items() if v.get("degenerate")}
    except (ValueError, OSError):
        pass
    print(f"degenerate methods on this reference: {sorted(degen) or 'none recorded'}")

    coefs = {}
    for m in methods:
        e = (est.loc[keep, m] - purity.loc[keep]).to_numpy(dtype="float64")
        for outcome, y in (("signed", e), ("abs", np.abs(e))):
            if not np.isfinite(y).all():
                continue
            beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
            coefs.setdefault(outcome, {})[m] = dict(zip(X.columns, beta[1:]))

    results = {}
    print(f"\n{'factor':20s} {'dir':>4s} {'median beta':>12s} {'signs OK':>9s} "
          f"{'sign p':>8s} {'Holm p':>8s}  verdict")
    print("-" * 78)
    raw = {}
    for f, (direction, outcome) in FACTORS.items():
        b = pd.Series({m: coefs[outcome][m][f] for m in coefs[outcome]})
        want_neg = direction == "neg"
        ok = int((b < 0).sum() if want_neg else (b > 0).sum())
        n = len(b)
        p = float(stats.binomtest(ok, n, 0.5).pvalue)
        raw[f] = p
        results[f] = {"direction": direction, "outcome": outcome,
                      "median_beta": round(float(b.median()), 5),
                      "n_methods": n, "n_in_predicted_direction": ok,
                      "sign_test_p": round(p, 5),
                      "per_method_beta": {k: round(float(v), 5) for k, v in b.items()}}
    order = sorted(raw, key=raw.get)
    for i, f in enumerate(order):
        holm = min(1.0, raw[f] * (len(order) - i))
        holm = max(holm, max((min(1.0, raw[g] * (len(order) - j))
                              for j, g in enumerate(order[:i])), default=0.0))
        r = results[f]
        r["holm_p"] = round(holm, 5)
        consistent = r["n_in_predicted_direction"] >= 0.75 * r["n_methods"]
        r["supported"] = bool(holm < ALPHA and consistent)
        v = ("SUPPORTED" if r["supported"] else
             ("wrong direction" if r["n_in_predicted_direction"] < r["n_methods"] / 2
              else "not significant"))
        print(f"{f:20s} {r['direction']:>4s} {r['median_beta']:12.5f} "
              f"{r['n_in_predicted_direction']:4d}/{r['n_methods']:<4d} "
              f"{r['sign_test_p']:8.4f} {holm:8.4f}  {v}")

    hits = [f for f, r in results.items() if r["supported"]]
    print(f"\n{len(hits)} of {len(FACTORS)} factors supported: {hits or 'NONE'}")
    if not hits:
        print("\nThe pre-specified falsification condition is met: deconvolution error in GBM "
              "is NOT predictable from the biology measured here. That is the finding.")

    out = {"prespecification": "prespecified/biological_failure_factors.md",
           "hypothesis": "Biological properties of GBM that alter the relationship between "
                         "malignant-cell abundance and transcriptomic composition will predict "
                         "systematic errors in bulk RNA deconvolution.",
           "n_samples_complete_cases": int(len(keep)), "n_methods": len(methods),
           "alpha": ALPHA, "correction": "Holm across seven primary tests",
           "degenerate_methods": sorted(degen),
           "factors": results, "supported": hits,
           "falsified": not hits}
    (config.RESULTS_DIR / "failure_factors_gbm.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/failure_factors_gbm.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
