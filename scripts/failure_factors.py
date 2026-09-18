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

import argparse
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

#: LGG's testable subset. Definitions, directions and model are those fixed in the
#: pre-specification; only the SET is smaller, and the two absentees are declared rather
#: than quietly dropped:
#:
#:   F5 (Verhaak Mesenchymal) does not exist for LGG -- lower-grade glioma is classified by
#:      IDH status and 1p/19q -- so it is carried by the SECONDARY continuous MES score
#:      (addendum 2) and reported separately, never as the pre-specified test.
#:   F7 (G-CIMP) is a methylation phenotype MC3 does not carry, and correlated with IDH1 at
#:      +0.868 in GBM, so it is not separately informative. Reported NOT-TESTABLE, not null.
#:
#: Holm therefore runs across FIVE primary tests in LGG rather than seven. Fewer tests is a
#: weaker correction, so this is stated rather than left for a reader to infer.
FACTORS_LGG = {
    "purity":             ("neg", "signed"),
    "ploidy":             ("pos", "abs"),
    "genome_doublings":   ("pos", "abs"),
    "subclonal_fraction": ("pos", "abs"),
    "idh1_mutant":        ("pos", "abs"),
}
SECONDARY_LGG = {"mes_score": ("neg", "signed")}


def k4(s): return "-".join(str(s).split("-")[:4])
def k3(s): return "-".join(str(s).split("-")[:3])


def k4s(s):
    """
    Barcode key with the ALIQUOT LETTER STRIPPED.

    Expression columns are `TCGA-FG-6692-01A`; the MC3 matrix uses `TCGA-CS-4938-01`. A naive
    four-field join shares ZERO samples between them and would have produced an empty IDH
    column silently rather than failing. Stripping the letter gives 511.
    """
    f = str(s).split("-")
    return "-".join(f[:3] + [f[3][:2]]) if len(f) > 3 else str(s)


def build_covariates_lgg(samples: list[str]) -> pd.DataFrame:
    """
    LGG covariates. Four factors come from ABSOLUTE exactly as in GBM; IDH1 comes from the
    MC3 gene-level mutation matrix; the mesenchymal variable is the SECONDARY continuous score
    (addendum 2), because the Verhaak call does not exist for LGG.

    F7 (G-CIMP) is NOT constructed. It is a methylation phenotype and MC3 does not carry it —
    and it is not separately informative anyway: in GBM it correlated with IDH1 at +0.868, so
    testing IDH1 properly here addresses that axis. Reported as not-testable rather than null.
    """
    a = pd.read_csv(ABS_T, sep="\t").dropna(subset=["Cancer DNA fraction"])
    a["k"] = a["sample"].map(k4s)
    a = a.drop_duplicates("k").set_index("k")
    mc3 = pd.read_csv(ROOT / "TCGA_LGG" / "mc3_gene_level_LGG_mc3_gene_level.txt",
                      sep="\t", index_col=0)
    mc3.columns = [k4s(c) for c in mc3.columns]
    idh = mc3.loc["IDH1"] if "IDH1" in mc3.index else None
    mes = pd.read_csv(config.RESULTS_DIR / "mes_score_lgg.csv", index_col=0)
    mes.index = [k4s(i) for i in mes.index]
    mes = mes[~mes.index.duplicated()]

    rows = {}
    for smp in samples:
        k = k4s(smp)
        if k not in a.index:
            continue
        r = {"purity": pd.to_numeric(a.loc[k, "purity"], errors="coerce"),
             "ploidy": pd.to_numeric(a.loc[k, "ploidy"], errors="coerce"),
             "genome_doublings": pd.to_numeric(a.loc[k, "Genome doublings"], errors="coerce"),
             "subclonal_fraction": pd.to_numeric(a.loc[k, "Subclonal genome fraction"],
                                                 errors="coerce"),
             "idh1_mutant": (float(idh[k]) if idh is not None and k in idh.index else np.nan),
             "mes_score": (float(mes.loc[k, "MES_minus_mean_other"])
                           if k in mes.index else np.nan)}
        rows[smp] = r
    return pd.DataFrame(rows).T


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
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["gbm", "lgg"], default="gbm")
    args = ap.parse_args()
    cohort = args.cohort
    global PER_SAMPLE, FACTORS
    if cohort == "lgg":
        PER_SAMPLE = config.RESULTS_DIR / "absolute_purity_per_sample_lgg.csv"
        FACTORS = {**FACTORS_LGG, **SECONDARY_LGG}
    if not PER_SAMPLE.exists():
        print(f"BLOCKED: {PER_SAMPLE.name} missing. Re-run "
              f"scripts/absolute_purity_yardstick.py (it now persists per-sample estimates).")
        return 2
    est = pd.read_csv(PER_SAMPLE, index_col=0)
    purity = est.pop("absolute_purity")
    methods = list(est.columns)
    print(f"per-sample estimates: {est.shape[0]} samples x {len(methods)} methods")

    cov = (build_covariates_lgg(list(est.index)) if cohort == "lgg"
           else build_covariates(list(est.index)))
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
    primary = set(FACTORS_LGG) if cohort == "lgg" else set(FACTORS)
    order = sorted([f for f in raw if f in primary], key=raw.get)
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

    # --- F1 IS PARTLY DEFINITIONAL, AND THE CORRECT STATISTIC IS THE RECOVERY FRACTION ---
    #
    # The outcome is e = estimate - purity, so regressing e on purity is mechanically negative:
    # a method carrying NO information about tumour content gives slope exactly -1.0, and only
    # a method tracking purity perfectly gives 0. The sign test therefore cannot discover that
    # "error depends on purity" -- that is arithmetic, not biology.
    #
    # What is informative is the MAGNITUDE. `recovered = 1 + slope` is the fraction of the true
    # purity variation the method actually captures. This is reported in place of F1's sign
    # test, and F1's sign test is marked as uninformative rather than deleted.
    #
    # Factors that are NOT subtracted from the outcome -- F2..F7 -- do not have this problem,
    # and because the model is joint they are already adjusted for purity.
    rec = {}
    for m in methods:
        e = (est.loc[keep, m] - purity.loc[keep]).to_numpy(dtype="float64")
        slope = float(np.polyfit(purity.loc[keep].to_numpy(dtype="float64"), e, 1)[0])
        rec[m] = {"slope": round(slope, 4), "recovered_fraction": round(1 + slope, 4),
                  "pearson_r_with_purity": round(float(np.corrcoef(
                      est.loc[keep, m], purity.loc[keep])[0, 1]), 4)}
    med = float(np.median([v["recovered_fraction"] for v in rec.values()]))
    print(f"\nF1 RECOVERY (the interpretable form): methods recover a median "
          f"{med:.1%} of the true purity variation")
    for m, v in sorted(rec.items(), key=lambda kv: kv[1]["recovered_fraction"]):
        print(f"  {m:16s} slope {v['slope']:+.3f}  recovers {v['recovered_fraction']:6.1%}"
              f"  r = {v['pearson_r_with_purity']:+.3f}")
    results["purity"]["sign_test_is_definitional"] = True
    results["purity"]["interpretable_statistic"] = {
        "what": "fraction of true purity variation recovered, = 1 + slope of (estimate - "
                "purity) on purity. 0 means no information, 1 means perfect.",
        "median_recovered_fraction": round(med, 4), "per_method": rec}

    # --- pre-specified repeat without degenerate methods ------------------------------
    nd = [m for m in methods if m not in degen]
    print(f"\nPRE-SPECIFIED REPEAT without the {len(degen)} degenerate methods "
          f"({len(nd)} remain):")
    for f, (direction, outcome) in FACTORS.items():
        b = pd.Series({m: coefs[outcome][m][f] for m in nd if m in coefs[outcome]})
        ok = int((b < 0).sum() if direction == "neg" else (b > 0).sum())
        p_nd = float(stats.binomtest(ok, len(b), 0.5).pvalue)
        results[f]["without_degenerate"] = {
            "n_methods": len(b), "n_in_predicted_direction": ok,
            "sign_test_p": round(p_nd, 5), "median_beta": round(float(b.median()), 5)}
        same = (ok >= 0.75 * len(b)) == (results[f]["n_in_predicted_direction"]
                                         >= 0.75 * results[f]["n_methods"])
        print(f"  {f:20s} {ok:2d}/{len(b):<2d} p={p_nd:.4f}  "
              f"{'consistent with the primary' if same else 'DIFFERS from the primary'}")

    for f in raw:
        if f not in primary:                      # secondary: reported, never Holm-corrected
            results[f]["holm_p"] = None
            results[f]["secondary_not_prespecified_variable"] = True
            results[f]["supported"] = False
            r = results[f]
            print(f"{f:20s} {r['direction']:>4s} {r['median_beta']:12.5f} "
                  f"{r['n_in_predicted_direction']:4d}/{r['n_methods']:<4d} "
                  f"{r['sign_test_p']:8.4f} {'--':>8s}  SECONDARY (not the pre-specified "
                  f"variable)")

    hits = [f for f, r in results.items() if r["supported"]]
    print(f"\n{len(hits)} of {len(FACTORS)} factors supported: {hits or 'NONE'}")
    if not hits:
        print("\nThe pre-specified falsification condition is met: deconvolution error in GBM "
              "is NOT predictable from the biology measured here. That is the finding.")

    out = {"cohort": cohort,
           "prespecification": "prespecified/biological_failure_factors.md",
           "hypothesis": "Biological properties of GBM that alter the relationship between "
                         "malignant-cell abundance and transcriptomic composition will predict "
                         "systematic errors in bulk RNA deconvolution.",
           "n_samples_complete_cases": int(len(keep)), "n_methods": len(methods),
           "alpha": ALPHA, "correction": "Holm across seven primary tests",
           "degenerate_methods": sorted(degen),
           "factors": results, "supported": hits,
           "f1_caveat": "F1's sign test is DEFINITIONAL -- e = estimate - purity, so the "
                        "slope is mechanically -1 for a method with no information. The "
                        "interpretable statistic is the recovery fraction, reported under "
                        "factors.purity.interpretable_statistic.",
           "falsified": not hits}
    (config.RESULTS_DIR / f"failure_factors_{cohort}.json").write_text(
        json.dumps(out, indent=2))
    print(f"\nwrote results/failure_factors_{cohort}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
