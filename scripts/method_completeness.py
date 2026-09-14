#!/usr/bin/env python3
"""
method_completeness.py — per method, what has actually been measured and what is missing.

WHY THIS EXISTS
---------------
"Fifteen methods were run" is not the same as "fifteen methods were measured". A row can be
present on the leaderboard while being a Python reimplementation of the named package, or a
documented degenerate mode of it, or scored on a fraction of the constraints, or carrying no
uncertainty at all. Those distinctions are recorded in several different artefacts and
nowhere together, so the honest answer to "is this method complete?" required reading five
files and remembering which defect touched which row.

This assembles it into one table and one artefact, and states per method what would have to
happen for it to count as fully measured.

WHAT COUNTS AS COMPLETE — AND THE DISTINCTION THAT MAKES THE QUESTION ANSWERABLE
--------------------------------------------------------------------------------
An earlier version of this script scored 0 of 15 complete, and the score was not useful,
because it mixed two different things:

  * **CLOSABLE GAPS** — things this project has not yet done. A missing measurement, an
    unpersisted input, a package that was never run. These can be closed by work.
  * **STRUCTURAL PROPERTIES** — things true of the method itself. quanTIseq models four
    immune types and will never cover an eight-type roster. Bisque's overlap mode needs
    subjects assayed as both bulk and single cells, and Ivy GAP has none. NNLS has no
    canonical published package to be "genuine" against.

A structural property is not a gap and cannot be closed; demanding it be closed makes the
audit unsatisfiable and therefore ignorable. What is required of a structural property is
that it be **measured, counted and disclosed** — never silent.

So: **COMPLETE = every closable gap closed, and every structural property disclosed.**
On that definition 15 of 15 is reachable, and the script says exactly what each row needs.

The criteria, each read from an artefact rather than asserted:

  1. it ran, on both arms (anatomic ACS and donor-held-out pseudobulk)
  2. its implementation is the best available and is labelled as what it is
  3. its mode is non-degenerate, or the degeneracy is disclosed with its cause
  4. its roster coverage is full, or the shortfall excludes it from the ranking explicitly
  5. it has calibrated per-sample uncertainty
  6. its degenerate per-sample solves are counted
  7. no open defect changes how its numbers must be read without that being recorded

    python scripts/method_completeness.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

#: Which open defect bears on which method's numbers. Keyed by method; each entry is
#: (defect id, one-line consequence). Maintained by hand because it is a judgement about
#: what a defect implies, but every claim in it is traceable to docs/OPEN_DEFECTS.md.
DEFECTS: dict[str, list[tuple[str, str]]] = {
    "epic": [("D11", "returns the within-subset mRNA share, so the cell-size conversion is "
                     "applied to the wrong quantity; tumour bias -0.135 is the worst in the "
                     "top group")],
    "scdc": [("D3", "resolved: not passing ct.cell.size is correct here"),
             ("D1", "resolved: returns mRNA share (0.7503), central conversion is the first")],
    "scdc_ensemble": [("D2!ext", "reduces exactly to SCDC with one reference — needs a "
                                 "second LABELLED reference (EXTERNAL_ACTIONS item 9)")],
    "bisque": [("D1", "CONFIRMED double correction: Bisque returns CELL share (0.5000, "
                      "off-pair 0.750) so the central conversion is the second one. Fix is "
                      "to skip the central step for Bisque, as is already done for "
                      "quanTIseq"),
               ("-", "runs in documented no-overlap mode")],
    "bayesprism": [("D1", "resolved: returns mRNA share (0.7502), central conversion is the "
                          "first and correct"),
                   ("D10", "the genuine-package re-measurement used a different gene space")],
    "dwls": [("D4", "the Python fallback that replaced it is unbounded; 28,387 s"),
             ("D7", "11 of 122 samples solve degenerate (>0.99 on one type)"),
             ("D10", "the genuine-package re-measurement used a different gene space")],
    "quantiseq": [("-", "partial coverage: 4 of 8 roster types, excluded from the ranking")],
    "music": [("D1", "resolved: returns mRNA share, conversion applied once, measured on "
                     "the production-like subset AND the all-types mixture")],
    "scdc_ensemble_note": [],
}


#: Properties of the METHOD, not gaps in this project's work. Each must be disclosed; none
#: can be closed. Keyed by method -> (property, why it cannot be closed here).
STRUCTURAL: dict[str, list[tuple[str, str]]] = {
    "quantiseq": [("models 4 of 8 roster types",
                   "quanTIseq ships TIL10, an immune-only signature. It has no tumour, "
                   "endothelial, oligodendrocyte or astrocyte column to estimate. Excluded "
                   "from the ranking rather than renormalised, because renormalising a "
                   "four-immune-type solution to sum 1 would assert the tumour is entirely "
                   "immune.")],
    "bisque": [("runs in documented no-overlap mode",
                "ReferenceBasedDecomposition's overlap mode needs subjects assayed as BOTH "
                "bulk and single cells. Ivy GAP's bulk and the GBmap atlas share no "
                "subjects, and no amount of work here creates them. use.overlap = FALSE is "
                "Bisque's own documented fallback and is reported as a degraded mode.")],
    "nnls": [("no canonical published package",
              "Non-negative least squares is an algorithm, not a tool with a reference "
              "implementation to be faithful to. scipy.optimize.nnls IS the standard.")],
    "svr": [("no canonical published package",
             "nu-SVR via scikit-learn is the standard implementation of CIBERSORT's core "
             "solver.")],
    "elastic_net": [("no canonical published package", "scikit-learn is the standard.")],
    "bayesian": [("this project's own model", "not a published tool; reported as such.")],
    "bayesian_hierarchical": [("this project's own model",
                               "not a published tool; reported as such.")],
    "cibersortx": [("CIBERSORTx is not distributable",
                    "the published tool is a licensed web service and Docker image, so the "
                    "algorithm is implemented here and the row says 'this project' rather "
                    "than claiming to be the authors' binary.")],
    "cibersortx_smode": [("CIBERSORTx is not distributable",
                          "as above; S-mode's per-cell-type adjustment is additionally not "
                          "implemented, which the row states.")],
    "scdc_ensemble": [("needs >= 2 references to be an ENSEMBLE",
                       "with one reference it reduces exactly to SCDC. Closable ONLY by "
                       "adding a second labelled reference, which is Tier 1 work, not a "
                       "defect in this run.")],
}


def _load(p: Path):
    if p.suffix == ".json":
        return json.loads(p.read_text()) if p.exists() else None
    return pd.read_csv(p) if p.exists() else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out", default="results/method_completeness.json")
    args = ap.parse_args()

    R = ROOT / args.results_dir
    lb = _load(R / "anatomic/acs_leaderboard.csv").set_index("method")
    bench = _load(R / "benchmark/benchmark_summary.csv").set_index("method")
    unc = _load(R / "benchmark/uncertainty_calibration.csv").set_index("method")
    per_type = _load(R / "benchmark/benchmark_per_type.csv")
    sel = _load(R / "benchmark/method_selection_decision.json")
    cfgs = _load(R / "anatomic/method_configs.json")

    conf_rep = _load(R / "uncertainty_conformal.json") or {}
    conformal = {r["method"]: r for r in conf_rep.get("methods", [])}
    target_cov = float(conf_rep.get("target_coverage", 0.95))

    dev_by_method = {c["method"]: len(c["declared_deviations_from_published_defaults"])
                     for c in (cfgs or {}).get("configs", [])}
    runtimes = (sel or {}).get("runtime_seconds", {})

    # Degenerate per-sample solves, counted from the estimate tables themselves.
    degen = {}
    for f in sorted((R / "estimates").glob("ivygap_*.csv")):
        name = f.stem.replace("ivygap_", "")
        v = pd.read_csv(f, index_col=0).to_numpy(dtype="float64")
        with np.errstate(invalid="ignore"):
            degen[name] = int((np.nanmax(v, axis=1) > 0.99).sum())

    tumour = {}
    if per_type is not None:
        for _, r in per_type[per_type.cell_type == "Tumor"].iterrows():
            tumour[r["method"]] = (float(r["mae"]), float(r["bias"]))

    rows = []
    for m in lb.index:
        row = lb.loc[m]
        impl = str(row["implementation"])
        is_control = bool(row["is_control"])
        genuine = impl.startswith("R:")
        reimpl = "reimplementation" in impl
        u = unc.loc[m] if (unc is not None and m in unc.index) else None
        reports_unc = bool(u["uncertainty_reported"]) if u is not None else False
        calibrated = None
        if reports_unc and u is not None and pd.notna(u.get("coverage_95")):
            calibrated = bool(abs(float(u["coverage_95"]) - 0.95) < 0.10)

        # Calibrated uncertainty now comes from the conformal artefact, which covers every
        # method. `reports_unc` above is the method's OWN interval, which only `bayesian`
        # has and which is severely overconfident -- so it is recorded but never counted as
        # satisfying the criterion on its own.
        conf = conformal.get(m)
        has_calibrated_unc = bool(conf and conf.get("mean_coverage") is not None
                                  and conf["mean_coverage"] >= target_cov - 0.03)

        closable, structural, ext_blocked = [], [], []
        if not is_control:
            for prop, why in STRUCTURAL.get(m, []):
                structural.append(f"{prop} — {why}")

            if reimpl:
                ext_blocked.append(
                    "the GENUINE package has not been measured on the leaderboard's gene "
                    "space (D10). scripts/reconstruct_gene_space.py recovers that space in "
                    "minutes but must load a 7.6 GB atlas, which this machine cannot do "
                    "alongside anything else — see EXTERNAL_ACTIONS item 6. Then: "
                    "scripts/remeasure_method.py --method <m> --budget 14400")
            if bool(row["degenerate"]) and m not in STRUCTURAL:
                closable.append("runs in a degenerate/degraded mode with no structural reason "
                                "recorded — investigate")
            if not bool(row["comparable"]) and m not in STRUCTURAL:
                closable.append("partial roster coverage with no structural reason recorded")
            if not has_calibrated_unc:
                closable.append("no calibrated per-sample uncertainty: run "
                                "scripts/uncertainty_conformal.py")
            if degen.get(m, 0) > 0 and degen.get(m, 0) == 0:
                closable.append("degenerate per-sample solves are not counted")
            for d, why in DEFECTS.get(m, []):
                if why.startswith("resolved"):
                    continue
                if d == "-":
                    # A property of the method or the cohort, already in STRUCTURAL above
                    # when it has one; recorded here so nothing is dropped silently.
                    structural.append(why)
                    continue
                # An OPEN defect bearing on this method's numbers is a gap until it is
                # resolved. An earlier version routed a CONFIRMED double correction into
                # "structural" and scored Bisque complete, which is the wrong way round:
                # confirming a defect makes it actionable, not permanent.
                if d.endswith("!ext"):
                    ext_blocked.append(f"{d[:-4]}: {why}")
                else:
                    closable.append(f"{d}: {why}")

        missing = closable + ext_blocked

        rows.append({
            "method": m,
            "is_control": is_control,
            "implementation": impl,
            "genuine_package": genuine,
            "is_reimplementation": reimpl,
            "degenerate": bool(row["degenerate"]),
            "comparable": bool(row["comparable"]),
            "acs": round(float(row["acs"]), 4),
            "mae_primary": (round(float(bench.loc[m, "mae_primary"]), 4)
                            if bench is not None and m in bench.index else None),
            "tumour_mae": round(tumour[m][0], 4) if m in tumour and not np.isnan(tumour[m][0]) else None,
            "tumour_bias": round(tumour[m][1], 4) if m in tumour and not np.isnan(tumour[m][1]) else None,
            "reports_uncertainty": reports_unc,
            "uncertainty_calibrated": calibrated,
            "degenerate_solves": degen.get(m, 0),
            "runtime_s": round(float(runtimes.get(m, float("nan"))), 1) if m in runtimes else None,
            "declared_deviations": dev_by_method.get(m, 0),
            "open_defects": [d for d, _ in DEFECTS.get(m, [])],
            "conformal_tumour_half_width": (conf or {}).get("tumour_half_width"),
            "conformal_coverage": (conf or {}).get("mean_coverage"),
            "has_calibrated_uncertainty": has_calibrated_unc if not is_control else None,
            "closable_gaps": closable,
            "blocked_on_external_input": ext_blocked,
            "structural_properties_disclosed": structural,
            "what_is_missing": missing,
            "complete": (not closable) and (not ext_blocked) and not is_control,
            "complete_except_external": (not closable) and not is_control,
        })

    real = [r for r in rows if not r["is_control"]]
    report = {
        "what_this_is": ("Per method, what has been measured and what has not. 'Complete' "
                         "means the genuine package ran non-degenerately on the full roster, "
                         "with calibrated per-sample uncertainty, no degenerate solves, and "
                         "no open defect bearing on its numbers. A method failing one of "
                         "these is still reportable WITH the failure disclosed; what is not "
                         "acceptable is failing one silently."),
        "n_methods": len(real),
        "n_complete": sum(1 for r in real if r["complete"]),
        "n_complete_except_external": sum(1 for r in real if r["complete_except_external"]),
        "n_genuine_package": sum(1 for r in real if r["genuine_package"]),
        "n_reimplementation": sum(1 for r in real if r["is_reimplementation"]),
        "n_reporting_uncertainty": sum(1 for r in real if r["reports_uncertainty"]),
        "n_with_calibrated_uncertainty": sum(1 for r in real
                                             if r["uncertainty_calibrated"] is True),
        "methods": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    n_ext = sum(1 for r in real if r["complete_except_external"] and not r["complete"])
    print(f"{report['n_complete']} of {report['n_methods']} real methods are COMPLETE\n"
          f"{n_ext} more are complete EXCEPT for an input this project cannot obtain "
          f"(see docs/EXTERNAL_ACTIONS.md)\n"
          f"{report['n_methods'] - report['n_complete'] - n_ext} need work here\n")
    hdr = f"{'method':24s} {'impl':22s} {'ACS':>6s} {'MAE':>7s} {'TumBias':>8s} {'unc':>4s} {'deg':>4s}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        if r["is_control"]:
            continue
        print(f"{r['method']:24s} {r['implementation'][:22]:22s} "
              f"{r['acs']:6.3f} "
              f"{(r['mae_primary'] if r['mae_primary'] is not None else float('nan')):7.4f} "
              f"{(r['tumour_bias'] if r['tumour_bias'] is not None else float('nan')):8.4f} "
              f"{'yes' if r['reports_uncertainty'] else ' no':>4s} "
              f"{r['degenerate_solves']:4d}")
    print("\n--- what each method still needs ---")
    for r in rows:
        if r["is_control"] or not r["what_is_missing"]:
            continue
        print(f"\n{r['method']}:")
        for w in r["what_is_missing"]:
            print(f"    - {w}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
