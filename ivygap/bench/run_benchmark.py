"""
run_benchmark.py — the method-selection stage.

WHAT THIS STAGE DECIDES, AND WHAT IT DOES NOT
---------------------------------------------
This is where the primary deconvolution method is chosen, and it is the ONLY place a
choice is made. The criterion is accuracy against known composition on donor-held-out
pseudobulk mixtures — an outcome-blind criterion. Survival is not consulted; the
anatomic claims are not consulted; nothing produced downstream can reach back and
change the selection.

The decision, its criterion, its evidence and the full ranking are written to
`method_selection_decision.json`, and `ivygap.survival` refuses to run against a
selection file whose criterion mentions an outcome. The rule is enforced by the code
rather than by the reader's memory of it.

WHAT IT PRODUCES
----------------
    benchmark_per_type.csv     per method, per cell type: error, bias, correlation,
                               false-positive rate on truly-absent types, LoD
    benchmark_summary.csv      one row per method, ranked
    method_selection_decision.json
    equal_footing_certificate.json
    implementation_report.json which of the four published tools ran as the genuine
                               R package and, when not, exactly why
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict

import pandas as pd

from ivygap import config
from ivygap.bench import calibration, metrics, pseudobulk as pb
from ivygap.bench.equal_footing import (build_certificate, check_estimates_aligned,
                                        require_comparable)
from ivygap.data.reference import (build_reference, build_synthetic,
                                   select_signature_genes)
from ivygap.deconv.base import DeconvolutionInput
from ivygap.deconv.controls import CONTROL_NAMES, build_controls
from ivygap.deconv.registry import build_methods, implementation_report


def run(expression: pd.DataFrame, meta: pd.DataFrame, prefer_r: bool = True,
        n_test: int = config.N_TEST_PSEUDOBULK,
        n_signature_genes: int = config.SIGNATURE_GENES_PER_TYPE,
        extra_references: list | None = None, verbose: bool = True) -> dict:
    """
    Run every method against donor-held-out pseudobulk and rank them.

    `expression`/`meta` are cell-level reference data. Donors are split first and the
    reference is rebuilt from training donors ONLY, so no test donor's cells are inside
    the reference any method sees.
    """
    config.ensure_dirs()
    t_start = time.time()

    # --- donor split, then reference from training donors only ---------------
    test_set, train_donors, test_donors = pb.build_train_test(expression, meta, n_test=n_test)
    train_cells = meta.index[meta["donor"].astype(str).isin(train_donors)]
    reference = build_reference(expression[train_cells], meta.loc[train_cells],
                                name=config.PRIMARY_REFERENCE)

    # DONOR LEAKAGE GUARD. The R-backed methods do not read `reference` — they read the
    # cell-level export, which r_bridge writes from whatever cell source is registered.
    # If that source still holds every donor, MuSiC, Bisque and SCDC see the HELD-OUT
    # donors' cells while NNLS, SVR and the rest see only the training signature. The
    # mixtures are pooled from exactly those held-out cells, so the R methods would be
    # scored against their own reference and would win for that reason alone.
    #
    # Registering the training-only cells here keeps the donor split intact and keeps
    # every method on the same evidence. The previous source is restored in the finally
    # block below, so a caller's registration is not clobbered by running a benchmark.
    from ivygap.deconv import r_bridge
    _prior_source = r_bridge.get_cell_source(config.PRIMARY_REFERENCE)
    r_bridge.set_cell_source(config.PRIMARY_REFERENCE,
                             expression[train_cells], meta.loc[train_cells])

    genes = [g for g in select_signature_genes(reference, n_per_type=n_signature_genes)
             if g in test_set.expression.index]
    if len(genes) < config.MIN_GENES_SHARED:
        raise ValueError(
            f"only {len(genes)} genes are shared between the signature and the "
            f"mixtures; {config.MIN_GENES_SHARED} required. Check gene identifiers."
        )

    refs = [reference.subset_genes(genes)]
    for extra in (extra_references or []):
        refs.append(extra.subset_genes([g for g in genes if g in extra.profile.index]))

    manifest = pd.DataFrame(
        {"patient_id": test_set.donors.astype(str),
         "structure": test_set.niche.astype(str)},
        index=test_set.expression.columns,
    )
    data = DeconvolutionInput(bulk=test_set.expression.loc[genes],
                              references=tuple(refs), manifest=manifest)

    if verbose:
        print(f"donors: {len(train_donors)} train / {len(test_donors)} held out")
        print(f"mixtures: {len(test_set)} | signature genes: {len(genes)} | "
              f"references: {[r.name for r in refs]}")

    # --- run every method on the identical input -----------------------------
    estimates: dict[str, pd.DataFrame] = {}
    implementations: dict[str, str] = {}
    timings: dict[str, float] = {}
    failures: dict[str, str] = {}
    fitted: list = []                      # kept so calibration can read posteriors

    # Controls run here as well as in the anatomic stage: the agreement test needs
    # both rankings over the same set, and a control's accuracy against KNOWN truth
    # is the calibration for what 'meaningless' scores on this benchmark.
    for method in build_methods(prefer_r=prefer_r) + build_controls():
        t0 = time.time()
        try:
            estimates[method.name] = method.fit_predict(data)
            fitted.append(method)
            timings[method.name] = time.time() - t0
            implementations[method.name] = getattr(
                method, "implementation_", "python") or "python"
            if verbose:
                note = ""
                reason = getattr(method, "fallback_reason_", None)
                if reason:
                    note = f"  [fell back to Python: {reason[:60]}]"
                print(f"  {method.name:24s} {timings[method.name]:7.1f}s{note}")
        except Exception as exc:                          # noqa: BLE001
            # A method that fails is recorded as failed, never silently dropped: a
            # ranking that quietly omits a method reads as though it lost.
            failures[method.name] = f"{type(exc).__name__}: {exc}"
            if verbose:
                print(f"  {method.name:24s} FAILED — {type(exc).__name__}: {exc}")

    if not estimates:
        raise RuntimeError(f"every method failed: {failures}")

    # --- prove they faced the same problem, then score ------------------------
    check_estimates_aligned(estimates)
    certificate = build_certificate(data, methods=list(estimates),
                                    implementations=implementations)
    require_comparable({name: certificate for name in estimates})

    summary = metrics.summarise_methods(estimates, test_set.truth, groups=test_set.donors)

    per_type = []
    for name, est in estimates.items():
        block = metrics.per_type_metrics(est, test_set.truth, groups=test_set.donors)
        block.insert(0, "method", name)
        per_type.append(block)
    per_type_df = pd.concat(per_type, ignore_index=True)

    # Per-niche breakdown: a method can win overall and still fail on leading-edge-like
    # mixtures, which are the ones this project most depends on.
    per_niche = []
    for name, est in estimates.items():
        for niche in sorted(test_set.niche.unique()):
            idx = test_set.niche.index[test_set.niche == niche]
            if len(idx) < 5:
                continue
            m = metrics.overall_metrics(est.loc[idx], test_set.truth.loc[idx],
                                        groups=test_set.donors.loc[idx])
            per_niche.append({"method": name, "niche": niche, **m})
    per_niche_df = pd.DataFrame(per_niche)

    # Uncertainty calibration: for every method that states a credible interval, how
    # often does the interval actually contain the truth? Measured, never assumed —
    # this is what caught the Bayesian sampler's original overconfident default.
    calib_df = calibration.calibration_table(fitted, test_set.truth)

    # Is the top row actually winning, or is the leaderboard reporting noise? A paired
    # bootstrap over donors answers that. The survival stage already refuses to rank an
    # underpowered comparison; applying a weaker standard to the criterion that SELECTS
    # the method would be backwards.
    ties_df = metrics.paired_bootstrap_ties(estimates, test_set.truth,
                                            groups=test_set.donors)

    selected = str(summary.index[0])

    tied = ([m for m in ties_df.index if bool(ties_df.loc[m, "tied_with_best"])]
            if not ties_df.empty else [selected])
    real_tied = [m for m in tied if m not in CONTROL_NAMES]
    if len(real_tied) > 1:
        # Selection still has to return one method — downstream stages need something
        # frozen — but the decision must record that it was a tie-break rather than a
        # measured win, and name who it was tied with.
        selection_basis = (
            f"TIE-BREAK: {len(real_tied)} methods are statistically indistinguishable "
            f"from the best on this benchmark ({', '.join(sorted(real_tied))}). "
            f"{selected} was taken as the lowest point estimate, which is a convention, "
            f"not evidence that it is better. Any downstream claim that {selected} is "
            f"the best method is unsupported."
        )
    else:
        selection_basis = (
            f"{selected} is distinguishable from every other method by a paired "
            f"bootstrap over donors."
        )

    decision = {
        "selected_method": selected,
        "selection_basis": selection_basis,
        "methods_tied_with_best": sorted(real_tied),
        "criterion": "mean absolute error on primary cell types, donor-held-out "
                     "pseudobulk, aggregated donor-equally",
        "criterion_is_outcome_blind": True,
        "explicitly_not_used": ["survival", "anatomic structure labels",
                                "the frozen biology claims"],
        "n_test_mixtures": int(len(test_set)),
        "train_donors": train_donors,
        "test_donors": test_donors,
        "n_signature_genes": len(genes),
        "ranking": summary.reset_index().to_dict(orient="records"),
        "failures": failures,
        "implementations": implementations,
        "runtime_seconds": {k: round(v, 2) for k, v in timings.items()},
        "total_runtime_seconds": round(time.time() - t_start, 1),
    }

    # --- write everything -----------------------------------------------------
    out = config.BENCH_DIR
    summary.to_csv(out / "benchmark_summary.csv")
    per_type_df.to_csv(out / "benchmark_per_type.csv", index=False)
    if not per_niche_df.empty:
        per_niche_df.to_csv(out / "benchmark_per_niche.csv", index=False)
    calib_df.to_csv(out / "uncertainty_calibration.csv")
    if not ties_df.empty:
        ties_df.to_csv(out / "method_ties.csv")
    (out / "method_selection_decision.json").write_text(json.dumps(decision, indent=2))
    (out / "equal_footing_certificate.json").write_text(certificate.to_json())
    (out / "implementation_report.json").write_text(
        json.dumps(implementation_report(config.PRIMARY_REFERENCE), indent=2))

    for name, est in estimates.items():
        est.to_csv(config.ESTIMATES_DIR / f"pseudobulk_{name}.csv")
    test_set.truth.to_csv(config.ESTIMATES_DIR / "pseudobulk_truth.csv")

    if verbose:
        print("\n--- ranking (lower MAE is better) ---")
        print(summary.round(4).to_string())
        reported = calib_df[calib_df.get("uncertainty_reported", False) == True]
        if not reported.empty:
            print("\n--- uncertainty calibration (nominal coverage 0.95) ---")
            print(reported[["coverage_95", "mean_ci_width", "interval_score",
                            "verdict"]].round(3).to_string())
        print(f"\nselected: {selected}")
        print(selection_basis)

    # Restore whatever cell source the caller had registered, so a benchmark run does
    # not silently leave the training-only subset in place for later stages.
    if _prior_source is None:
        r_bridge.clear_cell_source(config.PRIMARY_REFERENCE)
    else:
        r_bridge.set_cell_source(config.PRIMARY_REFERENCE, *_prior_source)

    return {"summary": summary, "per_type": per_type_df, "per_niche": per_niche_df,
            "calibration": calib_df, "ties": ties_df, "decision": decision,
            "certificate": certificate, "estimates": estimates,
            "truth": test_set.truth, "reference": reference, "genes": genes,
            # Returned at the top level, not only inside `decision`, because callers
            # need the split to keep later stages on the same reference the selection
            # was made on. Reaching into decision["train_donors"] for that is a trap.
            "train_donors": train_donors, "test_donors": test_donors}


def main(prefer_r: bool = True, synthetic: bool = False) -> int:
    """
    CLI entry point.

    `--synthetic` runs the whole benchmark against a generated reference, which is how
    the pipeline is validated end-to-end before any multi-gigabyte download. Results
    from that path are labelled synthetic everywhere they touch disk.
    """
    if synthetic:
        print("SYNTHETIC FIXTURE RUN — results are not real biological findings\n")
        _, expression, meta = build_synthetic(
            n_genes=1200, n_donors=10, cells_per_type_per_donor=30)
    else:
        from ivygap.data.reference import build_from_h5ad
        path = config.REFERENCE_DIR / "gbmap_core.h5ad"
        if not path.exists():
            print(f"ERROR: {path} not found.\n"
                  "Download a GBM single-cell atlas to that path, or run with "
                  "--synthetic to validate the pipeline without it.")
            return 1
        _, expression, meta = build_from_h5ad(path, export=False)
        from ivygap.deconv import r_bridge
        r_bridge.set_cell_source("gbmap", expression, meta)
    run(expression, meta, prefer_r=prefer_r)
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--synthetic", action="store_true",
                    help="run on a generated reference (fixture validation)")
    ap.add_argument("--no-r", action="store_true",
                    help="force the all-Python implementations")
    a = ap.parse_args()
    raise SystemExit(main(prefer_r=not a.no_r, synthetic=a.synthetic))
