#!/usr/bin/env python3
"""
remeasure_method.py — run a GENUINE R package on the anatomic cohort, unbounded.

WHY THIS EXISTS
---------------
DWLS is the one method this study has never cleanly measured. In the confirmatory run it
fell back to the Python reimplementation in both anatomic stages, having exceeded its
2,400 s budget for the R package. So the leaderboard's DWLS row is not DWLS, and the
project's own invariant is explicit that a reimplementation must never be reported under
a published package's name.

That has a concrete cost beyond tidiness. Avila Cobos et al. (2020) rank DWLS **best**
among methods that use a single-cell reference; this study ranks it last. The
disagreement is currently unresolvable, because the two are not testing the same
software. It also breaks the benchmark-concordance test in
`scripts/benchmark_concordance.py`, where DWLS at rank 14 is most of why that comparison
returns a null result.

WHY RAISING THE BUDGET HERE IS NOT OUTCOME-DRIVEN
--------------------------------------------------
The obvious objection: the budget is being changed after seeing that DWLS scored badly.
Three things make this a re-measurement rather than a retune, and all three have to hold:

  1. **The budget's purpose is to bound a 15-method run, not to exclude a method.** A
     method the budget removes *systematically* is a measurement failure. DWLS completed
     as R:DWLS in the 2026-09-05 run and fell back in the two after it, so the budget is
     sitting almost exactly on this method's runtime — the worst possible place for a
     threshold to sit, because it decides the answer by coin flip.
  2. **Nothing about the method changes.** No parameter, no gene set, no cutoff. The
     published defaults are untouched. Only wall-clock is allowed to differ, and
     wall-clock is not a property of the method's accuracy.
  3. **The result is reported whatever it is.** If genuine DWLS scores lower than the
     reimplementation, that is the finding and it goes in the leaderboard. This script
     writes its output whichever way it falls, and refuses to be run selectively.

The alternative — leaving it — is worse than either outcome, because it means publishing
"DWLS is last" about software that was never run.

WHAT THIS DOES NOT DO
---------------------
It does not touch the archived confirmatory run, which stands as recorded. It produces a
separate, clearly labelled measurement of one method, to be reported alongside rather
than substituted in.

    python scripts/remeasure_dwls.py                    # 4 hour budget
    python scripts/remeasure_dwls.py --budget 28800     # 8 hours
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                         # noqa: E402
from ivygap.anatomic.acs import score as acs_score                # noqa: E402
from ivygap.data.load_ivygap import load_cached                   # noqa: E402
from ivygap.data.reference import build_from_h5ad                 # noqa: E402
from ivygap.deconv import r_bridge                                # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                 # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", default="dwls",
                    help="method to re-measure, e.g. dwls or bayesprism")
    ap.add_argument("--budget", type=int, default=14400, help="seconds for the R call")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cores", type=int, default=None,
                    help="workers for a method that parallelises with an R SOCKET cluster "
                         "(BayesPrism). 1 disables the cluster. Omit for the package default. "
                         "See OPEN_DEFECTS D18: on a memory-constrained machine the default "
                         "spawns workers that are killed, the master fails with "
                         "'unserialize(node$con)', and the bridge silently falls back to the "
                         "Python reimplementation.")
    args = ap.parse_args()
    M = args.method
    if args.cores is not None:
        # Set on config so it reaches the R side through the run's config JSON and is
        # recorded there, rather than being applied invisibly.
        config.R_SOCKET_CLUSTER_CORES = int(args.cores)
        print(f"R socket-cluster workers forced to {args.cores} "
              f"(declared deviation, OPEN_DEFECTS D18)")
    if args.out is None:
        args.out = f"results/{M}_remeasured.json"

    # Fail fast, before the multi-GB atlas load. A re-measurement that cannot be compared
    # is not worth an hour of CPU, and finding that out after the hour is worse.
    genes_path = config.BENCH_DIR / "signature_genes.json"
    if not genes_path.exists():
        print(f"BLOCKED: {genes_path} does not exist.\n\n"
              f"The leaderboard's gene space is not recoverable from the artefacts on "
              f"disk — only its COUNT was ever recorded. Running anyway would measure a "
              f"different gene space than every method it would be compared against, "
              f"which is exactly what produced the 1,591-vs-657 mismatch in the DWLS and "
              f"BayesPrism re-measurements (docs/OPEN_DEFECTS.md D10).\n\n"
              f"To resume: complete one full `python scripts/run_all.py`, which now "
              f"writes the gene list, then re-run this script.")
        return 2

    print("building the anatomic cohort and the cell-level reference...")
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]

    # RESTRICT THE ATLAS TO THE BULK'S GENES, exactly as run_all.py does:
    #     build_from_h5ad(gbmap, restrict_to_genes=bulk.index, export=False)
    #
    # This is not an optimisation and it is not about the 657 gene space. In
    # `build_from_h5ad` the gene mask is applied to X BEFORE `raw_totals` is taken, so the
    # restriction decides two things that reach every method:
    #
    #   * `cell_size`, which is the per-type mean of raw_totals. Summed over 27,625 atlas
    #     genes it is a different vector than summed over the 16,758 the bulk also carries,
    #     and the RATIOS between cell types differ, so the central mRNA-to-cell conversion
    #     differs.
    #   * the normalisation base. Each cell is scaled to 1e6 across whatever genes survive
    #     the mask, so the reference PROFILE differs too.
    #
    # Without the restriction a re-measurement can use the leaderboard's exact 657 genes and
    # still not be input-equivalent to it. Found 2026-09-14 while checking the seven
    # equivalence conditions; two runs were killed and restarted because of it.
    ref, cells, cmeta = build_from_h5ad(
        config.REFERENCE_DIR / "gbmap_core.h5ad",
        restrict_to_genes=expr.index, export=False)

    # THE GENE SPACE MUST BE THE LEADERBOARD'S, OR THE NUMBER IS NOT COMPARABLE.
    #
    # This used to call `frozen_gene_space`, which returns the FALLBACK space the pipeline
    # uses only when it has no cell-level reference — 1,591 genes. The leaderboard's 16
    # methods all ran on the benchmark's 657-gene marker subset. Both re-measurements made
    # before 2026-09-12 therefore ran on 1,591 genes and were reported "alongside" a
    # leaderboard built on 657, with the DWLS report printing a direct comparison to the
    # reimplementation's 0.7385 as though only the implementation differed. It did not.
    #
    # The run now persists its exact gene list (`benchmark/signature_genes.json`). If that
    # file is absent this script ABORTS: a re-measurement on a gene space nobody can
    # reproduce is worse than no re-measurement, because it looks like evidence.
    sel_rec = json.loads(
        (config.BENCH_DIR / "method_selection_decision.json").read_text())
    gene_rec = json.loads(genes_path.read_text())
    genes = [g for g in gene_rec["genes"] if g in set(expr.index)]
    if len(genes) != len(gene_rec["genes"]):
        print(f"BLOCKED: {len(gene_rec['genes']) - len(genes)} of the recorded "
              f"{len(gene_rec['genes'])} genes are absent from the loaded bulk. The "
              f"recorded gene space does not match this bulk matrix, so the comparison "
              f"would not be like-for-like. Nothing is imputed and nothing is dropped "
              f"silently.")
        return 2
    prov = {"source": str(genes_path.relative_to(config.PROJECT_ROOT)),
            "n_genes": len(genes),
            "sha256": config.sha256_strings(genes),
            "matches_recorded_hash": config.sha256_strings(genes) == gene_rec.get("sha256"),
            "strategy": "the leaderboard's own gene space, read from the run artefact"}
    if not prov["matches_recorded_hash"]:
        print("BLOCKED: the gene list reproduces a different hash than the artefact "
              "records. Do not report a comparison built on it.")
        return 2
    print(f"  gene space: {len(genes)} genes, hash {prov['sha256'][:16]}… "
          f"(the leaderboard's own)")

    # THE REFERENCE MUST BE BUILT FROM THE 88 TRAINING DONORS ONLY.
    #
    # The leaderboard's anatomic stage received `bench["reference"]`, built from training
    # donors alone, and a cell source of `sc_expression.loc[genes, train_cells]`. This script
    # used to hand the method the reference and the cells for ALL 110 donors, which is
    # non-equivalent twice over: the reference profile is built from different cells, and the
    # method sees the 22 HELD-OUT donors it was never supposed to see. That is the donor
    # leakage `run_benchmark`'s own guard exists to prevent, reintroduced one level up.
    #
    # Found 2026-09-14 while checking the seven equivalence conditions. Two in-flight runs
    # were killed because of it.
    from ivygap.bench.pseudobulk import split_donors
    train_donors, test_donors = split_donors(cmeta, seed=config.RANDOM_SEED)
    rec_train = sorted(map(str, sel_rec["train_donors"]))
    rec_test = sorted(map(str, sel_rec["test_donors"]))
    if sorted(map(str, train_donors)) != rec_train or sorted(map(str, test_donors)) != rec_test:
        print(f"BLOCKED: the donor split does not reproduce the recorded one "
              f"({len(train_donors)}/{len(rec_train)} train, "
              f"{len(test_donors)}/{len(rec_test)} test). The atlas did not load to the "
              f"same cells, so no comparison against the leaderboard is valid.")
        return 2
    print(f"  donor split reproduces: {len(train_donors)} train / {len(test_donors)} "
          f"held out, by name")

    train_cells = cmeta.index[cmeta["donor"].astype(str).isin(train_donors)]
    from ivygap.data.reference import build_reference
    ref_train = build_reference(cells[train_cells], cmeta.loc[train_cells],
                                name=config.PRIMARY_REFERENCE)

    bulk = expr.loc[genes, anat]
    bulk = bulk / bulk.sum(axis=0) * 1e6
    data = DeconvolutionInput(bulk=bulk, references=(ref_train.subset_genes(genes),),
                              manifest=meta.loc[anat])
    # Training cells only, matching run_all stage 4 exactly.
    r_bridge.clear_cell_source(ref_train.name)
    r_bridge.set_cell_source(ref_train.name, cells.loc[genes, train_cells],
                             cmeta.loc[train_cells])
    print(f"  {bulk.shape[0]} genes x {bulk.shape[1]} samples | "
          f"reference {len(train_cells):,} TRAINING cells from {len(train_donors)} donors")

    # --- the seven equivalence conditions, recorded rather than assumed --------
    canon_path = config.ESTIMATES_DIR / f"ivygap_{M}.csv"
    canon = pd.read_csv(canon_path, index_col=0) if canon_path.exists() else None
    equivalence = {
        "1_gene_space": {"n_genes": len(genes), "sha256": prov["sha256"],
                         "matches_leaderboard": prov["matches_recorded_hash"]},
        "2_train_donors": {"n": len(train_donors),
                           "matches_recorded": sorted(map(str, train_donors)) == rec_train},
        "3_held_out_donors": {"n": len(test_donors),
                              "matches_recorded": sorted(map(str, test_donors)) == rec_test,
                              "reference_excludes_held_out_donors": True,
                              "note": "the reference and the cell export carry TRAINING "
                                      "donors only, as run_all stage 4 does, so none of "
                                      "the 22 held-out donors' cells reach the method"},
        "4_no_silent_sample_loss": {
            "n_anatomic_samples_selected": len(anat),
            "n_samples_in_bulk": int(bulk.shape[1]),
            "equal": len(anat) == int(bulk.shape[1]),
            "canonical_n_samples": None if canon is None else int(len(canon)),
        },
        "5_sample_ids": {
            "canonical_artefact": str(canon_path.relative_to(config.PROJECT_ROOT)),
            # Compared as STRINGS on purpose. Ivy GAP sample ids are numeric, so read_csv
            # returns them as int64 while the bulk carries them as object -- comparing the
            # lists directly reports a mismatch on dtype when the content is identical, and
            # the first version of this check did exactly that and blocked a valid run.
            "identical_and_in_order": (
                None if canon is None
                else [str(x) for x in canon.index] == [str(x) for x in bulk.columns]),
            "n_canonical": None if canon is None else int(len(canon)),
            "same_membership": (None if canon is None
                                else set(map(str, canon.index)) == set(map(str, bulk.columns))),
        },
        "6_cell_type_order": {
            "order": list(data.cell_types),
            "matches_config": list(data.cell_types) == list(config.CELL_TYPES),
            "matches_canonical": (None if canon is None
                                  else list(canon.columns) == list(data.cell_types)),
        },
        "7_normalisation_and_scale": {
            "atlas_restricted_to_bulk_genes": True,
            "n_atlas_genes_after_restriction": int(cells.shape[0]),
            "per_cell_normalisation": "1e6 across the restricted gene set, before profile "
                                      "averaging; raw totals captured first as cell_size",
            "bulk_scaling": "CPM over the 657-gene space (columns sum to 1e6)",
            "cell_size_sha256": config.sha256_strings(
                [f"{c}:{v:.6f}" for c, v in
                 ref_train.cell_size.reindex(list(data.cell_types)).items()]),
            # The VALUES, not just a hash. A hash proves two runs agree; it cannot tell a
            # reader whether the conversion these factors drive does anything at all, and
            # that turned out to be the question.
            "cell_size": {c: round(float(v), 4) for c, v in
                          ref_train.cell_size.reindex(list(data.cell_types)).items()},
            "cell_size_spread_max_over_min": round(
                float(ref_train.cell_size.max() / ref_train.cell_size.min()), 4),
            "apply_cell_size_correction": bool(data.apply_cell_size_correction),
        },
    }
    hard = [equivalence["1_gene_space"]["matches_leaderboard"],
            equivalence["2_train_donors"]["matches_recorded"],
            equivalence["3_held_out_donors"]["matches_recorded"],
            equivalence["3_held_out_donors"]["reference_excludes_held_out_donors"],
            equivalence["4_no_silent_sample_loss"]["equal"],
            equivalence["6_cell_type_order"]["matches_config"]]
    if equivalence["5_sample_ids"]["identical_and_in_order"] is not None:
        hard.append(equivalence["5_sample_ids"]["identical_and_in_order"])
    if not all(hard):
        print("BLOCKED: an equivalence condition failed.")
        print(json.dumps(equivalence, indent=2))
        return 2
    print(f"  equivalence: all {len(hard)} hard conditions PASS")

    print(f"\nrunning the genuine {M} R package with a {args.budget}s budget "
          f"(the pipeline's is {r_bridge.timeout_for(M)}s)...")
    t0 = time.perf_counter()
    failed = None
    est = None
    try:
        est = r_bridge.run_r_method(M, data, timeout=args.budget)
    except Exception as exc:                                       # noqa: BLE001
        failed = f"{type(exc).__name__}: {str(exc)[:600]}"
    elapsed = time.perf_counter() - t0
    print(f"  finished in {elapsed:.0f}s ({elapsed/3600:.2f} h)")

    report = {
        "method": M,
        "what_this_is": (f"A single-method re-measurement of the GENUINE {M} R package on "
                         "the anatomic cohort, after it fell back to the Python "
                         "reimplementation in the confirmatory run. Reported alongside "
                         "that run, never substituted into it."),
        "budget_seconds": args.budget,
        "pipeline_budget_seconds": r_bridge.timeout_for("dwls"),
        "elapsed_seconds": round(elapsed, 1),
        "n_samples": int(bulk.shape[1]), "n_genes": int(bulk.shape[0]),
        "gene_space": prov,
        "input_equivalence": equivalence,
        "central_cell_size_conversion_applied": M not in r_bridge.R_RETURNS_CELL_FRACTIONS,
        "failed": failed,
    }

    if est is not None:
        # Cell-size correction and simplex projection, exactly as fit_predict applies
        # them, so this number is on the same footing as the leaderboard's.
        from ivygap.deconv.registry import build_methods
        dwls = next(m for m in build_methods(prefer_r=False) if m.name == M)
        import numpy as np
        # Use the SHARED post-processing, not a copy of it.
        #
        # This block used to reimplement `fit_predict`'s three steps and applied
        # `to_cell_fractions` unconditionally. When the conversion was made conditional for
        # Bisque -- the one package measured to convert cell size itself -- the fix landed in
        # `fit_predict` and this copy kept double-correcting, producing an estimate table
        # BIT-IDENTICAL to the uncorrected one. The fix looked applied and was not, and it
        # took comparing the two tables to notice.
        from ivygap.deconv.base import finalize_estimates
        returns_cells = M in r_bridge.R_RETURNS_CELL_FRACTIONS
        if returns_cells:
            print(f"  {M} is measured to return CELL fractions already "
                  f"(scripts/verify_cell_size_semantics.py), so the central cell-size "
                  f"conversion is SKIPPED -- applying it would be the second one")
        frame = finalize_estimates(est.to_numpy(), data, covered=None,
                                   returns_cell_fractions=returns_cells)
        out_csv = Path(f"results/estimates/ivygap_{M}_genuine.csv")
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out_csv)

        res = acs_score(frame, meta.loc[anat], method=f"{M}_genuine",
                        n_permutations=10000, n_boot=2000)
        s = res.summary()
        report["acs"] = round(float(s["acs"]), 4)
        report["ci"] = [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)]
        report["null_p"] = float(s["null_p"])
        report["n_pairs"] = int(s["n_constraint_tumor_pairs"])
        report["n_tumors"] = int(s["n_tumors"])
        report["estimates"] = str(out_csv)
        print(f"\n  GENUINE {M}: ACS {s['acs']:.4f} "
              f"[{s['ci_low']:.4f}, {s['ci_high']:.4f}] "
              f"on {s['n_constraint_tumor_pairs']} pairs, null_p {s['null_p']:.4g}")
        print(f"  (compare against the confirmatory run's leaderboard row)")
    else:
        print(f"\n  FAILED: {failed}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
