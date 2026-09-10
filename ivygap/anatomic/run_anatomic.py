"""
run_anatomic.py — score every method (and both controls) with ACS, then run the
agreement test that is the study's primary result.

THE ORDER OF THE TWO QUESTIONS
------------------------------
1. **The leaderboard.** ACS per method, with a bootstrap CI over tumours and a p-value
   against the within-tumour permutation null. Plus the two negative controls, which are
   what make the leaderboard interpretable at all.
2. **The experiment.** Does the ACS ranking agree with the ranking from real ground
   truth? One number, judged against a threshold fixed before it was computed.

Question 1 alone is a tool. Question 2 is what makes this a study, and it is the one
whose answer is informative in both directions.

WHAT THIS MODULE REFUSES TO DO
------------------------------
It does not select a method. ACS is the *object under test* here — the whole point is to
find out whether it would be safe to select on, and using it to select before that
question is answered would assume the conclusion. The method the release freezes is still
chosen by the outcome-blind pseudobulk benchmark.

It also does not drop the controls from the leaderboard when they score well. If a
control scores highly, the protocol's instruction is to report it and state what a harder
constraint set would need — never to retune the constraints and re-report.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.anatomic import agreement, constraints as K
from ivygap.anatomic import control_calibration, coverage as coverage_mod
from ivygap.anatomic import registration
from ivygap.anatomic.acs import score as acs_score, tumor_structure_means
from ivygap.bench.equal_footing import (EqualFootingViolation, build_certificate,
                                        check_estimates_aligned, require_comparable)
from ivygap.deconv import r_bridge
from ivygap.deconv.base import DeconvolutionInput
from ivygap.deconv.controls import CONTROL_NAMES, build_controls
from ivygap.deconv.registry import build_methods


def _load_yardsticks() -> tuple[dict[str, dict[str, float]], dict]:
    """
    The protocol's three ground-truth yardsticks.

    Prefers this run's own benchmark when it exists (it covers every method), and
    otherwise falls back to the vendored TCGA evidence, which is real but covers only
    NNLS and SVR. Either way the provenance records which source was used and how many
    methods it reached — "ground truth exists and covers 2 of 10" is a measured
    limitation with a remedy, not the same thing as "nothing was found".
    """
    from ivygap.bench import real_yardsticks

    scores, prov = real_yardsticks.load_all()

    own = config.BENCH_DIR / "benchmark_summary.csv"
    if own.exists():
        df = pd.read_csv(own, index_col=0)
        if "mae_primary" in df.columns:
            local = df["mae_primary"].dropna().to_dict()
            if len(local) > len(scores.get("synthetic_mixtures", {})):
                scores["synthetic_mixtures"] = local
                prov["synthetic_mixtures"] = {
                    "available": True,
                    "source": "this run's own donor-held-out pseudobulk benchmark",
                    "is_real_ground_truth": True,
                    "methods_covered": sorted(local),
                }
    return scores, prov


def _atlas_label_counts() -> dict | None:
    """
    The atlas's own annotation counts, read straight from the .h5ad.

    Used to measure which populations the roster drops rather than assert it. Returns
    None when no atlas is present, in which case coverage falls back to the static
    statement in the frozen constraint file.
    """
    path = config.REFERENCE_DIR / "gbmap_core.h5ad"
    if not path.exists():
        return None
    try:
        import h5py
        import numpy as np

        with h5py.File(path, "r") as f:
            grp = f["obs"]["annotation_level_3"]
            cats = [c.decode() if isinstance(c, bytes) else str(c)
                    for c in grp["categories"][:]]
            codes = grp["codes"][:]
        counts = np.bincount(codes[codes >= 0], minlength=len(cats))
        return {cats[i]: int(counts[i]) for i in range(len(cats))}
    except Exception:                                       # noqa: BLE001
        return None


def _runtime_note(method, elapsed: float) -> str | None:
    """
    Flag a fallback that ran far past the budget its R path was held to.

    The R side is bounded by `r_bridge.timeout_for`; the Python reimplementation that
    replaces it on timeout is not bounded at all. So a method can be stopped at 2,400 s
    for being slow and then spend eight hours in its substitute, which is the opposite of
    what the budget was for. Killing the fallback would be worse — it would leave the
    method with no result at all — so this records the fact loudly instead of hiding it.
    """
    r_name = getattr(method, "r_method", method.name)
    budget = r_bridge.R_METHOD_TIMEOUTS.get(r_name)
    if budget is None or elapsed <= budget:
        return None
    return (f"ran {elapsed:.0f}s against a {budget}s R budget "
            f"({elapsed / budget:.1f}x). The R path is bounded; the Python fallback that "
            f"replaces it is not, so most of this is unbounded substitute time. Reported "
            f"rather than truncated, because a truncated method has no result at all.")


def run(bulk: pd.DataFrame, manifest: pd.DataFrame, references: tuple,
        prefer_r: bool = True, n_permutations: int = 10_000, n_boot: int = 2000,
        verbose: bool = True, deconvolve_all: bool = False,
        out_subdir: str | None = None,
        n_control_draws: int = control_calibration.DEFAULT_DRAWS,
        bulk_full: pd.DataFrame | None = None) -> dict:
    """
    Deconvolve real Ivy GAP tissue with every method plus both controls, score ACS, and
    run the agreement test.

    `manifest` must carry `patient_id` (the tumour) and `structure`. The structure label
    is used only for SCORING — no method receives it, and the equal-footing guard checks
    the manifest carries no outcome column.

    DECONVOLVED SET vs SCORED SET
    -----------------------------
    `deconvolve_all=False` (the default, and the pre-registered analysis) restricts both
    to the H&E-selected anatomic samples.

    `deconvolve_all=True` solves every sample handed in — including Ivy GAP's 148
    ISH-selected cluster samples — and still scores ACS on the anatomic subset alone.
    The two are deliberately separable because they answer different questions and
    because conflating them would break an invariant:

      * ACS may only ever be scored on the anatomic subset. The cluster samples had
        their structure assigned *from expression*, so an anatomic constraint scored on
        them is circular by construction. This is enforced below, not remembered.
      * Composition itself is perfectly well defined on all 270 samples, and the wider
        set is what makes a 29-tumour prognostic cohort possible instead of a 7-tumour
        one.

    The distinction matters for a second, less obvious reason: several methods use
    cross-sample statistics (Bisque's per-gene transformation, SCDC's weighting), so the
    same sample deconvolved inside a 270-sample cohort does not receive the same estimate
    as inside a 122-sample one. Running both and comparing is therefore a real
    sensitivity check on the leaderboard, not bookkeeping.

    `out_subdir` sends every artefact into a named subdirectory of the results tree
    instead of the top level. Two cohorts never share an output path — that is the same
    isolation rule the synthetic/real split already enforces, applied one level down,
    and it is why the pre-registered 122-sample analysis keeps the canonical filenames
    while the wider run cannot overwrite a single one of them.
    """
    config.ensure_dirs()

    out = config.ANATOMIC_DIR if out_subdir is None else config.ANATOMIC_DIR / out_subdir
    est_dir = (config.ESTIMATES_DIR if out_subdir is None
               else config.ESTIMATES_DIR / out_subdir)
    out.mkdir(parents=True, exist_ok=True)
    est_dir.mkdir(parents=True, exist_ok=True)

    primary = list(config.PRIMARY_STRUCTURES)
    scored = manifest.index[manifest["structure"].isin(primary)
                            & manifest["is_anatomic_study"].astype(bool)] \
        if "is_anatomic_study" in manifest.columns else \
        manifest.index[manifest["structure"].isin(primary)]
    scored = [s for s in bulk.columns if s in scored]
    if not scored:
        raise ValueError(
            f"no samples belong to the primary structures {primary}. "
            f"Structures present: {sorted(manifest['structure'].unique())}"
        )

    if deconvolve_all:
        keep = list(bulk.columns)
    else:
        keep = scored
    bulk, manifest = bulk[keep], manifest.loc[keep]
    if bulk_full is not None:
        bulk_full = bulk_full[keep]

    # The invariant, enforced rather than remembered: whatever was deconvolved, the ACS
    # scoring set is the H&E-selected anatomic samples and nothing else.
    if "is_anatomic_study" in manifest.columns:
        leaked = [s for s in scored if not bool(manifest.loc[s, "is_anatomic_study"])]
        if leaked:
            raise ValueError(
                f"{len(leaked)} ISH-cluster sample(s) reached the ACS scoring set "
                f"({leaked[:5]}). Their structure labels were assigned using expression, "
                f"so scoring anatomic constraints on them is circular."
            )
    scoring_manifest = manifest.loc[scored]

    if verbose:
        print(f"deconvolved: {len(keep)} samples | "
              f"{manifest['patient_id'].nunique()} tumours")
        print(f"ACS-scored : {len(scored)} samples | "
              f"{scoring_manifest['patient_id'].nunique()} tumours "
              f"(H&E anatomic study only)")
        print("per structure:", scoring_manifest["structure"].value_counts().to_dict())
        print(f"constraint file: {len(K.CONSTRAINTS)} constraints, "
              f"freeze hash {K.freeze_hash()[:16]}...\n")

    data = DeconvolutionInput(bulk=bulk, references=references, manifest=manifest,
                              bulk_full=bulk_full)

    # Controls run through the identical pipeline on the identical inputs. Anything less
    # and they would not be controls for this benchmark.
    all_methods = build_methods(prefer_r=prefer_r) + build_controls()

    estimates, implementations, failures, disclosure = {}, {}, {}, []
    for method in all_methods:
        try:
            _t0 = time.perf_counter()
            estimates[method.name] = method.fit_predict(data)
            _elapsed = time.perf_counter() - _t0
            implementations[method.name] = getattr(method, "implementation_", "python") or "python"
            # Two invariants are recorded here rather than in the benchmark stage,
            # because the benchmark does not run without a cell-level reference and
            # these disclosures must exist for every run that produces a leaderboard:
            #   * a Python reimplementation is never reported as the published package;
            #   * a method running in a degenerate mode says so under its own name.
            inner = getattr(method, "fallback", method)
            disclosure.append({
                "method": method.name,
                "implementation": implementations[method.name],
                "fallback_reason": (getattr(method, "fallback_reason_", None)
                                    or getattr(method, "r_path_disabled_reason_", None)),
                "r_path_disabled": bool(getattr(method, "r_path_disabled_reason_", None)),
                "degenerate": bool(getattr(inner, "degenerate_", False)
                                   or getattr(method, "degenerate_", False)),
                "degeneracy_reason": (getattr(inner, "degeneracy_reason_", None)
                                      or getattr(method, "degeneracy_reason_", None)),
                # Which roster types this method can estimate AT ALL, and which gene
                # space it actually received. Both are part of "what was this method
                # shown", which is half of any claim that the comparison was fair.
                "models_cell_types": (sorted(method.models_cell_types)
                                      if getattr(method, "models_cell_types", None)
                                      else None),
                "gene_space": r_bridge.LAST_GENE_SPACE.get(
                    getattr(method, "r_method", method.name),
                    {"n_genes": int(data.bulk.shape[0]),
                     "gene_space": "shared marker subset"}),
                # Wall-clock, recorded because it was previously printed to the console
                # and nowhere else. In the 2026-09-10 run DWLS took 28,387 s — 2,400 s of
                # bounded R plus ~26,000 s of UNBOUNDED Python fallback — and no artefact
                # said so. A reader comparing that run to the 3,471 s one before it had no
                # way to see the difference, and a fallback that runs for eight hours is
                # an operational fact about the result, not a detail of the console.
                "wall_clock_seconds": round(_elapsed, 1),
                "wall_clock_note": _runtime_note(method, _elapsed),
            })
            if verbose:
                tag = "  [CONTROL]" if method.name in CONTROL_NAMES else ""
                print(f"  {method.name:28s} ok{tag}")
        except Exception as exc:                          # noqa: BLE001
            failures[method.name] = f"{type(exc).__name__}: {exc}"
            if verbose:
                print(f"  {method.name:28s} FAILED — {exc}")

    if not estimates:
        raise RuntimeError(f"every method failed: {failures}")

    check_estimates_aligned(estimates)
    certificate = build_certificate(data, methods=list(estimates),
                                    implementations=implementations)
    require_comparable({n: certificate for n in estimates})

    # ---- 1. the leaderboard -------------------------------------------------
    if verbose:
        print(f"\nscoring ACS ({n_permutations:,} permutations per method)...")

    results, per_constraint, per_tumor = {}, [], []
    for name, est in estimates.items():
        # ACS is scored on the anatomic subset only, whatever was deconvolved.
        res = acs_score(est.loc[scored], scoring_manifest, method=name,
                        n_permutations=n_permutations, n_boot=n_boot)
        results[name] = res
        block = res.per_constraint.copy()
        block.insert(0, "method", name)
        per_constraint.append(block)
        tblock = res.per_tumor.copy()
        tblock.insert(0, "method", name)
        per_tumor.append(tblock)
        est.to_csv(est_dir / f"ivygap_{name}.csv")

    leaderboard = pd.DataFrame([r.summary() for r in results.values()]).set_index("method")
    leaderboard["is_control"] = [m in CONTROL_NAMES for m in leaderboard.index]

    # The invariant "never report a degenerate method under its own name" has to hold on
    # the HEADLINE artefact, not only in a sidecar JSON nobody opens next to the table.
    # Without these columns acs_leaderboard.csv shows ten independent real methods where
    # there are eight: MuSiC is NNLS on a collapsed reference, SCDC ENSEMBLE is SCDC
    # with one reference, and Bisque is running two degraded modes at once.
    disc = {d["method"]: d for d in disclosure}
    leaderboard["degenerate"] = [bool(disc.get(m, {}).get("degenerate", False))
                                 for m in leaderboard.index]
    leaderboard["implementation"] = [disc.get(m, {}).get("implementation", "python")
                                     for m in leaderboard.index]
    leaderboard["degeneracy_reason"] = [disc.get(m, {}).get("degeneracy_reason")
                                        for m in leaderboard.index]

    # COMPARABILITY. ACS is a proportion of SATISFIED constraint-tumour pairs, so two
    # methods scored on different numbers of pairs are not on one scale. A method that
    # structurally cannot estimate some roster types is scored only on the constraints
    # it can address — quanTIseq's TIL10 has no tumour, endothelial or glial population,
    # so it speaks to the two macrophage constraints and to nothing else. Ranking its
    # 2-constraint score against a 7-constraint score would be a category error, and
    # putting them in one sorted column invites exactly that reading.
    #
    # So the denominator decides: a method scored on the full pair count is comparable
    # and belongs in the leaderboard proper; anything less is reported, with its
    # denominator, as partial coverage. Reported, never dropped — "no method quietly
    # missing from a ranking" is the rule, and a method that cannot be ranked still has
    # to be visible along with the reason.
    real_mask = ~leaderboard["is_control"]
    pairs = leaderboard["n_constraint_tumor_pairs"].astype("float64")
    full_denominator = float(pairs[real_mask].max()) if real_mask.any() else float("nan")
    if not np.isfinite(full_denominator):
        # No real method produced an evaluable pair at all. That is a broken run, not a
        # leaderboard, and it must not be dressed up as one method being incomparable.
        raise RuntimeError(
            "no method was scored on a single constraint-tumour pair; the ACS scoring "
            "set is empty and there is nothing to rank"
        )
    leaderboard["n_constraint_tumor_pairs_full"] = full_denominator
    leaderboard["comparable"] = (pairs >= full_denominator) & pairs.notna()
    leaderboard["coverage_note"] = [
        "-" if c else
        f"partial coverage: scored on {int(n) if pd.notna(n) else 0} of "
        f"{int(full_denominator)} constraint-tumour pairs; not comparable to the "
        f"ranked methods and excluded from the ranking"
        for c, n in zip(leaderboard["comparable"], leaderboard["n_constraint_tumor_pairs"])
    ]

    # Methods whose ACS is identical to another's are marked as a group, so a reader
    # counting points for a rank correlation counts the right number. Computed among
    # the comparable methods only — a tie across two different denominators is not one.
    groups: dict[float, list[str]] = {}
    rankable = real_mask & leaderboard["comparable"]
    for m, v in leaderboard.loc[rankable, "acs"].items():
        groups.setdefault(round(float(v), 12), []).append(str(m))
    # "-" rather than "" for no tie: pandas writes an empty string and reads it back as
    # NaN, so a consumer doing .split(",") on the column gets an AttributeError on
    # exactly the rows that are fine. An explicit sentinel round-trips.
    leaderboard["acs_tie_group"] = [
        ",".join(groups[round(float(v), 12)])
        if (r and len(groups.get(round(float(v), 12), [])) > 1) else "-"
        for v, r in zip(leaderboard["acs"], rankable)
    ]

    # Comparable methods first, so the partial-coverage rows cannot be misread as the
    # bottom of the ranking. They are not last-placed; they are not placed.
    leaderboard = leaderboard.sort_values(["comparable", "acs"],
                                          ascending=[False, False])
    per_constraint_df = pd.concat(per_constraint, ignore_index=True)
    per_tumor_df = pd.concat(per_tumor, ignore_index=True)

    # ---- the control check, stated before anything is read as a finding -----
    control_rows = leaderboard[leaderboard["is_control"]]
    # The control check compares against the methods actually in the ranking. Letting a
    # partial-coverage method into this median would move the bar a control has to clear
    # using a score computed on a different denominator.
    real_rows = leaderboard[~leaderboard["is_control"] & leaderboard["comparable"]]
    best_control = float(control_rows["acs"].max()) if not control_rows.empty else float("nan")
    median_real = float(real_rows["acs"].median()) if not real_rows.empty else float("nan")

    # The headline verdict below compares the best control against the MEDIAN real
    # method. That is a weak test on its own: a control can tie the WORST real method,
    # and beat its own permutation null, while still sitting under the median. Both are
    # the protocol's third branch — "negative controls score highly, so the constraint
    # set is too permissive" — so both are measured and reported alongside it rather
    # than left for a reader to reconstruct. Nothing here changes a score, a control or
    # the constraint file.
    worst_real = float(real_rows["acs"].min()) if not real_rows.empty else float("nan")
    control_audit = {
        "best_control": (str(control_rows["acs"].idxmax()) if not control_rows.empty
                         else None),
        "best_control_acs": best_control,
        "worst_real_method": (str(real_rows["acs"].idxmin()) if not real_rows.empty
                              else None),
        "worst_real_acs": worst_real,
        "real_methods_tied_or_beaten_by_best_control": (
            sorted(real_rows.index[real_rows["acs"] <= best_control].tolist())
            if not control_rows.empty else []),
        "controls_beating_their_own_null": sorted(
            control_rows.index[(control_rows["null_p"] < 0.05)
                               & (control_rows["acs"] > control_rows["null_mean"])
                               ].tolist()) if not control_rows.empty else [],
    }
    if not control_rows.empty:
        bc = control_audit["best_control"]
        bc_rows = per_constraint_df[per_constraint_df["method"] == bc].copy()
        bc_rows["w_den"] = bc_rows["n_tumors_evaluable"] * bc_rows["weight"]
        perfect = bc_rows[bc_rows["n_satisfied"] == bc_rows["n_tumors_evaluable"]]
        den = float(bc_rows["w_den"].sum())
        control_audit["constraints_the_best_control_satisfies_perfectly"] = \
            sorted(perfect["constraint"].tolist())
        control_audit["weight_share_satisfied_perfectly_by_best_control"] = (
            float(perfect["w_den"].sum() / den) if den else float("nan"))

    # The control is a single permutation. Measure its sampling distribution so the
    # leaderboard row can be placed in it, and so "distinguishable from noise?" is
    # answered against a percentile rather than against one arbitrary draw.
    calibration = {}
    if n_control_draws > 0:
        if verbose:
            print(f"\ncalibrating the shuffled-signature control "
                  f"({n_control_draws} draws)...")
        calibration = control_calibration.calibrate(
            bulk[scored], scoring_manifest, references[0],
            observed={n: r.acs for n, r in results.items()},
            n_draws=n_control_draws)
        control_audit["calibration_verdict"] = calibration.get("verdict")
        control_audit["not_distinguishable_from_noise"] = calibration.get(
            "not_distinguishable_from_noise", [])
        if verbose:
            print(f"  {calibration.get('verdict')}")

    # The verdict used to compare the best control against the MEDIAN real method only.
    # That is too lenient to be worth printing: a control can tie the WORST real method
    # AND beat its own permutation null while still sitting under the median, and the
    # run would announce "CONTROLS BEHAVE". Both of those are the protocol's third
    # branch. The constraint file is not touched — only the honesty of the sentence.
    tied = control_audit["real_methods_tied_or_beaten_by_best_control"]
    beat_null = control_audit["controls_beating_their_own_null"]

    if control_rows.empty:
        control_verdict = "NO CONTROLS RAN — the leaderboard is not interpretable."
    elif best_control >= median_real:
        control_verdict = (
            f"CONSTRAINT SET TOO PERMISSIVE: the best negative control scores "
            f"{best_control:.3f}, at or above the median real method ({median_real:.3f}). "
            f"Per protocol this is reported as-is; the constraint file is NOT retuned. "
            f"A harder set would need constraints that a structured-but-meaningless "
            f"composition cannot satisfy."
        )
    elif tied or beat_null:
        bits = []
        if tied:
            bits.append(f"it ties or beats {', '.join(tied)}")
        if beat_null:
            bits.append(f"{', '.join(beat_null)} beats its own permutation null")
        control_verdict = (
            f"CONTROLS PARTLY BEHAVE: the best control ({best_control:.3f}) is below the "
            f"median real method ({median_real:.3f}), but {' and '.join(bits)}. The "
            f"constraint set separates most methods from noise and does NOT separate "
            f"all of them. Reported as-is; the constraint file is NOT retuned. Note the "
            f"control is a SINGLE permutation — see the control calibration artefact "
            f"before reading a tie as a property of the constraints."
        )
    else:
        control_verdict = (
            f"CONTROLS BEHAVE: best control {best_control:.3f} sits below the median "
            f"real method {median_real:.3f}, ties no real method, and does not beat its "
            f"own permutation null. The constraint set discriminates."
        )

    # ---- 2. the experiment ---------------------------------------------------
    # The agreement test asks whether ACS orders methods the way accuracy does, so it
    # may only see methods whose ACS is on one scale. A partial-coverage method's score
    # comes from a different denominator over a different subset of constraints; feeding
    # it into the rank correlation would move the study's headline number using a
    # quantity that is not the same quantity.
    comparable_methods = set(leaderboard.index[leaderboard["comparable"]])
    acs_scores = {n: r.acs for n, r in results.items() if n in comparable_methods}
    excluded_from_agreement = sorted(set(results) - comparable_methods)
    yardsticks, yardstick_prov = _load_yardsticks()
    agreement_table = agreement.run_all_yardsticks(
        acs_scores, yardsticks,
        # Every yardstick here is an error metric, so lower is better throughout.
        higher_is_better={k: False for k in yardsticks},
    )
    yardstick_prov["methods_excluded_for_partial_coverage"] = excluded_from_agreement
    (out / "yardstick_provenance.json").write_text(json.dumps(yardstick_prov, indent=2))

    # Where a yardstick covers too few methods to rank, a direct pairwise comparison is
    # still informative and is reported rather than discarded: on two methods a rank
    # correlation is meaningless, but "ACS orders these two the other way round" is not.
    pairwise = []
    for yname, ys in yardsticks.items():
        shared = [m for m in ys if m in acs_scores]
        if 2 <= len(shared) < agreement.MIN_METHODS_FOR_CORRELATION:
            for i, a in enumerate(shared):
                for b in shared[i + 1:]:
                    acs_says = a if acs_scores[a] > acs_scores[b] else b
                    truth_says = a if ys[a] < ys[b] else b     # lower error is better
                    pairwise.append({
                        "yardstick": yname, "method_a": a, "method_b": b,
                        "acs_prefers": acs_says, "truth_prefers": truth_says,
                        "agree": acs_says == truth_says,
                    })
    pairwise_df = pd.DataFrame(pairwise)

    # The agreement test in its legible form: each method's ACS rank beside its rank on
    # each yardstick that could rank it, with the gap. The correlation is one number and
    # this is where it becomes arguable — a large gap names a specific method the anatomy
    # likes and the truth does not.
    rank_rows = []
    for yname, ys in yardsticks.items():
        shared = [m for m in ys if m in acs_scores and m not in CONTROL_NAMES]
        if len(shared) < 2:
            continue
        block = pd.DataFrame({
            "acs": {m: acs_scores[m] for m in shared},
            "truth": {m: ys[m] for m in shared},
        })
        block["yardstick"] = yname
        block["acs_rank"] = block["acs"].rank(ascending=False)      # higher ACS is better
        block["truth_rank"] = block["truth"].rank(ascending=True)   # lower error is better
        block["rank_gap"] = (block["truth_rank"] - block["acs_rank"]).abs()
        rank_rows.append(block.reset_index(names="method"))
    ranks_df = (pd.concat(rank_rows, ignore_index=True) if rank_rows
                else pd.DataFrame(columns=["method", "yardstick", "acs", "truth",
                                           "acs_rank", "truth_rank", "rank_gap"]))

    # ---- write ---------------------------------------------------------------
    leaderboard.to_csv(out / "acs_leaderboard.csv")
    per_constraint_df.to_csv(out / "acs_per_constraint.csv", index=False)
    per_tumor_df.to_csv(out / "acs_per_tumor.csv", index=False)
    agreement_table.to_csv(out / "agreement_test.csv")
    if not pairwise_df.empty:
        pairwise_df.to_csv(out / "agreement_pairwise.csv", index=False)
    ranks_df.to_csv(out / "agreement_method_ranks.csv", index=False)
    agreement.write_report(agreement_table, out / "agreement_report.json")
    (out / "constraint_file.json").write_text(
        json.dumps(K.registration_payload(), indent=2))
    (out / "equal_footing_certificate.json").write_text(certificate.to_json())
    (out / "implementation_report.json").write_text(json.dumps(disclosure, indent=2))

    # The three risks Anatomy_Test.md names that nothing was checking.
    from ivygap.deconv import configs as configs_mod
    configs_mod.write(all_methods, out / "method_configs.json", implementations)

    reg = registration.status()
    (out / "registration_status.json").write_text(reg.to_json())

    cov = coverage_mod.reference_coverage(
        sc_meta=getattr(references[0], "cell_meta_", None),
        atlas_labels=_atlas_label_counts())
    (out / "reference_coverage.json").write_text(json.dumps(cov, indent=2, default=str))
    if calibration:
        (out / "control_calibration.json").write_text(json.dumps(calibration, indent=2))

    best_real = str(real_rows.index[0]) if not real_rows.empty else None
    if best_real:
        tumor_structure_means(estimates[best_real].loc[scored], scoring_manifest).to_csv(
            out / "composition_by_tumor_structure.csv")

    report = {
        "constraint_freeze_hash": K.freeze_hash(),
        "n_constraints": len(K.CONSTRAINTS),
        "coverage": K.coverage_report(),
        "out_subdir": out_subdir,
        "n_samples_deconvolved": len(keep),
        "n_tumors_deconvolved": int(manifest["patient_id"].nunique()),
        "n_samples": len(scored),
        "n_tumors": int(scoring_manifest["patient_id"].nunique()),
        "deconvolved_all_samples": bool(deconvolve_all),
        "acs_scoring_set": ("H&E-selected anatomic study only; ISH-cluster samples are "
                            "excluded from scoring by construction because their "
                            "structure labels were assigned using expression"),
        "structures": scoring_manifest["structure"].value_counts().to_dict(),
        "n_permutations": n_permutations,
        "best_acs_method": best_real,
        "best_acs": float(real_rows["acs"].max()) if not real_rows.empty else None,
        "control_verdict": control_verdict,
        "registration": json.loads(reg.to_json()),
        "reference_coverage": cov,
        "control_audit": control_audit,
        "control_calibration": calibration,
        "best_control_acs": best_control,
        "median_real_acs": median_real,
        "primary_result": agreement_table.to_dict(orient="index"),
        "yardstick_provenance": yardstick_prov,
        "pairwise_where_ranking_impossible": (
            pairwise_df.to_dict(orient="records") if not pairwise_df.empty else []),
        "implementation_disclosure": disclosure,
        "degenerate_methods": [d["method"] for d in disclosure if d["degenerate"]],
        "python_reimplementations": [d["method"] for d in disclosure
                                     if d["implementation"] == "python-reimplementation"],
        "failures": failures,
        "note": "ACS is the object under test here, not a selection criterion. The "
                "method frozen for release is chosen by the outcome-blind pseudobulk "
                "benchmark, never by ACS.",
    }
    (out / "anatomic_report.json").write_text(json.dumps(report, indent=2, default=str))

    if verbose:
        print("\n--- ACS leaderboard ---")
        print(leaderboard[["acs", "ci_low", "ci_high", "null_mean", "null_p",
                           "n_tumors", "is_control"]].round(3).to_string())
        print(f"\n{control_verdict}\n")
        print("--- agreement test (the primary result) ---")
        cols = [c for c in ["rho", "ci_low", "ci_high", "n_methods", "verdict"]
                if c in agreement_table.columns]
        print(agreement_table[cols].to_string())

    return {"leaderboard": leaderboard, "per_constraint": per_constraint_df,
            "per_tumor": per_tumor_df, "agreement": agreement_table,
            "results": results, "report": report, "estimates": estimates,
            "scored_index": scored, "scoring_manifest": scoring_manifest,
            "certificate": certificate}
