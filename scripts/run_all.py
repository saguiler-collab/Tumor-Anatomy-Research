#!/usr/bin/env python3
"""
run_all.py — the whole pipeline, in the only order that is scientifically valid.

    1. load        Ivy GAP expression + anatomic manifest + clinical table
    2. reference   build the signature matrices from single-cell data
    3. benchmark   score every method on donor-held-out pseudobulk, SELECT one
    4. anatomic    deconvolve real tissue, score the frozen biology claims  <- headline
    5. survival    out-of-fold prognosis, with the power verdict attached
    6. release     assemble the public bundle

Stage 3 must precede stage 4: the method is chosen against known truth, and the frozen
anatomic claims then validate that choice. Choosing the method that best reproduces the
claims would make the validation circular, so `run_anatomic` refuses to start without
the selection file stage 3 writes.

Stage 5 must follow stage 3 for the same reason in a stronger form: survival may never
influence method choice. `run_survival.assert_selection_frozen` checks the recorded
criterion and aborts if it mentions an outcome.

    python scripts/run_all.py --synthetic     validate end-to-end without downloads
    python scripts/run_all.py                 real Ivy GAP data, anatomic study only
    python scripts/run_all.py --full-database also deconvolve all 270 archive samples

WHAT --full-database DOES AND DOES NOT WIDEN
--------------------------------------------
It widens the DECONVOLUTION to every sample in the archive. It does not widen ACS
scoring, which stays on the 122 H&E-selected anatomic samples and is enforced inside
`run_anatomic.run` rather than left to this script to remember. Ivy GAP's other 148
samples had their anatomic structure assigned by in-situ hybridization — that is, from
expression — so an anatomic constraint scored on them would be circular, and the study
exists to eliminate exactly that circularity.

What the wider set actually buys:

  * composition for 37 tumours instead of 10, which is what takes the prognostic cohort
    from 7 tumours to 29;
  * a real sensitivity test of the leaderboard, because several methods use cross-sample
    statistics and therefore do not return the same estimate for a sample depending on
    what it was solved alongside.

The two runs write to separate directories and cannot overwrite each other.
"""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from ivygap import config
from ivygap.anatomic import run_anatomic
from ivygap.bench import run_benchmark
from ivygap.data import clinical as clinical_mod
from ivygap.deconv import r_bridge
from ivygap.data.reference import (build_reference, build_synthetic,
                                   select_signature_genes)
from ivygap.report import build_release
from ivygap.survival import run_survival


def _banner(n: int, title: str) -> None:
    print(f"\n{'=' * 72}\nSTAGE {n} — {title}\n{'=' * 72}")


def _synthetic_ivygap(expression: pd.DataFrame, meta: pd.DataFrame,
                      n_patients: int = 24, seed: int = 7):
    """
    Build an Ivy-GAP-shaped bulk cohort from the synthetic reference.

    Deliberately mimics the real archive's awkward shape: several anatomic structures
    per tumour, UNEVEN block counts, and a few tumours missing a structure entirely.
    A fixture with a tidy balanced design would validate the pipeline against a
    situation that does not occur and would never exercise the patient-equal
    aggregation that the real data demands.
    """
    from ivygap.bench.pseudobulk import NICHE_TEMPLATES
    import numpy as np

    rng = np.random.default_rng(seed)
    types = list(config.CELL_TYPES)
    donors = sorted(meta["donor"].astype(str).unique())

    pools = {}
    for d in donors:
        for t in types:
            sel = meta.index[(meta["donor"].astype(str) == d) & (meta["cell_type"] == t)]
            if len(sel):
                pools[(d, t)] = list(sel)

    cols, rows, n = {}, {}, 0
    for p in range(n_patients):
        pid = f"IVY{p:02d}"
        donor = donors[p % len(donors)]
        # Every patient contributes CT (the comparator); the rest are present with
        # high but not certain probability, and with a variable number of blocks.
        for structure in config.PRIMARY_STRUCTURES:
            if structure != config.REFERENCE_STRUCTURE and rng.random() < 0.15:
                continue
            for _ in range(int(rng.integers(1, 4))):
                target = np.array([NICHE_TEMPLATES[structure].get(t, 0.0) for t in types])
                target = rng.dirichlet(np.clip(target, 1e-4, None) * 80.0)
                counts = rng.multinomial(int(rng.integers(400, 900)), target)

                picked = []
                for k, t in enumerate(types):
                    pool = pools.get((donor, t))
                    if pool and counts[k]:
                        picked.extend(rng.choice(pool, size=int(counts[k]), replace=True))
                if not picked:
                    continue
                mix = expression[picked].sum(axis=1)
                sid = f"IVYS{n:04d}"
                cols[sid] = mix / mix.sum() * 1e6
                rows[sid] = {"patient_id": pid, "structure": structure,
                             "is_primary_structure": True}
                n += 1

    return pd.DataFrame(cols, index=expression.index), pd.DataFrame.from_dict(rows, orient="index")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--synthetic", action="store_true",
                    help="run end-to-end on a generated reference, no downloads")
    ap.add_argument("--no-r", action="store_true",
                    help="force the all-Python implementations")
    ap.add_argument("--n-test", type=int, default=config.N_TEST_PSEUDOBULK)
    ap.add_argument("--permutations", type=int, default=10_000,
                    help="within-tumour label permutations for the ACS null "
                         "(the protocol specifies 10,000)")
    ap.add_argument("--force-lock", action="store_true",
                    help="take the results-tree lock even if another run appears to "
                         "hold it. Only when you are certain that process is gone.")
    ap.add_argument("--control-draws", type=int, default=200,
                    help="independent permutations used to calibrate the "
                         "shuffled-signature control. The leaderboard's control row is "
                         "a SINGLE draw and its sd on real data is ~0.16 against a "
                         "real-method spread of 0.29, so one draw cannot support "
                         "'the controls behave' in either direction. 0 disables.")
    ap.add_argument("--full-database", action="store_true",
                    help="ALSO deconvolve all 270 archive samples, not just the 122 "
                         "H&E anatomic ones. ACS is still scored on the anatomic subset "
                         "alone — the ISH-cluster samples had their structure assigned "
                         "from expression, so scoring them would be circular. What the "
                         "wider set buys is composition for 37 tumours instead of 10, "
                         "and with it a prognostic cohort of 29 instead of 7.")
    ap.add_argument("--assume-recorded-survival-is-death", action="store_true",
                    help="Ivy GAP publishes survival_days and NO vital-status column. "
                         "Without this flag survival is BLOCKED, which is the default "
                         "and the invariant-preserving behaviour. With it, a recorded "
                         "time is read as an observed death and blank rows are dropped "
                         "(never censored, never imputed). The assumption is stamped "
                         "into every artefact it touches. Read "
                         "clinical.audit_missingness() before using it: the blanks are "
                         "NOT missing at random.")
    args = ap.parse_args()

    if args.synthetic:
        # Before anything else touches a path. A fixture run must never be able to
        # overwrite the output of a real one.
        config.use_synthetic_paths()
    config.ensure_dirs()

    # Exclusive from here on. Two runs writing one results tree is how this project lost
    # a completed run once, and nearly did so a second time when a background run
    # reported as killed had not actually died.
    try:
        release_lock = config.acquire_run_lock(force=args.force_lock)
    except config.ResultsTreeBusy as exc:
        print(f"ABORT: {exc}")
        return 1

    prefer_r = not args.no_r

    if args.synthetic:
        print("SYNTHETIC FIXTURE RUN")
        print("Everything below validates that the pipeline runs and that its scoring")
        print("behaves correctly. None of it is a biological finding.")
        print(f"Writing to {config.RESULTS_DIR.name}/ and {config.RELEASE_DIR.name}/ — "
              f"a real run's output is untouched.\n")

    # ---- stages 1-2 --------------------------------------------------------
    _banner(1, "load data")
    if args.synthetic:
        _, sc_expression, sc_meta = build_synthetic(
            n_genes=1500, n_donors=10, cells_per_type_per_donor=30)
        bulk, manifest = _synthetic_ivygap(sc_expression, sc_meta)
        clinical = None
        anatomic_index = list(bulk.columns)   # the fixture is entirely "anatomic"
        print(f"synthetic bulk: {bulk.shape[1]} samples, "
              f"{manifest['patient_id'].nunique()} patients")
        print("blocks per patient:",
              manifest.groupby('patient_id').size().describe()[['min', '50%', 'max']].to_dict())
    else:
        from ivygap.data.load_ivygap import load
        from ivygap.data import portal_metadata
        bulk, manifest, report = load(write=True)
        print(report.to_json())

        # --- protocol step 3's gate, closed against the portal rather than ourselves --
        # The loader's counts come from columns-samples.csv, which ships inside the same
        # archive it is describing. Checking them against the portal's live LIMS export
        # is what makes the gate an external check instead of a tautology.
        recon = portal_metadata.reconcile(manifest)
        (config.RESULTS_DIR / "data_reconciliation.json").write_text(recon.to_json())
        print("\nportal reconciliation:", recon.verdict)
        if recon.available and (recon.archive_only
                                or recon.study_assignment_disagreements):
            print("ABORT: the archive parse disagrees with the portal's own table.")
            return 1

        # --- clinical ---------------------------------------------------------------
        audit = clinical_mod.audit_missingness()
        (config.RESULTS_DIR / "clinical_missingness_audit.json").write_text(
            json.dumps(audit, indent=2))
        print("\nclinical missingness:", audit.get("verdict", audit.get("reason")))

        policy = ("observed-time-is-death" if args.assume_recorded_survival_is_death
                  else "require")
        try:
            clinical = clinical_mod.load_clinical(strict=False, event_policy=policy)
        except KeyError as exc:
            clinical = None
            print(f"\nclinical: NOT USABLE — {exc}")
        else:
            print("clinical:", json.dumps(clinical_mod.describe(clinical), default=str))

        # The full archive stays in `bulk`; `anatomic_index` names the H&E-selected
        # subset. Which of the two each stage receives is decided per stage, explicitly,
        # rather than by a filter applied here that a later stage has to remember.
        anatomic_index = [s for s in bulk.columns
                          if s in manifest.index[manifest["is_anatomic_study"]
                                                 & manifest["is_primary_structure"]]]
        print(f"\narchive        : {bulk.shape[1]} samples, "
              f"{manifest['patient_id'].nunique()} tumours")
        print(f"anatomic study : {len(anatomic_index)} samples, "
              f"{manifest.loc[anatomic_index, 'patient_id'].nunique()} tumours "
              f"(the rest are ISH-selected clusters — never ACS-scored)")

        if not args.full_database:
            bulk, manifest = bulk[anatomic_index], manifest.loc[anatomic_index]

        _banner(2, "reference")
        gbmap = config.REFERENCE_DIR / "gbmap_core.h5ad"
        if gbmap.exists():
            from ivygap.data.reference import build_from_h5ad
            # Restrict to genes the bulk actually carries BEFORE any cell is
            # materialised. On Core GBmap that turns a ~30,000-column read into a
            # ~16,000-column one, and a gene absent from the bulk could not have
            # contributed to a deconvolution anyway.
            # export=False: the cell-level export is written per gene set, on demand,
            # by the R bridge. Writing it for all 16,758 genes takes ~49 minutes and
            # 1.1 GB, and no method ever reads beyond its own gene space.
            _, sc_expression, sc_meta = build_from_h5ad(
                gbmap, restrict_to_genes=bulk.index, export=False)
            r_bridge.set_cell_source(config.PRIMARY_REFERENCE,
                                     sc_expression, sc_meta)
            print(f"cell-level atlas: {sc_expression.shape[1]:,} cells, "
                  f"{sc_meta['donor'].nunique()} donors")
        else:
            # No cell-level atlas. The frozen signature still supports the ACS
            # leaderboard, which is the stage that needs real tissue. It cannot support
            # pseudobulk with known truth, so that yardstick is reported UNAVAILABLE
            # rather than approximated.
            sc_expression = sc_meta = None
            print(f"no cell-level atlas at {gbmap.name}; using the vendored frozen "
                  f"signature in reference_frozen/.")
            print("  -> the pseudobulk benchmark cannot run without cells, so method")
            print("     selection and the synthetic_mixtures yardstick are skipped.")
            print("     The ACS leaderboard and the permutation null still run on real")
            print("     tissue; the agreement test will report its yardsticks as")
            print("     UNAVAILABLE, which is a reportable state, not a failure.")

    # ---- stage 3: benchmark and SELECT -------------------------------------
    if sc_expression is not None:
        _banner(3, "benchmark — score every method on known truth, select one")
        bench = run_benchmark.run(sc_expression, sc_meta, prefer_r=prefer_r,
                                  n_test=args.n_test, n_signature_genes=80)
        selected = bench["decision"]["selected_method"]
        genes = [g for g in bench["genes"] if g in bulk.index]
        reference = bench["reference"].subset_genes(genes)

        # EQUAL FOOTING ACROSS THE R BOUNDARY.
        #
        # `bench["reference"]` is built from TRAINING donors only. The Python methods
        # read that object. The R methods do not — they read the cell-level export,
        # written from whatever cell source is registered, and run_benchmark restores
        # the caller's FULL source when it finishes.
        #
        # Left alone, stage 4 would hand MuSiC, Bisque and SCDC all 110 donors while
        # NNLS, SVR, Elastic Net and the Bayesian methods worked from an 88-donor
        # signature. That is not a like-for-like leaderboard, and the certificate cannot
        # detect it because the difference is outside the DeconvolutionInput.
        #
        # The reference the benchmark selected on is also the reference the anatomy
        # validates, which keeps the two stages coherent; the cell source is matched to
        # it here rather than left at whatever the previous stage happened to leave.
        train_cells = sc_meta.index[
            sc_meta["donor"].astype(str).isin(bench["train_donors"])]
        r_bridge.set_cell_source(config.PRIMARY_REFERENCE,
                                 sc_expression.loc[genes, train_cells],
                                 sc_meta.loc[train_cells])
        print(f"cell source for stage 4: {len(train_cells):,} cells from "
              f"{len(bench['train_donors'])} training donors — the same reference the "
              f"selection was made on, for every method")
    else:
        _banner(3, "benchmark — SKIPPED (no cell-level reference)")
        print("Method selection needs known-truth mixtures, which need cells. Supply a\n"
              "single-cell atlas to enable it. No method is frozen by this run.")
        from ivygap.data.reference import frozen_gene_space, load_frozen_reference
        frozen = load_frozen_reference()
        genes, gene_prov = frozen_gene_space(frozen, bulk.index)
        print(f"gene space: {gene_prov['n_selected']:,} of {gene_prov['n_shared']:,} "
              f"shared genes — {gene_prov['strategy']}")
        (config.BENCH_DIR / "gene_space.json").write_text(
            json.dumps(gene_prov, indent=2))

        if len(genes) < config.MIN_GENES_SHARED:
            print(f"ERROR: gene selection left only {len(genes)} genes.")
            return 1
        reference = frozen.subset_genes(genes)
        selected = None

        # `assert_selection_frozen` needs a decision on disk to inspect. Writing "no
        # method was selected, and here is why" is a truthful record of a decision that
        # was not made; leaving the file absent would look identical to a run that
        # selected a method and forgot to say so.
        config.BENCH_DIR.mkdir(parents=True, exist_ok=True)
        (config.BENCH_DIR / "method_selection_decision.json").write_text(json.dumps({
            "selected_method": None,
            "criterion": "NO SELECTION MADE — the donor-held-out pseudobulk benchmark "
                         "requires a cell-level single-cell reference, which is not "
                         "present. No method is frozen by this run.",
            "outcome_used": None,
            "outcome_blind": True,
            "acs_used_for_selection": False,
        }, indent=2))

    # ---- stage 4: ACS leaderboard + the agreement test ----------------------
    # 4a is the pre-registered analysis: the 122 H&E-selected anatomic samples, both
    # deconvolved and scored. It writes the canonical paths and is the run whose numbers
    # are the study's leaderboard.
    _banner(4, "the Anatomy Test — ACS leaderboard, controls, and the agreement test")
    print("4a — PRE-REGISTERED COHORT: the H&E anatomic study, deconvolved and scored\n")
    anat = run_anatomic.run(bulk.loc[genes, anatomic_index], manifest.loc[anatomic_index],
                            references=(reference,), prefer_r=prefer_r,
                            n_permutations=args.permutations,
                            n_control_draws=args.control_draws)

    # 4b widens only the DECONVOLUTION. Scoring is still the anatomic subset — enforced
    # inside run_anatomic, not trusted to this caller — so the wider set can change a
    # method's estimates but can never contribute a circular constraint verdict.
    anat_full = None
    if args.full_database and len(anatomic_index) < bulk.shape[1]:
        _banner(4, "the Anatomy Test — FULL DATABASE (all archive samples deconvolved)")
        print("4b — every sample in the archive is deconvolved. ACS is still scored on\n"
              "     the anatomic subset alone. Several methods use cross-sample\n"
              "     statistics, so this is a genuine sensitivity test of 4a's ranking,\n"
              "     and it is what supplies composition for the prognostic cohort.\n")
        anat_full = run_anatomic.run(bulk.loc[genes], manifest,
                                     references=(reference,), prefer_r=prefer_r,
                                     n_permutations=args.permutations,
                                     deconvolve_all=True,
                                     out_subdir="full_database",
                                     n_control_draws=args.control_draws)

        # Does widening the input change the leaderboard? Reported as a table rather
        # than asserted either way.
        cmp = pd.DataFrame({
            "acs_anatomic_only": anat["leaderboard"]["acs"],
            "acs_full_database": anat_full["leaderboard"]["acs"],
        })
        cmp["delta"] = cmp["acs_full_database"] - cmp["acs_anatomic_only"]
        cmp["rank_anatomic_only"] = cmp["acs_anatomic_only"].rank(ascending=False)
        cmp["rank_full_database"] = cmp["acs_full_database"].rank(ascending=False)
        cmp["rank_change"] = cmp["rank_full_database"] - cmp["rank_anatomic_only"]
        cmp = cmp.sort_values("acs_anatomic_only", ascending=False)
        cmp.to_csv(config.ANATOMIC_DIR / "acs_cohort_sensitivity.csv")
        rho = cmp["acs_anatomic_only"].corr(cmp["acs_full_database"], method="spearman")
        print("\n--- does deconvolving the whole archive change the ACS ranking? ---")
        print(cmp.round(4).to_string())
        print(f"\nSpearman rho between the two ACS rankings: {rho:.4f}")

    # ---- stage 5: survival ---------------------------------------------------
    # Two things are written, always in this order and never merged:
    #
    #   results/survival/                       the strict result. Ivy GAP publishes no
    #                                           vital-status column, so this is BLOCKED,
    #                                           and BLOCKED is the honest answer.
    #   results/survival/declared_event_policy/ only when the run was asked for it: the
    #                                           same analysis under a named assumption,
    #                                           stamped, quarantined, and never promoted
    #                                           into the canonical tree.
    _banner(5, "survival — out-of-fold prognosis, with power stated first")

    strict_verdict = {
        "verdict": (
            "BLOCKED: tumor_details.csv is present and joins correctly, but it "
            "publishes survival_days and NO vital-status column. Who was censored is "
            "not knowable from the release, so a C-index cannot be computed. Nothing "
            "is imputed and no prognostic claim is made."
            if config.IVYGAP_TUMOR_DETAILS_PATH.exists() else
            "BLOCKED: no clinical/survival table available. Expected "
            f"{config.IVYGAP_TUMOR_DETAILS_PATH}. Survival was not computed."),
        "clinical_file_present": config.IVYGAP_TUMOR_DETAILS_PATH.exists(),
        "vital_status_column_published": False,
        "adequately_powered": False,
        "smallest_human_action": (
            "None available from Ivy GAP. The Allen Institute does not publish vital "
            "status. Either obtain follow-up status for the 42 tumours from the "
            "authors, or read the recorded times under an explicitly declared "
            "assumption (--assume-recorded-survival-is-death), which this pipeline "
            "will do only into a separate, labelled directory."),
    }
    if clinical is None or clinical.empty:
        print(strict_verdict["verdict"])
        (config.SURVIVAL_DIR / "survival_power.json").write_text(
            json.dumps(strict_verdict, indent=2))
    else:
        # A declared policy is in force. The strict verdict is still recorded, because a
        # reader must be able to see that the default answer was BLOCKED.
        (config.SURVIVAL_DIR / "survival_power.json").write_text(
            json.dumps(strict_verdict, indent=2))

        sub = config.SURVIVAL_DIR / "declared_event_policy"
        sub.mkdir(parents=True, exist_ok=True)
        print("STRICT RESULT:", strict_verdict["verdict"])
        print(f"\nProceeding under the declared policy "
              f"{clinical.attrs.get('event_policy')!r}; everything below is written to "
              f"{sub.relative_to(config.PROJECT_ROOT)}/ and is NOT the canonical result.\n")

        run_survival.assert_selection_frozen(
            config.BENCH_DIR / "method_selection_decision.json")
        contract = run_survival.assess_power(clinical, len(config.PRIMARY_CELL_TYPES) + 1)
        print(contract.verdict)

        power_doc = json.loads(contract.to_json())
        power_doc["event_policy"] = clinical.attrs.get("event_policy")
        power_doc["event_assumption"] = clinical.attrs.get("event_assumption")
        power_doc["clinical_provenance"] = clinical_mod.describe(clinical)
        power_doc["missingness_audit"] = clinical_mod.audit_missingness()
        power_doc["read_this_first"] = (
            "This directory exists only because the run was explicitly asked for it. "
            "The canonical survival result is BLOCKED — see ../survival_power.json. "
            "Every number here rests on reading a recorded survival time as an observed "
            "death, and on dropping the tumours with no recorded time. Those dropped "
            "tumours are NOT missing at random: they are strongly enriched for MGMT "
            "methylation and younger age, so this cohort is biased toward short "
            "survival.")
        (sub / "survival_power.json").write_text(json.dumps(power_doc, indent=2,
                                                            default=str))

        # Widest composition wins: the prognostic cohort is the binding constraint here
        # by a wide margin, and the full-database run covers every tumour with RNA-seq
        # rather than only the ten in the anatomic study.
        source = anat_full if anat_full is not None else anat
        source_name = ("full_database (all archive samples deconvolved)"
                       if anat_full is not None else
                       "anatomic study only (122 samples / 10 tumours)")
        src_manifest = manifest if anat_full is not None else manifest.loc[anatomic_index]
        print(f"composition source: {source_name}")

        results = {}
        from ivygap.deconv.controls import CONTROL_NAMES
        for name, est in source["estimates"].items():
            if name in CONTROL_NAMES:
                continue                       # a control has no prognostic claim to make
            per_patient = est.copy()
            per_patient["patient_id"] = src_manifest.reindex(est.index)["patient_id"].astype(str)
            per_patient = per_patient.groupby("patient_id").mean(numeric_only=True)
            results[name] = run_survival.evaluate_method(per_patient, clinical)
        table = run_survival.rank_methods(results, contract)
        table["composition_source"] = source_name
        table["event_policy"] = clinical.attrs.get("event_policy")
        table.to_csv(sub / "survival_metrics.csv")
        print(table.round(4).to_string())
        print(f"\n{contract.verdict}")
        if not contract.adequately_powered:
            print("No method is ranked on prognosis. That is the pre-specified "
                  "behaviour on this many events, not a failure of the run.")

    # ---- stage 6: release ----------------------------------------------------
    _banner(6, "release bundle")
    manifest_out = build_release.build(synthetic=args.synthetic)
    print(f"wrote {len(manifest_out['files'])} files to {config.RELEASE_DIR}")
    if manifest_out["missing_from_this_build"]:
        print(f"not produced this run: {manifest_out['missing_from_this_build']}")

    rep = anat["report"]
    print(f"\n{'=' * 72}")
    print(f"ACS cohort (pre-registered)   : {rep['n_samples']} samples / "
          f"{rep['n_tumors']} tumours")
    if anat_full is not None:
        rf = anat_full["report"]
        print(f"full-database deconvolution   : {rf['n_samples_deconvolved']} samples / "
              f"{rf['n_tumors_deconvolved']} tumours "
              f"(ACS still scored on {rf['n_samples']}/{rf['n_tumors']})")
    print(f"frozen method (outcome-blind) : {selected or 'none — benchmark skipped'}")
    print(f"best ACS method               : {rep['best_acs_method']} "
          f"({rep['best_acs']:.3f})" if rep['best_acs'] is not None else "")
    print(f"controls                      : {rep['control_verdict'][:70]}")
    primary = rep["primary_result"].get("synthetic_mixtures", {})
    if primary:
        print(f"PRIMARY RESULT                : rho = {primary.get('rho')} — "
              f"{str(primary.get('verdict'))[:60]}")
    if args.synthetic:
        print("run type                      : SYNTHETIC FIXTURE — not a finding")
    print(f"{'=' * 72}")

    # Archive before releasing the lock. results/ is by definition what the next run
    # overwrites; the archive is the citable snapshot, read-only, with a hash for every
    # file. This project has lost a completed run twice — once to a synthetic fixture
    # sharing paths, once to the test suite writing into the real tree — and both fixes
    # prevent a repeat of their own mechanism without making a finished run durable.
    # This does.
    if not args.synthetic:
        try:
            from scripts.archive_run import archive
        except ImportError:                                  # invoked as a script
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            from archive_run import archive
        try:
            dest = archive(label=f"auto · {len(anat['leaderboard'])} methods")
            print(f"archived to {dest.relative_to(config.PROJECT_ROOT)} "
                  f"(read-only; verify with scripts/archive_run.py --verify {dest.name})")
        except SystemExit as exc:
            print(f"NOT ARCHIVED: {exc}")
        except Exception as exc:                             # noqa: BLE001
            print(f"NOT ARCHIVED: {type(exc).__name__}: {exc}")

    release_lock()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise SystemExit(2)
