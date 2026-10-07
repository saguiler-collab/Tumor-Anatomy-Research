"""
verify_rerun.py -- re-run the project and show, file by file, whether every result reproduces.

The user asked (2026-10-06) to re-run the entire project to verify how correct it is. This harness
does that without trusting anything it re-runs:

  1. SNAPSHOT  copies results/ to results_snapshot_<label>/ before anything runs (once; never
               overwritten). The GIMiCC probe-drop INPUT matrices (regenerable, 280 MB) are skipped;
               every output is kept.
  2. RE-RUN    each step's command, by cost tier, logging stdout/stderr, exit code and wall time to
               results/verification/logs/<step>.log. Resumable: a finished step is skipped unless
               --force.
  3. COMPARE   every output the step declares against the snapshot copy:
                 JSON  leaf by leaf with scripts/compare_artefacts.py's own rules (timestamps,
                       durations and paths ignored; numbers within a relative 1e-9);
                 CSV   aligned on index and columns, largest numeric difference reported;
               and classifies it REPRODUCED, NUMERICAL NOISE (only numbers differ, all by <= 1e-6
               relative), REPRODUCED + NEW FIELDS (every original value reproduced; the re-run only adds
               fields the code gained after the original was written -- listed), DIFFERS (anything
               else -- listed), NEW (absent before) or NOT WRITTEN.
               Steps whose verdict is printed rather than written (independent_verification,
               extract_gimicc_cpgs --verify) are judged on their own words.
  4. REPORT    results/verification/verification_rerun.json and docs/VERIFICATION_RERUN.md.

Tiers (run in this order; the user's machine has 8 GB and 2 cores, so heavy tiers run detached):
  fast    analyses that read saved estimates and fits             minutes each
  medium  single-package refits and atlas streams                 tens of minutes
  heavy   full deconvolution re-runs (registered pipeline, TCGA)  hours each
  days    BayesPrism in its authors' configuration on TCGA        days; run only on request

    python3 scripts/verify_rerun.py --snapshot 20261006
    python3 scripts/verify_rerun.py --tier fast
    python3 scripts/verify_rerun.py --report
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
from compare_artefacts import DEFAULT_IGNORE, compare, flatten  # noqa: E402

PY = sys.executable
RES = config.RESULTS_DIR
VDIR = RES / "verification"
STATE = VDIR / "state.json"
REPORT_JSON = VDIR / "verification_rerun.json"
REPORT_MD = config.PROJECT_ROOT / "docs" / "VERIFICATION_RERUN.md"
LABEL_FILE = VDIR / "snapshot_label.txt"


def S(sid, tier, cmd, outputs=(), verdict=None, note=""):
    """One step: id, tier, command (list; 'PY' is this interpreter), outputs (relative to results/)."""
    return {"id": sid, "tier": tier, "cmd": [PY if c == "PY" else c for c in cmd], "outputs": list(outputs),
            "verdict": verdict, "note": note}


def sc(name):
    return str(ROOT / "scripts" / name)


STEPS = [
    # ------------------------------------------------------------------ fast: from saved estimates
    S("independent_verification", "fast", ["PY", sc("independent_verification.py")], verdict="no_disagrees",
      note="recomputes headline numbers from raw files with NO project code"),
    S("lymphoid_ordering_gbm", "fast", ["PY", sc("lymphoid_ordering.py"), "--cohort", "gbm", "--reference", "frozen"],
      ["lymphoid_ordering.json"]),
    S("lymphoid_ordering_lgg", "fast", ["PY", sc("lymphoid_ordering.py"), "--cohort", "lgg", "--reference", "frozen"],
      ["lymphoid_ordering_lgg.json"]),
    S("lymphoid_ordering_gbm_h5ad", "fast", ["PY", sc("lymphoid_ordering.py"), "--cohort", "gbm", "--reference", "h5ad"],
      ["lymphoid_ordering_h5ad.json"]),
    S("lymphoid_ordering_lgg_h5ad", "fast", ["PY", sc("lymphoid_ordering.py"), "--cohort", "lgg", "--reference", "h5ad"],
      ["lymphoid_ordering_lgg_h5ad.json"]),
    S("lymphoid_tracking", "fast", ["PY", sc("lymphoid_tracking.py")], ["lymphoid_tracking.json"]),
    S("immune_arm_gbm", "fast", ["PY", sc("immune_arm.py"), "--cohort", "gbm"], ["immune_arm.json"]),
    S("immune_arm_lgg", "fast", ["PY", sc("immune_arm.py"), "--cohort", "lgg"], ["immune_arm_lgg.json"]),
    S("failure_factors_gbm", "fast", ["PY", sc("failure_factors.py"), "--cohort", "gbm"], ["failure_factors_gbm.json"]),
    S("failure_factors_lgg", "fast", ["PY", sc("failure_factors.py"), "--cohort", "lgg"], ["failure_factors_lgg.json"]),
    S("cross_cohort_validation", "fast", ["PY", sc("validate_across_cohorts.py")], ["cross_cohort_validation.json"]),
    S("yardstick_agreement", "fast", ["PY", sc("yardstick_agreement.py"), "--reference", "frozen"],
      ["yardstick_agreement.json"], note="the registered orthogonal agreement statistic"),
    S("yardstick_agreement_h5ad", "fast", ["PY", sc("yardstick_agreement.py"), "--reference", "h5ad"],
      ["yardstick_agreement_h5ad.json"]),
    S("agreement_power", "fast", ["PY", sc("agreement_power.py")], ["agreement_power.json"]),
    S("anatomic_auc", "fast", ["PY", sc("anatomic_auc.py")], ["anatomic_auc.json"]),
    S("anatomy_vs_biology", "fast", ["PY", sc("anatomy_vs_biology.py"), "--write"], ["anatomy_vs_biology.json"]),
    S("method_composition", "fast", ["PY", sc("method_composition.py"), "--write"], ["method_composition.json"]),
    S("model_fit_residual", "fast", ["PY", sc("model_fit_residual.py")], ["model_fit_residual.json"]),
    S("identifiability_robustness", "fast", ["PY", sc("identifiability_robustness.py")], ["identifiability_robustness.json"]),
    S("identifiability_denominators", "fast", ["PY", sc("identifiability_denominators.py")],
      ["identifiability_denominators.json"]),
    S("agreement_selection_test", "fast", ["PY", sc("agreement_selection_test.py")], ["agreement_selection_test.json"]),
    S("extension_agreement", "fast", ["PY", sc("extension_agreement.py")], ["extension/agreement_extended.json"]),
    S("bayesprism_verdict", "fast", ["PY", sc("bayesprism_verdict.py")], ["extension/bayesprism_verdict.json"]),
    S("imc_anchored_test", "fast", ["PY", sc("imc_anchored_test.py")], ["imc_anchored_test.json"]),
    S("imc_anchored_robustness", "fast", ["PY", sc("imc_anchored_robustness.py")], ["imc_anchored_robustness.json"]),
    S("darmanis_constraint_check", "fast", ["PY", sc("darmanis_constraint_check.py")], ["darmanis_constraint_check.json"]),
    S("albiach_constraint_check", "fast", ["PY", sc("albiach_constraint_check.py")], ["albiach_constraint_check.json"]),
    S("ish_constraint_check", "fast", ["PY", sc("ish_constraint_check.py")], ["ish_constraint_check.json"]),
    S("benchmark_concordance", "fast", ["PY", sc("benchmark_concordance.py")], ["benchmark_concordance.json"]),
    S("gimicc_analyse", "fast", ["PY", sc("gimicc_truth.py"), "--stage", "analyse"], ["gimicc_truth_confirmation.json"]),
    S("gimicc_secondary", "fast", ["PY", sc("gimicc_secondary.py")], ["gimicc_secondary.json"]),
    S("gimicc_diagnostics", "fast", ["PY", sc("gimicc_diagnostics.py")], ["gimicc_diagnostics.json"]),
    S("gimicc_inputs_verify", "fast", ["PY", sc("extract_gimicc_cpgs.py"), "--verify"], verdict="exit0",
      note="GIMiCC's input matrices rebuilt from the raw 450K files and compared by sha256"),
    S("absolute_join_audit", "fast", ["PY", sc("absolute_join_audit.py")], ["absolute_join_audit.json"]),
    S("klemm_figure1f", "fast", ["PY", sc("klemm_figure1f.py")], ["klemm_t_vs_b.json"]),
    S("gbmap_t_vs_b", "fast", ["PY", sc("gbmap_t_vs_b.py")], ["gbmap_t_vs_b.json"]),
    S("spillover_lgg", "fast", ["PY", sc("spillover_test.py")], ["spillover_lgg.json"],
      note="recovered 2026-10-06 from a one-off command (provenance gap)"),
    S("zero_lymphoid_vs_purity", "fast", ["PY", sc("zero_lymphoid_vs_purity.py")], ["zero_lymphoid_vs_purity.json"],
      note="recovered 2026-10-06 from a one-off command (provenance gap)"),
    S("summarize_results_check", "fast", ["PY", sc("summarize_results.py"), "--check"], verdict="exit0",
      note="RESULTS.md agrees with the artefacts"),
    S("cptac_per_sample_analyse", "fast", ["PY", sc("cptac_per_sample.py"), "--stage", "analyse"],
      ["cptac_per_sample_truth.json"], note="added 2026-10-06 after the snapshot: NEW on the first pass"),
    S("identifiability_e3", "fast", ["PY", sc("identifiability_e3.py")], ["identifiability_e3.json"],
      note="added 2026-10-07 after the snapshot (NEW); recomputes E3's statistics from its saved subset fits"),
    S("evaluation_matrix", "fast", ["PY", sc("evaluation_matrix.py")], ["evaluation_matrix.json"],
      note="added 2026-10-07 after the snapshot: NEW on the first pass; exits 1 if a join drifts"),
    S("cptac_wgs_purity_analyse", "fast", ["PY", sc("cptac_wgs_purity.py"), "--stage", "analyse"],
      ["cptac_wgs_purity.json", "cptac/wgs_purity_per_sample.csv"],
      note="added 2026-10-06 after the snapshot: NEW on the first pass"),
    # ------------------------------------------------------------------ medium: refits and streams
    S("atlas_t_vs_b", "medium", ["PY", sc("atlas_t_vs_b.py")], ["atlas_t_vs_b.json"], note="streams the 2 GB matrix"),
    S("identifiability_probe", "medium", ["PY", sc("identifiability_probe.py")], ["identifiability_probe.json"]),
    S("equal_footing_ranking", "medium", ["PY", sc("equal_footing_ranking.py")], ["equal_footing_ranking.json"]),
    S("equal_footing_ranking_lgg", "medium", ["PY", sc("equal_footing_ranking.py"), "lgg"], ["equal_footing_ranking_lgg.json"]),
    S("methylation_celltypes_gbm", "medium", ["PY", sc("methylation_celltypes.py"), "--cohort", "gbm"],
      ["methylation_celltypes.json", "methylation_celltypes.csv"], note="EpiDISH re-run: the registered lymphoid truth itself"),
    S("methylation_celltypes_lgg", "medium", ["PY", sc("methylation_celltypes.py"), "--cohort", "lgg"],
      ["methylation_celltypes_lgg.json", "methylation_celltypes_lgg.csv"]),
    S("bisque_anchoring_synthetic", "medium", ["Rscript", sc("bisque_anchoring_synthetic.R")],
      ["bisque_anchoring_synthetic.json"], note="genuine BisqueRNA on planted truth"),
    S("nk_reannotation", "medium", ["PY", sc("nk_reannotation.py")], ["nk_reannotation.json"]),
    S("ribosomal_test", "medium", ["PY", sc("ribosomal_test.py")], ["ribosomal_test.json"]),
    S("b_column_diagnostics", "medium", ["PY", sc("b_column_diagnostics.py")], ["b_column_diagnostics.json"]),
    S("b_profile_tissue_likeness", "medium", ["PY", sc("b_profile_tissue_likeness.py")], ["b_profile_tissue_likeness.json"]),
    S("atlas_cell_fractions", "medium", ["PY", sc("atlas_cell_fractions.py")], ["atlas_cell_fractions.json"]),
    S("independent_atlas_test", "medium", ["PY", sc("independent_atlas_test.py")], ["independent_atlas_test.json"]),
    S("independent_atlas_test_ig", "medium", ["PY", sc("independent_atlas_test.py"), "--drop-ig"],
      ["independent_atlas_test_ig_removed.json"]),
    S("independent_atlas_seeds", "medium", ["PY", sc("independent_atlas_seeds.py")],
      ["independent_atlas_seed_sensitivity.json"]),
    S("bisque_anchoring", "medium", ["PY", sc("bisque_anchoring.py")], ["bisque_anchoring.json"]),
    S("unmix_loss_ablation_gbm", "medium", ["PY", sc("unmix_loss_ablation.py"), "--cohort", "gbm"],
      ["unmix_loss_ablation_gbm.json"]),
    S("identifiability_diagnostics", "medium", ["PY", sc("identifiability_diagnostics.py")],
      ["identifiability_diagnostics.json"], note="re-fits the pre-declared unmix setting; reads the other six"),
    S("linseed_anatomic", "medium", ["PY", sc("linseed_run.py"), "anatomic"], ["linseed_anatomic.json"]),
    S("linseed_tcga_gbm", "medium", ["PY", sc("linseed_run.py"), "tcga", "--cohort", "gbm"], ["linseed_tcga_gbm.json"]),
    S("linseed_tcga_lgg", "medium", ["PY", sc("linseed_run.py"), "tcga", "--cohort", "lgg"], ["linseed_tcga_lgg.json"]),
    # ------------------------------------------------------------------ heavy: full re-fits
    S("registered_pipeline", "heavy", ["PY", sc("run_all.py"), "--matrix", "raw/X"],
      ["anatomic/anatomic_report.json", "anatomic/implementation_report.json", "anatomic/acs_leaderboard.csv",
       "clinical_missingness_audit.json"],
      note="the registered pipeline end to end (reference, benchmark, anatomic ACS, controls)"),
    S("tcga_gbm_frozen", "heavy", ["PY", sc("absolute_purity_yardstick.py"), "--cohort", "gbm", "--reference", "frozen"],
      ["absolute_purity_yardstick.json", "estimates_full.csv"]),
    S("tcga_lgg_frozen", "heavy", ["PY", sc("absolute_purity_yardstick.py"), "--cohort", "lgg", "--reference", "frozen"],
      ["absolute_purity_yardstick_lgg.json", "estimates_full_lgg.csv"]),
    S("tcga_gbm_h5ad", "heavy", ["PY", sc("absolute_purity_yardstick.py"), "--cohort", "gbm", "--reference", "h5ad"],
      ["absolute_purity_yardstick_h5ad.json", "estimates_full_h5ad.csv"]),
    S("tcga_lgg_h5ad", "heavy", ["PY", sc("absolute_purity_yardstick.py"), "--cohort", "lgg", "--reference", "h5ad"],
      ["absolute_purity_yardstick_lgg_h5ad.json", "estimates_full_lgg_h5ad.csv"]),
    S("extension_tcga", "heavy", ["PY", sc("extension_tcga.py")],
      ["extension/tcga_gbm_frozen.json", "extension/tcga_lgg_frozen.json", "extension/tcga_gbm_h5ad.json",
       "extension/tcga_lgg_h5ad.json"]),
    S("cdseq_anatomic", "heavy", ["PY", sc("cdseq_anatomic.py")], ["cdseq_anatomic.json"]),
    S("reference_sensitivity_neftel", "heavy", ["PY", sc("reference_sensitivity.py"), "--baseline", "gbmap_linear",
      "--reference", "neftel"], ["reference_sensitivity_gbmap_linear_vs_neftel.json"]),
    S("reference_sensitivity_darmanis", "heavy", ["PY", sc("reference_sensitivity.py"), "--baseline", "gbmap_linear",
      "--reference", "darmanis"], ["reference_sensitivity_gbmap_linear_vs_darmanis.json"]),
    S("unmix_s2_anatomy", "heavy", ["PY", sc("remeasure_method.py"), "--method", "deseq2_unmix", "--e2-s2",
      "--out", "results/verification/rerun_outputs/ivygap_unmix_ablation.json"],
      ["verification/rerun_outputs/ivygap_unmix_ablation.json=>identifiability/ivygap_unmix_ablation.json"],
      note="written beside the original, which its producer refuses to overwrite"),
    S("recide_anatomic", "heavy", ["PY", sc("run_recide.py"),
      "--out", "results/verification/rerun_outputs/recide_anatomic.json"],
      ["verification/rerun_outputs/recide_anatomic.json=>extension/recide_anatomic.json"],
      note="ReCIDE on Ivy GAP (no time budget, as its default now is)"),
]
BY_ID = {s["id"]: s for s in STEPS}

#: Heavy re-fits of supplementary analyses, not re-run in the 2026-10-07 verification by the user's decision ("core steps
#: only"): their downstream analyses reproduced from the saved fits in the fast and medium tiers, and re-running them
#: mostly repeats the BayesPrism/DWLS 2,400 s timeouts (OPEN_DEFECTS D35). Run them with --include-supplementary.
SUPPLEMENTARY = {
    "cdseq_anatomic": "reference-free arm (supports R1)",
    "reference_sensitivity_neftel": "reference sensitivity (D14)",
    "reference_sensitivity_darmanis": "reference sensitivity (D14)",
    "unmix_s2_anatomy": "E2's secondary S2 arm",
    "recide_anatomic": "extension panel, anatomy arm",
}
TIERS = ["fast", "medium", "heavy", "days"]


# ------------------------------------------------------------------------------------ comparison
def compare_json(old: Path, new: Path) -> dict:
    a, b = json.loads(old.read_text()), json.loads(new.read_text())
    diffs = compare(a, b, 1e-9, DEFAULT_IGNORE)
    n = sum(1 for _ in flatten(a))
    if not diffs:
        return {"status": "REPRODUCED", "leaves": n}
    if all(d.startswith("extra in re-run:") for d in diffs):
        return {"status": "REPRODUCED + NEW FIELDS", "leaves": n, "n_differences": len(diffs), "examples": diffs}
    loose = compare(a, b, 1e-6, DEFAULT_IGNORE)
    status = "NUMERICAL NOISE" if not loose else "DIFFERS"
    return {"status": status, "leaves": n, "n_differences": len(diffs), "examples": diffs[:12]}


def _recompare_moved(sid: str, out: str, o: dict) -> dict:
    """A DIFFERS whose re-run copy was moved aside (keep_registered) is compared again with the current rules,
    so a rule added later (REPRODUCED + NEW ROWS) applies to records written before it existed."""
    if o.get("status") != "DIFFERS" or not out.endswith(".csv"):
        return o
    moved = VDIR / "rerun_outputs" / sid / out
    try:
        snap = snapshot_dir() / out
    except SystemExit:
        return o
    if not (moved.exists() and snap.exists()):
        return o
    r = compare_csv(snap, moved)
    return r if r["status"] != "DIFFERS" else o


def _reclassify(o: dict) -> dict:
    """A DIFFERS recorded before REPRODUCED + NEW FIELDS existed, judged by the same rule from its record:
    only when every difference was listed and every one is a field the re-run adds."""
    ex = o.get("examples", [])
    if (o.get("status") == "DIFFERS" and ex and o.get("n_differences") == len(ex)
            and all(e.startswith("extra in re-run:") for e in ex)):
        return {**o, "status": "REPRODUCED + NEW FIELDS"}
    return o


def _keyed_rows(a: pd.DataFrame, b: pd.DataFrame) -> dict | None:
    """When a re-run CSV ADDS rows (a method that now runs) or renames a key column, compare the registered rows
    by key: the non-numeric columns, aligned by position. Returns a verdict only if every registered row is present
    and identical (relative 1e-9); otherwise None, and the caller reports DIFFERS."""
    na, nb = a.select_dtypes("number").columns, b.select_dtypes("number").columns
    ka, kb = [c for c in a.columns if c not in na], [c for c in b.columns if c not in nb]
    if list(na) != list(nb) or not ka or len(ka) != len(kb):
        return None
    renamed = {y: x for x, y in zip(ka, kb) if x != y}
    b = b.rename(columns=renamed)
    if a.duplicated(ka).any() or b.duplicated(ka).any():
        return None
    m = a.merge(b, on=ka, how="left", suffixes=("_a", "_b"), indicator=True)
    if (m["_merge"] != "both").any():
        return None
    x = m[[f"{c}_a" for c in na]].to_numpy(float)
    y = m[[f"{c}_b" for c in na]].to_numpy(float)
    if (np.isnan(x) != np.isnan(y)).any():
        return None
    d = np.where(np.isnan(x), 0.0, np.abs(x - y)) / np.maximum(1.0, np.abs(np.nan_to_num(x)))
    if d.size and float(d.max()) > 1e-9:
        return None
    extra = b.merge(a[ka], on=ka, how="left", indicator=True)
    extra = extra[extra["_merge"] == "left_only"]
    return {"status": "REPRODUCED + NEW ROWS", "registered_rows": int(len(a)), "new_rows": int(len(extra)),
            "new_row_keys": sorted({str(v) for c in ka for v in extra[c].unique()} - {str(v) for c in ka for v in a[c].unique()})[:20],
            "renamed_key_columns": {v: k for k, v in renamed.items()}}


def compare_csv(old: Path, new: Path) -> dict:
    a, b = pd.read_csv(old), pd.read_csv(new)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        keyed = _keyed_rows(a, b)
        return keyed or {"status": "DIFFERS", "why": f"shape/columns differ: {a.shape} vs {b.shape}"}
    num = a.select_dtypes("number").columns
    txt = [c for c in a.columns if c not in num]
    if txt and not a[txt].astype(str).equals(b[txt].astype(str)):
        return {"status": "DIFFERS", "why": "non-numeric columns differ"}
    x, y = a[num].to_numpy(float), b[num].to_numpy(float)
    both_nan = np.isnan(x) & np.isnan(y)
    d = np.where(both_nan, 0.0, np.abs(x - y))
    if np.isnan(d).any():
        return {"status": "DIFFERS", "why": "NaN pattern differs", "n_nan_mismatch": int(np.isnan(d).sum())}
    rel = d / np.maximum(1.0, np.abs(x))
    m = float(rel.max()) if rel.size else 0.0
    status = "REPRODUCED" if m <= 1e-9 else "NUMERICAL NOISE" if m <= 1e-6 else "DIFFERS"
    return {"status": status, "cells": int(x.size), "max_relative_difference": m,
            "n_cells_over_1e-6": int((rel > 1e-6).sum())}


def compare_output(rel: str, snap: Path) -> dict:
    """`rel` is a path under results/; "written=>original" compares a file a step had to write elsewhere
    (because its producer refuses to overwrite) with the snapshot's copy of the original."""
    new_rel, _, old_rel = rel.partition("=>")
    new, old = RES / new_rel, snap / (old_rel or new_rel)
    if not new.exists():
        return {"status": "NOT WRITTEN"}
    if not old.exists():
        return {"status": "NEW"}
    try:
        if rel.endswith(".json"):
            return compare_json(old, new)
        if rel.endswith(".csv"):
            return compare_csv(old, new)
        same = old.read_bytes() == new.read_bytes()
        return {"status": "REPRODUCED" if same else "DIFFERS"}
    except Exception as e:                                     # noqa: BLE001
        return {"status": "UNREADABLE", "error": str(e)[:300]}


def keep_registered(st: dict, rec: dict, snap: Path) -> list[str]:
    """A re-run writes in place. Where an output DIFFERS from the registered one, move the re-run's copy to
    results/verification/rerun_outputs/<step>/ and put the registered artefact back, so results/ always holds what
    the documents report while the difference stays on record (state.json and the moved file). Added 2026-10-07,
    after the registered pipeline's re-run differed in two timing-dependent methods (OPEN_DEFECTS D35)."""
    moved = []
    for o, res in rec["outputs"].items():
        if res.get("status") != "DIFFERS" or "=>" in o:
            continue
        new, old = RES / o, snap / o
        if not (new.exists() and old.exists()):
            continue
        dest = VDIR / "rerun_outputs" / st["id"] / o
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(new), str(dest))
        shutil.copy2(old, new)
        moved.append(o)
    return moved


# ------------------------------------------------------------------------------------ running
def load_state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def snapshot_dir() -> Path:
    if not LABEL_FILE.exists():
        raise SystemExit("BLOCKED: no snapshot yet -- run with --snapshot LABEL first")
    return config.PROJECT_ROOT / f"results_snapshot_{LABEL_FILE.read_text().strip()}"


def take_snapshot(label: str) -> None:
    dst = config.PROJECT_ROOT / f"results_snapshot_{label}"
    if dst.exists():
        print(f"snapshot {dst.name} already exists; not overwritten"); return
    skip = lambda d, names: [n for n in names if n.endswith(".tsv") and "robust" in n]  # noqa: E731
    shutil.copytree(RES, dst, ignore=skip)
    VDIR.mkdir(parents=True, exist_ok=True)
    LABEL_FILE.write_text(label)
    print(f"snapshot -> {dst.name}")


def run_step(st: dict, snap: Path) -> dict:
    VDIR.joinpath("logs").mkdir(parents=True, exist_ok=True)
    log = VDIR / "logs" / f"{st['id']}.log"
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write(" ".join(st["cmd"]) + "\n\n"); fh.flush()
        r = subprocess.run(st["cmd"], cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    secs = round(time.time() - t0, 1)
    rec = {"tier": st["tier"], "exit": r.returncode, "seconds": secs, "log": str(log.relative_to(config.PROJECT_ROOT)),
           "outputs": {o: compare_output(o, snap) for o in st["outputs"]}}
    rec["kept_registered"] = keep_registered(st, rec, snap)
    text = log.read_text(errors="ignore")
    if st["verdict"] == "no_disagrees":
        rec["verdict"] = "AGREES" if (r.returncode == 0 and "DISAGREES" not in text) else "DISAGREES"
    elif st["verdict"] == "exit0":
        rec["verdict"] = "PASSED" if r.returncode == 0 else "FAILED"
    return rec


def run_tier(tier: str, force: bool, only: list[str] | None, include_supplementary: bool = False) -> None:
    snap = snapshot_dir()
    state = load_state()
    for st in STEPS:
        if st["tier"] != tier or (only and st["id"] not in only):
            continue
        if st["id"] in state and not force:
            print(f"skip {st['id']} (done)"); continue
        if st["id"] in SUPPLEMENTARY and not include_supplementary:
            print(f"skip {st['id']} (supplementary; --include-supplementary to run)"); continue
        print(f"run  {st['id']} ...", flush=True)
        state[st["id"]] = run_step(st, snap)
        STATE.write_text(json.dumps(state, indent=2))
        o = state[st["id"]]
        stat = ", ".join(f"{k}: {v['status']}" for k, v in o["outputs"].items()) or o.get("verdict", "")
        print(f"     exit {o['exit']} in {o['seconds']} s | {stat}", flush=True)


# ------------------------------------------------------------------------------------ report
def provenance_gaps() -> list[str]:
    """Artefacts the manuscript reads that no verification step produces."""
    import re  # noqa: PLC0415
    src = (ROOT / "scripts" / "build_manuscript.py").read_text()
    read = sorted(set(re.findall(r'J\("([^"]+)"\)', src)))
    produced = {o.partition("=>")[2] or o for s in STEPS for o in s["outputs"]}
    return [a for a in read if a not in produced]


#: For each explained step, the pattern every one of its differences must match (checked in the report against the
#: full difference list recomputed from the moved-aside re-run copies, not just the stored examples).
EXPLAINED_PATTERNS = {
    "registered_pipeline": r"bayesprism|dwls|\[8\]|\[14\]|median_real_acs|control_verdict|primary_result\.|"
                           r"methods_excluded_for_partial_coverage",
    "tcga_gbm_frozen": r"bayesian_hierarchical|equal_footing|reference = 'frozen'",
    "tcga_lgg_frozen": r"bayesian_hierarchical|equal_footing|reference = 'frozen'",
    "anatomy_vs_biology": r"registered_status",
}
#: Where a step's re-run copies were moved by hand rather than by keep_registered.
MOVED_DIRS = {"registered_pipeline": "registered_pipeline_20261007"}

#: A difference with a diagnosed cause, stated where the difference is listed. Each names the defect record.
EXPLAINED = {
    "registered_pipeline": (
        "the per-method differences are confined to `bayesprism` and `dwls`: in the registered run (2026-09-21) both "
        "overran the genuine packages' 2,400 s budget by 5 and 16 s and fell back to the Python versions; in this "
        "re-run, with identical code, both genuine packages finished (BayesPrism ACS 0.815; DWLS 0.70 on 26 of 57 "
        "pairs, so it is excluded as partially covered). Every other method and both controls reproduce exactly. "
        "Through those two rows the run's own summary statistics move: the ACS-versus-ABSOLUTE agreement on this "
        "run's panel becomes +0.165 on 11 methods (the archived file carries the pre-D21 sign, -0.081, on 12; the "
        "registered value is +0.081 from `yardstick_agreement.py`), the synthetic arm 0.637 -> 0.828, and the median "
        "real ACS 0.854 -> 0.862. No reading changes. Registered artefacts were restored; the re-run's are in "
        "`results/verification/rerun_outputs/registered_pipeline_20261007/` (OPEN_DEFECTS D35, D21)."),
    "tcga_gbm_frozen": (
        "the only differences are additions. `bayesian_hierarchical` now runs: the registered frozen arm predates the "
        "`patient_id` fix and recorded `KeyError: 'patient_id'` (OPEN_DEFECTS D20); it now scores tumour rho 0.750. "
        "There are also two new provenance fields (`equal_footing`; `reference` = 'frozen'), and the estimates file names its key column `sample` "
        "instead of `sample_id`. Every registered method's statistics and all 1,848 registered estimate rows "
        "reproduce exactly (largest difference 0.0, checked 2026-10-07)."),
    "tcga_lgg_frozen": (
        "the same as tcga_gbm_frozen: `bayesian_hierarchical` now runs (it recorded `KeyError: 'patient_id'` in the "
        "registered frozen arm, OPEN_DEFECTS D20) and scores tumour rho 0.541; two provenance fields are new "
        "(`equal_footing`; `reference` = 'frozen'); the key column is `sample`. Every registered method's statistics and all 6,120 registered estimate "
        "rows reproduce exactly (largest difference 0.0, checked 2026-10-07)."),
    "anatomy_vs_biology": (
        "a deliberate text correction: the `registered_status` sentence now says the prediction was committed 94 "
        "minutes (2026-09-17 23:37:27 to 2026-09-18 01:11:52) before the LGG methylation existed, where it said 95, a "
        "rounding up (corrected 2026-10-07 in the 'verify ts' pass). No number the analysis computes changes."),
}


def explanation_check(sid: str, outputs: list[str]) -> str:
    """Recompute every difference of a step's DIFFERS JSON outputs and confirm each matches EXPLAINED_PATTERNS."""
    import re  # noqa: PLC0415
    from compare_artefacts import compare  # noqa: PLC0415
    pat = EXPLAINED_PATTERNS.get(sid)
    if not pat:
        return "no pattern declared"
    base = VDIR / "rerun_outputs" / MOVED_DIRS.get(sid, sid)
    n_all, outside = 0, []
    for o in outputs:
        if not o.endswith(".json"):
            continue
        new_p, old_p = base / o, snapshot_dir() / o
        if not (new_p.exists() and old_p.exists()):
            return f"cannot check: re-run copy of {o} not found"
        d = compare(json.loads(old_p.read_text()), json.loads(new_p.read_text()), 1e-9, DEFAULT_IGNORE)
        n_all += len(d)
        outside += [x for x in d if not re.search(pat, x, re.I)]
    if outside:
        return f"UNEXPLAINED: {len(outside)} of {n_all} differences fall outside the stated cause, e.g. {outside[0][:160]}"
    return f"checked: all {n_all} differences match the stated cause"


def report() -> int:
    state = load_state()
    rows, counts = [], {}
    for st in STEPS:
        rec = state.get(st["id"])
        if rec is None:
            rows.append((st, None)); continue
        rec = {**rec, "outputs": {k: _reclassify(_recompare_moved(st["id"], k, v)) for k, v in rec["outputs"].items()}}
        for o in rec["outputs"].values():
            counts[o["status"]] = counts.get(o["status"], 0) + 1
        if "verdict" in rec:
            counts[rec["verdict"]] = counts.get(rec["verdict"], 0) + 1
        rows.append((st, rec))
    gaps = provenance_gaps()
    out = {"snapshot": LABEL_FILE.read_text().strip() if LABEL_FILE.exists() else None, "counts": counts,
           "steps": state, "not_yet_run": [s["id"] for s, r in rows if r is None],
           "manuscript_artefacts_without_a_producing_step": gaps}
    REPORT_JSON.write_text(json.dumps(out, indent=2))
    L = ["# Verification re-run", "",
         f"Generated by `scripts/verify_rerun.py` from `{REPORT_JSON.relative_to(config.PROJECT_ROOT)}`. Every output is "
         f"compared with the snapshot `results_snapshot_{out['snapshot']}/` taken before anything re-ran.", "",
         "**Status key:**",
         "- REPRODUCED: identical within a relative 1e-9.",
         "- NUMERICAL NOISE: only numbers differ, all by at most 1e-6.",
         "- REPRODUCED + NEW FIELDS: every value the original recorded is reproduced; the re-run also writes "
         "fields the code gained after the original was written (listed below the table).",
         "- REPRODUCED + NEW ROWS: every registered row is present and identical; the re-run adds rows (a method "
         "that now runs) or renames a key column.",
         "- DIFFERS: anything else; the differences are listed below the table.",
         "- NEW: the file did not exist before the re-run.",
         "- NOT WRITTEN: the step should have written the file and did not.", "",
         "**Totals:** " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), "",
         "| step | tier | exit | time | outputs / verdict |", "|---|---|---|---|---|"]
    for st, rec in rows:
        if rec is None:
            why = (f"not re-run: supplementary ({SUPPLEMENTARY[st['id']]}), by decision 2026-10-07"
                   if st["id"] in SUPPLEMENTARY else "not yet run")
            L.append(f"| `{st['id']}` | {st['tier']} | -- | -- | {why} |"); continue
        outs = "; ".join(f"`{k}` {v['status']}" for k, v in rec["outputs"].items())
        if "verdict" in rec:
            outs = (outs + "; " if outs else "") + f"**{rec['verdict']}**"
        L.append(f"| `{st['id']}` | {st['tier']} | {rec['exit']} | {rec['seconds']:.0f} s | {outs} |")
    diffs = [(st, k, v) for st, rec in rows if rec for k, v in rec["outputs"].items() if v["status"] == "DIFFERS"]
    if diffs:
        L += ["", "## What differs", ""]
        for sid in sorted({st["id"] for st, _, _ in diffs} & set(EXPLAINED)):
            outs = [k for st, k, _ in diffs if st["id"] == sid]
            L += [f"> **Step `{sid}`, explained:** {EXPLAINED[sid]} *[{explanation_check(sid, outs)}]*", ""]
        for st, k, v in diffs:
            L.append(f"**`{k}`** (step `{st['id']}`):")
            for e in v.get("examples", [])[:8]:
                L.append(f"- `{e[:220]}`")
            if "why" in v:
                L.append(f"- {v['why']}")
            if "max_relative_difference" in v:
                L.append(f"- largest relative difference {v['max_relative_difference']:.3g}, "
                         f"{v['n_cells_over_1e-6']} of {v['cells']} cells over 1e-6")
            L.append("")
    added = [(st, k, v) for st, rec in rows if rec for k, v in rec["outputs"].items()
             if v["status"] == "REPRODUCED + NEW FIELDS"]
    if added:
        L += ["", "## Fields added since the original was written", ""]
        for st, k, v in added:
            L.append(f"**`{k}`** (step `{st['id']}`), {v['n_differences']} field(s):")
            L += [f"- `{e[:220]}`" for e in v.get("examples", [])]
            L.append("")
    failed = [st["id"] for st, rec in rows if rec and rec["exit"] != 0]
    if failed:
        L += ["", "## Steps that exited with an error", ""] + [f"- `{f}` (see its log)" for f in failed]
    if gaps:
        L += ["", "## Manuscript artefacts no step re-produces", "",
              "The manuscript reads these, but no script in this harness writes them. Each is a provenance gap to "
              "close (or an artefact a heavier step writes that is not yet listed):", ""] + [f"- `{g}`" for g in gaps]
    REPORT_MD.write_text("\n".join(L) + "\n")
    print("\n".join(L[:12]))
    print(f"wrote {REPORT_JSON.relative_to(config.PROJECT_ROOT)} and {REPORT_MD.relative_to(config.PROJECT_ROOT)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--snapshot", metavar="LABEL")
    ap.add_argument("--tier", choices=TIERS)
    ap.add_argument("--only", nargs="+")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--include-supplementary", action="store_true",
                    help="also re-run the supplementary heavy re-fits (SUPPLEMENTARY), skipped by default")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for s in STEPS:
            print(f"{s['tier']:6s} {s['id']:30s} -> {', '.join(s['outputs']) or s['verdict']}")
        return 0
    if a.snapshot:
        take_snapshot(a.snapshot)
    if a.tier:
        run_tier(a.tier, a.force, a.only, a.include_supplementary)
    if a.report or a.tier:
        report()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
