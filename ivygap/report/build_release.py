"""
build_release.py — assemble the stable, public-facing snapshot.

WHAT A RELEASE BUNDLE IS FOR
----------------------------
Anything consuming this project — a website, a figure script, a collaborator — should
read `release/` and nothing else. The pipeline directories underneath it accumulate
intermediates, superseded runs and sensitivity analyses, and asking a consumer to work
out which of those is current is how stale numbers end up in a talk.

Every file in the bundle is copied from a live source path and recorded with that
path and its SHA-256, so any number in the release can be traced back to the artefact
that produced it. The bundle is always rebuilt from source; it never copies from
itself.

PROVENANCE THIS BUNDLE REFUSES TO LOSE
--------------------------------------
Three things are carried into `manifest.json` because losing them would make the
numbers misleading rather than merely unsourced:

  * whether each of the four published tools ran as the genuine R package or as this
    project's Python reimplementation,
  * whether the run used real Ivy GAP data or the synthetic fixture reference,
  * the survival power verdict, so a C-index can never be quoted without the note
    that the cohort could not distinguish methods on it.
"""

from __future__ import annotations

import json
import shutil
from datetime import date
from pathlib import Path

from ivygap import config
from ivygap.anatomic import constraints as K

#: (source relative to results/, destination name, description)
BUNDLE = [
    ("anatomic/acs_leaderboard.csv", "acs_leaderboard.csv",
     "THE LEADERBOARD. Anatomic Concordance Score per method, with a bootstrap CI over "
     "tumours and a p-value against the within-tumour label-permutation null. The two "
     "rows flagged is_control are negative controls and are what make the rest "
     "interpretable - read the control verdict in anatomic_report.json first."),
    ("anatomic/agreement_test.csv", "agreement_test.csv",
     "THE PRIMARY RESULT. Spearman rank correlation between the ACS ranking and the "
     "ranking from real ground truth, per yardstick, against a threshold fixed before "
     "it was computed. A yardstick marked UNAVAILABLE produced no scores in this run - "
     "that is not the same as disagreeing."),
    ("anatomic/agreement_report.json", "agreement_report.json",
     "The primary result with its pre-registered threshold and the requirement it was "
     "judged against (rho >= 0.60 AND a bootstrap CI excluding zero)."),
    ("anatomic/agreement_method_ranks.csv", "agreement_method_ranks.csv",
     "Each method's ACS rank beside its ground-truth rank, with the gap. This is where "
     "the agreement test is actually legible: a large gap means a method the anatomy "
     "likes and the truth does not, or vice versa."),
    ("benchmark/simulated_yardstick_accuracy.csv", "simulated_yardstick_accuracy.csv",
     "Per-method accuracy on the SIMULATED yardstick. Not the protocol's yardstick 1 - "
     "mixtures come from donor-perturbed copies of the frozen signature, not from real "
     "cells. Read simulated_yardstick_provenance.json before using it."),
    ("benchmark/simulated_yardstick_provenance.json",
     "simulated_yardstick_provenance.json",
     "What the simulated yardstick is, how hard it was made, and the explicit warning "
     "that it must not be reported as the protocol's ground truth."),
    ("anatomic/acs_per_constraint.csv", "acs_per_constraint.csv",
     "Every constraint x every method: how many tumours could evaluate it and how many "
     "satisfied it. A constraint no method satisfies is evidence about the constraint, "
     "or about the reference - read it before reading the leaderboard."),
    ("anatomic/constraint_file.json", "constraint_file.json",
     "The frozen, pre-registered constraint file: C1-C7 with directions, weights, "
     "histological citations, the declared T-cell exclusion, and the facts this "
     "signature is structurally unable to test. Its SHA-256 is the registration hash."),
    ("anatomic/anatomic_report.json", "anatomic_report.json",
     "Machine-readable summary, including the control verdict and the note that ACS "
     "was never used to select a method."),
    ("anatomic/composition_by_tumor_structure.csv", "composition.csv",
     "Composition per (tumour, anatomic structure) from the best-scoring real method. "
     "One row per tumour per structure - already collapsed, so uneven block counts "
     "cannot skew anything downstream."),
    ("benchmark/benchmark_summary.csv", "benchmark_summary.csv",
     "Per-method accuracy on donor-held-out pseudobulk with known composition. This is "
     "the ground-truth ranking the agreement test correlates ACS against, and the only "
     "evidence any method is selected on."),
    ("benchmark/benchmark_per_type.csv", "benchmark_per_type.csv",
     "Per method, per cell type: MAE, RMSE, bias, correlation, the fraction invented "
     "where the type is truly absent, and the limit of detection."),
    ("benchmark/benchmark_per_niche.csv", "benchmark_per_niche.csv",
     "The same accuracy broken down by which anatomic niche each synthetic mixture was "
     "shaped after. A method can win overall and still fail at the leading edge."),
    ("benchmark/uncertainty_calibration.csv", "uncertainty_calibration.csv",
     "For every method reporting credible intervals: how often the interval actually "
     "contained the truth, against a nominal 0.95. Read before quoting any interval."),
    ("benchmark/method_selection_decision.json", "method_selection_decision.json",
     "Which method was selected, on what criterion, with the full ranking - and the "
     "explicit record that survival and the anatomic labels were not used."),
    ("benchmark/gene_space.json", "gene_space.json",
     "Which genes every method was given, how many were shared with the frozen "
     "signature, and why the subset was taken. Part of the equal-footing claim: the "
     "gene set is chosen once, for everyone, before any method runs."),
    ("anatomic/equal_footing_certificate.json", "equal_footing_certificate.json",
     "Content hashes proving every method and control was given identical inputs, "
     "identical tumours and identical folds."),
    ("anatomic/implementation_report.json", "implementation_report.json",
     "For each of MuSiC, DWLS, Bisque and SCDC: whether the genuine R package ran, and "
     "if not, exactly why the Python reimplementation was used instead."),
    ("anatomic/control_calibration.json", "control_calibration.json",
     "Each negative control scored over many independent draws, not once. The "
     "leaderboard's control rows are single draws; on real data the shuffled-signature "
     "control has a standard deviation of about 0.16 against a real-method spread of "
     "0.29, so one draw cannot establish that the controls behave in either direction. "
     "Read this before reading a control row."),
    ("anatomic/acs_per_tumor.csv", "acs_per_tumor.csv",
     "Every method x every tumour: which constraints that tumour could evaluate and "
     "which it satisfied. The protocol asks for per-tumour results explicitly, because "
     "with 9 evaluable tumours a pooled ACS can be carried by one or two of them and "
     "the pooled number alone cannot show that."),
    ("anatomic/acs_cohort_sensitivity.csv", "acs_cohort_sensitivity.csv",
     "ACS for each method when only the 122 anatomic samples are deconvolved, beside "
     "ACS when all 270 archive samples are, with the rank change. Several methods use "
     "cross-sample statistics, so this measures how much of the leaderboard is a "
     "property of the method rather than of the cohort it was handed."),
    ("data_reconciliation.json", "data_reconciliation.json",
     "Protocol step 3's gate: our parse of the frozen 2014-11-25 archive checked "
     "sample-by-sample and structure-by-structure against the Allen Institute portal's "
     "own live metadata export, which is produced independently of this pipeline."),
    ("clinical_missingness_audit.json", "clinical_missingness_audit.json",
     "Whether a blank survival_days is missing at random. It is not: the blanks are "
     "strongly enriched for MGMT methylation and younger age, the two strongest "
     "favourable prognostic factors in GBM. Read before any survival number."),
    ("anatomic/full_database/acs_leaderboard.csv", "full_database_acs_leaderboard.csv",
     "The leaderboard recomputed when every archive sample is deconvolved rather than "
     "only the anatomic ones. ACS is still scored on the anatomic subset alone - the "
     "ISH-cluster samples had their structure assigned from expression, so scoring "
     "them would be circular."),
    ("anatomic/full_database/anatomic_report.json", "full_database_anatomic_report.json",
     "The full-database run's machine-readable summary, including how many samples "
     "were deconvolved versus how many were scored."),
    ("survival/survival_power.json", "survival_power.json",
     "The canonical survival verdict. Ivy GAP publishes survival_days and no "
     "vital-status column, so who was censored is unknowable and this is BLOCKED. "
     "Nothing is imputed."),
    ("survival/declared_event_policy/survival_power.json",
     "survival_power_declared_policy.json",
     "Present only if the run was explicitly asked to read recorded survival times as "
     "observed deaths. Carries the assumption, the missingness audit, and what the "
     "event count can actually detect - computed before any model was fitted."),
    ("survival/declared_event_policy/survival_metrics.csv",
     "survival_metrics_declared_policy.csv",
     "Out-of-fold prognostic evaluation per method under that declared assumption. "
     "NOT the canonical result. Read survival_power_declared_policy.json first; on this "
     "many events every row is INCONCLUSIVE by design."),
]


def _copy(src: Path, dst: Path) -> dict | None:
    if not src.exists():
        return None
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return {"sha256": config.sha256_file(dst), "bytes": dst.stat().st_size}


def build(synthetic: bool = False, notes: str = "") -> dict:
    """Rebuild `release/` from the live result files. Returns the manifest."""
    config.ensure_dirs()
    out = config.RELEASE_DIR

    files, missing = [], []
    for rel_src, dst_name, description in BUNDLE:
        src = config.RESULTS_DIR / rel_src
        info = _copy(src, out / dst_name)
        if info is None:
            missing.append(rel_src)
            continue
        files.append({
            "public_path": f"release/{dst_name}",
            "source_path": str(src.relative_to(config.PROJECT_ROOT)),
            "description": description,
            **info,
        })

    decision_path = config.BENCH_DIR / "method_selection_decision.json"
    decision = json.loads(decision_path.read_text()) if decision_path.exists() else {}

    impl_path = config.BENCH_DIR / "implementation_report.json"
    implementations = json.loads(impl_path.read_text()) if impl_path.exists() else []

    power_path = config.SURVIVAL_DIR / "survival_power.json"
    power = json.loads(power_path.read_text()) if power_path.exists() else {}

    agree_path = config.ANATOMIC_DIR / "agreement_report.json"
    agree = json.loads(agree_path.read_text()) if agree_path.exists() else {}

    manifest = {
        "project": "Ivy GAP GBM anatomic cell-type deconvolution",
        "built": date.today().isoformat(),
        "data_source": {
            "bulk": "Ivy Glioblastoma Atlas Project (Allen Institute), "
                    "anatomic-structure RNA-seq",
            "url": "https://glioblastoma.alleninstitute.org/",
            "citation": "Puchalski RB et al. Science 360:660-663 (2018). "
                        "doi:10.1126/science.aaf2666",
            "replaces": "TCGA-GBM bulk RNA-seq (the predecessor project's cohort)",
        },
        "run_type": "SYNTHETIC FIXTURE" if synthetic else "real data",
        "synthetic_warning": (
            "This bundle was produced from a GENERATED reference, not from real "
            "single-cell data. The numbers validate that the pipeline runs and that "
            "its scoring behaves correctly. They are NOT biological findings and must "
            "never be presented as such."
        ) if synthetic else None,
        "selected_method": decision.get("selected_method"),
        "selection_criterion": decision.get("criterion"),
        "selection_is_outcome_blind": decision.get("criterion_is_outcome_blind"),
        "published_tool_implementations": implementations,
        "survival_power": power.get("verdict"),
        "primary_result_status": agree.get("primary_result_status"),
        "yardsticks_with_real_ground_truth": agree.get(
            "yardsticks_with_real_ground_truth", []),
        "constraint_freeze_hash": K.freeze_hash(),
        "constraint_file": K.registration_payload(),
        "constraint_coverage": K.coverage_report(),
        "cell_types_primary": config.PRIMARY_CELL_TYPES,
        "cell_types_sidecar": config.SIDECAR_CELL_TYPES,
        "anatomic_structures": config.PRIMARY_STRUCTURES,
        "files": files,
        "missing_from_this_build": missing,
        "notes": notes,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    _write_index(out, manifest)
    return manifest


def _write_index(out: Path, manifest: dict) -> None:
    """A short human-readable front page for the bundle."""
    lines = [
        "# Release bundle — Ivy GAP GBM anatomic deconvolution",
        "",
        f"Built {manifest['built']} · run type: **{manifest['run_type']}**",
        "",
    ]
    if manifest.get("synthetic_warning"):
        lines += ["> **" + manifest["synthetic_warning"] + "**", ""]

    if manifest.get("selected_method"):
        lines += [
            f"Selected method: **{manifest['selected_method']}**",
            f"Selected on: {manifest.get('selection_criterion')}",
        ]
    else:
        # "None" on its own reads like a bug. Say what actually happened: selection needs
        # known-truth mixtures, which need cells, and no method is frozen without them.
        lines += [
            "**No method is frozen by this run.** Method selection requires known-truth "
            "mixtures, which require cell-level single-cell data; none was available, so "
            "the benchmark stage did not run. The ACS leaderboard below still measures "
            "real tissue, but ACS is the object under test and is never used to select.",
        ]
    lines += ["", "## Files", ""]
    for f in manifest["files"]:
        lines.append(f"- `{Path(f['public_path']).name}` — {f['description']}")

    if manifest["missing_from_this_build"]:
        lines += ["", "## Not in this build", "",
                  "These stages had not produced output when the bundle was built:", ""]
        lines += [f"- `{m}`" for m in manifest["missing_from_this_build"]]

    if manifest.get("survival_power"):
        lines += ["", "## Survival", "", manifest["survival_power"]]

    primary = manifest.get("primary_result_status")
    if primary:
        lines += ["", "## The primary result", "", primary]

    lines += [
        "", "## How the published tools were run", "",
        "| tool | genuine R package | reason |", "|---|---|---|",
    ]
    for impl in manifest.get("published_tool_implementations", []):
        lines.append(f"| {impl['method']} | {'yes' if impl['r_available'] else 'no'} "
                     f"| {impl['reason']} |")

    (out / "RELEASE_INDEX.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--synthetic", action="store_true")
    a = ap.parse_args()
    m = build(synthetic=a.synthetic)
    print(f"wrote {len(m['files'])} files to {config.RELEASE_DIR}")
    if m["missing_from_this_build"]:
        print(f"missing: {m['missing_from_this_build']}")
