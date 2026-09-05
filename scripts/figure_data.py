#!/usr/bin/env python3
"""
figure_data.py — collect everything the figures need into one JSON.

The figures are rendered from artefacts, never from numbers typed by hand. This script
is the single seam between the two: it reads the results tree and emits one object, so a
figure that shows a number nobody computed is impossible rather than merely unlikely.

    python scripts/figure_data.py --out results/figures/figure_data.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from ivygap import config


def _csv(path: Path, **kw):
    return pd.read_csv(path, **kw) if path.exists() else None


def _json(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:                                       # noqa: BLE001
        return None


def collect(results: Path) -> dict:
    anat = results / "anatomic"
    bench = results / "benchmark"
    surv = results / "survival"

    out: dict = {"generated_from": str(results), "figures_available": []}

    def have(name: str, ok: bool) -> None:
        if ok:
            out["figures_available"].append(name)

    # --- provenance -----------------------------------------------------------
    recon = _json(results / "data_reconciliation.json")
    load = _json(config.PROCESSED_DIR / "load_report.json")
    sampling = _json(config.REFERENCE_DIR / "reference_sampling_gbmap.json")
    out["provenance"] = {"reconciliation": recon, "load": load, "atlas_sampling": sampling}
    have("provenance_flow", bool(recon and load))

    # --- the constraint set, drawn --------------------------------------------
    cfile = _json(anat / "constraint_file.json")
    out["constraints"] = cfile
    have("constraint_map", bool(cfile))

    # --- leaderboard + controls -----------------------------------------------
    lb = _csv(anat / "acs_leaderboard.csv", index_col=0)
    cal = _json(anat / "control_calibration.json")
    if lb is not None:
        out["leaderboard"] = [
            {"method": m, **{k: (None if pd.isna(v) else v) for k, v in r.items()}}
            for m, r in lb.iterrows()]
        out["control_calibration"] = cal
    have("leaderboard", lb is not None)
    have("control_band", bool(cal))

    # --- resolution -----------------------------------------------------------
    pc = _csv(anat / "acs_per_constraint.csv")
    if lb is not None and pc is not None:
        one = pc[pc["method"] == pc["method"].iloc[0]]
        den = float((one["n_tumors_evaluable"] * one["weight"]).sum())
        real = lb[~lb["is_control"]]
        out["resolution"] = {
            "weighted_denominator": den,
            "n_possible_values": int(den) + 1,
            "granularity": (1.0 / den) if den else None,
            "n_methods": int(len(lb)),
            "n_distinct_all": int(lb["acs"].nunique()),
            "n_real_methods": int(len(real)),
            "n_distinct_real": int(real["acs"].nunique()),
            "ties": {f"{v:.4f}": list(g.index)
                     for v, g in lb.groupby(lb["acs"].round(9)) if len(g) > 1},
        }
    have("resolution_lattice", "resolution" in out)

    # --- per constraint and per tumour ----------------------------------------
    out["per_constraint"] = pc.to_dict(orient="records") if pc is not None else None
    comp = _csv(anat / "composition_by_tumor_structure.csv")
    out["composition"] = comp.to_dict(orient="records") if comp is not None else None
    pt = _csv(anat / "acs_per_tumor.csv")
    out["per_tumor"] = pt.to_dict(orient="records") if pt is not None else None
    have("evidence_matrix", pt is not None)

    # --- the primary result ---------------------------------------------------
    agree = _json(anat / "agreement_report.json")
    ranks = _csv(anat / "agreement_method_ranks.csv")
    out["agreement"] = agree
    out["agreement_ranks"] = ranks.to_dict(orient="records") if ranks is not None else None
    have("agreement_scatter", ranks is not None and len(ranks) > 0)

    # --- benchmark accuracy ---------------------------------------------------
    bs = _csv(bench / "benchmark_summary.csv", index_col=0)
    out["benchmark"] = ([{"method": m, **{k: (None if pd.isna(v) else v)
                                          for k, v in r.items()}}
                         for m, r in bs.iterrows()] if bs is not None else None)
    out["selection"] = _json(bench / "method_selection_decision.json")
    have("benchmark_accuracy", bs is not None)

    # --- cohort sensitivity ---------------------------------------------------
    sens = _csv(anat / "acs_cohort_sensitivity.csv", index_col=0)
    out["cohort_sensitivity"] = ([{"method": m, **{k: (None if pd.isna(v) else v)
                                                   for k, v in r.items()}}
                                  for m, r in sens.iterrows()]
                                 if sens is not None else None)
    have("cohort_sensitivity", sens is not None)

    # --- survival -------------------------------------------------------------
    out["survival_power"] = _json(surv / "survival_power.json")
    out["survival_declared"] = _json(surv / "declared_event_policy" / "survival_power.json")
    sm = _csv(surv / "declared_event_policy" / "survival_metrics.csv", index_col=0)
    out["survival_metrics"] = ([{"method": m, **{k: (None if pd.isna(v) else v)
                                                 for k, v in r.items()}}
                                for m, r in sm.iterrows()] if sm is not None else None)
    out["missingness"] = _json(results / "clinical_missingness_audit.json")
    have("power_wall", bool(out["survival_declared"]))
    have("missingness_mosaic", bool(out["missingness"]))

    # --- integrity ------------------------------------------------------------
    out["registration"] = _json(anat / "registration_status.json")
    out["reference_coverage"] = _json(anat / "reference_coverage.json")
    out["method_configs"] = _json(anat / "method_configs.json")
    out["implementations"] = _json(anat / "implementation_report.json")
    have("coverage", bool(out["reference_coverage"]))

    return out


def trim_for_figures(d: dict) -> dict:
    """
    The compact object the figure page embeds.

    figure_data.json carries everything, including per-method x per-tumour detail that
    the plates summarise rather than plot point-by-point. This keeps only what a figure
    actually draws, so the page stays small and its script stays legible.
    """
    lb = d.get("leaderboard") or []
    real = [r for r in lb if not r.get("is_control")]
    ctrl = [r for r in lb if r.get("is_control")]

    cal = d.get("control_calibration") or {}
    hardest = cal.get("hardest_control")
    band = (cal.get("controls", {}).get(hardest, {}) or {}).get("acs") if hardest else None

    out = {
        "cohort": {
            "n_samples_scored": (d.get("provenance", {}).get("load") or {}).get(
                "n_anatomic_study_samples"),
            "n_tumors_scored": (d.get("provenance", {}).get("load") or {}).get(
                "n_anatomic_study_tumors"),
            "n_samples_archive": (d.get("provenance", {}).get("load") or {}).get("n_samples"),
            "structures": (d.get("provenance", {}).get("load") or {}).get("structures_found"),
            "atlas": {
                "cells_kept": (d.get("provenance", {}).get("atlas_sampling") or {}).get("n_cells_kept"),
                "donors": (d.get("provenance", {}).get("atlas_sampling") or {}).get("n_donors_kept"),
                "genes": (d.get("provenance", {}).get("atlas_sampling") or {}).get("n_genes_kept"),
            },
        },
        "constraints": [
            {"id": c["id"], "kind": c["kind"], "cell_type": c["cell_type"],
             "structures": c["structures"], "among": c.get("among", []),
             "weight": c["weight"], "rationale": c["rationale"]}
            for c in (d.get("constraints") or {}).get("constraints", [])
        ],
        "leaderboard": [
            {"method": r["method"], "acs": r["acs"], "ci_low": r["ci_low"],
             "ci_high": r["ci_high"], "null_mean": r.get("null_mean"),
             "null_p": r.get("null_p"), "is_control": bool(r.get("is_control")),
             "degenerate": bool(r.get("degenerate")),
             "implementation": r.get("implementation"),
             "tie_group": r.get("acs_tie_group") or ""}
            for r in lb
        ],
        "control_band": band,
        "hardest_control": hardest,
        "not_distinguishable": cal.get("not_distinguishable_from_noise", []),
        "resolution": d.get("resolution"),
        "per_constraint": d.get("per_constraint"),
        "composition": d.get("composition"),
        "per_tumor": d.get("per_tumor"),
        "agreement": d.get("agreement"),
        "agreement_ranks": d.get("agreement_ranks"),
        "benchmark": d.get("benchmark"),
        "selection": {k: (d.get("selection") or {}).get(k)
                      for k in ("selected_method", "selection_basis", "criterion",
                                "methods_tied_with_best", "n_test_mixtures",
                                "train_donors", "test_donors")},
        "cohort_sensitivity": d.get("cohort_sensitivity"),
        "survival": {
            "strict": d.get("survival_power"),
            "declared": {k: (d.get("survival_declared") or {}).get(k)
                         for k in ("n_patients", "n_events", "detectable_delta_c_index",
                                   "adequately_powered", "verdict")},
            "metrics": d.get("survival_metrics"),
        },
        "missingness": d.get("missingness"),
        "registration": d.get("registration"),
        "provenance_recon": d.get("provenance", {}).get("reconciliation"),
        "coverage": {k: (d.get("reference_coverage") or {}).get(k)
                     for k in ("atlas_labels_dropped", "n_cells_dropped",
                               "fraction_of_atlas_dropped", "what_this_means",
                               "correction_to_the_frozen_constraint_file")},
        "implementations": d.get("implementations"),
        "deviations": (d.get("method_configs") or {}).get(
            "methods_with_declared_deviations", []),
        "n_real_methods": len(real),
        "n_controls": len(ctrl),
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=str(config.RESULTS_DIR))
    ap.add_argument("--out", default=str(config.FIGURES_DIR / "figure_data.json"))
    a = ap.parse_args()

    data = collect(Path(a.results))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2, default=str))

    trimmed = trim_for_figures(data)
    tpath = out.with_name("figure_data_compact.json")
    tpath.write_text(json.dumps(trimmed, separators=(",", ":"), default=str))
    print(f"wrote {tpath} ({tpath.stat().st_size / 1024:.0f} KB)")

    print(f"wrote {out}")
    print(f"figures with data: {len(data['figures_available'])}")
    for f in data["figures_available"]:
        print(f"  - {f}")
    missing = [k for k in ("leaderboard", "agreement_ranks", "per_tumor",
                           "benchmark", "survival_metrics") if not data.get(k)]
    if missing:
        print(f"\nNOT YET AVAILABLE: {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
