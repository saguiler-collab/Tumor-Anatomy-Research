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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=str(config.RESULTS_DIR))
    ap.add_argument("--out", default=str(config.FIGURES_DIR / "figure_data.json"))
    a = ap.parse_args()

    data = collect(Path(a.results))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2, default=str))

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
