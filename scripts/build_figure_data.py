#!/usr/bin/env python3
"""
build_figure_data.py — emit the JSON blob the figure plates are drawn from.

WHY
---
The published plates carry a `const D = {...}` blob that every SVG is rendered from.
It was assembled by hand once, and it went stale the moment the pipeline was rerun: the
plates still showed ten real methods, DWLS at 0.877 as `R:DWLS`, and 200 test mixtures,
long after the run had thirteen comparable methods, DWLS at 0.738 as a reimplementation,
and 500 mixtures. A figure that disagrees with the artefacts is worse than no figure,
because it is quoted with confidence.

So the blob is generated, exactly as RESULTS.md is. Every key below names the artefact it
reads. Nothing is transcribed.

    python scripts/build_figure_data.py                     # to stdout
    python scripts/build_figure_data.py -o figure_data.json
    python scripts/build_figure_data.py --results-dir results_archive/<stamp>/results
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def _json(path: Path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _csv(path: Path):
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def _records(path: Path) -> list:
    """
    Rows of a CSV, or [] if it is absent.

    Not `_csv(p) or []` — a DataFrame has no truth value, and pandas raises rather than
    silently picking one.
    """
    df = _csv(path)
    return [] if df is None else _clean(df.to_dict("records"))


def _clean(o):
    """NaN is not JSON. Convert to None so the plates get null, not a syntax error."""
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
        return None
    return o


def build(results: Path) -> dict:
    anat = results / "anatomic"
    bench = results / "benchmark"
    surv = results / "survival"

    rep = _json(anat / "anatomic_report.json") or {}
    cert = _json(anat / "equal_footing_certificate.json") or {}
    cov = _json(anat / "reference_coverage.json") or {}
    recon = _json(results / "data_reconciliation.json") or {}
    cal = _json(anat / "control_calibration.json") or {}

    lb = _csv(anat / "acs_leaderboard.csv")
    if lb is None:
        raise SystemExit(f"no acs_leaderboard.csv under {anat}")
    lb = lb.copy()
    if "acs_tie_group" in lb.columns:
        # the sentinel "-" means "no tie group"; the plates expect an empty string
        lb["tie_group"] = lb["acs_tie_group"].fillna("").replace("-", "")
    real = lb[~lb["is_control"].astype(bool)]
    ctrl = lb[lb["is_control"].astype(bool)]
    comparable = real[real["comparable"].astype(bool)] if "comparable" in real else real

    # --- the control band: the HARDEST control's calibrated distribution ------------
    hardest = cal.get("hardest_control")
    controls = cal.get("controls") or {}
    # `controls` is keyed by control name, and each entry holds its ACS distribution
    # under "acs". Iterating it as a list of dicts yields bare strings and silently
    # produces an empty band, which the plate would then draw as a zero-width line.
    entry = controls.get(hardest) if isinstance(controls, dict) else None
    band = dict((entry or {}).get("acs") or {})
    if band:
        band["n_draws"] = (entry or {}).get("n_draws") or cal.get("n_draws")
    # Calibration measures the control as a DISTRIBUTION. A band built from one draw
    # would be a point, and the whole reason the plate shows a band is that a single
    # draw of this control has a standard deviation of ~0.066.

    # --- resolution: how coarse ACS actually is on this cohort ----------------------
    denom = None
    pc = _csv(anat / "acs_per_constraint.csv")
    if pc is not None and len(pc):
        one = pc[pc["method"] == real.iloc[0]["method"]]
        denom = float((one["weight"] * one["n_tumors_evaluable"]).sum())
    ties = defaultdict(list)
    for _, r in comparable.iterrows():
        ties[f"{r['acs']:.4f}"].append(r["method"])
    resolution = {
        "weighted_denominator": denom,
        "n_possible_values": int(denom) + 1 if denom else None,
        "granularity": 1.0 / denom if denom else None,
        "n_methods": int(len(lb)),
        "n_distinct_all": int(lb["acs"].round(6).nunique()),
        "n_real_methods": int(len(real)),
        "n_distinct_real": int(real["acs"].round(6).nunique()),
        "ties": {k: v for k, v in ties.items() if len(v) > 1},
    }

    # --- which real methods the hardest control cannot be told apart from -----------
    not_dist = []
    if band.get("p95") is not None:
        not_dist = sorted(comparable[comparable["acs"] <= band["p95"]]["method"])

    agree_rep = _json(anat / "agreement_report.json") or {}
    at = _csv(anat / "agreement_test.csv")
    if at is not None and "by_yardstick" not in agree_rep:
        agree_rep["by_yardstick"] = _clean(at.to_dict("records"))

    surv_block = {
        "strict": _json(surv / "survival_power.json") or {},
        "declared": _json(surv / "declared_event_policy" / "survival_power.json") or {},
    }
    dm = _csv(surv / "declared_event_policy" / "survival_metrics.csv")
    if dm is not None:
        surv_block["declared_metrics"] = _clean(dm.to_dict("records"))

    cfgs = _json(anat / "method_configs.json") or {}

    return _clean({
        "cohort": {
            "n_samples_scored": rep.get("n_samples"),
            "n_tumors_scored": rep.get("n_tumors"),
            "n_samples_archive": rep.get("n_samples_deconvolved") or recon.get("n_archive_samples"),
            "structures": rep.get("structures"),
            "atlas": {
                "cells_kept": cov.get("n_cells_kept"),
                "donors": cert.get("n_patients"),
                "genes": cert.get("input_hashes", {}).get("n_genes"),
            },
        },
        "constraints": (_json(anat / "constraint_file.json") or {}).get("constraints"),
        "leaderboard": _clean(lb.to_dict("records")),
        "control_band": band,
        "hardest_control": hardest,
        "not_distinguishable": not_dist,
        "resolution": resolution,
        "per_constraint": _clean(pc.to_dict("records")) if pc is not None else [],
        "composition": _records(anat / "composition_by_tumor_structure.csv"),
        "per_tumor": _records(anat / "acs_per_tumor.csv"),
        "agreement": agree_rep,
        "agreement_ranks": _records(anat / "agreement_method_ranks.csv"),
        "benchmark": _records(bench / "benchmark_summary.csv"),
        "selection": _json(bench / "method_selection_decision.json") or {},
        "cohort_sensitivity": _records(anat / "acs_cohort_sensitivity.csv"),
        "survival": surv_block,
        "missingness": _json(results / "clinical_missingness_audit.json") or {},
        "registration": _json(anat / "registration_status.json") or {},
        "provenance_recon": recon,
        "coverage": cov,
        "implementations": _json(anat / "implementation_report.json") or [],
        "deviations": cfgs.get("methods_with_declared_deviations", []),
        "n_real_methods": int(len(real)),
        "n_comparable_methods": int(len(comparable)),
        "n_controls": int(len(ctrl)),
    })


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("-o", "--out", help="write here instead of stdout")
    args = ap.parse_args()

    data = build(Path(args.results_dir))
    text = json.dumps(data, separators=(",", ":"), sort_keys=False)
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {args.out}  ({len(text) / 1024:.1f} KB)")
        print(f"  {data['n_real_methods']} real methods "
              f"({data['n_comparable_methods']} comparable), "
              f"{data['n_controls']} controls")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
