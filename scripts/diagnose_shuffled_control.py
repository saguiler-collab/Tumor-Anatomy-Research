#!/usr/bin/env python3
"""
diagnose_shuffled_control.py — run the control calibration on its own.

The calibration itself lives in `ivygap/anatomic/control_calibration.py` and runs as
part of every pipeline run. This script exists so it can be re-run against an existing
results tree — with more draws, say — without re-deconvolving anything.

It deliberately holds no arithmetic of its own. An earlier version of this script
reproduced the control's solve and its normalisation inline, which is exactly how a
diagnostic ends up disagreeing with the thing it is diagnosing.

    python scripts/diagnose_shuffled_control.py --draws 500
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from ivygap import config
from ivygap.anatomic import control_calibration
from ivygap.data.reference import frozen_gene_space, load_frozen_reference


def _load_scoring_cohort():
    """The ACS scoring set: H&E anatomic samples only, on the run's gene space."""
    bulk = pd.read_csv(config.BULK_EXPRESSION_PATH, sep="\t", index_col=0)
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", dtype=str).set_index("sample_id")
    man["is_anatomic_study"] = man["is_anatomic_study"] == "True"
    man["structure"] = man["structure"].map(config.canonical_structure)

    keep = [s for s in bulk.columns
            if s in man.index
            and man.loc[s, "is_anatomic_study"]
            and man.loc[s, "structure"] in config.PRIMARY_STRUCTURES]

    frozen = load_frozen_reference()
    genes, prov = frozen_gene_space(frozen, bulk.index)
    return bulk.loc[genes, keep], man.loc[keep], frozen.subset_genes(genes), prov


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draws", type=int, default=control_calibration.DEFAULT_DRAWS)
    ap.add_argument("--out", default="results/anatomic/control_calibration.json")
    args = ap.parse_args()

    bulk, man, ref, gene_prov = _load_scoring_cohort()
    print(f"cohort: {bulk.shape[1]} samples, {man['patient_id'].nunique()} tumours, "
          f"{bulk.shape[0]:,} genes")

    observed = None
    lb = Path("results/anatomic/acs_leaderboard.csv")
    if lb.exists():
        observed = pd.read_csv(lb, index_col=0)["acs"].to_dict()

    report = control_calibration.calibrate(bulk, man, ref, observed=observed,
                                           n_draws=args.draws)
    report["gene_space"] = gene_prov

    a = report["acs"]
    print(f"\n--- shuffled-signature control, {args.draws} independent permutations ---")
    print(f"ACS  mean {a['mean']:.3f}  sd {a['sd']:.3f}  "
          f"range {a['min']:.3f}-{a['max']:.3f}  "
          f"(5-95%: {a['p05']:.3f}-{a['p95']:.3f})")
    if "as_run_draw" in report:
        d = report["as_run_draw"]
        print(f"the leaderboard's single draw scored {d['acs']:.3f}, at the "
              f"{d['percentile_among_draws']:.0%} percentile")
    print("\nfraction of evaluable pairs a MEANINGLESS signature satisfies:")
    for cid, rate in sorted(report["per_constraint_satisfied_rate"].items(),
                            key=lambda kv: -kv[1]):
        print(f"  {cid}: {rate:.2f}")
    print(f"\n{report['verdict']}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
