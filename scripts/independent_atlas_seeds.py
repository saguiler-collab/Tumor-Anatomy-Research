"""
independent_atlas_seeds.py -- apply Addendum 2 of prespecified/abdelfattah_cluster_mapping.md.

A cohort-method cell is SEED-STABLE if its cohort-level lymphoid ordering under every seed build
matches the registered build (seed 0), SEED-DEPENDENT otherwise. Written before any seed result
existed. Declared after the main result was seen, so it can only qualify the main reading.

    python3 scripts/independent_atlas_seeds.py     # writes results/independent_atlas_seed_sensitivity.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402

SEEDS = (1, 2, 3)


def classify(registered: str, others: list[str | None]) -> str:
    if any(o is None for o in others) or not others:
        return "BLOCKED"
    return "SEED-STABLE" if all(o == registered for o in others) else "SEED-DEPENDENT"


def orderings(path: Path) -> dict:
    d = json.loads(path.read_text())
    return {(c, m): v.get("ordering") for c, cc in d["cohorts"].items()
            for m, v in cc["methods"].items()}


def main() -> int:
    base = config.RESULTS_DIR / "independent_atlas_test.json"
    if not base.exists():
        print("BLOCKED: no registered independent-atlas result"); return 2
    reg = orderings(base)
    runs = {k: (config.RESULTS_DIR / f"independent_atlas_test_seed{k}.json") for k in SEEDS}
    got = {k: orderings(p) if p.exists() else {} for k, p in runs.items()}
    cells = {}
    for key, o in reg.items():
        others = [got[k].get(key) for k in SEEDS]
        cells[f"{key[0]}_{key[1]}"] = {"registered": o,
                                       "by_seed": dict(zip(map(str, SEEDS), others)),
                                       "reading": classify(o, others)}
    out = {"rule": "prespecified/abdelfattah_cluster_mapping.md, Addendum 2",
           "seeds": list(SEEDS), "missing_runs": [k for k in SEEDS if not runs[k].exists()],
           "cells": cells}
    (config.RESULTS_DIR / "independent_atlas_seed_sensitivity.json").write_text(json.dumps(out, indent=2))
    for k, v in cells.items():
        print(f"{k:10s} registered {v['registered']}  seeds {v['by_seed']}  -> {v['reading']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
