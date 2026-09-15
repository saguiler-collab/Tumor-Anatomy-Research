#!/usr/bin/env python3
"""
darmanis_constraint_check.py — C1's direction across FOUR patients, stratified by sorting gate.

WHY THIS DESIGN AND NOT THE OBVIOUS ONE
---------------------------------------
Darmanis et al. (2017), GSE84465, dissected four glioblastomas into **tumour core** and
**periphery** and annotated 3,589 cells. That reads like a direct replacement for the
composition test in `scripts/albiach_constraint_check.py`, whose crippling limit is n = 1
patient. It is not, and the reason matters:

**3,589 of those cells are mostly FACS-panned.** Only **665 are `Unpanned`**; the rest were
sorted on CD45, HEPACAM, Thy1, GC or BSC. A composition computed over sorted cells measures
the sort, not the tissue. And among the unpanned cells the **periphery is 13 cells** — twelve
from `BT_S4`, one from `BT_S6`, none at all from `BT_S1` or `BT_S2`. So the unbiased-composition
route gives one patient's worth of periphery, which is no better than Albiach and far thinner.

**What IS available is a stratified test of the same claim.** Within a fixed sorting gate,
compare the fraction of cells that are **neoplastic** in the core against the periphery. The
gate is held constant, so the selection bias is held constant, and what remains is the thing
C1 asserts: tumour cells are denser in cellular tumour than at the margin. Four patients, and
hundreds of cells per gate.

WHAT THIS CAN AND CANNOT TEST
-----------------------------
  * **C1 (Tumor: CT > LE) — YES**, as a direction, within gate, across four patients. This is
    the only cross-patient test of any constraint this project has.
  * **C2 (Oligodendrocyte: LE > CT) — NO.** It needs oligodendrocyte *abundance* per region,
    which needs unbiased composition, which the 13 unpanned periphery cells cannot give.
  * **C3, C4, C6 — NO.** No microvascular-proliferation region exists in this dissection.
  * **C5, C7 — NO.** No peri-necrotic region, and only two regions rather than three.

So this is one constraint, tested one way, on four patients. Narrow and real, and it is
reported as such rather than dressed up as a composition replication.

EXPLORATORY. The constraint file is frozen and registered and does not move on this.

    python scripts/darmanis_constraint_check.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap.anatomic import constraints as K                       # noqa: E402

META = ROOT / "data/raw/darmanis_2017/cell_metadata.csv"

#: Darmanis's `tissue` onto Ivy GAP structures. Only two regions exist.
TISSUE_TO_STRUCTURE = {
    "Tumor": ("CT", "STRONG", "Darmanis's 'Tumor' is the resected tumour core, which is what "
                              "Ivy GAP's cellular tumour is."),
    "Periphery": ("LE", "STRONG", "Darmanis dissected the periphery specifically to capture "
                                  "the migrating front — the same region Ivy GAP calls the "
                                  "leading edge, and the paper's stated purpose."),
}

#: A gate needs this many cells in BOTH regions before it is scored. Below it the fraction is
#: too noisy to mean anything, and reporting it would invite reading noise as evidence.
MIN_PER_CELL = 30


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=20000)
    ap.add_argument("--out", default="results/darmanis_constraint_check.json")
    args = ap.parse_args()

    if not META.exists():
        print(f"BLOCKED: {META} not found. Fetch GSE84465 first "
              f"(see docs/EXTERNAL_ACTIONS.md item 15).")
        return 2

    m = pd.read_csv(META)
    m = m.rename(columns={"patient id": "patient", "cell type": "cell_type"})
    print(f"{len(m)} cells, {m['patient'].nunique()} patients, "
          f"{m['tissue'].nunique()} regions, {m['selection'].nunique()} sorting gates")

    unp = m[m["selection"] == "Unpanned"]
    per_unp = unp[unp["tissue"] == "Periphery"]
    print(f"\nunpanned cells: {len(unp)}; unpanned PERIPHERY: {len(per_unp)} "
          f"from {per_unp['patient'].nunique()} patient(s)")
    print("  -> the unbiased-composition route is NOT viable; using the stratified test")

    rows = []
    for gate, g in m.groupby("selection"):
        t = g[g["tissue"] == "Tumor"]
        p = g[g["tissue"] == "Periphery"]
        if len(t) < MIN_PER_CELL or len(p) < MIN_PER_CELL:
            rows.append({"gate": gate, "n_tumour": int(len(t)), "n_periphery": int(len(p)),
                         "scored": False,
                         "why": f"fewer than {MIN_PER_CELL} cells in one region"})
            continue
        ft = float((t["cell_type"] == "Neoplastic").mean())
        fp = float((p["cell_type"] == "Neoplastic").mean())

        # Permutation over the region label WITHIN this gate, which is the null that the
        # region carries no information once the gate is fixed.
        lab = (g["cell_type"] == "Neoplastic").to_numpy()
        is_t = (g["tissue"] == "Tumor").to_numpy()
        obs = ft - fp
        rng = np.random.default_rng(0)
        hits = 0
        for _ in range(args.permutations):
            perm = rng.permutation(is_t)
            if (lab[perm].mean() - lab[~perm].mean()) >= obs:
                hits += 1
        pval = (hits + 1) / (args.permutations + 1)

        # Per patient, so a single patient cannot carry the gate.
        per_pat = {}
        for pat, gp in g.groupby("patient"):
            tt = gp[gp["tissue"] == "Tumor"]; pp = gp[gp["tissue"] == "Periphery"]
            if len(tt) >= 10 and len(pp) >= 10:
                per_pat[pat] = {
                    "n_tumour": int(len(tt)), "n_periphery": int(len(pp)),
                    "neoplastic_tumour": round(float((tt["cell_type"] == "Neoplastic").mean()), 4),
                    "neoplastic_periphery": round(float((pp["cell_type"] == "Neoplastic").mean()), 4),
                    "direction_supports_C1": bool(
                        (tt["cell_type"] == "Neoplastic").mean() > (pp["cell_type"] == "Neoplastic").mean()),
                }
        rows.append({
            "gate": gate, "scored": True,
            "n_tumour": int(len(t)), "n_periphery": int(len(p)),
            "neoplastic_fraction_tumour": round(ft, 4),
            "neoplastic_fraction_periphery": round(fp, 4),
            "difference": round(obs, 4),
            "null_p_one_sided": round(pval, 5),
            "supports_C1": bool(obs > 0),
            "per_patient": per_pat,
            "n_patients_scored": len(per_pat),
            "n_patients_supporting": sum(1 for v in per_pat.values()
                                         if v["direction_supports_C1"]),
        })

    scored = [r for r in rows if r.get("scored")]
    print(f"\n=== C1 (Tumor: CT > LE) within sorting gate, {len(scored)} gates scored ===\n")
    print(f"{'gate':24s} {'core':>7s} {'periph':>7s} {'diff':>8s} {'p':>8s}  patients")
    print("-" * 70)
    for r in sorted(scored, key=lambda x: -x["difference"]):
        print(f"{r['gate']:24s} {r['neoplastic_fraction_tumour']:7.3f} "
              f"{r['neoplastic_fraction_periphery']:7.3f} {r['difference']:+8.3f} "
              f"{r['null_p_one_sided']:8.4f}  "
              f"{r['n_patients_supporting']}/{r['n_patients_scored']} support")
    for r in rows:
        if not r.get("scored"):
            print(f"{r['gate']:24s} not scored: {r['why']} "
                  f"(core {r['n_tumour']}, periphery {r['n_periphery']})")

    support = [r for r in scored if r["supports_C1"]]
    sig = [r for r in scored if r["supports_C1"] and r["null_p_one_sided"] < 0.05]
    report = {
        "what_this_is": ("C1's direction tested across FOUR patients by comparing the "
                         "neoplastic fraction within a fixed FACS sorting gate, core versus "
                         "periphery. The only cross-patient test of any constraint in this "
                         "project. EXPLORATORY; the constraint file is frozen and does not "
                         "move on this."),
        "why_not_composition": (
            f"Only {len(unp)} of {len(m)} cells are unpanned, and the unpanned PERIPHERY is "
            f"{len(per_unp)} cells from {per_unp['patient'].nunique()} patient(s). A "
            "composition over sorted cells measures the sort, not the tissue, so the "
            "unbiased-composition route is not viable and the stratified test is used "
            "instead."),
        "cannot_test": {
            "C2": "needs oligodendrocyte abundance per region, which needs unbiased "
                  "composition",
            "C3/C4/C6": "no microvascular-proliferation region in this dissection",
            "C5": "no peri-necrotic region",
            "C7": "needs three ordered regions; this dissection has two",
        },
        "constraint_freeze_hash": K.freeze_hash(),
        "n_cells": int(len(m)), "n_patients": int(m["patient"].nunique()),
        "min_cells_per_region_to_score": MIN_PER_CELL,
        "n_permutations": args.permutations,
        "n_gates_scored": len(scored),
        "n_gates_supporting_C1": len(support),
        "n_gates_supporting_at_p_lt_0.05": len(sig),
        "per_gate": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\n{len(support)} of {len(scored)} gates support C1's direction; "
          f"{len(sig)} at p < 0.05")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
