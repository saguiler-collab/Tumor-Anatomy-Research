#!/usr/bin/env python3
"""
albiach_constraint_check.py — the constraints against MEASURED single-cell composition.

WHAT THIS IS
------------
Mossa Albiach et al. (2023) published a glioblastoma resection dissected into 27 samples
across 12 anatomically labelled locations, with per-cell type annotations: 135,482 cells,
15 cell types, and a `Zone` label of High fluorescence / Low fluorescence / Periphery /
Necrotic core. The file has been sitting in `pipeline packages/ repos/SCDC/albiach data/`
unused.

It is the only asset available to this project that gives **measured cell-type composition
per anatomic region of a glioblastoma**. Not inferred, not deconvolved — counted from
annotated single cells. So it can test the pre-registered constraints with no deconvolution
method anywhere in the chain, which is a stronger form of the same argument
`scripts/ish_constraint_check.py` makes from ISH energy.

WHY IT IS EXPLORATORY AND CANNOT MOVE THE CONSTRAINTS
-----------------------------------------------------
The constraint file is frozen and registered. Nothing here may edit it. Three further
limits, all structural:

1. **One donor.** Every one of the 135,482 cells is from `SL040`. The 27 samples give
   within-tumour replication and no cross-patient inference at all. This is n = 1 patient.
2. **The zone mapping is an interpretation, not an identity.** Ivy GAP's structures were
   assigned by a neuropathologist on H&E; Albiach's zones come from fluorescence-guided
   resection. The mapping used here is declared in `ZONE_TO_STRUCTURE` below and argued in
   the docstring beside it. One mapping is weak and is reported separately: **Albiach's
   "Necrotic core" is NOT Ivy GAP's PAN.** PAN is the hypercellular pseudopalisading rim
   *around* necrosis; the necrotic core is the dying centre. A constraint about PAN tested on
   the necrotic core is testing a different region.
3. **Single-cell composition is not tissue composition.** Dissociation under-represents
   tumour cells, which lyse more readily than leukocytes, so absolute fractions here are
   immune-inflated — macrophages run 0.33 to 0.65, far above any bulk estimate. This is why
   only ORDINAL comparisons of the SAME cell type ACROSS regions are scored: a bias that
   multiplies one cell type by a constant everywhere cannot reverse that type's ranking
   between two regions. Absolute fractions are reported for context and scored on nothing.

    python scripts/albiach_constraint_check.py
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

H5AD = ROOT / "pipeline packages / repos/SCDC/albiach data/d45b4ce6-9725-4d79-b97a-70a44158bdbf.h5ad"

#: Albiach's 15 labels onto this project's 8-type roster. Declared from the label names
#: alone. Types with no roster counterpart are EXCLUDED rather than forced into the nearest
#: column -- OPCs, pericytes, VSMC and fibroblasts have no roster home, and putting OPCs into
#: Oligodendrocyte would inflate exactly the type C2 is about.
ROSTER = {
    "Tumor": "Tumor",
    "Microglia": "Macrophage_Microglia", "MgTAM": "Macrophage_Microglia",
    "MoTAM": "Macrophage_Microglia", "Macrophages": "Macrophage_Microglia",
    "Monocytes": "Macrophage_Microglia",
    "Oligodendrocytes": "Oligodendrocyte",
    "Endothelial": "Endothelial",
    "Astrocytes": "Astrocyte",
    "B-cells": "B_cell",
    "NK, T-cells": "T_cell",
}

#: Albiach Zone -> Ivy GAP structure. Each entry carries the argument for it and a
#: confidence, because the constraints are scored on this mapping and a reader must be able
#: to reject a specific row rather than the whole analysis.
ZONE_TO_STRUCTURE = {
    "Periphery": ("LE", "STRONG",
                  "Ivy GAP's leading edge is the outermost region, mostly brain infiltrated "
                  "by few tumour cells. Albiach's periphery measures Tumor at 0.034 and "
                  "Oligodendrocyte at 0.166 — the highest oligodendrocyte and lowest tumour "
                  "of any zone, which is what LE is defined by."),
    "High": ("CT", "STRONG",
             "Fluorescence-high tissue is the dense tumour bulk; Albiach measures Tumor at "
             "0.536 there, the highest of the non-necrotic zones."),
    "Low": ("IT", "MODERATE",
            "Fluorescence-low sits between the bulk and the periphery, which is what "
            "infiltrating tumour is, and its Tumor fraction (0.419) sits between them. But "
            "low fluorescence is an optical property, not a histological one."),
    "Necrotic core": ("PAN", "WEAK — do not score without this caveat",
                      "Ivy GAP's PAN is the HYPERCELLULAR pseudopalisading rim AROUND "
                      "necrosis. Albiach's necrotic core is the dying centre. These are "
                      "adjacent and different, and the pseudopalisading rim is precisely "
                      "where macrophages accumulate. A PAN constraint tested here is tested "
                      "on the wrong side of the boundary."),
}


def _cat(f, col: str) -> pd.Series:
    import h5py
    g = f["obs"][col]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]]
        return pd.Series([cats[i] if i >= 0 else "NA" for i in g["codes"][:]])
    return pd.Series([x.decode() if isinstance(x, bytes) else str(x) for x in g[:]])


def load_composition() -> tuple[pd.DataFrame, dict]:
    """Per-SAMPLE composition over the roster, with the sample's zone. No expression read."""
    import h5py
    with h5py.File(H5AD, "r") as f:
        obs = pd.DataFrame({c: _cat(f, c)
                            for c in ("CellType", "Zone", "Location", "sample_id", "donor_id")})
    prov = {
        "file": str(H5AD.relative_to(ROOT)),
        "n_cells_in_file": int(len(obs)),
        "n_donors": int(obs["donor_id"].nunique()),
        "donors": sorted(obs["donor_id"].unique().tolist()),
        "n_samples": int(obs["sample_id"].nunique()),
        "n_locations": int(obs["Location"].nunique()),
        "labels_excluded_no_roster_home": sorted(
            set(obs.loc[~obs["CellType"].isin(ROSTER), "CellType"].unique())),
    }
    obs["roster"] = obs["CellType"].map(ROSTER)
    prov["n_cells_mapped"] = int(obs["roster"].notna().sum())
    m = obs[obs["roster"].notna()]

    # Per sample, so the unit of analysis is a dissected piece of tissue rather than a cell.
    # Pooling cells would let one deeply sequenced sample dominate its zone, the same nesting
    # rule the protocol states for tumours.
    counts = pd.crosstab(m["sample_id"], m["roster"])
    frac = counts.div(counts.sum(axis=1), axis=0)
    zone = m.groupby("sample_id")["Zone"].agg(lambda s: s.mode().iat[0])
    frac["Zone"] = zone
    frac["n_cells"] = counts.sum(axis=1)
    return frac, prov


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=20000)
    ap.add_argument("--out", default="results/albiach_constraint_check.json")
    args = ap.parse_args()

    if not H5AD.exists():
        print(f"BLOCKED: {H5AD} not found.")
        return 2

    print("reading annotations only (no expression matrix) ...")
    frac, prov = load_composition()
    print(f"  {prov['n_cells_mapped']:,} of {prov['n_cells_in_file']:,} cells mapped onto the "
          f"roster; {prov['n_samples']} samples, {prov['n_donors']} donor")
    print(f"  excluded, no roster home: {prov['labels_excluded_no_roster_home']}")

    struct_of = {z: v[0] for z, v in ZONE_TO_STRUCTURE.items()}
    frac["structure"] = frac["Zone"].map(struct_of)

    rng = np.random.default_rng(0)
    rows = []
    for c in K.CONSTRAINTS:
        ct = c.cell_type
        if ct not in frac.columns:
            rows.append({"constraint": c.id, "cell_type": ct, "verdict": "NOT TESTABLE",
                         "why": f"{ct} has no Albiach label mapped onto it"})
            continue

        if c.kind == "pairwise":
            hi, lo = c.structures
            a = frac.loc[frac["structure"] == hi, ct]
            b = frac.loc[frac["structure"] == lo, ct]
            if len(a) == 0 or len(b) == 0:
                rows.append({"constraint": c.id, "cell_type": ct, "verdict": "NOT TESTABLE",
                             "why": f"no Albiach zone maps to {hi if len(a)==0 else lo}"})
                continue
            obs_diff = float(a.mean() - b.mean())
            pool = np.concatenate([a.to_numpy(), b.to_numpy()])
            n_a = len(a)
            hits = 0
            for _ in range(args.permutations):
                p = rng.permutation(pool)
                if (p[:n_a].mean() - p[n_a:].mean()) >= obs_diff:
                    hits += 1
            pval = (hits + 1) / (args.permutations + 1)
            conf = [v[1] for z, v in ZONE_TO_STRUCTURE.items() if v[0] in (hi, lo)]
            rows.append({
                "constraint": c.id, "cell_type": ct, "claim": c.describe(),
                "kind": "pairwise", "hi": hi, "lo": lo,
                "mean_hi": round(float(a.mean()), 4), "mean_lo": round(float(b.mean()), 4),
                "n_samples_hi": int(n_a), "n_samples_lo": int(len(b)),
                "satisfied": bool(obs_diff > 0),
                "fold": round(float(a.mean() / b.mean()), 2) if b.mean() > 0 else None,
                "null_p_one_sided": round(pval, 5),
                "mapping_confidence": conf,
                "verdict": ("SATISFIED" if obs_diff > 0 else "VIOLATED"),
            })
        elif c.kind == "monotone":
            seq = c.structures
            means = [frac.loc[frac["structure"] == s, ct].mean() for s in seq]
            if any(pd.isna(m) for m in means):
                rows.append({"constraint": c.id, "cell_type": ct, "verdict": "NOT TESTABLE",
                             "why": f"not every structure in {seq} has an Albiach zone"})
                continue
            ok = all(means[i] < means[i + 1] for i in range(len(means) - 1))
            rows.append({"constraint": c.id, "cell_type": ct, "claim": c.describe(),
                         "kind": "monotone", "sequence": list(seq),
                         "means": [round(float(m), 4) for m in means],
                         "satisfied": bool(ok),
                         "verdict": "SATISFIED" if ok else "VIOLATED"})
        else:                                     # maximum
            target = c.structures[0]
            if target not in set(struct_of.values()):
                rows.append({"constraint": c.id, "cell_type": ct, "verdict": "NOT TESTABLE",
                             "why": f"no Albiach zone maps to {target}"})
                continue
            means = {s: frac.loc[frac["structure"] == s, ct].mean()
                     for s in set(frac["structure"].dropna())}
            ok = all(means[target] > v for s, v in means.items() if s != target)
            rows.append({"constraint": c.id, "cell_type": ct, "claim": c.describe(),
                         "kind": "maximum", "target": target,
                         "means": {k: round(float(v), 4) for k, v in means.items()},
                         "satisfied": bool(ok),
                         "verdict": "SATISFIED" if ok else "VIOLATED"})

    print("\n=== the pre-registered constraints against MEASURED composition ===\n")
    for r in rows:
        v = r["verdict"]
        extra = ""
        if v in ("SATISFIED", "VIOLATED") and r.get("kind") == "pairwise":
            extra = (f"  {r['hi']} {r['mean_hi']:.4f} vs {r['lo']} {r['mean_lo']:.4f}"
                     f"  p={r['null_p_one_sided']:.4f}")
        elif r.get("kind") == "monotone":
            extra = "  " + " < ".join(f"{s}:{m:.4f}" for s, m in zip(r["sequence"], r["means"]))
        print(f"  {r['constraint']}  {v:13s} {r.get('claim', r.get('why', ''))[:52]:52s}{extra}")

    testable = [r for r in rows if r["verdict"] in ("SATISFIED", "VIOLATED")]
    sat = [r for r in testable if r["verdict"] == "SATISFIED"]
    report = {
        "what_this_is": ("The pre-registered anatomic constraints checked against MEASURED "
                         "single-cell composition per anatomic region (Mossa Albiach et al. "
                         "2023). No deconvolution method is involved. EXPLORATORY: the "
                         "constraint file is frozen and registered and does not move on this."),
        "limitations": [
            f"ONE donor ({prov['donors']}). {prov['n_samples']} samples give within-tumour "
            "replication and no cross-patient inference. n = 1 patient.",
            "The zone-to-structure mapping is an interpretation; see zone_mapping, which "
            "carries a confidence and an argument per row.",
            "Albiach's 'Necrotic core' is NOT Ivy GAP's PAN — PAN is the hypercellular "
            "pseudopalisading rim AROUND necrosis, and macrophages accumulate precisely "
            "there. Any PAN result here is tested on the wrong side of that boundary.",
            "Single-cell dissociation under-represents tumour cells, so absolute fractions "
            "are immune-inflated (macrophages 0.33-0.65). Only ordinal comparisons of the "
            "same cell type across regions are scored, which a per-type constant bias "
            "cannot reverse.",
        ],
        "provenance": prov,
        "roster_mapping": ROSTER,
        "zone_mapping": {z: {"structure": v[0], "confidence": v[1], "argument": v[2]}
                         for z, v in ZONE_TO_STRUCTURE.items()},
        "constraint_freeze_hash": K.freeze_hash(),
        "n_testable": len(testable), "n_satisfied": len(sat),
        "n_permutations": args.permutations,
        # dict.fromkeys, not a list comprehension: ROSTER maps five Albiach labels onto
        # Macrophage_Microglia, so the values are not unique and the duplicate column names
        # made to_json(orient="index") fail.
        "composition_by_zone": json.loads(
            frac.groupby("Zone")[[c for c in dict.fromkeys(ROSTER.values())
                                  if c in frac.columns]]
            .mean().round(4).to_json(orient="index")),
        "per_constraint": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\n{len(sat)} of {len(testable)} testable constraints SATISFIED")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
