#!/usr/bin/env python3
"""
reference_sensitivity.py — is the leaderboard a property of the METHODS, or of GBmap?

THE QUESTION
------------
Every number in this study rests on one single-cell atlas. Nothing so far distinguishes
"this is how these methods rank" from "this is how these methods rank **on GBmap**", and a
reviewer will ask. `docs/OPEN_DEFECTS.md` D8 sharpens it: the constraint carrying most of the
ranking sits on the roster column that absorbs the 7.0% of the atlas the roster discards.

Darmanis 2017 gives a second reference — 5 roster types, 4 donors
(`scripts/build_darmanis_reference.py`). It is the first chance to vary the reference and hold
everything else fixed.

THE DESIGN, AND WHY BOTH ARMS ARE NEEDED
----------------------------------------
Swapping GBmap for Darmanis changes two things at once: the reference AND the roster, since
Darmanis cannot resolve T_cell, NK_cell or B_cell. Comparing a 5-type Darmanis run against the
published 8-type GBmap run would confound them, and the confound is not small — dropping three
columns changes what every solver has to explain.

So both arms run on the **same declared 5-type sub-roster** and the **same gene space**:

    arm A   GBmap reference,    5 types, shared genes
    arm B   Darmanis reference, 5 types, shared genes

The only difference is which cells built the reference. Arm A is also reported against the
published 8-type leaderboard, so the cost of dropping the three immune columns is visible
separately rather than absorbed.

**All seven constraints remain scoreable.** The constraint file names only Tumor,
Macrophage_Microglia, Endothelial and Oligodendrocyte, and Darmanis has all four. So ACS is
computed on the same 7 constraints and the same 9 tumours in both arms, and the orderings are
directly comparable.

WHAT THIS IS NOT
----------------
Not a selection procedure. Neither arm may reorder the published leaderboard, which stands as
the confirmatory run recorded it. A method that does better under one reference has shown
reference sensitivity, not merit.

Python path only by default: the R packages read a cell-level export keyed to a registered cell
source, and registering Darmanis's cells is separate work. The Python solvers are enough to
answer whether the ORDERING is reference-dependent.

    python scripts/reference_sensitivity.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402
from ivygap.data.reference import ReferenceBundle, build_from_h5ad  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                  # noqa: E402
from ivygap.deconv.registry import build_methods                   # noqa: E402

#: Each alternative reference, with the sub-roster it can support. The sub-roster is the
#: reference's OWN type list -- it is not a choice, it is what the reference has.
ALTERNATIVES = {
    "darmanis": (ROOT / "data/reference/darmanis_2017",
                 ["Tumor", "Macrophage_Microglia", "Endothelial", "Oligodendrocyte",
                  "Astrocyte"]),
    "neftel": (ROOT / "data/reference/neftel_2019",
               ["Tumor", "Macrophage_Microglia", "T_cell", "Oligodendrocyte"]),
}


def load_alt(name: str, sub_roster: list[str]) -> ReferenceBundle:
    d, _ = ALTERNATIVES[name]
    prof = pd.read_csv(d / "profile.csv", index_col=0)[sub_roster]
    sig = pd.read_csv(d / "sigma.csv", index_col=0)[sub_roster]
    cs = pd.read_csv(d / "cell_size.csv", index_col=0)["cell_size"].reindex(sub_roster)
    prov = json.loads((d / "provenance.json").read_text())
    return ReferenceBundle(name=name, profile=prof, sigma=sig, cell_size=cs,
                           n_donors=len(prov["donors"]), donor_profiles=None,
                           has_cross_donor_variance=True)


def to_sub_roster(ref: ReferenceBundle, SUB_ROSTER: list[str]) -> ReferenceBundle:
    """Restrict a bundle to the 5 shared types, keeping every other field aligned."""
    return ReferenceBundle(
        name=ref.name + "_5type",
        profile=ref.profile[SUB_ROSTER], sigma=ref.sigma[SUB_ROSTER],
        cell_size=ref.cell_size.reindex(SUB_ROSTER),
        n_donors=ref.n_donors,
        donor_profiles=({d: p[SUB_ROSTER] for d, p in ref.donor_profiles.items()}
                        if ref.donor_profiles else None),
        has_cross_donor_variance=ref.has_cross_donor_variance,
    )


def score_arm(label: str, ref: ReferenceBundle, bulk: pd.DataFrame,
              manifest: pd.DataFrame, n_perm: int, n_boot: int,
              SUB_ROSTER: list[str] | None = None) -> dict:
    data = DeconvolutionInput(bulk=bulk, references=(ref,), manifest=manifest,
                              cell_types=tuple(ref.cell_types))
    out = {}
    for m in build_methods(prefer_r=False):
        if m.name.startswith("control_"):
            continue
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                   # noqa: BLE001
            out[m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:160]}"}
            print(f"    {m.name:24s} FAILED {type(exc).__name__}")
            continue
        res = acs_score(est, manifest, method=f"{m.name}__{label}",
                        n_permutations=n_perm, n_boot=n_boot)
        s = res.summary()
        out[m.name] = {"acs": round(float(s["acs"]), 4),
                       "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
                       "null_p": float(s["null_p"]),
                       "n_pairs": int(s["n_constraint_tumor_pairs"]),
                       "n_tumors": int(s["n_tumors"])}
        print(f"    {m.name:24s} ACS {s['acs']:.4f}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--reference", choices=sorted(ALTERNATIVES), default="darmanis")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    ALT_DIR, SUB_ROSTER = ALTERNATIVES[args.reference]
    out_path = Path(args.out) if args.out else Path(
        f"results/reference_sensitivity_{args.reference}.json")

    if not (ALT_DIR / "profile.csv").exists():
        print(f"BLOCKED: {ALT_DIR} not built. "
              f"Run scripts/build_{args.reference}_reference.py first.")
        return 2
    print(f"alternative reference: {args.reference}; sub-roster {SUB_ROSTER}")

    print("loading the Ivy GAP anatomic cohort ...")
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    print(f"  {len(anat)} anatomic samples")

    dar = load_alt(args.reference, SUB_ROSTER)
    print(f"  {args.reference} reference: {dar.profile.shape[0]:,} genes x {len(SUB_ROSTER)} types")

    print("\nloading the atlas for arm A (restricted to the bulk's genes) ...")
    gb_full, _, _ = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad",
                                    restrict_to_genes=expr.index, export=False)
    gb = to_sub_roster(gb_full, SUB_ROSTER)
    print(f"  GBmap reference: {gb.profile.shape[0]:,} genes x {len(SUB_ROSTER)} types")

    # ONE shared gene space, so the only difference between arms is the reference.
    #
    # The LEADERBOARD's 657-gene marker space, not the full shared transcriptome. Two
    # reasons, and the second is the honest caveat.
    #
    #   * Speed, and it is not marginal: nu-SVR is superlinear in genes, which is why the
    #     pipeline uses a marker subset at all. A first version used all 16,356 shared genes
    #     and SVR alone had not finished after five minutes.
    #   * Production match: this is the space the published leaderboard was computed on, so
    #     arm A is directly comparable to it.
    #
    # THE CAVEAT: those 657 genes were selected ON GBMAP, so the gene space is chosen by
    # arm A's own reference. That favours arm A. It is therefore CONSERVATIVE for any finding
    # that Darmanis does as well or better, and it would be the wrong space to use for a
    # finding that Darmanis does worse. Stated in the report.
    #
    # Sets built ONCE: written inline the first time, which rebuilt both on every one of
    # 16,758 iterations -- roughly 800 million hash insertions.
    genes_path = config.BENCH_DIR / "signature_genes.json"
    if not genes_path.exists():
        print(f"BLOCKED: {genes_path} not found. Run scripts/reconstruct_gene_space.py.")
        return 2
    lb_genes = json.loads(genes_path.read_text())["genes"]
    dar_genes, bulk_genes = set(dar.profile.index), set(expr.index)
    gb_genes = set(gb.profile.index)
    shared = [g for g in lb_genes if g in dar_genes and g in bulk_genes and g in gb_genes]
    print(f"  leaderboard gene space: {len(lb_genes)}; usable in BOTH references: "
          f"{len(shared)}")
    print(f"\nshared gene space: {len(shared):,} genes "
          f"(GBmap {gb.profile.shape[0]:,}, Darmanis {dar.profile.shape[0]:,})")
    if len(shared) < config.MIN_GENES_SHARED:
        print(f"BLOCKED: only {len(shared)} shared genes; {config.MIN_GENES_SHARED} required.")
        return 2

    bulk = expr.loc[shared, anat]
    bulk = bulk / bulk.sum(axis=0) * 1e6
    man = meta.loc[anat]

    print(f"\n=== arm A: GBmap reference, {len(SUB_ROSTER)}-type sub-roster ===")
    arm_a = score_arm(f"gbmap_{len(SUB_ROSTER)}", gb.subset_genes(shared), bulk, man,
                      args.permutations, args.boot, SUB_ROSTER)
    print(f"\n=== arm B: {args.reference} reference, {len(SUB_ROSTER)}-type "
          f"sub-roster ===")
    arm_b = score_arm(f"{args.reference}_{len(SUB_ROSTER)}", dar.subset_genes(shared),
                      bulk, man, args.permutations, args.boot, SUB_ROSTER)

    common = [m for m in arm_a if "acs" in arm_a[m] and m in arm_b and "acs" in arm_b[m]]
    a = pd.Series({m: arm_a[m]["acs"] for m in common})
    b = pd.Series({m: arm_b[m]["acs"] for m in common})
    rho = float(stats.spearmanr(a, b).statistic) if len(common) > 2 else float("nan")
    ra = a.rank(ascending=False, method="average")
    rb = b.rank(ascending=False, method="average")
    moved = (rb - ra).abs()

    print(f"\n=== reference sensitivity across {len(common)} methods ===\n")
    print(f"{'method':24s} {'GBmap':>8s} {'Darmanis':>9s} {'delta':>8s} {'rank move':>10s}")
    print("-" * 64)
    for m in a.sort_values(ascending=False).index:
        print(f"{m:24s} {a[m]:8.4f} {b[m]:9.4f} {b[m]-a[m]:+8.4f} "
              f"{rb[m]-ra[m]:+10.1f}")
    print(f"\nSpearman between the two ACS orderings: {rho:.4f}")
    print(f"methods whose rank moves: {int((moved > 0).sum())} of {len(common)}; "
          f"largest move {moved.max():.1f}")

    report = {
        "what_this_is": ("Does the ACS ordering depend on which single-cell atlas built the "
                         "reference? Both arms use the SAME declared 5-type sub-roster and "
                         "the SAME gene space, so the only difference is the reference. NOT a "
                         "selection procedure: neither arm may reorder the published "
                         "leaderboard."),
        "why_both_arms": ("Swapping GBmap for Darmanis changes the reference AND the roster, "
                          "since Darmanis cannot resolve T_cell, NK_cell or B_cell. Running "
                          "GBmap on the same 5 types isolates the reference."),
        "alternative_reference": args.reference,
        "sub_roster": SUB_ROSTER,
        "n_types": len(SUB_ROSTER),
        "constraints_use_only": sorted({c.cell_type for c in __import__(
            "ivygap.anatomic.constraints", fromlist=["CONSTRAINTS"]).CONSTRAINTS}),
        "n_shared_genes": len(shared),
        "gene_space": ("the leaderboard's own 657-gene marker space, intersected with both "
                       "references. Those genes were selected ON GBmap, so the space favours "
                       "arm A -- conservative for a finding that Darmanis does as well or "
                       "better, and the wrong space for a finding that it does worse."),
        "n_samples": len(anat),
        "implementation": "python path only; the R packages need a cell-level export keyed to "
                          "a registered cell source, which Darmanis does not yet have",
        "spearman_between_orderings": None if np.isnan(rho) else round(rho, 4),
        "n_methods": len(common),
        "n_ranks_moved": int((moved > 0).sum()),
        "largest_rank_move": round(float(moved.max()), 1) if len(common) else None,
        "arm_a_gbmap": arm_a,
        "arm_b_alternative": arm_b,
    }
    out = out_path
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
