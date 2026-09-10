#!/usr/bin/env python3
"""
ish_constraint_check.py — check the pre-registered constraints against Ivy GAP's own ISH.

WHAT THIS ASKS
--------------
The seven constraints in `ivygap/anatomic/constraints.py` are ordinal claims about which
cell types are enriched in which anatomic structures, justified by neuropathology. Ivy GAP
independently quantified in-situ hybridization signal for 480 genes across the same five
structures, in the same tumours. So the claims can be checked against a measurement that
involves no deconvolution and no cell-type model at all.

If ISH agrees, a constraint rests on measurement rather than on expert expectation. If it
disagrees, that is a finding about the constraint and is reported as one.

WHAT IT IS NOT
--------------
**Exploratory. The constraints do not move on this result.** They are frozen, hashed
(`2d1fb47c…`) and registered at <https://osf.io/dm2t8>. This belongs under the
registration's "Other planned analysis" and must be labelled that way wherever reported.

**Not a measurement of composition.** ISH expression energy is transcript signal in a
tissue neighbourhood. A marker higher in one structure is consistent with more cells of
that type, more transcript per cell, or both. It supports an ordinal claim; it cannot
supply a fraction, and it is not a substitute for deconvolution.

WHY THIS REPLACED THE FIRST ATTEMPT
-----------------------------------
An earlier version pulled per-experiment values from the Allen API and was wrong in three
ways that all pointed the same direction — toward a confident answer built on very little:

  * it used 72 of 632 experiments (11.4%), whichever the API returned first;
  * it pooled expression energy across DIFFERENT sub-block specimens, conflating
    tumour-to-tumour variation with structure differences — the exact error `acs.py`
    avoids by collapsing to (tumour, structure) before aggregating;
  * it pooled markers that contradict each other, then reported the pooled verdict.

The correct input was already on disk. `gene_expression_details.csv` carries one row per
(gene x sub-block) with **all five structures as columns on that row**, so a comparison
between two structures is paired inside one physical block by construction.

THE DESIGN, MIRRORING acs.py
----------------------------
1. A constraint is evaluated in a sub-block only where BOTH structures it names have a
   value. Missing structures are excluded from numerator and denominator alike — a
   sub-block is not penalised for anatomy it does not contain.
2. Sub-blocks nest within donors, so satisfaction is averaged within donor first and
   across donors second. A donor contributing many sub-blocks does not dominate.
3. Every marker is reported separately. Markers are never pooled, because when ESM1 and
   CD34 disagree that disagreement is the result.
4. Significance comes from permuting structure labels WITHIN each sub-block, which
   preserves that block's overall expression level and destroys only the assignment.

CIRCULARITY
-----------
ISH is how Ivy GAP labelled the 148 "cluster" RNA-seq samples this project excludes from
scoring. Using ISH here is legitimate because the constraints are scored on the **H&E**
samples, whose structure labels come from histology. Different samples, different
labelling basis. Any write-up must say so, or a reader will assume it was overlooked.

    python scripts/ish_constraint_check.py
    python scripts/ish_constraint_check.py --permutations 10000
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                        # noqa: E402

SOURCE = config.RAW_DIR / "ivygap" / "gene_expression_details.csv"

#: ISH column -> this project's constraint vocabulary.
STRUCT = {
    "expression_energy_le":    "LE",
    "expression_energy_it":    "IT",
    "expression_energy_ct":    "CT",
    "expression_energy_ctmvp": "MVP",
    "expression_energy_ctpan": "PAN",
}

#: Markers, declared before any value is read, and chosen ONLY for being canonical for
#: their cell type and present in Ivy GAP's panel. Every one is reported whether it
#: agrees or not; none is dropped for being inconvenient.
MARKERS: dict[str, list[str]] = {
    "Endothelial":          ["ESM1", "CD34", "KDR", "CAV1"],
    "Macrophage_Microglia": ["CD163", "LAPTM5"],
    "Tumor":                ["SOX2", "PTPRZ1", "EGFR", "CD44", "BIRC5", "TOP2A"],
}

#: Constraints this panel can speak to, and how.
CHECKS = [
    ("C1", "Tumor",                "CT",  "LE",  "pairwise"),
    ("C3", "Endothelial",          "MVP", "CT",  "pairwise"),
    ("C4", "Endothelial",          "MVP", None,  "maximum"),
    ("C5", "Macrophage_Microglia", "PAN", "LE",  "pairwise"),
    ("C6", "Macrophage_Microglia", "MVP", "CT",  "pairwise"),
    ("C7", "Tumor",                None,  None,  "monotone"),   # LE < IT < CT
]

#: Declared before running, so the silence is visible rather than mistaken for agreement.
NOT_CHECKABLE = {
    "C2": ("Oligodendrocyte: LE > CT cannot be checked. Ivy GAP's ISH panel contains no "
           "myelin or oligodendrocyte-lineage marker — MBP, PLP1, MOG, MAG, CNP, SOX10, "
           "MOBP and CLDN11 are all absent. OLIG2 IS present and is deliberately not "
           "used: in glioma it is expressed by the tumour cells themselves and by OPCs, "
           "so scoring an oligodendrocyte claim on it would measure tumour content and "
           "report it as oligodendrocyte content."),
}


def load() -> pd.DataFrame:
    d = pd.read_csv(SOURCE)
    keep = ["donor_id", "tumor_name", "sub_block_id", "gene_symbol", *STRUCT]
    d = d[keep].rename(columns=STRUCT)
    return d


def satisfied_mask(sub: pd.DataFrame, kind: str, hi: str | None, lo: str | None):
    """
    Per sub-block: True satisfied, False violated, NaN not evaluable.

    A sub-block that lacks either structure the claim names is NOT evaluable — never a
    violation. This is the same rule the ACS scorer applies to tumours.
    """
    if kind == "pairwise":
        ok = sub[hi].notna() & sub[lo].notna()
        out = pd.Series(np.nan, index=sub.index, dtype="float64")
        out[ok] = (sub.loc[ok, hi] > sub.loc[ok, lo]).astype(float)
        return out
    if kind == "maximum":
        others = [s for s in STRUCT.values() if s != hi]
        ok = sub[hi].notna() & sub[others].notna().any(axis=1)
        out = pd.Series(np.nan, index=sub.index, dtype="float64")
        vals = sub.loc[ok, others].max(axis=1)
        out[ok] = (sub.loc[ok, hi] > vals).astype(float)
        return out
    # monotone LE < IT < CT
    ok = sub[["LE", "IT", "CT"]].notna().all(axis=1)
    out = pd.Series(np.nan, index=sub.index, dtype="float64")
    out[ok] = ((sub.loc[ok, "LE"] < sub.loc[ok, "IT"]) &
               (sub.loc[ok, "IT"] < sub.loc[ok, "CT"])).astype(float)
    return out


def donor_equal_rate(sub: pd.DataFrame, mask: pd.Series) -> tuple[float | None, int, int]:
    """
    Average within donor, then across donors.

    Sub-blocks nest in donors. Pooling them would let a donor contributing many blocks
    dominate, which is the nesting rule the protocol states for tumours.
    """
    ok = mask.notna()
    if not ok.any():
        return None, 0, 0
    df = pd.DataFrame({"donor": sub.loc[ok, "donor_id"], "sat": mask[ok]})
    per_donor = df.groupby("donor")["sat"].mean()
    return float(per_donor.mean()), int(ok.sum()), int(per_donor.size)


def permutation_p(sub: pd.DataFrame, kind: str, hi: str | None, lo: str | None,
                  observed: float, n_perm: int, rng) -> float:
    """
    Null: shuffle the structure labels WITHIN each sub-block.

    That preserves the block's own expression level and its set of measured structures,
    and destroys only which value belongs to which structure — the same construction the
    ACS null uses on tumours.

    Written in numpy rather than pandas. The obvious version rebuilds a DataFrame and
    runs a groupby on every draw, which for CD44 (613 sub-blocks) never finishes.
    """
    cols = list(STRUCT.values())
    idx = {c: i for i, c in enumerate(cols)}
    vals = sub[cols].to_numpy(dtype="float64")
    donors = pd.factorize(sub["donor_id"])[0]
    n_donor = donors.max() + 1

    def rate_of(v: np.ndarray) -> float | None:
        if kind == "pairwise":
            a, b = v[:, idx[hi]], v[:, idx[lo]]
            ok = ~np.isnan(a) & ~np.isnan(b)
            sat = np.where(ok, a > b, np.nan)
        elif kind == "maximum":
            others = [idx[s] for s in cols if s != hi]
            a = v[:, idx[hi]]
            rest = v[:, others]
            ok = ~np.isnan(a) & (~np.isnan(rest)).any(axis=1)
            with np.errstate(invalid="ignore"):
                mx = np.nanmax(np.where(np.isnan(rest), -np.inf, rest), axis=1)
            sat = np.where(ok, a > mx, np.nan)
        else:
            le, it, ct = v[:, idx["LE"]], v[:, idx["IT"]], v[:, idx["CT"]]
            ok = ~np.isnan(le) & ~np.isnan(it) & ~np.isnan(ct)
            sat = np.where(ok, (le < it) & (it < ct), np.nan)
        if not np.any(ok):
            return None
        # donor-equal: mean within donor, then across donors that contributed
        sums = np.zeros(n_donor); cnts = np.zeros(n_donor)
        np.add.at(sums, donors[ok], sat[ok])
        np.add.at(cnts, donors[ok], 1.0)
        live = cnts > 0
        return float((sums[live] / cnts[live]).mean())

    hits = 0
    present = [np.flatnonzero(~np.isnan(vals[i])) for i in range(vals.shape[0])]
    for _ in range(n_perm):
        v = vals.copy()
        for i, pres in enumerate(present):
            if pres.size > 1:
                v[i, pres] = v[i, rng.permutation(pres)]
        r = rate_of(v)
        if r is not None and r >= observed:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--out", default="results/ish_constraint_check.json")
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED)
    args = ap.parse_args()

    if not SOURCE.exists():
        print(f"missing {SOURCE}")
        return 1
    d = load()
    rng = np.random.default_rng(args.seed)

    try:
        from ivygap.anatomic import constraints as K
        freeze = K.freeze_hash()
    except Exception:
        freeze = None

    report = {
        "what_this_is": ("Ivy GAP ISH expression energy checked against the "
                         "pre-registered constraints. EXPLORATORY — the constraints are "
                         "frozen and registered and do not move on this result."),
        "source": str(SOURCE.relative_to(config.PROJECT_ROOT)),
        "constraint_freeze_hash": freeze,
        "n_rows": int(len(d)), "n_genes": int(d.gene_symbol.nunique()),
        "n_sub_blocks": int(d.sub_block_id.nunique()),
        "n_donors": int(d.donor_id.nunique()),
        "markers_declared_before_reading_values": MARKERS,
        "design": ("paired within sub-block; averaged within donor then across donors; "
                   "per marker, never pooled; null permutes structure labels within "
                   "sub-block"),
        "not_checkable": NOT_CHECKABLE,
        "per_marker": [],
    }

    print(f"{len(d):,} rows | {d.gene_symbol.nunique()} genes | "
          f"{d.sub_block_id.nunique()} sub-blocks | {d.donor_id.nunique()} donors\n")
    hdr = f"{'constraint':>4}  {'marker':<8} {'rate':>6} {'blocks':>7} {'donors':>7} {'p':>8}"
    print(hdr); print("-" * len(hdr))

    for cid, cell_type, hi, lo, kind in CHECKS:
        for gene in MARKERS[cell_type]:
            sub = d[d.gene_symbol == gene]
            if sub.empty:
                continue
            mask = satisfied_mask(sub, kind, hi, lo)
            rate, n_blocks, n_donors = donor_equal_rate(sub, mask)
            row = {"constraint": cid, "cell_type": cell_type, "marker": gene,
                   "claim": (f"{cell_type}: {hi} > {lo}" if kind == "pairwise"
                             else f"{cell_type}: {hi} is the maximum" if kind == "maximum"
                             else f"{cell_type}: LE < IT < CT"),
                   "donor_equal_rate": None if rate is None else round(rate, 4),
                   "n_sub_blocks_evaluable": n_blocks, "n_donors": n_donors}
            if rate is not None and n_blocks >= 3:
                row["null_p"] = round(permutation_p(sub, kind, hi, lo, rate,
                                                    args.permutations, rng), 5)
            report["per_marker"].append(row)
            r = "  n/a " if rate is None else f"{rate:6.3f}"
            p = row.get("null_p")
            print(f"{cid:>4}  {gene:<8} {r} {n_blocks:>7} {n_donors:>7} "
                  f"{('%8.4f' % p) if p is not None else '       -'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")
    print("\nEXPLORATORY. Markers are reported separately and never pooled; where markers")
    print("for one cell type disagree, that disagreement is the result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
