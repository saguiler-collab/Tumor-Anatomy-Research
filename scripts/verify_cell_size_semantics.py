#!/usr/bin/env python3
"""
verify_cell_size_semantics.py — does each package already divide out cell size?

WHY THIS EXISTS
---------------
`OPEN_DEFECTS.md` D1 and `ROAD_TO_PAPER.md` Tier 0.1 both turn on a five-row table:
MuSiC returns 0.500, Bisque 0.500, EPIC offers both, BayesPrism 0.750, SCDC 0.653. That
table decides which packages get this project's cell-size factors passed *into* them and
which keep the central conversion — a blocking design decision for the whole leaderboard.

Until this script existed, those five numbers lived **only as text typed into two
markdown files**. The probe that produced them was run once and never saved, and no
artefact recorded it. That breaks the project's own rule — *measure and record, never
validate against hard-coded literals* — and it means neither a reviewer nor this project
could check the claim that justifies the fix.

THE PROBE
---------
Two cell types, one carrying 3x the mRNA of the other, mixed **50/50 by cell count**:

    true CELL fraction of the big type   = 0.500
    true mRNA fraction of the big type   = 0.750     (3 / (3 + 1))

Those two numbers are far apart and there is nothing in between them to confuse, so what
a package returns identifies what it thinks it is reporting:

    returns ~0.500  ->  the package already divided out cell size.
                        This project's central conversion would be the SECOND one.
    returns ~0.750  ->  the package reports mRNA share.
                        This project's central conversion is the FIRST, and correct.

Every method runs with `apply_cell_size_correction=False`, so what is measured is the
package's own behaviour and nothing of this project's.

WHAT IT IS NOT
--------------
Not a benchmark, and not a ranking. A package is not better for returning either number —
they are different, equally defensible conventions. What is indefensible is applying a
conversion twice, which is what this measures.

    python scripts/verify_cell_size_semantics.py                  # all available methods
    python scripts/verify_cell_size_semantics.py --methods music,scdc
    python scripts/verify_cell_size_semantics.py --no-r           # Python path only
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap.data.reference import build_reference                  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                  # noqa: E402

#: The mRNA-content ratio between the two types. Chosen so the cell answer (0.500) and
#: the mRNA answer (0.750) are far enough apart that no plausible solver error spans them.
SIZE_RATIO = 3.0
BIG, SMALL = "Tumor", "T_cell"          # names from the real roster, so R scripts map them
TRUE_CELL_FRACTION = 0.5
TRUE_MRNA_FRACTION = SIZE_RATIO / (SIZE_RATIO + 1.0)

#: How close to an answer counts as that answer. The midpoint is 0.625; this is well
#: inside either side of it, and a value that lands between the bands is reported as
#: BETWEEN rather than forced into one — SCDC genuinely lands there.
TOL = 0.08


def build_probe(n_genes: int = 240, n_donors: int = 6, cells_per_type_per_donor: int = 30,
                n_samples: int = 12, seed: int = 0):
    """
    A reference on the real 8-type roster whose BIG type carries 3x the mRNA per cell,
    and bulk mixed 50/50 by CELL count between BIG and SMALL only.

    The full roster is used rather than a two-type toy for two reasons: `build_reference`
    builds against `config.CELL_TYPES` by construction, and the R scripts map cell-type
    names, so a toy roster would exercise a path production never takes. Every type gets
    real cells, so no method meets a zero-variance column it was never designed for.

    The bulk contains only BIG and SMALL. The other six are absent from the mixture, so
    any mass a method assigns them is leakage, which is measured and reported rather than
    hidden by the renormalisation.

    The bulk is built by literally summing single-cell count vectors, so the truth is
    exact by construction rather than asserted: 30 cells of each type in means the cell
    fraction IS 0.5, and the mRNA fraction IS whatever those cells carry.
    """
    from ivygap import config

    rng = np.random.default_rng(seed)
    roster = list(config.CELL_TYPES)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    per = n_genes // len(roster)

    # Each type owns its own block of marker genes, so the mixture is easy to solve. A
    # method that lands between 0.500 and 0.750 is telling us about its cell-size
    # convention, not about its ability to deconvolve.
    shape = {}
    for i, ctype in enumerate(roster):
        v = np.full(n_genes, 0.5)
        v[i * per:(i + 1) * per] = 10.0
        shape[ctype] = v

    # Only BIG differs in size; every other type matches SMALL, so the single ratio under
    # test is BIG:SMALL and nothing else can confound it.
    total = {ctype: 1000.0 for ctype in roster}
    total[BIG] = 1000.0 * SIZE_RATIO

    cols, meta_rows, cells = [], [], {}
    for d in range(n_donors):
        for ctype in roster:
            for c in range(cells_per_type_per_donor):
                cid = f"D{d:02d}_{ctype}_{c:03d}"
                p = shape[ctype] * rng.lognormal(0.0, 0.15, n_genes)
                cells[cid] = p / p.sum() * total[ctype]
                cols.append(cid)
                meta_rows.append({"donor": f"D{d:02d}", "cell_type": ctype})

    expr = pd.DataFrame(cells, index=genes)[cols]          # RAW per-cell mRNA
    meta = pd.DataFrame(meta_rows, index=cols)

    # Reproduce `build_from_h5ad` exactly, because the answer depends on it. Production
    # normalises every cell to 1e6 before averaging into the profile and captures the raw
    # library sizes separately as the cell-size factors. That makes every profile column
    # sum to the same total, so a least-squares solve against it returns mRNA share and
    # the central conversion is the first one.
    #
    # An earlier version of this probe passed the RAW cells to build_reference. That
    # smuggled cell size into the profile columns, and plain NNLS duly came out at 0.500 —
    # "already cell share" — which is an artefact of the probe, not a property of NNLS.
    # Getting this wrong here would have produced a false table with more authority than
    # the prose one it was written to replace.
    raw_totals = expr.sum(axis=0)
    norm = expr.div(raw_totals.replace(0.0, np.nan), axis=1).mul(1e6).fillna(0.0)
    ref = build_reference(norm, meta, name="cellsize_probe", cell_totals=raw_totals)

    # --- bulk: equal CELL counts of BIG and SMALL, summed --------------------
    big_ids = [c for c in cols if f"_{BIG}_" in c]
    small_ids = [c for c in cols if f"_{SMALL}_" in c]
    k = 30
    bulk_cols, truth = {}, []
    for s in range(n_samples):
        a = list(rng.choice(big_ids, k, replace=False))
        b = list(rng.choice(small_ids, k, replace=False))
        v = expr[a].sum(axis=1) + expr[b].sum(axis=1)
        bulk_cols[f"S{s:02d}"] = v / v.sum() * 1e6
        ma, mb = expr[a].to_numpy().sum(), expr[b].to_numpy().sum()
        truth.append({"sample": f"S{s:02d}",
                      "cell_fraction_big": k / (2 * k),
                      "mrna_fraction_big": float(ma / (ma + mb))})
    bulk = pd.DataFrame(bulk_cols, index=genes)
    manifest = pd.DataFrame({"patient_id": [f"P{i}" for i in range(n_samples)],
                             "structure": ["CT"] * n_samples}, index=bulk.columns)
    return ref, expr, meta, bulk, manifest, pd.DataFrame(truth)


def classify(value: float) -> tuple[str, str]:
    """Which convention the returned number is consistent with."""
    if abs(value - TRUE_CELL_FRACTION) <= TOL:
        return ("CELL share",
                "the package already divided out cell size; this project's central "
                "conversion would be the SECOND one — a defect")
    if abs(value - TRUE_MRNA_FRACTION) <= TOL:
        return ("mRNA share",
                "the package reports mRNA share; this project's central conversion is "
                "the FIRST and is correct")
    return ("BETWEEN",
            f"lands between the cell answer ({TRUE_CELL_FRACTION:.3f}) and the mRNA "
            f"answer ({TRUE_MRNA_FRACTION:.3f}); the package does partial size handling, "
            "so neither applying nor skipping the central conversion is right for it")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--methods", default="",
                    help="comma-separated subset; default is every available method")
    ap.add_argument("--no-r", action="store_true", help="force the all-Python path")
    ap.add_argument("--out", default="results/cell_size_semantics.json")
    ap.add_argument("--budget", type=int, default=1800, help="seconds per R method")
    args = ap.parse_args()

    print("building the probe: two types, "
          f"{SIZE_RATIO:g}x mRNA ratio, mixed 50/50 by cell count")
    ref, expr, meta, bulk, manifest, truth = build_probe()
    cs = ref.cell_size
    print(f"  reference cell_size: {BIG}={cs[BIG]:.1f}  {SMALL}={cs[SMALL]:.1f}  "
          f"ratio {cs[BIG] / cs[SMALL]:.3f} (target {SIZE_RATIO:g})")
    print(f"  bulk truth: cell fraction {truth['cell_fraction_big'].mean():.4f}, "
          f"mRNA fraction {truth['mrna_fraction_big'].mean():.4f}")

    data = DeconvolutionInput(
        bulk=bulk, references=(ref,), manifest=manifest,
        apply_cell_size_correction=False,      # the whole point: measure the package alone
    )

    from ivygap.deconv import r_bridge
    from ivygap.deconv.registry import build_methods
    # The same normalised cells production exports, not the raw ones -- otherwise the R
    # packages would build their own basis from a different scaling than the Python ones.
    raw_totals = expr.sum(axis=0)
    norm = expr.div(raw_totals.replace(0.0, np.nan), axis=1).mul(1e6).fillna(0.0)
    r_bridge.set_cell_source(ref.name, norm, meta)

    wanted = {m.strip() for m in args.methods.split(",") if m.strip()}
    methods = [m for m in build_methods(prefer_r=not args.no_r)
               if not wanted or m.name in wanted]

    rows = []
    for m in methods:
        t0 = time.perf_counter()
        try:
            est = m.fit_predict(data)
            # Read the ratio BETWEEN the two types that are actually in the mixture. Any
            # mass on the other six is leakage: it is reported, never renormalised away,
            # because a method that puts 30% on absent types is not answering the
            # question this probe asks and must not be read as if it were.
            big, small = est[BIG].to_numpy(), est[SMALL].to_numpy()
            pair = big + small
            leak = float(1.0 - np.mean(pair))
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = np.where(pair > 0, big / np.maximum(pair, 1e-12), np.nan)
            val = float(np.nanmean(ratio))
            conv, meaning = classify(val)
            if leak > 0.25:
                conv, meaning = "UNREADABLE", (
                    f"assigned {leak:.1%} of its mass to the six cell types absent from "
                    "the mixture, so the BIG:SMALL ratio does not measure a cell-size "
                    "convention here")
            impl = getattr(m, "last_implementation", None) or (
                "R" if getattr(m, "requires_r", False) else "python")
            rows.append({"method": m.name, "returned_big_fraction": round(val, 4),
                         "sd_across_samples": round(float(np.nanstd(ratio)), 4),
                         "leakage_to_absent_types": round(leak, 4),
                         "convention": conv, "meaning": meaning,
                         "implementation": str(impl),
                         "elapsed_s": round(time.perf_counter() - t0, 1),
                         "failed": None})
            print(f"  {m.name:24s} {val:.4f}  {conv:11s}  "
                  f"leak {leak:+.3f}  ({impl})")
        except Exception as exc:                                   # noqa: BLE001
            rows.append({"method": m.name, "returned_big_fraction": None,
                         "convention": "FAILED", "meaning": None,
                         "implementation": None,
                         "elapsed_s": round(time.perf_counter() - t0, 1),
                         "failed": f"{type(exc).__name__}: {str(exc)[:300]}"})
            print(f"  {m.name:24s} FAILED: {type(exc).__name__}: {str(exc)[:120]}")

    report = {
        "what_this_is": ("Measures, per package, whether it already converts mRNA share to "
                         "CELL share. Replaces five numbers that previously existed only "
                         "as prose in OPEN_DEFECTS.md and ROAD_TO_PAPER.md."),
        "probe": {
            "size_ratio": SIZE_RATIO,
            "big_type": BIG, "small_type": SMALL,
            "true_cell_fraction_big": TRUE_CELL_FRACTION,
            "true_mrna_fraction_big": round(float(truth["mrna_fraction_big"].mean()), 4),
            "measured_cell_size_ratio": round(float(cs[BIG] / cs[SMALL]), 4),
            "n_samples": int(bulk.shape[1]), "n_genes": int(bulk.shape[0]),
            "tolerance": TOL,
            "central_correction_applied": False,
        },
        "methods": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {out}")

    n_fail = sum(1 for r in rows if r["convention"] == "FAILED")
    print(f"{len(rows) - n_fail} measured, {n_fail} failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
