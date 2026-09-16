#!/usr/bin/env python3
"""
cdseq_matched_comparison.py — is CDSeq's score remarkable, or is its constraint subset easy?

CDSeq is reference-free and its estimated cell types mapped onto only part of the roster, so it
was scored on fewer constraint-tumour pairs than the leaderboard methods. Comparing those
numbers directly would be meaningless, and it would flatter whichever arm happened to draw the
easier constraints — D8 already established that the constraints are nowhere near equally
discriminating (C3 is unanimous across methods and separates nothing; C6 carries the ranking).

So every leaderboard method is re-scored on **exactly the cell types CDSeq covered**, by setting
the others to NaN and letting the scorer drop those pairs — the same partial-coverage path
quanTIseq goes through. The constraint file is untouched.

Two outcomes, and they mean opposite things:

  * every method scores near 1.0 on this subset  -> the subset is easy, CDSeq's score is
    unremarkable, and the comparison says nothing about reference-free deconvolution;
  * the methods spread out and CDSeq sits high    -> a method that never saw an atlas recovers
    the anatomy as well as methods built on one, which is direct evidence for C4/C13.

    python scripts/cdseq_matched_comparison.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

EST = config.RESULTS_DIR / "estimates"
CD = config.RESULTS_DIR / "cdseq"


def main() -> int:
    rep = json.loads((config.RESULTS_DIR / "cdseq_anatomic.json").read_text())
    labellings = rep["labellings"]
    covered = sorted(set.intersection(*[set(v["roster_types_covered"])
                                        for v in labellings.values()]))
    print(f"cell types CDSeq covered under BOTH labellings: {covered}")

    _, meta = load_cached()
    expr, _ = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    man = meta.loc[anat]

    def mask_and_score(est: pd.DataFrame, tag: str) -> dict | None:
        e = est.copy()
        e.index = e.index.astype(str)
        for c in config.CELL_TYPES:
            if c not in e.columns:
                e[c] = np.nan
        e = e[list(config.CELL_TYPES)]
        for c in config.CELL_TYPES:
            if c not in covered:
                e[c] = np.nan
        keep = [s for s in e.index if s in man.index]
        if len(keep) < 50:
            print(f"  {tag}: only {len(keep)} samples match the manifest; skipped")
            return None
        r = acs_score(e.loc[keep], man, method=tag, n_permutations=10000, n_boot=2000)
        s = r.summary()
        return {"acs": round(float(s["acs"]), 4),
                "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
                "null_p": float(s["null_p"]),
                "n_pairs": int(s["n_constraint_tumor_pairs"]),
                "beats_null": bool(s["beats_null"])}

    out: dict[str, dict] = {}
    print(f"\nre-scoring every leaderboard method on {covered}:")
    for f in sorted(EST.glob("ivygap_*.csv")):
        tag = f.stem.replace("ivygap_", "")
        if tag.endswith("_genuine"):
            continue
        est = pd.read_csv(f, index_col=0)
        r = mask_and_score(est, tag)
        if r:
            out[tag] = r
            print(f"  {tag:28s} ACS {r['acs']:.4f}  CI [{r['ci'][0]:.4f}, {r['ci'][1]:.4f}]  "
                  f"p={r['null_p']:.4f}  pairs={r['n_pairs']}")

    # CDSeq is RE-SCORED here, not copied from its own report. Its profile labelling covered
    # Tumor and the marker labelling did not, so its published ACS is over 57 pairs while every
    # row above is over 41. Copying that number across would be the exact apples-to-oranges
    # this script exists to prevent.
    print("\nCDSeq, reference-free, RE-SCORED on the same subset:")
    for name in labellings:
        f = CD / f"estimates_{name}.csv"
        if not f.exists():
            print(f"  cdseq [{name}]: {f.name} missing; skipped"); continue
        r = mask_and_score(pd.read_csv(f, index_col=0), f"cdseq__{name}")
        if r:
            out[f"cdseq__{name}"] = r
            print(f"  cdseq [{name:26s}] ACS {r['acs']:.4f}  CI "
                  f"[{r['ci'][0]:.4f}, {r['ci'][1]:.4f}]  p={r['null_p']:.4f}  "
                  f"pairs={r['n_pairs']}")

    npairs = {v["n_pairs"] for k, v in out.items() if k != "quantiseq"}
    if len(npairs) != 1:
        print(f"\nWARNING: rows are not on the same denominator: {sorted(npairs)} pairs. "
              f"The comparison below is not valid.")
    real = {k: v["acs"] for k, v in out.items() if not k.startswith("control_")}
    ctrl = {k: v["acs"] for k, v in out.items() if k.startswith("control_")}
    spread = max(real.values()) - min(real.values()) if real else float("nan")
    print(f"\nspread across real methods on this subset: {spread:.4f} "
          f"({min(real, key=real.get)} {min(real.values()):.4f} -> "
          f"{max(real, key=real.get)} {max(real.values()):.4f})")
    print(f"controls: {ctrl}")
    verdict = ("the subset does NOT discriminate — every real method is near the top, so "
               "CDSeq's score is unremarkable and says little about reference-free "
               "deconvolution" if spread < 0.10 else
               "the subset DOES discriminate, so CDSeq's position on it is informative")
    print(f"\nVERDICT: {verdict}")

    (config.RESULTS_DIR / "cdseq_matched_comparison.json").write_text(json.dumps({
        "covered_cell_types": covered,
        "why": "CDSeq's estimated types mapped onto part of the roster, so it was scored on "
               "fewer pairs; every method is re-scored on the same subset to make the "
               "comparison mean anything.",
        "spread_across_real_methods": round(float(spread), 4),
        "verdict": verdict, "methods": out}, indent=2))
    print("wrote results/cdseq_matched_comparison.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
