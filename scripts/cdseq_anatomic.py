#!/usr/bin/env python3
"""
cdseq_anatomic.py — the reference-free arm. Does the anatomy hold up without an atlas?

THE QUESTION THIS ANSWERS AND NOTHING ELSE DOES
-----------------------------------------------
`CORRECTIONS_REGISTRATION.md` C4 and C13: the ACS ordering is substantially a property of
GBmap, and every method in the panel is handed GBmap, so none of them can separate "the methods
recover the anatomy" from "GBmap recovers the anatomy". CDSeq estimates cell types **de novo
from the bulk alone** and is the only instrument here that can.

Blocked until 2026-09-15, because CDSeq is a multinomial read-count model and Ivy GAP publishes
FPKM. The retraction in C9 unblocked it: the Allen portal serves per-sample RSEM counts.

TWO INDEPENDENT LABELLINGS, COMPARED
------------------------------------
CDSeq returns anonymous cell types. Mapping them onto the roster needs *something*, and the
choice could smuggle the reference back in. So it is done twice:

  (a) correlation with the GBmap profile  — uses the atlas, for LABELLING ONLY
  (b) marker enrichment using the GBM-specific marker sets from GBMdeconvoluteR — independent
      of GBmap entirely

If ACS agrees under both, the result does not depend on the atlas even for labelling. If it
disagrees, that is itself the finding and is reported as one.

GENE SELECTION IS REFERENCE-FREE, which is the point
----------------------------------------------------
Genes are chosen by variance in the BULK. Using the leaderboard's 657 markers would reintroduce
exactly the dependence this arm exists to remove — those genes were selected on GBmap.

    python scripts/cdseq_anatomic.py [--genes 1500] [--types 8] [--iters 700]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

WORK = config.RESULTS_DIR / "cdseq"
COUNTS = config.PROCESSED_DIR / "ivygap_expected_counts.csv.gz"


def label_by_profile(gep: pd.DataFrame) -> tuple[dict[int, str], pd.DataFrame]:
    """Correlate each estimated GEP with the GBmap profile. LABELLING ONLY."""
    prof = pd.read_csv(ROOT / "data/reference/gbmap_linear/profile.csv", index_col=0)
    shared = gep.index.intersection(prof.index)
    if len(shared) < 100:
        return {}, pd.DataFrame()
    a = np.log1p(gep.loc[shared].to_numpy(float))
    b = np.log1p(prof.loc[shared].to_numpy(float))
    C = pd.DataFrame(np.corrcoef(a.T, b.T)[:a.shape[1], a.shape[1]:],
                     index=gep.columns, columns=prof.columns)
    return {k: C.loc[k].idxmax() for k in C.index}, C


def label_by_markers(gep: pd.DataFrame) -> tuple[dict[int, str], pd.DataFrame]:
    """
    Marker enrichment with GBM-specific markers. Uses NO atlas.

    Score = mean rank-percentile of a type's markers within an estimated GEP, so a GEP is
    assigned the roster type whose markers sit highest in it.
    """
    mk = pd.read_csv(WORK.parent / "imc_anchored" / "markers_MCP_GBM_moreno.tsv", sep="\t")
    pct = gep.rank(axis=0, pct=True)
    rows = {}
    for ctype, g in mk.groupby("Cell population"):
        genes = [x for x in g["HUGO symbols"] if x in pct.index]
        if len(genes) >= 5:
            rows[ctype] = pct.loc[genes].mean(axis=0)
    C = pd.DataFrame(rows)                                  # estimated types x roster types
    return {k: C.loc[k].idxmax() for k in C.index}, C


def to_roster(prop: pd.DataFrame, labels: dict) -> pd.DataFrame:
    """
    Estimated types -> roster columns. Several estimated types can map to one roster type;
    their proportions are SUMMED, which is what a partition of the same tissue means.
    """
    est = pd.DataFrame(0.0, index=prop.columns, columns=list(config.CELL_TYPES))
    used = set()
    for k, ctype in labels.items():
        if ctype in est.columns:
            est[ctype] = est[ctype] + prop.loc[k].to_numpy()
            used.add(ctype)
    for c in est.columns:                     # a roster type no estimated type claimed is
        if c not in used:                     # ABSENT, not zero: zero would be an assertion
            est[c] = np.nan
    est.index.name = "sample_id"
    return est


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=int, default=1500)
    ap.add_argument("--types", type=int, default=len(config.CELL_TYPES))
    ap.add_argument("--iters", type=int, default=700)
    ap.add_argument("--dilution", type=float, default=100.0)
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)

    if not COUNTS.exists():
        print(f"BLOCKED: {COUNTS} missing. Run scripts/build_ivygap_counts_matrix.py.")
        return 2
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    man = meta.loc[anat]

    counts = pd.read_csv(COUNTS, index_col=0)
    counts.columns = counts.columns.astype(str)
    counts.index = counts.index.astype(str)
    keep = [s for s in anat if s in counts.columns]
    print(f"anatomic samples with counts: {len(keep)} of {len(anat)}")
    if len(keep) != len(anat):
        print("BLOCKED: not every anatomic sample has counts."); return 2
    sub = counts[keep]

    # Ivy GAP gene_id -> HUGO, so the labelling steps can speak the same vocabulary.
    genes = pd.read_csv(config.IVYGAP_GENES_PATH)
    i2s = dict(zip(genes["gene_id"].astype(str), genes["gene_symbol"].astype(str)))
    sub.index = [i2s.get(g, g) for g in sub.index]
    sub = sub[~sub.index.duplicated(keep="first")]

    # REFERENCE-FREE gene selection: variance of log-CPM in the bulk itself.
    cpm = np.log1p(sub / sub.sum(axis=0) * 1e6)
    top = cpm.var(axis=1).sort_values(ascending=False).head(args.genes).index
    sub = sub.loc[top]
    print(f"genes: {len(sub)} chosen by bulk variance (NOT the leaderboard's markers, which "
          f"were selected on GBmap)")

    cin = WORK / "cdseq_input.csv"
    sub.to_csv(cin)
    pref = str(WORK / f"cdseq_T{args.types}")
    print(f"\nrunning CDSeq: T={args.types}, {args.iters} iterations, "
          f"dilution_factor={args.dilution:g}")
    r = subprocess.run(["Rscript", str(ROOT / "scripts/run_cdseq_anatomic.R"), str(cin),
                        str(args.types), str(args.iters), str(args.dilution), pref],
                       capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print("CDSeq FAILED:\n" + r.stderr.strip()[-1500:]); return 2

    prop = pd.read_csv(f"{pref}_prop.csv", index_col=0)
    gep = pd.read_csv(f"{pref}_gep.csv", index_col=0)
    prop.columns = prop.columns.astype(str)

    report = {"n_samples": len(keep), "n_genes": int(len(sub)), "n_types": args.types,
              "mcmc_iterations": args.iters, "dilution_factor": args.dilution,
              "gene_selection": "top variance of log-CPM in the Ivy GAP bulk; reference-free "
                                "by design, because the leaderboard's 657 markers were "
                                "selected on GBmap",
              "reference_gep_passed_to_cdseq": False, "labellings": {}}

    for name, fn in [("gbmap_profile_correlation", label_by_profile),
                     ("gbm_marker_enrichment", label_by_markers)]:
        labels, C = fn(gep)
        if not labels:
            print(f"\n[{name}] SKIPPED — too few shared genes"); continue
        est = to_roster(prop, labels)
        est.to_csv(WORK / f"estimates_{name}.csv")
        covered = sorted({v for v in labels.values() if v in config.CELL_TYPES})
        print(f"\n=== labelling: {name} ===")
        print(f"    assignments: {dict(sorted(labels.items(), key=lambda x: str(x[0])))}")
        print(f"    roster types covered: {covered}")
        res = acs_score(est, man, method=f"cdseq__{name}", n_permutations=10000, n_boot=2000)
        s = res.summary()
        print(f"    ACS {s['acs']:.4f}  CI [{s['ci_low']:.4f}, {s['ci_high']:.4f}]  "
              f"p={s['null_p']:.4f}  pairs={s['n_constraint_tumor_pairs']}  "
              f"tumours={s['n_tumors']}")
        report["labellings"][name] = {
            "assignments": {str(k): v for k, v in labels.items()},
            "roster_types_covered": covered,
            "acs": round(float(s["acs"]), 4),
            "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
            "null_p": float(s["null_p"]),
            "n_pairs": int(s["n_constraint_tumor_pairs"]),
            "beats_null": bool(s["beats_null"])}

    (config.RESULTS_DIR / "cdseq_anatomic.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote results/cdseq_anatomic.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
