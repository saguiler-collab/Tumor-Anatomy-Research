#!/usr/bin/env python3
"""
build_ivygap_counts_matrix.py — assemble the per-sample RSEM files into one counts matrix,
and prove it is the source of the published FPKM matrix before writing anything.

Why
---
`docs/CORRECTIONS_REGISTRATION.md` C9 (a retraction) established that Ivy GAP does publish read
counts, per sample, as RSEM `genes.results`. This turns those 270 files into the matrix the
count-based methods need — CDSeq above all, which is the one instrument that can test C4's
"is the ordering a property of the reference atlas?" without using a reference at all.

THE GATE, and it is the point of this script
--------------------------------------------
Two id spaces are in play and both are numeric, so a position-based join silently appears to
work: it shares 7,741 of 25,873 ids by coincidence and yields `corr = 0.12` against the
published matrix. Ivy GAP's internal `gene_id` is NOT RSEM's Entrez id; `rows-genes.csv` is the
bridge. This is the same alignment trap that once put 23 of 24 marker genes in the wrong
cell-type column in this project's reference builder.

So the mapping is not trusted, it is TESTED, on every sample:

  * `corr(published FPKM, this file's FPKM)` must be 1.0 to 1e-9;
  * `published / raw` must be ONE CONSTANT across genes within a sample, since the authors'
    normalisation is a per-sample scale factor.

A sample failing either is reported and excluded; the matrix is written only if every sample
passes. That check also recovers the authors' per-sample normalisation factors as a by-product,
and confirms the normalisation is a no-op for this pipeline, which rescales each bulk column to
sum to 1e6 regardless.

    python scripts/build_ivygap_counts_matrix.py

Writes `data/processed/ivygap_expected_counts.csv.gz` (genes x samples, Ivy GAP `gene_id`
space, matching the published matrix exactly) plus a provenance json.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

SRC = config.RAW_DIR / "ivygap_counts"
OUT = config.PROCESSED_DIR / "ivygap_expected_counts.csv.gz"
CORR_MIN = 1.0 - 1e-9
FACTOR_CV_MAX = 1e-6


def main() -> int:
    files = sorted(SRC.glob("*.genes.results"))
    if not files:
        print(f"BLOCKED: no files in {SRC}. Run scripts/fetch_ivygap_counts.py first.")
        return 2
    genes = pd.read_csv(config.IVYGAP_GENES_PATH)
    e2i = dict(zip(genes["gene_entrez_id"].astype(str), genes["gene_id"].astype(str)))
    pub = pd.read_csv(config.IVYGAP_FPKM_PATH, index_col=0)
    pub.index = pub.index.astype(str)
    print(f"{len(files)} RSEM files; published matrix {pub.shape[0]} genes x {pub.shape[1]} "
          f"samples")

    cols, lib, factors, failures, effl = {}, {}, {}, [], {}
    for f in files:
        well = f.name.split(".")[0]
        r = pd.read_csv(f, sep="\t")
        r["ivy"] = r["gene_id"].astype(str).map(e2i)
        n_unmapped = int(r["ivy"].isna().sum())
        r = r.dropna(subset=["ivy"]).set_index("ivy")
        if r.index.duplicated().any():
            r = r[~r.index.duplicated(keep="first")]

        if well not in pub.columns:
            failures.append((well, "not a column of the published matrix"))
            continue
        shared = pub.index.intersection(r.index)
        if len(shared) != pub.shape[0]:
            failures.append((well, f"only {len(shared)} of {pub.shape[0]} genes mapped "
                                   f"({n_unmapped} entrez ids unmapped)"))
            continue
        p = pub.loc[shared, well].to_numpy(float)
        raw = r.loc[shared, "FPKM"].to_numpy(float)
        m = (raw > 1.0) & (p > 0)          # above 2-dp rounding noise
        if m.sum() < 1000:
            failures.append((well, f"only {m.sum()} genes above the rounding floor"))
            continue
        c = float(np.corrcoef(p[m], raw[m])[0, 1])
        ratio = p[m] / raw[m]
        cv = float(np.std(ratio) / np.mean(ratio))
        if c < CORR_MIN or cv > FACTOR_CV_MAX:
            failures.append((well, f"corr {c:.8f}, factor CV {cv:.2e} — MAPPING IS WRONG"))
            continue

        ec = r.loc[shared, "expected_count"].to_numpy(float)
        cols[well] = ec
        lib[well] = float(ec.sum())
        factors[well] = float(np.median(ratio))
        effl[well] = float(np.median(r.loc[shared, "effective_length"]))

    if failures:
        print(f"\n{len(failures)} SAMPLE(S) FAILED THE GATE:")
        for w, why in failures[:10]:
            print(f"  {w}: {why}")
        print("\nBLOCKED: nothing written. The gate exists because a wrong gene mapping "
              "produces a full-looking matrix.")
        return 2

    mat = pd.DataFrame(cols, index=pub.index.copy())
    mat = mat[[c for c in pub.columns if c in mat.columns]]     # published column order
    mat.index.name = "gene_id"
    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt") as fh:
        mat.to_csv(fh, float_format="%.4f")

    fac = pd.Series(factors)
    libs = pd.Series(lib)
    print(f"\nALL {mat.shape[1]} samples passed. corr = 1.0 and a single scale factor each.")
    print(f"  matrix: {mat.shape[0]:,} genes x {mat.shape[1]} samples -> "
          f"{OUT.relative_to(ROOT)}")
    print(f"  library size (sum of expected_count): median {libs.median():,.0f}  "
          f"min {libs.min():,.0f}  max {libs.max():,.0f}")
    print(f"  authors' per-sample normalisation factor: median {fac.median():.4f}  "
          f"range {fac.min():.4f}-{fac.max():.4f}")
    print(f"  integer-valued as delivered: {bool(np.all(mat.to_numpy() == np.round(mat.to_numpy())))} "
          f"(RSEM expected counts are posterior expectations, so fractional)")

    (SRC / "matrix_provenance.json").write_text(json.dumps({
        "what_this_is": "Ivy GAP gene-level RSEM expected counts, assembled from the portal's "
                        "per-sample files into the published matrix's gene and column space.",
        "output": str(OUT.relative_to(ROOT)),
        "n_genes": int(mat.shape[0]), "n_samples": int(mat.shape[1]),
        "gate": {
            "what_it_checks": "that RSEM Entrez ids were mapped onto Ivy GAP gene_id correctly, "
                              "by requiring corr(published FPKM, raw FPKM) = 1 and a single "
                              "per-sample scale factor between them",
            "why": "both id spaces are numeric, so a position-based join shares 7,741 of "
                   "25,873 ids by coincidence and gives corr 0.12 while producing a "
                   "full-looking matrix",
            "corr_min_required": CORR_MIN, "factor_cv_max_allowed": FACTOR_CV_MAX,
            "samples_passed": int(mat.shape[1]), "samples_failed": 0,
        },
        "expected_count_is_fractional": True,
        "library_size": {"median": float(libs.median()), "min": float(libs.min()),
                         "max": float(libs.max())},
        "authors_normalisation_factor_recovered": {
            "median": float(fac.median()), "min": float(fac.min()), "max": float(fac.max()),
            "note": "published matrix = raw FPKM * this factor, per sample. It is a no-op for "
                    "this pipeline, which rescales every bulk column to sum to 1e6.",
            "per_sample": {k: round(v, 6) for k, v in factors.items()},
        },
        "median_effective_length_per_sample": {k: round(v, 2) for k, v in effl.items()},
    }, indent=2))
    print(f"  wrote {(SRC / 'matrix_provenance.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
