#!/usr/bin/env python3
"""
mes_score.py — a continuous mesenchymal score, computed identically in GBM and LGG.

WHY THIS IS A SECONDARY ANALYSIS AND NOT THE PRE-SPECIFIED F5
--------------------------------------------------------------
F5 was pre-specified as the **Verhaak categorical call** `EXPRESSION_SUBTYPE == Mesenchymal`,
taken from cBioPortal. That classification exists for GBM and **not** for TCGA-LGG: lower-grade
glioma is classified by IDH status and 1p/19q codeletion instead. So F5 cannot be validated in
LGG with the same variable, and substituting a different variable and calling it the same test
would be exactly the kind of quiet swap this project refuses.

Instead: a **continuous MES program score** from the Neftel four-state neoplastic markers already
vendored in the GBMdeconvoluteR drop, computed the same way in both cohorts. Reported as a
**labelled secondary analysis**.

The HYPOTHESIS is unchanged — mesenchymal character should depress the tumour estimate, because
MES tumours are transcriptionally myeloid-like and malignant signal leaks into the macrophage
column. The pre-specified DIRECTION is unchanged. Only the measurement changes, and that change
is declared here rather than absorbed.

**Score definition, fixed before it is related to any error:** per sample, the mean
rank-percentile of the MES marker genes within that sample's own expression profile. Percentile
within-sample makes it independent of library size and comparable across cohorts without any
cross-cohort normalisation, which is the property that matters here.

    python scripts/mes_score.py --cohort gbm
    python scripts/mes_score.py --cohort lgg
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

GBMD = (ROOT / "pipeline_packages " / " repos" / "GBMDeconvoluteR"
        / "GBMDeconvoluteR-main" / "data")
RDS = "Neftel_et_al_2019_four_state_neoplastic_markers.rds"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["gbm", "lgg"], required=True)
    args = ap.parse_args()
    bulk_p = config.PROCESSED_DIR / f"tcga_{args.cohort}_bulk_cpm.csv.gz"
    if not bulk_p.exists():
        print(f"BLOCKED: {bulk_p} missing."); return 2

    out_tsv = config.RESULTS_DIR / f"_neftel_states_{args.cohort}.tsv"
    r = subprocess.run(["Rscript", "-e", f'''
        x <- readRDS("{GBMD / RDS}")
        write.table(data.frame(gene = unlist(x, use.names = FALSE),
                               state = rep(names(x), lengths(x))),
                    "{out_tsv}", sep = "\\t", row.names = FALSE, quote = FALSE)
    '''], capture_output=True, text=True)
    if r.returncode != 0:
        print("BLOCKED: could not read the Neftel marker object:\n" + r.stderr[-400:])
        return 2
    mk = pd.read_csv(out_tsv, sep="\t")
    print(f"Neftel four-state markers: {len(mk)} genes over "
          f"{mk['state'].nunique()} states {sorted(mk['state'].unique())}")

    with gzip.open(bulk_p, "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    pct = bulk.rank(axis=0, pct=True)          # within-sample percentile, library-size free
    scores = {}
    for state, g in mk.groupby("state"):
        genes = [x for x in g["gene"] if x in pct.index]
        if len(genes) >= 5:
            scores[state] = pct.loc[genes].mean(axis=0)
        print(f"  {state:4s} {len(genes):3d} of {len(g):3d} markers present in the bulk")
    s = pd.DataFrame(scores)
    s.index.name = "sample"
    s["MES_minus_mean_other"] = s["MES"] - s[[c for c in s.columns if c != "MES"]].mean(axis=1)
    out = config.RESULTS_DIR / f"mes_score_{args.cohort}.csv"
    s.to_csv(out)
    print(f"\n{args.cohort.upper()} MES score: median {s['MES'].median():.4f}  "
          f"range {s['MES'].min():.4f}-{s['MES'].max():.4f}")
    print(f"  MES minus mean(AC,NPC,OPC): median {s['MES_minus_mean_other'].median():+.4f}")
    print(f"wrote {out.relative_to(ROOT)}")
    out_tsv.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
