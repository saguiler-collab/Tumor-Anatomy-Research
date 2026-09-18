"""Extract the EpiDISH reference CpGs from a Xena HM450 matrix.

WHY THIS EXISTS AS A SCRIPT. The LGG probe file was produced by an ad-hoc step, which makes the
GBM arm unreproducible and the LGG arm unverifiable. This does it the same way for any cohort,
and it records what it dropped instead of quietly dropping it.

The probe list is taken from EpiDISH itself (`rownames(centDHSbloodDMC.m)`), never typed out
here -- a hard-coded list would silently diverge if the package version changed, and
`CLAUDE.md` forbids validating against hard-coded literals.

INPUT: a Xena `TCGA.<COHORT>.sampleMap_HumanMethylation450[.gz]` matrix -- probes x samples,
tab-separated, first column the probe id.

Usage:
    python scripts/extract_epidish_probes.py --cohort gbm \
        --source TCGA_GBM/TCGA.GBM.sampleMap_HumanMethylation450.gz
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import pandas as pd

from ivygap import config


def reference_probes() -> list[str]:
    """The CpGs EpiDISH's blood reference is defined on, read from the package."""
    r = subprocess.run(
        ["Rscript", "-e",
         'suppressMessages(library(EpiDISH)); data(centDHSbloodDMC.m); '
         'cat(rownames(centDHSbloodDMC.m), sep="\\n")'],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"could not read centDHSbloodDMC.m from EpiDISH:\n{r.stderr[-800:]}")
    probes = [x.strip() for x in r.stdout.splitlines() if x.strip()]
    if not probes:
        raise RuntimeError("EpiDISH returned an empty probe list")
    return probes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=["gbm", "lgg"])
    ap.add_argument("--source", required=True, help="Xena HumanMethylation450 matrix")
    a = ap.parse_args()

    src = pathlib.Path(a.source)
    if not src.is_absolute():
        src = config.PROJECT_ROOT / src
    if not src.exists():
        print(f"BLOCKED: {src} does not exist.")
        return 2
    probes = reference_probes()
    print(f"EpiDISH reference: {len(probes)} CpGs")

    # Read only the reference probes. The full Xena matrix is ~485k probes x ~hundreds of
    # samples; loading all of it to keep 333 rows is what makes this step feel like it needs a
    # better machine, and it does not.
    keep = set(probes)
    chunks = []
    for ch in pd.read_csv(src, sep="\t", index_col=0, chunksize=50_000,
                          compression="infer", low_memory=False):
        hit = ch.loc[[i for i in ch.index if i in keep]]
        if len(hit):
            chunks.append(hit)
    if not chunks:
        print(f"BLOCKED: none of the {len(probes)} reference CpGs are in {a.source}. "
              f"Wrong platform, or probe ids in a different namespace.")
        return 2
    beta = pd.concat(chunks)
    beta = beta[~beta.index.duplicated()]

    found = len(beta)
    missing = sorted(keep - set(beta.index))
    n_all_na = int(beta.isna().all(axis=1).sum())
    complete = beta.dropna()
    print(f"  found {found} of {len(probes)} CpGs in the matrix"
          f"{'' if not missing else f'; {len(missing)} absent, e.g. {missing[:4]}'}")
    print(f"  {n_all_na} CpGs are all-NA across samples")
    print(f"  {len(complete)} CpGs complete across all {beta.shape[1]} samples "
          f"(EpiDISH uses complete cases)")
    if len(complete) < 0.8 * len(probes):
        print(f"WARNING: only {len(complete) / len(probes):.0%} of the reference survives "
              f"complete-case filtering. Deconvolution on a truncated reference is not the "
              f"published method; report this rather than proceeding silently.")

    out = config.PROCESSED_DIR / f"{a.cohort}_methylation_epidish_probes.csv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    beta.to_csv(out, compression="gzip")
    print(f"\nwrote {out.relative_to(config.PROJECT_ROOT)} "
          f"({beta.shape[0]} probes x {beta.shape[1]} samples)")
    print(f"next: python scripts/methylation_celltypes.py --cohort {a.cohort}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
