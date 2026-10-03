"""
extract_gimicc_cpgs.py -- the input step of prespecified/gimicc_truth_confirmation.md, as a script.

Until 2026-10-03 13:35 these two steps existed only as one-off commands in the session transcript
(found when the data inventory asked which script writes the files every GIMiCC run reads). They are
reproduced here unchanged in logic, with paths from ivygap.config, and `--verify` shows the committed
inputs are byte-identical to what this script writes.

  1. The union of every GIMiCC library's CpGs (all layers, all four tumour-type Layer-0 libraries),
     from ExperimentHub EH9483, sorted -> results/gimicc/gimicc_library_cpgs.txt
  2. Those rows of the Xena 450K matrices (D07 GBM, D08 LGG), copied verbatim with the header ->
     results/gimicc/tcga_{gbm,lgg}_gimicc_cpgs.tsv, plus results/gimicc/extraction_summary.json
     (rows, samples, CpGs found, CpGs with any missing value, missing cells). Nothing is imputed
     here or downstream: R/run_gimicc.R drops incomplete CpGs.

    python3 scripts/extract_gimicc_cpgs.py            # (re)write the inputs
    python3 scripts/extract_gimicc_cpgs.py --verify   # write to a temp dir, compare sha256 to the inputs
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

G = config.RESULTS_DIR / "gimicc"
NA = ("", "NA", "NaN", "nan")


def library_cpgs(out: Path) -> int:
    code = ('suppressPackageStartupMessages({library(ExperimentHub)});'
            'lib <- query(ExperimentHub(), "GIMiCC")[["EH9483"]];'
            'cpgs <- sort(unique(unlist(lapply(lib, rownames))));'
            f'writeLines(cpgs, "{out}"); cat(length(cpgs), "\\n")')
    r = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True)
    if r.returncode != 0 or not out.exists():
        raise SystemExit(f"BLOCKED: could not read the GIMiCC library from ExperimentHub: {r.stderr[-600:]}")
    return int(r.stdout.split()[-1])


def extract(keep: set[str], out_dir: Path) -> dict:
    summary = {}
    for h, p in config.TCGA_HM450.items():
        op = gzip.open if p.suffix == ".gz" else open
        n_rows = n_kept = n_na_cells = n_rows_with_na = 0
        with op(p, "rt") as fh, open(out_dir / f"tcga_{h}_gimicc_cpgs.tsv", "w") as fo:
            header = fh.readline(); fo.write(header)
            n_samples = len(header.rstrip("\n").split("\t")) - 1
            for line in fh:
                n_rows += 1
                if line.split("\t", 1)[0] in keep:
                    n_kept += 1
                    fo.write(line)
                    na = sum(1 for v in line.rstrip("\n").split("\t")[1:] if v in NA)
                    n_na_cells += na; n_rows_with_na += na > 0
        summary[h] = {"source_rows": n_rows, "samples": n_samples, "library_cpgs_found": n_kept,
                      "library_cpgs_total": len(keep), "cpgs_with_any_missing": n_rows_with_na,
                      "missing_cells": n_na_cells}
        print(h, summary[h], flush=True)
    (out_dir / "extraction_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    if a.verify:
        with tempfile.TemporaryDirectory(prefix="gimicc_extract_") as tmp:
            tmp = Path(tmp)
            n = library_cpgs(tmp / "gimicc_library_cpgs.txt")
            extract(set((tmp / "gimicc_library_cpgs.txt").read_text().split()), tmp)
            names = ["gimicc_library_cpgs.txt", "tcga_gbm_gimicc_cpgs.tsv", "tcga_lgg_gimicc_cpgs.tsv",
                     "extraction_summary.json"]
            res = {f: (sha(tmp / f), sha(G / f) if (G / f).exists() else None) for f in names}
        ok = all(x == y for x, y in res.values())
        for f, (x, y) in res.items():
            print(f"{'IDENTICAL' if x == y else 'DIFFERS  '} {f}  {x[:16]}  vs  {(y or 'missing')[:16]}")
        print(f"library CpGs: {n}; verify {'PASSED' if ok else 'FAILED'}")
        return 0 if ok else 1
    G.mkdir(parents=True, exist_ok=True)
    n = library_cpgs(G / "gimicc_library_cpgs.txt")
    print("library CpGs:", n)
    extract(set((G / "gimicc_library_cpgs.txt").read_text().split()), G)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
