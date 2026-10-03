"""
recide_tcga.py -- ReCIDE on TCGA's raw/X arm: prespecified/recide_tcga.md.

POST-REGISTRATION EXTENSION PANEL. ReCIDE runs through `extension_tcga.run` UNCHANGED -- the same
bulk, gene space, samples, purity, metrics and per-sample file as every other raw/X-arm method --
by handing `run` an adapter in place of `build_extension_methods`. The adapter builds ReCIDE's
own reference exactly as its anatomic run does (`run_recide.py`: the pre-fixed donor draw,
uncapped cells, counts checked integral) on the arm's gene space, and calls R/run_recide.R.

Controls, each aborting on failure:
  1. the arm reproduces an existing row: deseq2_unmix in tcga_{cohort}_h5ad.json, every field;
  2. the donor draw reproduces the anatomic run's donors.

ReCIDE's raw proportions are saved the moment R returns (results/extension/recide_tcga_{cohort}_raw.csv),
so a failure afterwards cannot cost the hours of compute before it.

    python3 scripts/recide_tcga.py --cohorts gbm [lgg] [--n-cores 1]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
from ivygap.deconv import r_bridge  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput, finalize_estimates  # noqa: E402
import extension_tcga  # noqa: E402
import run_recide  # noqa: E402

OUT = config.RESULTS_DIR / "extension"
ARM_MATRIX = "raw/X"          # extension_tcga's "h5ad" arm builds the atlas from raw/X


class ReCIDEAdapter:
    """What `extension_tcga.run` needs from a method: a name, an implementation label, fit_predict."""

    name = "recide"
    implementation_ = run_recide.IMPLEMENTATION

    def __init__(self, cohort: str, donors: list[str], n_cores: int, budget: int | None):
        self.cohort, self.donors, self.n_cores, self.budget = cohort, donors, n_cores, budget
        self.report: dict = {}

    def fit_predict(self, data: DeconvolutionInput):
        genes = list(data.bulk.index)
        # ReCIDE never uses the bridge's cell source; free it before an hours-long R call.
        r_bridge.clear_cell_source(config.PRIMARY_REFERENCE)
        ref, cells, cmeta, counts = run_recide.recide_reference(self.donors, genes, ARM_MATRIX)
        print(f"    ReCIDE reference: {counts.shape[1]:,} cells from {len(self.donors)} donors, "
              f"{len(genes)} genes; {data.bulk.shape[1]} samples", flush=True)
        est, elapsed, tail, failed = run_recide.call_recide(counts, genes, cells, cmeta, data.bulk,
                                                            self.n_cores, self.budget)
        self.report = {"elapsed_seconds": round(elapsed, 1), "n_reference_cells": int(counts.shape[1]),
                       "r_stdout_tail": tail[-1500:], "failed": failed}
        if failed:
            raise RuntimeError(failed)
        est.to_csv(OUT / f"recide_tcga_{self.cohort}_raw.csv")          # keep the expensive part
        self.report["types_dropped_by_recide"] = [c for c in config.CELL_TYPES if c not in est.columns]
        est = est.reindex(index=[str(c) for c in data.bulk.columns], columns=list(config.CELL_TYPES))
        d2 = DeconvolutionInput(bulk=data.bulk, references=(ref.subset_genes(genes),),
                                manifest=data.manifest)
        return finalize_estimates(est.to_numpy(), d2, covered=None, returns_cell_fractions=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohorts", nargs="+", default=["gbm"], choices=["gbm", "lgg"])
    ap.add_argument("--n-cores", type=int, default=1)
    ap.add_argument("--budget", type=int, default=None, help="seconds; default: no limit")
    a = ap.parse_args()

    prov = json.loads((config.RESULTS_DIR / "run_provenance.json").read_text())["matrix"]
    if prov != ARM_MATRIX:
        print(f"BLOCKED: the leaderboard's layer is {prov!r}, the arm's is {ARM_MATRIX!r}"); return 2
    donors, k, n_elig = run_recide.donor_draw()
    rec = json.loads((OUT / "recide_anatomic.json").read_text())["donor_draw"]
    if not (donors == rec["donors"] and k == rec["k"] and n_elig == rec["n_eligible"]):
        print("BLOCKED (control 2): the donor draw does not reproduce the anatomic run's"); return 2
    print(f"control 2: donor draw reproduces the anatomic run's {k} donors", flush=True)

    for cohort in a.cohorts:
        path = OUT / f"tcga_{cohort}_h5ad.json"
        stored = json.loads(path.read_text())["methods"]["deseq2_unmix"]
        got = extension_tcga.run(cohort, "h5ad", ["deseq2_unmix"])["methods"]["deseq2_unmix"]
        diffs = {k2: (v, got.get(k2)) for k2, v in stored.items() if got.get(k2) != v}
        if diffs:
            print(f"BLOCKED (control 1, {cohort}): the arm does not reproduce deseq2_unmix's row: {diffs}")
            return 2
        print(f"control 1 ({cohort}): the arm reproduces deseq2_unmix's row, every field", flush=True)

        adapter = ReCIDEAdapter(cohort, donors, a.n_cores, a.budget)
        original = extension_tcga.build_extension_methods
        extension_tcga.build_extension_methods = lambda names: [adapter]
        try:
            res = extension_tcga.run(cohort, "h5ad", ["recide"])
        finally:
            extension_tcga.build_extension_methods = original
        row = res["methods"]["recide"]
        row.update({k2: v for k2, v in adapter.report.items() if k2 != "failed" or v})
        row["donor_draw"] = {"rule": run_recide.DRAW_RULE, "k": k, "n_eligible": n_elig, "donors": donors}
        row["controls"] = {"arm_reproduces_deseq2_unmix_row": True, "donor_draw_reproduces_anatomic": True}
        summary = json.loads(path.read_text())
        summary.setdefault("methods", {})["recide"] = row
        path.write_text(json.dumps(summary, indent=2))
        print(f"  wrote {path.relative_to(config.PROJECT_ROOT)} (recide)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
