"""ReCIDE on the anatomic cohort -- POST-REGISTRATION EXTENSION PANEL.

ReCIDE builds one signature per reference DONOR and integrates the per-donor deconvolutions. Its
own code keeps only donors with more than 500 cells and cell types with more than 20 cells, so
the leaderboard's capped reference (at most 50 cells per donor and type, <= 400 per donor) cannot
feed it -- every donor would be dropped. It is therefore given UNCAPPED cells from a seeded random
draw of the leaderboard's own TRAINING donors, as its README directs for limited resources.

THE DRAW RULE WAS FIXED BEFORE ANY RESULT EXISTED (2026-10-01): among training donors with more
than 500 mapped cells, the smallest k in (10, 15, 20) for which a draw with the project seed gives
every roster type more than ReCIDE's 20-cell threshold. Coverage of the roster, never an outcome.

Everything else is the leaderboard's: atlas layer (from run_provenance.json), gene space (hash-
verified), anatomic samples and order, CPM bulk on that gene space, `finalize_estimates`, and the
ACS scorer with 10,000 permutations. Writes results/extension/recide_anatomic.json.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402
from ivygap.data.reference import build_from_h5ad                  # noqa: E402
from ivygap.deconv.base import DeconvolutionInput, finalize_estimates  # noqa: E402


IMPLEMENTATION = "R:ReCIDE vendored @31bbd6b (+ disclosed summarize_each shim; v3 assays)"
DRAW_RULE = ("smallest k in (10,15,20) covering every roster type >20 cells; "
             "eligible = training donors with >500 mapped cells; project seed")


def donor_draw() -> tuple[list[str], int | None, int]:
    """The draw rule fixed before any result existed (2026-10-01), from the atlas obs only.
    Returns (donors, k, n_eligible); donors is empty when no k covers the roster."""
    import h5py                                                    # noqa: PLC0415
    from anndata.io import read_elem                               # noqa: PLC0415
    from ivygap.data.reference import _collapse_to_roster          # noqa: PLC0415
    sel = json.loads((config.BENCH_DIR / "method_selection_decision.json").read_text())
    with h5py.File(config.REFERENCE_DIR / "gbmap_core.h5ad", "r") as f:
        obs = read_elem(f["obs"])
    ct = _collapse_to_roster(obs["annotation_level_3"].astype(str), config.GBMAP_CELL_TYPE_MAP)
    d = pd.DataFrame({"donor": obs["donor_id"].astype(str).to_numpy(), "ct": ct.to_numpy()}).dropna()
    train = set(map(str, sel["train_donors"]))
    tab = d[d.donor.isin(train)].groupby(["donor", "ct"]).size().unstack(fill_value=0)
    elig = sorted(tab.index[tab.sum(1) > 500])
    for k in (10, 15, 20):
        cand = list(np.random.default_rng(config.RANDOM_SEED).choice(elig, size=k, replace=False))
        if (tab.loc[cand].reindex(columns=config.CELL_TYPES, fill_value=0).sum() > 20).all():
            return sorted(cand), k, len(elig)
    return [], None, len(elig)


def recide_reference(pick: list[str], genes: list[str], matrix: str):
    """Uncapped cells from the drawn donors on `genes`, with counts recovered and checked integral.
    Returns (ref, cells, cmeta, counts) or raises ValueError."""
    ref, cells, cmeta = build_from_h5ad(
        config.REFERENCE_DIR / "gbmap_core.h5ad", matrix=matrix, restrict_to_genes=genes,
        cell_filter=lambda o: o["donor_id"].astype(str).isin(pick).to_numpy(),
        max_cells_per_donor_type=10**9, max_total_cells=10**9, export=False)
    counts = cells.loc[genes].to_numpy(dtype="float64") * (
        cmeta.loc[cells.columns, "library_size"].to_numpy(dtype="float64") / 1e6)
    from ivygap.deconv.extension import counts_are_integral       # noqa: PLC0415
    ok, rel = counts_are_integral(counts)
    if not ok:
        raise ValueError(f"recovered counts are not integers (max relative deviation {rel:.3g})")
    return ref, cells, cmeta, sparse.csc_matrix(np.round(counts))


def call_recide(counts, genes: list[str], cells, cmeta, bulk: pd.DataFrame, n_cores: int,
                budget: int | None) -> tuple[pd.DataFrame | None, float, str, str | None]:
    """Run R/run_recide.R. Returns (proportions or None, seconds, stdout tail, failure or None)."""
    with tempfile.TemporaryDirectory(prefix="ivygap_recide_") as tmp:
        tmp = Path(tmp)
        sio.mmwrite(tmp / "counts.mtx", counts)
        (tmp / "genes.txt").write_text("\n".join(genes) + "\n")
        pd.DataFrame({"cell": list(cells.columns), "donor": cmeta.loc[cells.columns, "donor"].astype(str).to_numpy(),
                      "cell_type": cmeta.loc[cells.columns, "cell_type"].astype(str).to_numpy()}
                     ).to_csv(tmp / "meta.csv", index=False)
        bulk.to_csv(tmp / "bulk.csv")
        cfg = {"counts": str(tmp / "counts.mtx"), "genes": str(tmp / "genes.txt"),
               "meta": str(tmp / "meta.csv"), "bulk": str(tmp / "bulk.csv"),
               "out": str(tmp / "props.csv"), "seed": config.RANDOM_SEED, "n_cores": n_cores}
        (tmp / "args.json").write_text(json.dumps(cfg))
        t0 = time.perf_counter()
        proc = subprocess.run(["Rscript", str(config.PROJECT_ROOT / "R" / "run_recide.R"),
                               str(tmp / "args.json")], capture_output=True, text=True,
                              timeout=budget)
        elapsed = time.perf_counter() - t0
        (config.DIAGNOSTICS_DIR / "recide_r_stdout.log").write_text(proc.stdout[-20000:] + "\n--- stderr ---\n" + proc.stderr[-8000:])
        if proc.returncode != 0 or not (tmp / "props.csv").exists():
            return None, elapsed, proc.stdout[-3000:], f"Rscript exited {proc.returncode}: {proc.stderr[-800:]}"
        est = pd.read_csv(tmp / "props.csv", index_col=0)
    est.index = est.index.astype(str)
    return est, elapsed, proc.stdout[-3000:], None


DEFAULT_OUT = str(config.RESULTS_DIR / "extension" / "recide_anatomic.json")


def estimates_path(out: Path) -> Path:
    """Only the default run writes the shared estimates file; a run that names its own report (--out) keeps its
    estimates beside it, so a verification re-run cannot overwrite the registered estimates (OPEN_DEFECTS D38's rule,
    applied here 2026-10-07)."""
    return (config.RESULTS_DIR / "extension" / "ivygap_recide.csv" if str(out) == DEFAULT_OUT
            else out.with_suffix(".csv"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-cores", type=int, default=2)
    # No budget by default (user directive 2026-10-01: "I don't want no budget"). The 6-h default
    # that stood here killed the first real run at 06:05 on 2026-10-02 after ~9 h of work.
    ap.add_argument("--budget", type=int, default=None, help="seconds; default: no limit")
    ap.add_argument("--out", default=DEFAULT_OUT)
    a = ap.parse_args()
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)

    matrix = json.loads((config.RESULTS_DIR / "run_provenance.json").read_text())["matrix"]
    gene_rec = json.loads((config.BENCH_DIR / "signature_genes.json").read_text())
    expr, meta = load_cached()
    anat = [s for s in expr.columns if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    genes = [g for g in gene_rec["genes"] if g in set(expr.index)]
    if config.sha256_strings(genes) != gene_rec.get("sha256"):
        print("BLOCKED: gene space does not reproduce the leaderboard's hash"); return 2

    pick, k_used, n_elig = donor_draw()
    if not pick:
        print("BLOCKED: no k in (10, 15, 20) covers every roster type above ReCIDE's threshold")
        return 2
    print(f"donor draw: k = {k_used} of {n_elig} eligible training donors")
    try:
        ref, cells, cmeta, counts = recide_reference(pick, genes, matrix)
    except ValueError as exc:
        print(f"BLOCKED: {exc}"); return 2
    bulk = expr.loc[genes, anat]
    bulk = bulk / bulk.sum(axis=0) * 1e6
    print(f"reference: {counts.shape[1]:,} cells from {len(pick)} donors, {len(genes)} genes")

    est, elapsed, tail, failed = call_recide(counts, genes, cells, cmeta, bulk, a.n_cores, a.budget)
    report = {"panel": "extension (post-registration)", "method": "recide",
              "implementation": IMPLEMENTATION,
              "matrix": matrix, "n_genes": len(genes), "gene_sha256": gene_rec.get("sha256"),
              "donor_draw": {"rule": DRAW_RULE, "k": k_used, "n_eligible": n_elig, "donors": pick},
              "n_reference_cells": int(counts.shape[1]), "elapsed_seconds": round(elapsed, 1),
              "r_stdout_tail": tail, "failed": failed}
    if failed:
        out.write_text(json.dumps(report, indent=2)); print(failed); return 1
    missing = [c for c in config.CELL_TYPES if c not in est.columns]
    report["types_dropped_by_recide"] = missing
    est = est.reindex(index=[str(c) for c in bulk.columns], columns=list(config.CELL_TYPES))
    data = DeconvolutionInput(bulk=bulk, references=(ref.subset_genes(genes),), manifest=meta.loc[anat])
    # A sample ReCIDE could not estimate stays NaN -- never zero-filled into a fake composition.
    frame = finalize_estimates(est.to_numpy(), data, covered=None, returns_cell_fractions=False)
    est_csv = estimates_path(out)
    frame.to_csv(est_csv)
    res = acs_score(frame, meta.loc[anat], method="recide", n_permutations=10000, n_boot=2000).summary()
    report.update({"acs": round(float(res["acs"]), 4),
                   "ci": [round(float(res["ci_low"]), 4), round(float(res["ci_high"]), 4)],
                   "null_p": float(res["null_p"]), "n_pairs": int(res["n_constraint_tumor_pairs"]),
                   "n_tumors": int(res["n_tumors"]),
                   "n_samples_unestimated": int(frame.isna().all(axis=1).sum()),
                   "estimates": str(est_csv.relative_to(config.PROJECT_ROOT))})
    out.write_text(json.dumps(report, indent=2))
    print(f"ReCIDE: ACS {report['acs']} {report['ci']} on {report['n_pairs']} pairs, "
          f"null_p {report['null_p']:.4g}; {elapsed:.0f}s; dropped types {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
