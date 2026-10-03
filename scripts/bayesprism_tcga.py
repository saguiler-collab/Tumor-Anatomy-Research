"""GENUINE BayesPrism on TCGA GBM and LGG -- unbudgeted, in two configurations.

USER DIRECTIVE (2026-10-01): BayesPrism must run as the real package on the GBM and LGG data, with
no time budget. In the registered TCGA arms it never did -- every row labelled `bayesprism` there
is the Python reimplementation (it exceeded its 2,400 s budget or lost its socket workers, D18).

Two configurations, on the registered raw/X (h5ad) arm's exact inputs -- same ABSOLUTE samples,
same reference cells (build_from_h5ad, raw/X, all donors, capped sampler, seed), same
reference-only marker space, CPM re-normalised on it:

  pipeline  this project's driver (R/run_bayesprism.R): key = NULL, <= 300 cells per type,
            no gene cleanup. What the registered h5ad arm would have recorded had the package
            finished.
  authors   the package authors' tutorial (R/run_bayesprism_authors.R): key = "Tumor" (BayesPrism
            learns each tumour's own malignant expression), cleanup.genes (ribosomal,
            mitochondrial-ribosomal, ribosomal pseudogenes, chrM, MALAT1, chrX, chrY),
            protein-coding genes only, all reference cells as exact counts.

PREDICTION, FIXED BEFORE EITHER RAN (2026-10-01; verdict rule fixed before any output in
prespecified/bayesprism_authors_prediction.md, which supersedes the section reference that stood
here): the B-column diagnostics
showed B acting as an overflow for tumour signal (SVR, GBM: 92% of B's estimate is re-absorbed by
Tumor when B is removed). If inter-tumour expression differences are that overflow, the authors'
configuration -- which models each tumour's own expression -- should LOWER the B share of the
lymphoid compartment relative to the pipeline configuration.

Writes results/extension/bayesprism_tcga_{cohort}_{variant}.json and the estimates CSV.
Selects nothing; enters no registered statistic.
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402
from lymphoid_ordering import COMPARE, _truth, k4s, ordering_of, renormalise  # noqa: E402

from ivygap import config                                      # noqa: E402
from ivygap.data.reference import build_from_h5ad, select_signature_genes  # noqa: E402
from ivygap.deconv import r_bridge                              # noqa: E402
from ivygap.deconv.base import DeconvolutionInput, finalize_estimates  # noqa: E402

OUT = config.RESULTS_DIR / "extension"
UNBOUNDED = 10 ** 9          # seconds: no budget, by directive


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True, choices=["gbm", "lgg"])
    ap.add_argument("--variant", required=True, choices=["authors", "pipeline"])
    ap.add_argument("--n-cores", type=int, default=1)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"bayesprism_tcga_{a.cohort}_{a.variant}"
    base = "" if a.cohort == "gbm" else "_lgg"

    with gzip.open(config.PROCESSED_DIR / f"tcga_{a.cohort}_bulk_cpm.csv.gz", "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(c): c for c in bulk.columns}
    shared_s = sorted(set(absol.index) & set(bmap))
    cols = [bmap[k] for k in shared_s]
    purity = absol.loc[shared_s, PURITY_COL].to_numpy(dtype="float64")

    ref, sc_expr, sc_meta = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad",
                                            matrix="raw/X", restrict_to_genes=bulk.index,
                                            export=False)
    shared_g = [g for g in ref.profile.index if g in bulk.index]
    genes = [g for g in select_signature_genes(ref.subset_genes(shared_g),
                                               n_per_type=config.SIGNATURE_GENES_PER_TYPE)
             if g in bulk.index]
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6
    manifest = pd.DataFrame({"patient_id": ["-".join(str(c).split("-")[:3]) for c in sub.columns]},
                            index=sub.columns)
    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),), manifest=manifest,
                              cell_types=tuple(ref.cell_types), bulk_full=bulk.loc[:, cols])
    report = {"panel": "extension (post-registration): genuine BayesPrism, unbudgeted",
              "cohort": a.cohort, "variant": a.variant, "matrix": "raw/X",
              "n_samples": len(cols), "n_marker_genes": len(genes),
              "n_reference_cells": int(sc_expr.shape[1]), "n_cores": a.n_cores, "failed": None}
    print(f"{tag}: {len(cols)} samples x {len(genes)} marker genes; {sc_expr.shape[1]:,} cells",
          flush=True)

    t0 = time.perf_counter()
    try:
        if a.variant == "pipeline":
            r_bridge.set_cell_source(config.PRIMARY_REFERENCE, sc_expr, sc_meta)
            config.R_SOCKET_CLUSTER_CORES = a.n_cores
            est = r_bridge.run_r_method("bayesprism", data, timeout=UNBOUNDED)
            report["r_stdout_tail"] = r_bridge.LAST_R_STDOUT.get("bayesprism", "")[-3000:]
        else:
            counts = sc_expr.loc[genes].to_numpy(dtype="float64") * (
                sc_meta.loc[sc_expr.columns, "library_size"].to_numpy(dtype="float64") / 1e6)
            from ivygap.deconv.extension import counts_are_integral
            ok, rel = counts_are_integral(counts)
            if not ok:
                raise RuntimeError(f"recovered counts are not integers (max relative deviation {rel:.3g})")
            with tempfile.TemporaryDirectory(prefix="ivygap_bp_authors_") as tmp:
                tmp = Path(tmp)
                sio.mmwrite(tmp / "counts.mtx", sparse.csc_matrix(np.round(counts)))
                (tmp / "genes.txt").write_text("\n".join(genes) + "\n")
                pd.DataFrame({"cell": list(sc_expr.columns),
                              "cell_type": sc_meta.loc[sc_expr.columns, "cell_type"].astype(str).to_numpy()}
                             ).to_csv(tmp / "meta.csv", index=False)
                sub.T.to_csv(tmp / "bulk.csv")                       # samples x genes
                cfg = {"counts": str(tmp / "counts.mtx"), "genes": str(tmp / "genes.txt"),
                       "meta": str(tmp / "meta.csv"), "bulk": str(tmp / "bulk.csv"),
                       "out": str(tmp / "theta.csv"), "seed": config.RANDOM_SEED,
                       "n_cores": a.n_cores, "key": "Tumor", "cleanup": True}
                (tmp / "args.json").write_text(json.dumps(cfg))
                # FREE THE ATLAS BEFORE THE LONG R CALL. Everything R needs is now on disk; holding the
                # ~2 GB dense cell matrix in this process for the hours BayesPrism samples starved the
                # 8 GB machine (2026-10-01: 2.2 GB idle in the parent during the GBM run).
                import gc                                                   # noqa: PLC0415
                n_ref_cells = int(sc_expr.shape[1])
                del sc_expr, sc_meta, counts
                gc.collect()
                log = OUT / f"{tag}_r.log"
                with open(log, "w") as fh:
                    proc = subprocess.run(["Rscript", str(config.PROJECT_ROOT / "R" / "run_bayesprism_authors.R"),
                                           str(tmp / "args.json")], stdout=fh, stderr=subprocess.STDOUT)
                report["r_log"] = str(log.relative_to(config.PROJECT_ROOT))
                if proc.returncode != 0 or not (tmp / "theta.csv").exists():
                    raise RuntimeError(f"Rscript exited {proc.returncode}; see {log.name}")
                est = pd.read_csv(tmp / "theta.csv", index_col=0)
                used = (tmp / "theta_genes_used.txt").read_text().split()
                report["n_genes_after_authors_cleanup"] = len(used)
                report["genes_removed_by_cleanup"] = sorted(set(genes) - set(used))
                first = pd.read_csv(tmp / "theta_first.csv", index_col=0)
                first.to_csv(OUT / f"{tag}_theta_first.csv")
            est.index = est.index.astype(str)
            est = est.reindex(index=[str(c) for c in sub.columns], columns=list(ref.cell_types))
    except Exception as exc:                                       # noqa: BLE001
        report["failed"] = f"{type(exc).__name__}: {str(exc)[:600]}"
    report["elapsed_seconds"] = round(time.perf_counter() - t0, 1)
    if report["failed"]:
        (OUT / f"{tag}.json").write_text(json.dumps(report, indent=2))
        print("FAILED:", report["failed"]); return 1

    frame = finalize_estimates(est.to_numpy(), data, covered=None, returns_cell_fractions=False)
    frame.index = list(sub.columns)
    frame.to_csv(OUT / f"{tag}.csv")
    t = frame["Tumor"].to_numpy(dtype="float64"); ok = np.isfinite(t)
    report.update({"spearman_vs_purity": round(float(stats.spearmanr(t[ok], purity[ok]).statistic), 4),
                   "bias_tumor_minus_purity": round(float(np.mean(t[ok] - purity[ok])), 4)})
    truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
    g = frame.copy(); g.index = [k4s(i) for i in g.index]; g = g[~g.index.duplicated()]
    idx = [k for k in g.index if k in truth.index]
    ren = renormalise(g.loc[idx]); mu = ren.mean(); ok2 = ren.dropna()
    report.update({"n_lymphoid_matched": len(idx),
                   "lymphoid_mean": {c: round(float(mu[c]), 4) for c in COMPARE},
                   "lymphoid_ordering": ordering_of(mu), "T_exceeds_B": bool(mu["T_cell"] > mu["B_cell"]),
                   "frac_samples_B_over_T": round(float((ok2.B_cell > ok2.T_cell).mean()), 4) if len(ok2) else None,
                   "n_zero_lymphoid": int(ren.isna().all(axis=1).sum()),
                   "methylation_ordering": ordering_of(truth.mean())})
    anch = config.RESULTS_DIR / "bisque_anchoring.json"
    if anch.exists():
        prior = pd.Series(json.loads(anch.read_text())["reference_prior_all_types"])
        report["L1_cohort_mean_to_reference_prior"] = round(float((frame[list(prior.index)].mean() - prior).abs().sum()), 4)
    (OUT / f"{tag}.json").write_text(json.dumps(report, indent=2))
    print(f"{tag}: purity rho {report['spearman_vs_purity']:+.4f}; lymphoid {report['lymphoid_ordering']} "
          f"(B>T in {report['frac_samples_B_over_T']}); {report['elapsed_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
