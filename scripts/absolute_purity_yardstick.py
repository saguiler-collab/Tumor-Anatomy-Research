#!/usr/bin/env python3
"""
absolute_purity_yardstick.py — YARDSTICK 2: rank the methods against DNA, not RNA.

WHY THIS IS WORTH MORE THAN ANOTHER RNA BENCHMARK
--------------------------------------------------
The registered primary outcome correlates the ACS ranking with an accuracy ranking from
synthetic pseudobulk — and **both arms use GBmap** (correction C4). So "anatomy agrees with
truth" has meant "two GBmap-based rankings agree with each other".

ABSOLUTE purity is measured from **DNA** (SNP arrays, PanCanAtlas). It shares no failure mode
with any RNA deconvolution method under test, and it measures the same quantity the Tumor
column estimates. Nguyen et al. (2024) criticise simulation-based benchmarking for rewarding
methods that share the simulator's assumptions; a DNA yardstick is immune to that in a way no
RNA yardstick can be.

WHAT IS AND IS NOT INDEPENDENT — stated plainly, because it is easy to overclaim
-------------------------------------------------------------------------------
INDEPENDENT: the yardstick. ABSOLUTE purity comes from DNA copy number, not from GBmap, not
from this pipeline, not from any simulation.

NOT INDEPENDENT: the reference. Every method here solves against the vendored frozen signature,
which `reference_frozen/PROVENANCE.json` records as "derived from GBmap single-cell atlas,
annotation_level_3, 100 donors". So this measures *methods scored against a GBmap reference on
a DNA yardstick* — it removes the shared dependence on the **ground truth**, not on the
reference. That is the improvement it offers and the limit of it.

TWO DECLARED CONSEQUENCES OF THE FROZEN SIGNATURE
--------------------------------------------------
  * It carries **no cross-donor variance**, so MuSiC degenerates to residual-weighted NNLS.
    The pipeline detects this and labels it; the label is carried into the output here.
  * Its `cell_size` factors have real spread (L_median: Tumor 5,897, T_cell 3,483), unlike the
    leaderboard's uniform vector (D12). So on THIS cohort the mRNA-to-cell conversion is NOT
    the identity, which makes the Tumor column a cell fraction — the right quantity to compare
    against a DNA-derived cell purity.

    python scripts/absolute_purity_yardstick.py
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.deconv.base import DeconvolutionInput                  # noqa: E402
from ivygap.deconv.registry import build_methods                   # noqa: E402
from ivygap.data.reference import load_frozen_reference            # noqa: E402

BULK = config.PROCESSED_DIR / "tcga_gbm_bulk_cpm.csv.gz"
ABS_T = config.RAW_DIR / "tcga" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
PURITY_COL = "Cancer DNA fraction"


def key4(s: str) -> str:
    return "-".join(str(s).split("-")[:4])


def main() -> int:
    for p in (BULK, ABS_T):
        if not p.exists():
            print(f"BLOCKED: {p} missing."); return 2

    with gzip.open(BULK, "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    print(f"bulk: {bulk.shape[0]:,} genes x {bulk.shape[1]} samples")

    absol = pd.read_csv(ABS_T, sep="\t")
    absol = absol.dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(c): c for c in bulk.columns}
    shared = sorted(set(absol.index) & set(bmap))
    print(f"samples with a DNA purity call: {len(shared)}")
    if len(shared) < 50:
        print("BLOCKED: too few samples to rank anything."); return 2
    cols = [bmap[k] for k in shared]
    purity = absol.loc[shared, PURITY_COL].to_numpy(dtype="float64")
    print(f"  purity: median {np.median(purity):.3f}  range "
          f"{purity.min():.3f}-{purity.max():.3f}")

    ref = load_frozen_reference()
    shared_genes = [g for g in ref.profile.index if g in bulk.index]

    # MARKER SUBSET, and why. The full shared space is 17,611 genes and nu-SVR is superlinear
    # in genes: a first run of this script spent minutes inside SVR alone with three SVR-based
    # methods still queued. Every other arm of this project scores on a marker subset (the
    # leaderboard's is 657), so using one here is consistent rather than special.
    #
    # The rule is `select_signature_genes` on the frozen reference itself — the same rule used
    # everywhere else, and it reads ONLY the reference. It does not see purity, so it cannot
    # be selecting for the outcome this script measures.
    from ivygap.data.reference import select_signature_genes      # noqa: PLC0415
    picked = select_signature_genes(ref.subset_genes(shared_genes),
                                    n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    genes = [g for g in picked if g in bulk.index]
    print(f"reference: frozen GBmap-derived signature, {len(shared_genes):,} genes shared "
          f"with the bulk")
    print(f"  scored on {len(genes):,} markers "
          f"({config.SIGNATURE_GENES_PER_TYPE}/type, chosen from the reference alone)")
    if len(genes) < config.MIN_GENES_SHARED:
        print("BLOCKED: gene space too small."); return 2
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6                 # re-normalise on the marker space
    manifest = pd.DataFrame(index=sub.columns)

    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),),
                              manifest=manifest, cell_types=tuple(ref.cell_types),
                              bulk_full=bulk.loc[:, cols])
    print(f"\nrunning methods on {sub.shape[1]} samples x {sub.shape[0]:,} genes:")
    out: dict[str, dict] = {}
    per_sample: dict[str, pd.Series] = {}
    for m in build_methods(prefer_r=True):
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  {m.name:26s} FAILED {type(exc).__name__}: {str(exc)[:90]}")
            out[m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:160]}"}
            continue
        if "Tumor" not in est.columns:
            print(f"  {m.name:26s} SKIPPED — models no Tumor column")
            out[m.name] = {"skipped": "does not model Tumor"}
            continue
        t = est.loc[sub.columns, "Tumor"].to_numpy(dtype="float64")
        ok = np.isfinite(t)
        if ok.sum() < 50:
            print(f"  {m.name:26s} SKIPPED — only {ok.sum()} finite estimates")
            out[m.name] = {"skipped": f"only {int(ok.sum())} finite estimates"}
            continue
        # PERSIST THE PER-SAMPLE ESTIMATES. The aggregate correlation cannot support the
        # failure-factor study (prespecified/biological_failure_factors.md), whose outcome is
        # the per-sample signed error. Writing them here avoids a second 30-minute run.
        per_sample[m.name] = pd.Series(t, index=list(sub.columns))
        rho = float(stats.spearmanr(t[ok], purity[ok]).statistic)
        pr = float(stats.pearsonr(t[ok], purity[ok]).statistic)
        out[m.name] = {
            "spearman_vs_purity": round(rho, 4), "pearson_vs_purity": round(pr, 4),
            "n": int(ok.sum()), "mean_tumor": round(float(np.mean(t[ok])), 4),
            "mean_purity": round(float(np.mean(purity[ok])), 4),
            "bias_tumor_minus_purity": round(float(np.mean(t[ok] - purity[ok])), 4),
            "implementation": getattr(m, "implementation_", "python"),
            "degenerate": bool(getattr(m, "degenerate_", False)),
            "degeneracy_reason": getattr(m, "degeneracy_reason_", "") or "",
            "is_control": m.name.startswith("control_")}
        flag = "  [CONTROL]" if out[m.name]["is_control"] else (
            "  [DEGENERATE]" if out[m.name]["degenerate"] else "")
        print(f"  {m.name:26s} rho {rho:+.4f}  bias {out[m.name]['bias_tumor_minus_purity']:+.4f}"
              f"  n={ok.sum()}{flag}")

    real = {k: v["spearman_vs_purity"] for k, v in out.items()
            if "spearman_vs_purity" in v and not v["is_control"]}
    ctrl = {k: v["spearman_vs_purity"] for k, v in out.items()
            if "spearman_vs_purity" in v and v["is_control"]}
    print(f"\nreal methods: {len(real)}; best {max(real, key=real.get)} "
          f"{max(real.values()):+.4f}; worst {min(real, key=real.get)} {min(real.values()):+.4f}")
    print(f"controls: {ctrl}")

    report = {
        "what_this_is": "Yardstick 2: each method's Tumor fraction correlated against "
                        "ABSOLUTE tumour purity, measured from DNA.",
        "independent": "the YARDSTICK only. ABSOLUTE purity comes from DNA copy number and "
                       "shares no failure mode with any RNA method under test.",
        "not_independent": "the REFERENCE. Every method solves against the vendored frozen "
                           "signature, which is GBmap-derived (reference_frozen/"
                           "PROVENANCE.json). This removes the shared dependence on the "
                           "ground truth, not on the reference.",
        "purity_column": PURITY_COL,
        "n_samples": len(shared), "n_genes": len(genes),
        "n_genes_shared_before_marker_selection": len(shared_genes),
        "gene_selection": f"select_signature_genes on the frozen reference, "
                          f"{config.SIGNATURE_GENES_PER_TYPE} per type. Reads the reference "
                          f"only and never the purity calls, so it cannot select for the "
                          f"outcome. Chosen because nu-SVR is superlinear in genes and the "
                          f"full 17,611-gene space did not finish.",
        "bulk_provenance": "data/processed/tcga_gbm_bulk_provenance.json — the input was "
                           "log2(x+1) despite being named counts, and was inverted first",
        "frozen_signature_caveats": [
            "no cross-donor variance, so MuSiC degenerates to residual-weighted NNLS",
            "cell_size factors have real spread here (unlike D12's uniform vector), so the "
            "mRNA-to-cell conversion is NOT the identity on this cohort",
        ],
        "methods": out,
    }
    (config.RESULTS_DIR / "absolute_purity_yardstick.json").write_text(
        json.dumps(report, indent=2))
    if per_sample:
        est = pd.DataFrame(per_sample)
        est.insert(0, "absolute_purity", purity)
        est.index.name = "sample"
        est.to_csv(config.RESULTS_DIR / "absolute_purity_per_sample.csv")
        print(f"wrote results/absolute_purity_per_sample.csv "
              f"({est.shape[0]} samples x {est.shape[1] - 1} methods + purity)")
    print("\nwrote results/absolute_purity_yardstick.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
