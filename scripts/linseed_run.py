"""
linseed_run.py -- the reference-free arm's second instrument (prespecified/linseed_config.md).

POST-REGISTRATION, EXPLORATORY. Genuine linseed 0.99.3 (vendor/linseed) through R/run_linseed.R:
the package README's tutorial step for step, k fixed a priori at the roster size (as CDSeq's T),
topGenes 10,000. Components are named by CDSeq's two labellings, imported unchanged from
scripts/cdseq_anatomic.py: (a) GBmap profile correlation, for labelling only; (b) GBM-specific
marker enrichment, with no atlas. Components mapped to one roster type are summed; an unclaimed
type is absent (NaN), never zero.

  anatomic   Ivy GAP's 122 anatomic samples, CPM from the RSEM expected counts (CDSeq's source
             data); ACS per labelling with the registered scorer (10,000 permutations).
  tcga       TCGA-GBM / TCGA-LGG bulk CPM, on the samples with ABSOLUTE purity (as
             extension_tcga); the tumour-labelled component against purity (Spearman), and the
             lymphoid ordering against methylation where lymphoid types are claimed.

Needs the memory of an otherwise idle 8 GB machine (topGenes 10,000): run only when no other
heavy job is running.

    python3 scripts/linseed_run.py anatomic
    python3 scripts/linseed_run.py tcga --cohort gbm
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
from ivygap.anatomic.acs import score as acs_score  # noqa: E402
from ivygap.data.load_ivygap import load_cached  # noqa: E402
from cdseq_anatomic import COUNTS, label_by_markers, label_by_profile, to_roster  # noqa: E402

OUT = config.RESULTS_DIR / "linseed"
K = len(config.CELL_TYPES)
LABELLINGS = (("gbmap_profile_correlation", label_by_profile), ("gbm_marker_enrichment", label_by_markers))
PIPELINE = {"top_genes": 10000, "sig_iters": 100, "pval": 0.01, "k": K, "seed": config.RANDOM_SEED}


def run_linseed(bulk: pd.DataFrame, tag: str) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Proportions (samples x k) and signatures (genes x k, gene-expression scale)."""
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ivygap_linseed_") as tmp:
        tmp = Path(tmp)
        bulk.to_csv(tmp / "bulk.csv")
        args = dict(PIPELINE, bulk=str(tmp / "bulk.csv"), out_prop=str(OUT / f"{tag}_prop.csv"),
                    out_sig=str(OUT / f"{tag}_sig.csv"))
        (tmp / "args.json").write_text(json.dumps(args))
        r = subprocess.run(["Rscript", str(ROOT / "R" / "run_linseed.R"), str(tmp / "args.json")],
                           capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"run_linseed.R exited {r.returncode}: {r.stderr[-1500:]}")
    prop = pd.read_csv(OUT / f"{tag}_prop.csv", index_col=0)
    prop.index = prop.index.astype(str)
    sig = pd.read_csv(OUT / f"{tag}_sig.csv", index_col=0)
    return prop, sig, r.stdout.strip().splitlines()[-1]


def labelled(prop: pd.DataFrame, sig: pd.DataFrame) -> dict[str, tuple[dict, pd.DataFrame]]:
    """CDSeq's labelling and mapping, on Linseed's outputs (components as columns of `sig`)."""
    out = {}
    for name, fn in LABELLINGS:
        labels, _ = fn(sig)
        if labels:
            out[name] = (labels, to_roster(prop.T, labels))       # to_roster wants components x samples
    return out


def anatomic() -> int:
    expr, meta = load_cached()
    anat = [s for s in expr.columns if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    counts = pd.read_csv(COUNTS, index_col=0)
    counts.columns = counts.columns.astype(str)
    counts.index = counts.index.astype(str)
    if not all(s in counts.columns for s in anat):
        print("BLOCKED: not every anatomic sample has counts"); return 2
    sub = counts[anat]
    genes = pd.read_csv(config.IVYGAP_GENES_PATH)
    i2s = dict(zip(genes["gene_id"].astype(str), genes["gene_symbol"].astype(str)))
    sub.index = [i2s.get(g, g) for g in sub.index]
    sub = sub[~sub.index.duplicated(keep="first")]
    cpm = sub / sub.sum(axis=0) * 1e6
    prop, sig, log = run_linseed(cpm, "ivygap_anatomic")
    print(log)
    report = {"rule": "prespecified/linseed_config.md", "implementation": "R:linseed 0.99.3 (vendored @2435a8ce)",
              "arm": "Ivy GAP anatomic", "n_samples": len(anat), "pipeline": PIPELINE, "linseed_log": log,
              "reference_passed_to_linseed": False, "labellings": {}}
    for name, (labels, est) in labelled(prop, sig).items():
        est.to_csv(OUT / f"ivygap_estimates_{name}.csv")
        res = acs_score(est, meta.loc[anat], method=f"linseed__{name}", n_permutations=10000, n_boot=2000)
        s = res.summary()
        report["labellings"][name] = {
            "assignments": {str(k): v for k, v in labels.items()},
            "roster_types_covered": sorted({v for v in labels.values() if v in config.CELL_TYPES}),
            "acs": round(float(s["acs"]), 4), "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
            "null_p": float(s["null_p"]), "n_pairs": int(s["n_constraint_tumor_pairs"]),
            "beats_null": bool(s["beats_null"])}
        print(f"  {name}: ACS {s['acs']:.4f} [{s['ci_low']:.4f}, {s['ci_high']:.4f}] p {s['null_p']:.4g} "
              f"pairs {s['n_constraint_tumor_pairs']} covered {report['labellings'][name]['roster_types_covered']}")
    (config.RESULTS_DIR / "linseed_anatomic.json").write_text(json.dumps(report, indent=2))
    print("wrote results/linseed_anatomic.json")
    return 0


def tcga(cohort: str) -> int:
    from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: PLC0415
    from lymphoid_ordering import COMPARE, _truth, k4s, ordering_of, renormalise  # noqa: PLC0415
    with gzip.open(config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz", "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(c): c for c in bulk.columns}
    shared = sorted(set(absol.index) & set(bmap))
    cols = [bmap[k] for k in shared]
    purity = pd.Series(absol.loc[shared, PURITY_COL].to_numpy(float), index=cols)
    prop, sig, log = run_linseed(bulk[cols], f"tcga_{cohort}")
    print(log)
    base = "" if cohort == "gbm" else f"_{cohort}"
    truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
    report = {"rule": "prespecified/linseed_config.md", "implementation": "R:linseed 0.99.3 (vendored @2435a8ce)",
              "arm": f"TCGA-{cohort.upper()} bulk only", "n_samples": len(cols), "pipeline": PIPELINE,
              "linseed_log": log, "methylation_ordering": ordering_of(truth.mean()), "labellings": {}}
    for name, (labels, est) in labelled(prop, sig).items():
        est.to_csv(OUT / f"tcga_{cohort}_estimates_{name}.csv")
        covered = sorted({v for v in labels.values() if v in config.CELL_TYPES})
        rec = {"assignments": {str(k): v for k, v in labels.items()}, "roster_types_covered": covered}
        if "Tumor" in covered:
            t = est.loc[cols, "Tumor"].to_numpy(float)
            ok = np.isfinite(t)
            rec["tumour_spearman_vs_purity"] = round(float(stats.spearmanr(t[ok], purity.to_numpy()[ok]).statistic), 4)
            rec["n_purity"] = int(ok.sum())
        claimed = [c for c in COMPARE if c in covered]
        rec["lymphoid_types_claimed"] = claimed
        if len(claimed) == len(COMPARE):
            # only when T, NK and B are all claimed: an unclaimed type is ABSENT, and zero-filling
            # it would manufacture an ordering (B claimed, T absent -> "B > T")
            g = est.copy()
            g.index = [k4s(i) for i in g.index]
            g = g[~g.index.duplicated()]
            idx = [k for k in g.index if k in truth.index]
            mu = renormalise(g.loc[idx]).mean()
            rec.update({"lymphoid_ordering": ordering_of(mu), "T_exceeds_B": bool(mu["T_cell"] > mu["B_cell"]),
                        "n_lymphoid": len(idx)})
        else:
            rec["lymphoid_ordering"] = None
            rec["lymphoid_note"] = "not every lymphoid type is claimed by a component; no ordering is reported"
        report["labellings"][name] = rec
        print(f"  {name}: covered {covered} | {rec}")
    (config.RESULTS_DIR / f"linseed_tcga_{cohort}.json").write_text(json.dumps(report, indent=2, default=float))
    print(f"wrote results/linseed_tcga_{cohort}.json")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("arm", choices=["anatomic", "tcga"])
    ap.add_argument("--cohort", choices=["gbm", "lgg"], default="gbm")
    a = ap.parse_args()
    return anatomic() if a.arm == "anatomic" else tcga(a.cohort)


if __name__ == "__main__":
    raise SystemExit(main())
