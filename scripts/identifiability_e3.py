"""
identifiability_e3.py -- prespecified/identifiability_e3_gene_resampling.md (Extension E3, exploratory): does a second
truth-free perturbation -- which genes inform the fit -- order the compartments the way agreement with DNA does, as
E2's loss-scale stability did?

Fits genuine DESeq2 `unmix` (E2's pre-declared setting) on 10 random halves of the marker genes, per cohort, and
measures D3 stability per compartment. E2's truth agreement is read, not refit. Resumable: each fit is saved under
results/identifiability/e3/ and reused.

    python3 scripts/identifiability_e3.py
"""
from __future__ import annotations

import gzip
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import extension_tcga  # noqa: E402
import identifiability_diagnostics as idd  # noqa: E402

OUT_DIR = config.RESULTS_DIR / "identifiability" / "e3"
OUT = config.RESULTS_DIR / "identifiability_e3.json"
K = 10
SEED = config.RANDOM_SEED + 3


def marker_genes(cohort: str) -> list[str]:
    """extension_tcga.run's marker list for the frozen arm, in its stored order (reference-only selection)."""
    from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: PLC0415
    with gzip.open(config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz", "rt") as fh:
        genes_in_bulk = set(pd.read_csv(fh, index_col=0, usecols=[0]).index)
    ref = load_frozen_reference()
    shared = [g for g in ref.profile.index if g in genes_in_bulk]
    picked = select_signature_genes(ref.subset_genes(shared), n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    return [g for g in picked if g in genes_in_bulk]


def subsets(genes: list[str]) -> list[list[str]]:
    rng = np.random.default_rng(SEED)
    return [[genes[i] for i in sorted(rng.choice(len(genes), size=len(genes) // 2, replace=False))] for _ in range(K)]


def fit(cohort: str, k: int, keep: list[str]) -> pd.DataFrame:
    path = OUT_DIR / f"unmix_{cohort}_sub{k:02d}.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0)
    for v in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):         # E2's pre-declared setting
        os.environ.pop(v, None)
    res = extension_tcga.run(cohort, "frozen", ["deseq2_unmix"], return_estimates=True, genes_keep=keep)
    est = res["_estimates"]["deseq2_unmix"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    est.to_csv(path)
    return est


def stability(fits: list[pd.DataFrame], comp: str) -> tuple[float, int, int]:
    series = [idd.comp(f, comp) for f in fits]
    usable = [s for s in series if s.nunique() > 1]
    if len(usable) < 3:
        return float("nan"), len(usable), len(series) - len(usable)
    idx = usable[0].index
    for s in usable[1:]:
        idx = idx.intersection(s.index)
    r = [stats.spearmanr(a[idx], b[idx]).statistic for a, b in itertools.combinations(usable, 2)]
    return float(np.mean(r)), len(usable), len(series) - len(usable)


def main() -> int:
    e2 = json.loads((config.RESULTS_DIR / "identifiability_diagnostics.json").read_text())
    units_e2 = e2["units"]
    controls = {"e2_tumour_truth_read": {h: units_e2[f"Tumor|{h}"]["truth_rho_unmix"] for h in ("gbm", "lgg")}}
    controls["e2_tumour_truth_equals_reported"] = (round(controls["e2_tumour_truth_read"]["gbm"], 4) == 0.7801
                                                   and round(controls["e2_tumour_truth_read"]["lgg"], 4) == 0.5916)
    gbm_genes = marker_genes("gbm")
    subs = subsets(gbm_genes)
    out: dict = {"rule": "prespecified/identifiability_e3_gene_resampling.md", "k": K, "seed": SEED,
                 "n_marker_genes_gbm": len(gbm_genes), "subset_size": len(subs[0]), "controls": controls, "units": {}}
    tumour_rhos: dict = {}
    for cohort in ("gbm", "lgg"):
        genes = set(marker_genes(cohort))
        tr = idd.truths(cohort)
        fits = []
        for k, s in enumerate(subs):
            keep = [g for g in s if g in genes]
            print(f"{cohort}: subset {k + 1}/{K} ({len(keep)} genes) ...", flush=True)
            fits.append(fit(cohort, k, keep))
        tumour_rhos[cohort] = [round(idd.truth_rho(idd.comp(f, "Tumor"), "Tumor", tr)[0], 4) for f in fits]
        for comp in idd.PRIMARY + idd.SECONDARY:
            d3, n_ok, n_const = stability(fits, comp)
            u = units_e2.get(f"{comp}|{cohort}", {})
            out["units"][f"{comp}|{cohort}"] = {"d3_stability": None if not np.isfinite(d3) else round(d3, 4),
                                                "fits_used": n_ok, "fits_constant": n_const,
                                                "d1_stability_e2": u.get("d1_stability"),
                                                "truth_rho_unmix_e2": u.get("truth_rho_unmix")}
    controls["subset_fit_tumour_rho"] = {h: {"values": v, "median": float(np.median(v))} for h, v in tumour_rhos.items()}
    controls["subsets_preserve_signal"] = all(np.median(v) >= 0.40 for v in tumour_rhos.values())

    def test(names, xk, yk):
        pts = [(out["units"][f"{c}|{h}"][xk], out["units"][f"{c}|{h}"][yk]) for c in names for h in ("gbm", "lgg")]
        pts = [p for p in pts if p[0] is not None and p[1] is not None and all(np.isfinite(p))]
        if len(pts) < 4:
            return {"n_units": len(pts), "note": "too few estimable units"}
        rho, p = idd.exact_one_sided([q[0] for q in pts], [q[1] for q in pts])
        return {"n_units": len(pts), "rho": round(rho, 4), "p_one_sided": round(p, 4)}

    h1 = test(idd.PRIMARY, "d3_stability", "truth_rho_unmix_e2")
    h1["reading"] = ("uninformative: the subsets destroyed the signal" if not controls["subsets_preserve_signal"] else
                     "SUPPORTED" if (h1.get("rho", 0) > 0 and h1.get("p_one_sided", 1) < 0.05) else
                     "NOT SUPPORTED" if h1.get("rho", 0) <= 0 else "INCONCLUSIVE")
    out["H1_E3_primary"] = h1
    out["H1s_E3_secondary"] = test(idd.PRIMARY + idd.SECONDARY, "d3_stability", "truth_rho_unmix_e2")
    c1 = test(idd.PRIMARY, "d3_stability", "d1_stability_e2")
    out["C1_d3_vs_d1_descriptive"] = c1
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print("controls:", json.dumps(controls)[:400])
    print("H1-E3:", h1, "| H1s:", out["H1s_E3_secondary"], "| C1 (D3 vs D1):", c1)
    for k_, v in out["units"].items():
        print(f"  {k_:20s} D3 {v['d3_stability']}  D1 {v['d1_stability_e2']}  truth {v['truth_rho_unmix_e2']}")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
