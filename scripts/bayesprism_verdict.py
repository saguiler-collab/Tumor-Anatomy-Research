"""
bayesprism_verdict.py -- apply the rule in prespecified/bayesprism_authors_prediction.md.

Written before either BayesPrism configuration produced output, so the rule is code before it
is a number. Reads results/extension/bayesprism_tcga_{cohort}_{authors,pipeline}.json/.csv and
writes results/extension/bayesprism_verdict.json. Reports BLOCKED, not a guess, for a cohort
whose two runs are not both complete.

    python3 scripts/bayesprism_verdict.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
from lymphoid_ordering import COMPARE, _truth, k4s, renormalise  # noqa: E402

EXT = config.RESULTS_DIR / "extension"
MIN_PAIRED = 10


def verdict_from(authors_b: float | None, pipeline_b: float | None) -> str:
    """The pre-declared reading: lower B share under the authors' configuration."""
    if authors_b is None or pipeline_b is None or not (np.isfinite(authors_b) and np.isfinite(pipeline_b)):
        return "BLOCKED"
    return "SUPPORTED" if authors_b < pipeline_b else "NOT SUPPORTED"


def paired(a: pd.DataFrame, p: pd.DataFrame, truth_index) -> dict:
    """Per-sample B share, both runs, on samples with lymphoid signal in both."""
    def share(f):
        g = f.copy()
        g.index = [k4s(i) for i in g.index]
        g = g[~g.index.duplicated()]
        g = g.loc[[k for k in g.index if k in truth_index]]
        return renormalise(g)["B_cell"]
    sa, sp = share(a), share(p)
    both = sa.dropna().index.intersection(sp.dropna().index)
    out = {"n_samples_with_signal_in_both": int(len(both))}
    if len(both) < MIN_PAIRED:
        out["per_sample"] = "INCONCLUSIVE"
        return out
    d = (sa[both] - sp[both]).to_numpy()
    nz = d[d != 0]
    out.update({"median_difference_authors_minus_pipeline": round(float(np.median(d)), 4),
                "n_lower_under_authors": int((d < 0).sum()),
                "n_higher_under_authors": int((d > 0).sum()),
                "wilcoxon_p_two_sided": (round(float(stats.wilcoxon(nz).pvalue), 6)
                                         if nz.size >= MIN_PAIRED else None)})
    return out


def main() -> int:
    report = {"rule": "prespecified/bayesprism_authors_prediction.md", "cohorts": {}}
    for cohort in ("gbm", "lgg"):
        runs = {v: EXT / f"bayesprism_tcga_{cohort}_{v}" for v in ("authors", "pipeline")}
        rec = {v: json.loads(r.with_suffix(".json").read_text())
               if r.with_suffix(".json").exists() else None for v, r in runs.items()}
        missing = [v for v, r in rec.items() if r is None or r.get("failed")]
        if missing:
            report["cohorts"][cohort] = {"verdict": "BLOCKED",
                                         "why": f"no completed run for: {', '.join(missing)}"}
            continue
        ab, pb = (rec[v]["lymphoid_mean"]["B_cell"] for v in ("authors", "pipeline"))
        base = "" if cohort == "gbm" else "_lgg"
        truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
        est = {v: pd.read_csv(r.with_suffix(".csv"), index_col=0) for v, r in runs.items()}
        report["cohorts"][cohort] = {
            "verdict": verdict_from(ab, pb),
            "B_share_authors": ab, "B_share_pipeline": pb,
            "per_sample": paired(est["authors"], est["pipeline"], truth.index),
            "alongside": {v: {k: rec[v].get(k) for k in ("lymphoid_ordering", "T_exceeds_B",
                                                          "frac_samples_B_over_T", "n_zero_lymphoid",
                                                          "spearman_vs_purity")}
                          for v in ("authors", "pipeline")},
        }
    (EXT / "bayesprism_verdict.json").write_text(json.dumps(report, indent=2))
    for c, r in report["cohorts"].items():
        print(f"{c}: {r['verdict']}" + (f" -- {r['why']}" if "why" in r else
                                       f" (B share authors {r['B_share_authors']} vs pipeline "
                                       f"{r['B_share_pipeline']})"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
