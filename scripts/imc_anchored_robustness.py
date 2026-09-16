#!/usr/bin/env python3
"""
imc_anchored_robustness.py — is the IMC-anchored result an artefact of marker COUNT?

The pre-specified test (prespecified/imc_anchored_prediction.md) predicted
ACS(MCP_GBM) > ACS(MCP_GBmap) and it came out the other way. Before that is written up as
evidence against the anatomy test, the obvious alternative explanation has to be killed:

    MCP_GBmap carries 200 markers per population; MCP_GBM carries about 79. MCPcounter's
    estimate is a MEAN over marker genes, so more markers means a smoother, less noisy score,
    and a smoother score could satisfy ordinal constraints more often for reasons that have
    nothing to do with whether the markers are biologically right.

One observation already argues against it — MCP_default won with the FEWEST markers, about 15
per population — but that is an aside, not a test. This is the test: rebuild the GBmap marker
set at MATCHED density and re-run the primary comparison. If GBmap still wins at 79 markers per
population, count is not the explanation.

It also reports a PAIRED comparison. The two arms' bootstrap CIs overlap, and overlapping CIs
do not mean the difference is not significant — the arms are scored on the same nine tumours, so
the difference should be bootstrapped paired, over tumours, rather than eyeballed from two
marginal intervals.

    python scripts/imc_anchored_robustness.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from imc_anchored_test import gbmap_markers, run_mcp, WORK         # noqa: E402

SHARED = {"Tumor", "Macrophage_Microglia"}          # the primary comparison's shared types
DENSITIES = [200, 79, 15]                           # GBmap's own, Ajaib's, MCPcounter default's


def masked(est: pd.DataFrame, types: set[str]) -> pd.DataFrame:
    e = est.copy()
    for c in config.CELL_TYPES:
        if c not in types:
            e[c] = np.nan
    return e


def main() -> int:
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    man = meta.loc[anat]
    bulk_csv = WORK / "bulk_anatomic.csv"
    if not bulk_csv.exists():
        print("BLOCKED: run scripts/imc_anchored_test.py first.")
        return 2

    # INDEX DTYPE. Ivy GAP sample ids are numeric strings, so a CSV round-trip turns the
    # index into int64 while the manifest's index is object. They then share nothing, every
    # (tumour, structure) mean is empty, and the ACS comes back NaN rather than failing —
    # the same int64-vs-object trap this project already hit once on a sample-id check.
    gbm = pd.read_csv(WORK / "estimates_MCP_GBM.csv", index_col=0)
    gbm.index = gbm.index.astype(str)
    missing = [s for s in gbm.index if s not in man.index]
    if missing:
        print(f"BLOCKED: {len(missing)} estimate rows are not in the manifest, e.g. "
              f"{missing[:3]}"); return 2
    ref = acs_score(masked(gbm, SHARED), man, method="MCP_GBM", n_permutations=10000,
                    n_boot=2000)
    rs = ref.summary()
    if not np.isfinite(rs["acs"]):
        print("BLOCKED: the reference arm scored NaN; every comparison against it would be "
              "a NaN comparison silently reported as a loss."); return 2
    print(f"MCP_GBM (Ajaib, ~79 markers/pop): ACS {rs['acs']:.4f}  p={rs['null_p']:.4f}\n")

    out = {"reference_arm": {"tag": "MCP_GBM", "acs": round(float(rs["acs"]), 4),
                             "null_p": float(rs["null_p"])},
           "gbmap_at_matched_density": {}}
    per_tumor = {"MCP_GBM": ref.per_tumor.set_index("tumor_id")["acs"]}

    for n in DENSITIES:
        tag = f"MCP_GBmap_n{n}"
        mk = gbmap_markers(n_per_type=n)
        est = run_mcp(tag, mk, bulk_csv)
        if est.empty:
            print(f"  {tag}: FAILED"); continue
        est.to_csv(WORK / f"estimates_{tag}.csv")
        r = acs_score(masked(est, SHARED), man, method=tag, n_permutations=10000, n_boot=2000)
        s = r.summary()
        per_tumor[tag] = r.per_tumor.set_index("tumor_id")["acs"]
        print(f"  {tag:18s} {len(mk):5d} markers  ACS {s['acs']:.4f}  "
              f"CI [{s['ci_low']:.4f}, {s['ci_high']:.4f}]  p={s['null_p']:.4f}  "
              f"-> GBmap {'WINS' if s['acs'] > rs['acs'] else 'loses'}")
        out["gbmap_at_matched_density"][tag] = {
            "markers_per_type": n, "n_markers": int(len(mk)),
            "acs": round(float(s["acs"]), 4),
            "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
            "null_p": float(s["null_p"]),
            "beats_MCP_GBM": bool(s["acs"] > rs["acs"])}

    # --- PAIRED difference, bootstrapped over tumours --------------------------------
    print("\npaired difference, bootstrapped over the same tumours:")
    rng = np.random.default_rng(config.RANDOM_SEED)
    base = per_tumor["MCP_GBM"]
    for tag, other in per_tumor.items():
        if tag == "MCP_GBM":
            continue
        idx = base.index.intersection(other.index)
        d = (other.loc[idx] - base.loc[idx]).dropna()
        if len(d) < 3:
            print(f"  {tag}: only {len(d)} shared tumours; not testable"); continue
        boots = np.array([np.mean(rng.choice(d.to_numpy(), len(d), replace=True))
                          for _ in range(10000)])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        p = 2 * min((boots <= 0).mean(), (boots >= 0).mean())
        print(f"  {tag:18s} mean delta {d.mean():+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]  "
              f"p={p:.4f}  n={len(d)} tumours  "
              f"{'SIGNIFICANT' if lo > 0 or hi < 0 else 'not distinguishable'}")
        out.setdefault("paired_vs_MCP_GBM", {})[tag] = {
            "mean_delta": round(float(d.mean()), 4),
            "ci": [round(float(lo), 4), round(float(hi), 4)],
            "p": round(float(p), 4), "n_tumors": int(len(d)),
            "significant": bool(lo > 0 or hi < 0)}

    (config.RESULTS_DIR / "imc_anchored_robustness.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote results/imc_anchored_robustness.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
