"""EXPLORATORY: does ACS predict accuracy once the post-registration methods are added?

The registered agreement test (ACS on Ivy GAP vs Spearman with ABSOLUTE purity on TCGA-GBM, frozen
signature) is fixed at the registered panel: n = 12, rho = +0.0810, INCONCLUSIVE. That result
stands. This script asks only what the same statistic looks like when the extension methods
(ivygap/deconv/extension.py) are added, using the identical pairing: each method's anatomic ACS
(raw/X) against its TCGA-GBM frozen-signature purity correlation.

CONTROL FIRST: the registered subset is recomputed from the shipped per-method values and must
reproduce the registered rho exactly; if it does not, nothing else here is reported.
Selects nothing; enters no registered statistic. Writes results/extension/agreement_extended.json.
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import spearmanr

from ivygap import config

R = config.RESULTS_DIR


def rho_with_ci(acs, acc, seed, n_perm=10_000, n_boot=5_000):
    acs, acc = np.asarray(acs, float), np.asarray(acc, float)
    rho = float(spearmanr(acs, acc).statistic)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(acs, rng.permutation(acc)).statistic for _ in range(n_perm)])
    p = (1 + int(np.sum(np.abs(null) >= abs(rho) - 1e-12))) / (n_perm + 1)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(acs), len(acs))
        if len(set(acs[i])) > 1 and len(set(acc[i])) > 1:
            boots.append(spearmanr(acs[i], acc[i]).statistic)
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"rho": round(rho, 4), "p_permutation": round(p, 6), "ci95": [round(float(lo), 4), round(float(hi), 4)],
            "n": int(len(acs))}


def main() -> int:
    ya = json.loads((R / "yardstick_agreement.json").read_text())
    reg = {m: v for m, v in ya["per_method"].items() if v.get("rho_purity") is not None}
    registered_rho = ya["arm_2_absolute_purity"]["all_methods"]
    reg_res = rho_with_ci([v["acs"] for v in reg.values()], [v["rho_purity"] for v in reg.values()], config.RANDOM_SEED)
    reported = registered_rho.get("spearman") if isinstance(registered_rho, dict) else None
    control_ok = reported is not None and abs(reg_res["rho"] - reported) < 1e-3
    out = {"what_this_is": __doc__.split("\n")[0], "registered_reproduced": control_ok,
           "registered": reg_res, "registered_reported_rho": reported}
    if not control_ok:
        out["blocked"] = "the registered subset did not reproduce the registered rho; nothing else reported"
        (R / "extension" / "agreement_extended.json").write_text(json.dumps(out, indent=2))
        print(out["blocked"], reg_res, reported); return 2
    tg = json.loads((R / "extension" / "tcga_gbm_frozen.json").read_text()).get("methods", {})
    ext = {}
    for f in sorted((R / "extension").glob("*_anatomic.json")):
        m = f.name.replace("_anatomic.json", ""); a = json.loads(f.read_text())
        rp = (tg.get(m) or {}).get("spearman_vs_purity")
        if a.get("acs") is not None and rp is not None:
            ext[m] = {"acs": a["acs"], "rho_purity": rp}
    out["extension_methods"] = ext
    if ext:
        allm = {**{m: v for m, v in reg.items()}, **ext}
        out["registered_plus_extension"] = rho_with_ci([v["acs"] for v in allm.values()],
                                                       [v["rho_purity"] for v in allm.values()], config.RANDOM_SEED)
    (R / "extension" / "agreement_extended.json").write_text(json.dumps(out, indent=2))
    print(f"registered (control): {reg_res} | reported {reported} | reproduced {control_ok}")
    print(f"extension methods with both values: {list(ext)}")
    if ext:
        print(f"registered + extension (EXPLORATORY): {out['registered_plus_extension']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
