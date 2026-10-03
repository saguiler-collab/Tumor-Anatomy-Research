"""
unmix_loss_ablation.py -- prespecified/unmix_loss_scale_ablation.md (post hoc; rule fixed first).

Re-runs the genuine DESeq2 `unmix` on the registered TCGA frozen arm with `shift` walked from the
pre-declared value toward the near-linear end, and with power 2 at shift 1, through environment
overrides read by R/run_deseq2_unmix.R. Control: (shift 1, power 1) must reproduce the reported row.

    python3 scripts/unmix_loss_ablation.py [--cohort gbm]   # writes results/unmix_loss_ablation_<cohort>.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ivygap import config  # noqa: E402
import extension_tcga  # noqa: E402

#: The first setting is the CONTROL: no override, so the pre-declared per-cohort rule picks shift
#: exactly as in the reported row. (A forced shift = 1 control was right for GBM, whose rule picks 1,
#: and wrong for LGG -- found 2026-10-02 when the LGG control did not reproduce.)
SETTINGS = [(None, None), (1, 1), (10, 1), (100, 1), (1000, 1), (10000, 1), (1, 2)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="gbm", choices=["gbm", "lgg"])
    a = ap.parse_args()
    rows = []
    for shift, power in SETTINGS:
        for k in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):
            os.environ.pop(k, None)
        if shift is not None:
            os.environ["IVYGAP_UNMIX_SHIFT"], os.environ["IVYGAP_UNMIX_POWER"] = str(shift), str(power)
        res = extension_tcga.run(a.cohort, "frozen", ["deseq2_unmix"])
        r = res["methods"]["deseq2_unmix"]
        r.update({"shift": shift, "power": power})
        rows.append(r)
        print(f"shift {str(shift if shift is not None else 'pre-declared'):>12} power {power}: "
              f"{r.get('lymphoid_ordering')} B>T {r.get('frac_samples_B_over_T')} "
              f"purity rho {r.get('spearman_vs_purity')}", flush=True)
    for k in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):
        os.environ.pop(k, None)
    main_row = json.loads((config.RESULTS_DIR / "extension" / f"tcga_{a.cohort}_frozen.json").read_text())[
        "methods"].get("deseq2_unmix", {})
    ctrl = rows[0]
    reproduced = all(ctrl.get(k) == main_row.get(k) for k in ("spearman_vs_purity", "frac_samples_B_over_T",
                                                              "lymphoid_ordering"))
    s1 = next(r for r in rows if r["shift"] == 1 and r["power"] == 1)
    lin = [r for r in rows if r["power"] == 1 and (r["shift"] or 0) >= 1000]
    reading = ("LOSS SCALE IMPLICATED" if s1.get("T_exceeds_B") and all(not r.get("T_exceeds_B") for r in lin)
               else "NOT IMPLICATED" if all(r.get("T_exceeds_B") for r in rows if r["power"] == 1)
               else "MIXED")
    out = {"rule": "prespecified/unmix_loss_scale_ablation.md", "cohort": a.cohort, "reference": "frozen",
           "control_reproduces_reported_row": reproduced, "reading": reading, "rows": rows}
    (config.RESULTS_DIR / f"unmix_loss_ablation_{a.cohort}.json").write_text(json.dumps(out, indent=2))
    print(f"control reproduced: {reproduced} | reading: {reading}")
    return 0 if reproduced else 2


if __name__ == "__main__":
    raise SystemExit(main())
