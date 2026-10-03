"""
absolute_join_audit.py -- OPEN_DEFECTS D27: the ABSOLUTE join silently drops samples whose `sample`
field is a Broad-internal ID ("GBM-TCGA-06-5416-Tumor-SM-1QETM") rather than a TCGA barcode.

Found 2026-10-03 by an independent recomputation of GIMiCC's tumour control (tests/
test_gimicc_independent.py), which matched ABSOLUTE through its `array` column and found more samples.
`key4(sample)` turns a Broad ID into a key that matches nothing, so those rows vanish from every
ABSOLUTE comparison keyed on `sample`.

This AUDIT measures the effect; it changes no registered statistic.
  1. Reproduction control: the registered join reproduces every method's recorded purity rho exactly.
  2. Corrected join: standard rows keyed as registered (key4, vial-specific); Broad-ID rows keyed by
     their `array` barcode (TCGA-xx-xxxx-01), matched to the RNA sample with that barcode.
  3. Per-method purity rho under both joins, both cohorts; the registered agreement statistic (ACS vs
     GBM purity rho, 12 methods) recomputed with the corrected per-method values.
  Also reported: a vial-insensitive join (every row by `array`), as a sensitivity only -- the
  registered join requires the RNA and DNA vial letters to agree, which is a defensible convention.

    python3 scripts/absolute_join_audit.py
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402
from lymphoid_ordering import k4s  # noqa: E402

OUT = config.RESULTS_DIR / "absolute_join_audit.json"


def joins(cols: list[str], a: pd.DataFrame) -> dict[str, dict[str, float]]:
    """RNA column -> purity, under the registered, corrected and vial-insensitive joins."""
    std = a[a["sample"].astype(str).str.startswith("TCGA-")]
    broad = a[~a["sample"].astype(str).str.startswith("TCGA-")]
    reg = std.assign(k=std["sample"].map(key4)).drop_duplicates("k").set_index("k")[PURITY_COL]
    # the registered code keys ALL rows by key4(sample); Broad rows produce keys that match nothing,
    # so keying only the standard rows reproduces it exactly
    by_arr_broad = broad.drop_duplicates("array").set_index("array")[PURITY_COL]
    by_arr_all = a.drop_duplicates("array").set_index("array")[PURITY_COL]
    out = {"registered": {}, "corrected": {}, "vial_insensitive": {}}
    for c in cols:
        if key4(c) in reg.index:
            out["registered"][c] = float(reg[key4(c)])
            out["corrected"][c] = float(reg[key4(c)])
        elif k4s(c) in by_arr_broad.index:
            out["corrected"][c] = float(by_arr_broad[k4s(c)])
        if k4s(c) in by_arr_all.index:
            out["vial_insensitive"][c] = float(by_arr_all[k4s(c)])
    return out


def main() -> int:
    a = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    n_broad = int((~a["sample"].astype(str).str.startswith("TCGA-")).sum())
    out = {"what": __doc__.split("\n")[1], "absolute_rows_with_purity": int(len(a)),
           "rows_with_broad_internal_sample_id": n_broad, "cohorts": {}}
    ok_all = True
    for c, base in (("gbm", ""), ("lgg", "_lgg")):
        with gzip.open(config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz", "rt") as fh:
            cols = list(pd.read_csv(fh, index_col=0, nrows=1).columns)
        J = joins(cols, a)
        restored = sorted(set(J["corrected"]) - set(J["registered"]))
        full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{base}.csv")
        sid = "sample_id" if "sample_id" in full.columns else "sample"
        recorded = json.loads((config.RESULTS_DIR / f"absolute_purity_yardstick{base}.json").read_text())["methods"]
        per = {}
        for m, g in full.groupby("method"):
            t = g.set_index(sid)["Tumor"]
            row = {}
            for name, mp in J.items():
                idx = [s for s in t.index if s in mp]
                row[name] = {"n": len(idx),
                             "rho": round(float(stats.spearmanr(t[idx], [mp[s] for s in idx]).statistic), 4)}
            rec = (recorded.get(m) or {}).get("spearman_vs_purity")
            row["recorded"] = rec
            row["reproduces_recorded"] = rec is not None and abs(row["registered"]["rho"] - rec) < 1e-4
            ok_all &= bool(row["reproduces_recorded"]) if rec is not None else True
            per[m] = row
        out["cohorts"][c] = {"n_rna": len(cols), "n_registered": len(J["registered"]),
                             "n_corrected": len(J["corrected"]), "n_vial_insensitive": len(J["vial_insensitive"]),
                             "restored_samples": restored, "methods": per}
        print(f"{c}: registered n {len(J['registered'])} | corrected n {len(J['corrected'])} "
              f"(restored {restored}) | vial-insensitive n {len(J['vial_insensitive'])}")
        for m, r in per.items():
            print(f"   {m:18s} recorded {r['recorded']} | registered {r['registered']['rho']} "
                  f"| corrected {r['corrected']['rho']} | vial-insens {r['vial_insensitive']['rho']}"
                  f"{'' if r['reproduces_recorded'] else '  <-- DOES NOT REPRODUCE'}")
    out["reproduction_control_passes"] = ok_all
    # the registered agreement statistic, recomputed with corrected per-method GBM purity rho
    ya = json.loads((config.RESULTS_DIR / "yardstick_agreement.json").read_text())
    pm = ya["per_method"]
    ms = [m for m in pm if m in out["cohorts"]["gbm"]["methods"]]
    acs = [pm[m]["acs"] for m in ms]
    for name in ("registered", "corrected", "vial_insensitive"):
        rh = [out["cohorts"]["gbm"]["methods"][m][name]["rho"] for m in ms]
        r = stats.spearmanr(acs, rh)
        out[f"agreement_all_methods_{name}"] = {"n_methods": len(ms), "spearman": round(float(r.statistic), 4),
                                                "p": round(float(r.pvalue), 4)}
    out["agreement_recorded"] = ya["arm_2_absolute_purity"]["all_methods"]
    OUT.write_text(json.dumps(out, indent=2))
    print("reproduction control:", "PASSED" if ok_all else "FAILED")
    for name in ("recorded", "all_methods_registered", "all_methods_corrected", "all_methods_vial_insensitive"):
        print(f"agreement {name}: {out['agreement_' + name]}")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
