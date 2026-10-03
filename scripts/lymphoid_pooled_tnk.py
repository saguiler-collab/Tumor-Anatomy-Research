"""Does B exceed T and NK POOLED? A test of the "NK column carries T markers" suspect.

`docs/WHY_B_OVER_T.md` 4a names a structural suspect: the reference's NK column holds the
defining T-cell markers (CD3D/E/G, LCK), so real T-cell signal could be attributed to NK. If that
leak were the whole story, T + NK pooled would still exceed B even where T alone loses. If methods
place B above T + NK, the leak cannot explain the inversion: the mass lands on B, not merely
between T and NK. Summing two columns is the lenient direction -- any T<->NK split is undone.

Reads existing estimates only (registered arms and, where present, the extension panel).
Methylation is used at the cohort level (see results/lymphoid_tracking.json). Selects nothing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lymphoid_ordering import _truth, k4s, renormalise  # noqa: E402

from ivygap import config  # noqa: E402

R = config.RESULTS_DIR


def _load(p: Path) -> pd.DataFrame:
    f = pd.read_csv(p)
    scol = next(c for c in ("sample", "sample_id") if c in f.columns)
    return f.rename(columns={scol: "sample"})


def main() -> int:
    out = {"what_this_is": __doc__.split("\n")[0], "arms": {}}
    for c in ("gbm", "lgg"):
        base = "" if c == "gbm" else "_lgg"
        tr = _truth(R / f"methylation_celltypes{base}.csv").dropna()
        for ref in ("frozen", "h5ad"):
            tag = base + ("_h5ad" if ref == "h5ad" else "")
            frames = [_load(R / f"estimates_full{tag}.csv")]
            ext = R / "extension" / f"estimates_full{tag}_extension.csv"
            if ext.exists():
                frames.append(_load(ext))
            full = pd.concat(frames, ignore_index=True)
            full["k"] = full["sample"].map(k4s)
            rows = {}
            for m, g in full.groupby("method"):
                g = g.drop_duplicates("k").set_index("k")
                idx = [k for k in g.index if k in tr.index]
                ren = renormalise(g.loc[idx]).dropna()
                if len(ren) < 20:
                    continue
                rows[m] = {"n": int(len(ren)),
                           "frac_B_over_T": round(float((ren.B_cell > ren.T_cell).mean()), 4),
                           "frac_B_over_T_plus_NK": round(float((ren.B_cell > ren.T_cell + ren.NK_cell).mean()), 4),
                           "mean_B_share": round(float(ren.B_cell.mean()), 4),
                           "extension_panel": ext.exists() and m in set(_load(ext)["method"])}
            maj = [m for m, v in rows.items() if v["frac_B_over_T_plus_NK"] > 0.5 and not v["extension_panel"]]
            out["arms"][f"{c}_{ref}"] = {
                "methylation_frac_B_over_T_plus_NK": round(float((tr.B_cell > tr.T_cell + tr.NK_cell).mean()), 4),
                "methylation_mean_B_share": round(float(tr.B_cell.mean()), 4),
                "registered_methods_with_B_over_T_plus_NK_in_majority": sorted(maj),
                "n_registered_methods": sum(1 for v in rows.values() if not v["extension_panel"]),
                "methods": dict(sorted(rows.items(), key=lambda kv: -kv[1]["frac_B_over_T_plus_NK"]))}
    (R / "lymphoid_pooled_tnk.json").write_text(json.dumps(out, indent=2))
    for arm, a in out["arms"].items():
        print(f"\n=== {arm}: methylation B > T+NK in {a['methylation_frac_B_over_T_plus_NK']:.1%} of samples | "
              f"registered methods with B > T+NK in most samples: {len(a['registered_methods_with_B_over_T_plus_NK_in_majority'])} "
              f"of {a['n_registered_methods']}")
        for m, v in a["methods"].items():
            print(f"   {m:<22}{v['n']:>5}  B>T {v['frac_B_over_T']:>6.1%}  B>T+NK {v['frac_B_over_T_plus_NK']:>6.1%}  "
                  f"mean B {v['mean_B_share']:.3f}{'  [extension]' if v['extension_panel'] else ''}")
    print("\nwrote results/lymphoid_pooled_tnk.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
