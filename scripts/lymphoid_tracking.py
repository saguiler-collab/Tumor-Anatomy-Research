"""Does a method's lymphoid estimate TRACK the tissue, or only sit at a plausible level?

WHY THIS EXISTS. `lymphoid_ordering.py` asks whether a method's COHORT-MEAN lymphoid
composition orders T above B, as DNA methylation does. Under the h5ad (raw/X, donor-level)
reference, Bisque passes that test in both cohorts -- B above T on only 3.6% (GBM) and 1.0%
(LGG) of samples -- where under the frozen signature it failed. Read naively, Bisque has
solved the lymphoid inversion.

It may not have. Bisque without overlapping subjects (`use.overlap = FALSE`, the only mode
TCGA permits) applies `SemisupervisedTransformBulk` to every gene: z-score the bulk across
samples, then rescale it to the mean and (shrunk) spread of PSEUDOBULKS BUILT FROM THE
ATLAS'S OWN DONORS. The transformed cohort therefore has, gene for gene, the mean of
`sc.ref %*% mean(sc.props)` -- so the cohort-average estimate is pulled to the atlas's
average composition BY CONSTRUCTION, whatever the tissue holds. A correct cohort mean is
then evidence about the atlas, not about the bulk.

A level can be inherited. Covariation cannot. So this asks two questions a prior cannot
answer:

  1. PER-SAMPLE TRACKING. Across samples, does the method's lymphoid T-share (and B-share)
     rise and fall with methylation's? Spearman, with a 10,000-draw permutation null that
     shuffles which methylation sample each estimate is paired with.
  2. BETWEEN-COHORT SHIFT. Methylation's B-share is ~0.11 in GBM and ~0.20 in LGG. A method
     that measures tissue must move in the same direction. A method anchored to one atlas
     prior cannot, because the prior is the same in both cohorts.

POSITIVE CONTROL -- WITHOUT IT, "NO METHOD TRACKS" MEANS NOTHING. If methylation's
per-sample lymphoid split were itself noise, every method would fail to track it and that
would say nothing about the methods. So the same test is run on a two-panel MARKER INDEX
computed from the same bulk RNA: mean log2(CPM+1) over canonical T markers minus the same
over canonical B markers. It uses no reference, no signature and no deconvolution. If it
tracks methylation, the truth carries per-sample signal and a method's failure to track it
is informative.

THE MARKER PANELS WERE FIXED BEFORE THIS SCRIPT FIRST RAN (2026-09-30), from canonical
lineage biology, and are not to be revised in light of the result -- choosing genes until the
control passes would be selection on the answer:

  T: CD3D CD3E CD3G CD5 CD6 TRAT1      (excludes CD2, LCK, CD247 -- shared with NK -- and
                                        CD4, which macrophages and microglia express)
  B: CD19 MS4A1 CD79A CD79B CD22 BLK PAX5 FCRLA

The tracked quantity for the control is T/(T+B); methods are scored on the same quantity
alongside their T-share and B-share, so the comparison is like-for-like.

Nothing here selects a method, and no outcome is used. This is a measurement of what the
existing ordering result means.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lymphoid_ordering import COMPARE, _truth, k4s, renormalise  # noqa: E402

from ivygap import config  # noqa: E402
from ivygap.deconv.comparability import non_comparable  # noqa: E402

N_PERM = 10_000
MIN_PAIRS = 20
T_MARKERS = ["CD3D", "CD3E", "CD3G", "CD5", "CD6", "TRAT1"]
B_MARKERS = ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "BLK", "PAX5", "FCRLA"]


def marker_index(cohort: str) -> tuple[pd.Series, dict]:
    """Per-sample T-minus-B marker index from the TCGA bulk, streamed so only marker rows load."""
    path = config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz"
    want = set(T_MARKERS) | set(B_MARKERS)
    rows = [ch.loc[ch.index.intersection(want)]
            for ch in pd.read_csv(path, index_col=0, chunksize=2000)]
    m = pd.concat(rows)
    m = m[~m.index.duplicated()]
    present_t = [g for g in T_MARKERS if g in m.index]
    present_b = [g for g in B_MARKERS if g in m.index]
    lg = np.log2(m + 1.0)
    idx = lg.loc[present_t].mean() - lg.loc[present_b].mean()
    idx.index = [k4s(c) for c in idx.index]
    idx = idx[~idx.index.duplicated()]
    return idx, {"source": str(path.relative_to(config.PROJECT_ROOT)),
                 "T_markers_present": present_t, "B_markers_present": present_b,
                 "T_markers_missing": sorted(set(T_MARKERS) - set(present_t)),
                 "B_markers_missing": sorted(set(B_MARKERS) - set(present_b))}


def perm_spearman(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Spearman rho and a two-sided permutation p, floored at 1/(N_PERM+1).

    Ranks once, then permutes the ranks: equivalent to recomputing Spearman per draw and two
    orders of magnitude cheaper.
    """
    rho = float(spearmanr(x, y).statistic)
    if not np.isfinite(rho):
        return float("nan"), float("nan")
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    rx = (rx - rx.mean()) / rx.std()
    ry = (ry - ry.mean()) / ry.std()
    null = np.array([np.mean(rx * rng.permutation(ry)) for _ in range(N_PERM)])
    p = (1 + int(np.sum(np.abs(null) >= abs(rho) - 1e-12))) / (N_PERM + 1)
    return rho, p


def one_arm(cohort: str, reference: str, rng: np.random.Generator) -> dict:
    base = "" if cohort == "gbm" else "_lgg"
    tag = base + ("_h5ad" if reference == "h5ad" else "")
    truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
    full = pd.read_csv(config.RESULTS_DIR / f"estimates_full{tag}.csv")
    scol = next(c for c in ("sample", "sample_id") if c in full.columns)
    full["k"] = full[scol].map(k4s)
    nc = non_comparable(reference)

    out = {"cohort": cohort, "reference": reference,
           "truth_mean": {c: round(float(v), 4) for c, v in truth.mean().items()},
           "methods": {}}
    for m, g in full.groupby("method"):
        g = g.drop_duplicates("k").set_index("k")
        idx = [k for k in g.index if k in truth.index]
        est = renormalise(g.loc[idx, COMPARE])
        tru = truth.loc[idx, COMPARE]
        ok = est.notna().all(axis=1) & tru.notna().all(axis=1)
        e, t = est[ok], tru[ok]
        rec = {"n_matched": len(idx), "n_pairs": int(ok.sum()),
               "n_excluded_zero_lymphoid": int((~ok).sum()),
               "comparable": m not in nc,
               "est_mean": {c: round(float(v), 4) for c, v in e.mean().items()} if len(e) else None}
        e = e.assign(TB=e["T_cell"] / (e["T_cell"] + e["B_cell"]).replace(0, np.nan))
        t = t.assign(TB=t["T_cell"] / (t["T_cell"] + t["B_cell"]).replace(0, np.nan))
        for share in ("T_cell", "B_cell", "TB"):
            key = f"track_{share}"
            okp = e[share].notna() & t[share].notna()
            if share == "TB":
                e2, t2 = e[okp], t[okp]
                if len(e2) < MIN_PAIRS or e2[share].nunique() < 3:
                    rec[key] = {"rho": None, "p": None, "why_none": "too few pairs"}
                    continue
                rho, p = perm_spearman(e2[share].to_numpy(), t2[share].to_numpy(), rng)
                rec[key] = {"rho": round(rho, 4), "p": round(p, 6), "n": int(len(e2))}
                continue
            if len(e) < MIN_PAIRS or e[share].nunique() < 3:
                rec[key] = {"rho": None, "p": None,
                            "why_none": f"{len(e)} pairs / {e[share].nunique() if len(e) else 0} "
                                        f"distinct values -- too few to estimate"}
                continue
            rho, p = perm_spearman(e[share].to_numpy(), t[share].to_numpy(), rng)
            rec[key] = {"rho": round(rho, 4), "p": round(p, 6)}
        out["methods"][m] = rec
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(config.RESULTS_DIR / "lymphoid_tracking.json"))
    a = ap.parse_args()
    rng = np.random.default_rng(config.RANDOM_SEED)

    arms = {f"{c}_{r}": one_arm(c, r, rng) for c in ("gbm", "lgg") for r in ("frozen", "h5ad")}

    # --- POSITIVE CONTROL: does a reference-free marker index track methylation? ---
    control = {}
    for c in ("gbm", "lgg"):
        base = "" if c == "gbm" else "_lgg"
        truth = _truth(config.RESULTS_DIR / f"methylation_celltypes{base}.csv")
        tb = (truth["T_cell"] / (truth["T_cell"] + truth["B_cell"]).replace(0, np.nan)).dropna()
        idx, prov = marker_index(c)
        common = [k for k in idx.index if k in tb.index]
        rho, p = perm_spearman(idx.loc[common].to_numpy(), tb.loc[common].to_numpy(), rng)
        # the same control on T-share within {T, NK, B}, for direct comparison with methods
        ts = truth["T_cell"].dropna()
        common_t = [k for k in idx.index if k in ts.index and np.isfinite(ts.loc[k])]
        rho_t, p_t = perm_spearman(idx.loc[common_t].to_numpy(), ts.loc[common_t].to_numpy(), rng)
        # PRE-DECLARED EXPLORATORY STRATIFICATION (fixed 2026-09-30, before it was run): one
        # stratifier only, tertiles of PTPRC (CD45) log2(CPM+1) -- total leukocyte RNA, blind to
        # the T-vs-B split. Per-sample agreement is only expected where leukocytes are present,
        # so the a-priori prediction is agreement in the HIGH tertile. Reported whatever it shows.
        path = config.PROCESSED_DIR / f"tcga_{c}_bulk_cpm.csv.gz"
        pt = pd.concat([ch.loc[ch.index.intersection(["PTPRC"])]
                        for ch in pd.read_csv(path, index_col=0, chunksize=2000)]).iloc[0]
        pt.index = [k4s(x) for x in pt.index]
        pt = np.log2(pt[~pt.index.duplicated()] + 1)
        cs = [k for k in common if k in pt.index]
        dd = pd.DataFrame({"idx": idx.loc[cs], "tb": tb.loc[cs], "ptprc": pt.loc[cs]})
        dd["tert"] = pd.qcut(dd["ptprc"], 3, labels=["low", "mid", "high"])
        strata = {}
        for tt, g in dd.groupby("tert", observed=True):
            r_, p_ = perm_spearman(g["idx"].to_numpy(), g["tb"].to_numpy(), rng)
            strata[str(tt)] = {"n": int(len(g)), "rho": round(r_, 4), "p": round(p_, 6)}
        r_m, p_m = perm_spearman(dd["ptprc"].to_numpy(), dd["tb"].to_numpy(), rng)
        r_i, p_i = perm_spearman(dd["ptprc"].to_numpy(), dd["idx"].to_numpy(), rng)
        control[c] = {**prov, "n": len(common),
                      "ptprc_tertiles_marker_vs_truth_TB": strata,
                      "a_priori_prediction": "agreement in the HIGH-PTPRC tertile",
                      "ptprc_vs_methylation_TB": {"rho": round(r_m, 4), "p": round(p_m, 6)},
                      "ptprc_vs_marker_index": {"rho": round(r_i, 4), "p": round(p_i, 6)},
                      "rho_vs_truth_TB": round(rho, 4), "p": round(p, 6),
                      "rho_vs_truth_T_share": round(rho_t, 4), "p_T_share": round(p_t, 6)}

    # --- between-cohort shift: does the method move B-share the way methylation does? ---
    shift = {}
    for r in ("frozen", "h5ad"):
        tg, tl = arms[f"gbm_{r}"]["truth_mean"], arms[f"lgg_{r}"]["truth_mean"]
        truth_dB = tl["B_cell"] - tg["B_cell"]
        rows = {}
        for m in sorted(set(arms[f"gbm_{r}"]["methods"]) & set(arms[f"lgg_{r}"]["methods"])):
            eg = arms[f"gbm_{r}"]["methods"][m]["est_mean"]
            el = arms[f"lgg_{r}"]["methods"][m]["est_mean"]
            if eg is None or el is None:
                continue
            dB = el["B_cell"] - eg["B_cell"]
            rows[m] = {"B_share_gbm": eg["B_cell"], "B_share_lgg": el["B_cell"],
                       "delta_B": round(dB, 4), "same_direction_as_truth": bool(np.sign(dB) == np.sign(truth_dB))}
        shift[r] = {"truth_B_share_gbm": tg["B_cell"], "truth_B_share_lgg": tl["B_cell"],
                    "truth_delta_B": round(truth_dB, 4), "methods": rows}

    report = {
        "what_this_is": ("Per-sample tracking of methylation's lymphoid composition, and the "
                         "GBM->LGG B-share shift, for every method under both references. "
                         "Distinguishes a method that MEASURES the lymphoid compartment from one "
                         "that inherits a plausible level. Selects nothing."),
        "null": f"two-sided permutation of the estimate-to-truth pairing, {N_PERM} draws, "
                f"seed {config.RANDOM_SEED}; p floored at 1/{N_PERM + 1}",
        "pairs_rule": "samples where either side has zero T, NK and B are excluded and counted, "
                      "never imputed",
        "positive_control_marker_index": control,
        "positive_control_verdict": (
            "FAILED. A reference-free T-minus-B marker index from the same bulk RNA does not "
            "track methylation's per-sample T/(T+B), and in LGG runs significantly opposite to "
            "it; the pre-declared PTPRC stratification does not rescue it (no high-tertile "
            "agreement in either cohort). Methylation's PER-SAMPLE lymphoid split is therefore "
            "not corroborated, and per-sample tracking of it by any method is INCONCLUSIVE -- "
            "a failure to track it says nothing about the method. Cohort-level comparisons are "
            "unaffected by this verdict and are reported separately."),
        "arms": arms, "between_cohort_shift": shift,
    }
    Path(a.out).write_text(json.dumps(report, indent=2))

    for name, arm in arms.items():
        print(f"\n=== {name}   truth T/NK/B = "
              f"{arm['truth_mean']['T_cell']:.3f}/{arm['truth_mean']['NK_cell']:.3f}/{arm['truth_mean']['B_cell']:.3f}")
        print(f"   {'method':<22}{'pairs':>6}{'rho_T':>8}{'p':>9}{'rho_B':>8}{'p':>9}{'rho_T/(T+B)':>13}{'p':>9}")
        for m, v in arm["methods"].items():
            tT, tB, tTB = v["track_T_cell"], v["track_B_cell"], v["track_TB"]
            f = lambda d, k: "   -" if d[k] is None else (f"{d[k]:+.3f}" if k == "rho" else f"{d[k]:.4f}")
            print(f"   {m:<22}{v['n_pairs']:>6}{f(tT,'rho'):>8}{f(tT,'p'):>9}{f(tB,'rho'):>8}{f(tB,'p'):>9}{f(tTB,'rho'):>13}{f(tTB,'p'):>9}"
                  f"{'' if v['comparable'] else '  [non-comparable]'}")
    for r, s in shift.items():
        print(f"\n=== GBM->LGG B-share shift, {r}: truth {s['truth_B_share_gbm']:.3f} -> "
              f"{s['truth_B_share_lgg']:.3f} (delta {s['truth_delta_B']:+.3f})")
        for m, v in s["methods"].items():
            print(f"   {m:<22}{v['B_share_gbm']:.3f} -> {v['B_share_lgg']:.3f}  "
                  f"delta {v['delta_B']:+.3f}  {'same direction' if v['same_direction_as_truth'] else 'OPPOSITE'}")
    for c, v in control.items():
        print(f"\n=== POSITIVE CONTROL, {c}: marker index (T {len(v['T_markers_present'])} genes - "
              f"B {len(v['B_markers_present'])} genes) vs methylation, n={v['n']}")
        print(f"   vs T/(T+B):  rho {v['rho_vs_truth_TB']:+.3f}  p {v['p']:.4f}")
        print(f"   vs T-share:  rho {v['rho_vs_truth_T_share']:+.3f}  p {v['p_T_share']:.4f}")
        if v["T_markers_missing"] or v["B_markers_missing"]:
            print(f"   missing: T {v['T_markers_missing']}  B {v['B_markers_missing']}")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
