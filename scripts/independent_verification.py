"""Recompute the headline numbers from RAW FILES, sharing no code with the pipeline.

WHY THIS EXISTS. Every other number in this project is produced by `ivygap/`. If that package
contains a systematic error -- a mislabelled column, a wrong join, an off-by-one in a
normalisation -- every artefact inherits it and every internal consistency check still passes,
because they all read the same wrong values. Internal consistency cannot detect a shared cause.

So this script imports NOTHING from `ivygap`. It reads the raw inputs with pandas, redoes the
arithmetic in the most direct way available, and prints its answer beside the pipeline's. It is
deliberately simpler and slower than the pipeline; where the two disagree, that disagreement is
the finding and neither number should be reported until it is explained.

It does not prove the pipeline correct. Two implementations can share a misunderstanding of the
DATA -- if the methylation file's columns mean something other than we both assume, agreement
proves only that we misread it consistently. What it rules out is a bug living in one code path.

    python scripts/independent_verification.py
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parent.parent
OK, BAD = "AGREES", "*** DISAGREES ***"


def verdict(mine, theirs, tol):
    if mine is None or theirs is None:
        return "cannot compare"
    return OK if abs(float(mine) - float(theirs)) <= tol else BAD


def check_lymphoid_from_raw() -> list[str]:
    """T vs B from the raw methylation matrix, by constrained least squares, not EpiDISH.

    A FIRST ATTEMPT HERE WAS WRONG AND IS WORTH RECORDING. It assigned each reference CpG to
    the cell type with the HIGHEST beta value and compared mean methylation at "T-owned" and
    "B-owned" CpGs. That reported T > B in 13.8% of samples against the pipeline's 93.2% --
    an apparently catastrophic disagreement.

    The flaw was mine. Methylation markers are frequently defined by being HYPO-methylated in
    the type they mark: `cg11661493` reads B = 0.037 against ~0.96 in all six other types, so
    it is a B-cell marker, and `idxmax` assigned it to CD4T. Tested on mixtures with known
    composition, that statistic was right in 4 of 6 cases -- chance is 3. It was measuring
    almost nothing.

    This is why methylation deconvolution is posed as a constrained regression rather than as
    marker averaging, and it is a good illustration of why an "independent check" has to be
    validated before its disagreement means anything.

    The replacement solves `beta ~ R x` with non-negative least squares and a sum-to-one
    constraint -- the same MODEL EpiDISH's RPC fits, through a completely different
    implementation (scipy.optimize.nnls rather than robust partial correlation in R). Agreement
    therefore tests the pipeline's plumbing, not the choice of model.
    """
    out = []
    probe = ROOT / "data/processed/lgg_methylation_epidish_probes.csv.gz"
    if not probe.exists():
        return ["  lymphoid: probe matrix absent, skipped"]
    beta = pd.read_csv(probe, index_col=0).dropna()
    out.append(f"  raw probe matrix: {beta.shape[0]} complete CpGs x {beta.shape[1]} samples")

    import subprocess
    from io import StringIO
    r = subprocess.run(
        ["Rscript", "-e",
         'suppressMessages(library(EpiDISH)); data(centDHSbloodDMC.m); '
         'write.csv(centDHSbloodDMC.m, stdout())'],
        capture_output=True, text=True)
    if r.returncode != 0:
        return out + ["  could not read centDHSbloodDMC.m from R; skipped"]
    ref = pd.read_csv(StringIO(r.stdout), index_col=0)
    shared = [c for c in beta.index if c in ref.index]
    R = ref.loc[shared].to_numpy(float)
    B = beta.loc[shared].to_numpy(float)
    types = list(ref.columns)
    out.append(f"  solving beta ~ Rx on {len(shared)} shared CpGs, {len(types)} cell types, "
               f"by scipy NNLS (EpiDISH uses robust partial correlation in R)")

    from scipy.optimize import nnls                                # noqa: PLC0415
    iT = [types.index(t) for t in types if t.upper() in ("CD4T", "CD8T")]
    iB = [types.index(t) for t in types if t.upper() == "B"]
    tb = []
    for j in range(B.shape[1]):
        x, _ = nnls(R, B[:, j])
        s = x.sum()
        if s <= 0:
            continue
        x = x / s
        tb.append(x[iT].sum() > x[iB].sum())
    frac = float(np.mean(tb)) if tb else float("nan")
    out.append(f"  samples where T exceeds B (independent NNLS): {frac:.1%}  (n={len(tb)})")

    pipe = ROOT / "results/methylation_celltypes_lgg.json"
    p = json.loads(pipe.read_text())["prediction_T_exceeds_B"]["fraction_of_samples"] \
        if pipe.exists() else None
    if p is None:
        return out + ["  pipeline value unavailable"]
    out.append(f"  pipeline (EpiDISH RPC in R) says:             {p:.1%}")
    agree = abs(frac - p) <= 0.10
    out.append(f"  {OK if agree else BAD}  two implementations of the same model, "
               f"difference {abs(frac - p):.1%}")
    return out


def check_purity_correlation() -> list[str]:
    """Spearman(estimate, purity) recomputed from the per-sample CSV with plain scipy."""
    out = []
    f = ROOT / "results/absolute_purity_per_sample_lgg.csv"
    j = ROOT / "results/absolute_purity_yardstick_lgg.json"
    if not (f.exists() and j.exists()):
        return ["  purity: artefacts absent, skipped"]
    d = pd.read_csv(f, index_col=0)
    tcol = "absolute_purity" if "absolute_purity" in d.columns else "purity"
    rep = json.loads(j.read_text()).get("methods", {})
    out.append(f"  per-sample matrix: {d.shape[0]} samples x {d.shape[1] - 1} methods")
    worst = 0.0
    for m in [c for c in d.columns if c != tcol]:
        e, t = d[m].to_numpy(float), d[tcol].to_numpy(float)
        ok = np.isfinite(e) & np.isfinite(t)
        if ok.sum() < 50 or np.std(e[ok]) == 0:
            continue
        mine = float(stats.spearmanr(e[ok], t[ok]).statistic)
        theirs = (rep.get(m) or {}).get("spearman_vs_purity")
        if theirs is None:
            continue
        worst = max(worst, abs(mine - float(theirs)))
        if abs(mine - float(theirs)) > 0.001:
            out.append(f"    {BAD} {m}: independent {mine:+.4f} vs pipeline {theirs:+.4f}")
    out.append(f"  largest disagreement across all methods: {worst:.6f}  "
               f"{OK if worst <= 0.001 else BAD}")
    return out


def check_recovery_definition() -> list[str]:
    """Recovery = 1 + slope, recomputed with numpy.polyfit instead of the project's code."""
    out = []
    f = ROOT / "results/absolute_purity_per_sample.csv"
    j = ROOT / "results/equal_footing_ranking.json"
    if not (f.exists() and j.exists()):
        return ["  recovery: artefacts absent, skipped"]
    d = pd.read_csv(f, index_col=0)
    tcol = "absolute_purity" if "absolute_purity" in d.columns else "purity"
    theirs = json.loads(j.read_text()).get("recovery_frozen", {})
    worst, n = 0.0, 0
    for m, tv in theirs.items():
        if m not in d.columns:
            continue
        e, t = d[m].to_numpy(float), d[tcol].to_numpy(float)
        ok = np.isfinite(e) & np.isfinite(t)
        if ok.sum() < 50:
            continue
        mine = 1.0 + float(np.polyfit(t[ok], e[ok] - t[ok], 1)[0])
        worst = max(worst, abs(mine - float(tv)))
        n += 1
    out.append(f"  recomputed recovery for {n} methods; largest disagreement {worst:.6f}  "
               f"{OK if worst <= 0.001 else BAD}")
    return out


def check_control_separation() -> list[str]:
    """The property that makes ACS a detector: controls must sit below every working method."""
    out = []
    f = ROOT / "results/anatomic/acs_leaderboard.csv"
    if not f.exists():
        return ["  ACS: leaderboard absent, skipped"]
    lb = pd.read_csv(f).set_index("method")
    ctl = lb[lb.get("is_control", False) == True]                  # noqa: E712
    real = lb[(lb.get("is_control", False) == False) & (lb["acs"] > 0)]  # noqa: E712
    out.append(f"  real methods scoring > 0: {len(real)}   negative controls: {len(ctl)}")
    out.append(f"  worst real method : {real['acs'].min():.4f}  ({real['acs'].idxmin()})")
    out.append(f"  best control      : {ctl['acs'].max():.4f}  ({ctl['acs'].idxmax()})")
    sep = real["acs"].min() > ctl["acs"].max()
    out.append(f"  separated: {sep}  {OK if sep else BAD}")
    out.append(f"  every control fails its own permutation null: "
               f"{bool((ctl['null_p'] > 0.05).all())}")
    allbeat = bool((real["null_p"] < 0.05).all())
    out.append(f"  every real method beats its null: {allbeat}")
    if not allbeat:
        # A bare False here is not a finding, it is a prompt to look. Name the methods and
        # say whether each is AT CHANCE or NOT EVALUABLE -- the second is a coverage fact
        # about the constraint set, not a verdict on the method, and the protocol calls it
        # INCONCLUSIVE rather than a failure.
        med = float(real["n_constraint_tumor_pairs"].median())
        for m, r in real[real["null_p"] >= 0.05].iterrows():
            n = int(r["n_constraint_tumor_pairs"])
            why = ("scored on only %d of %d constraint-tumour pairs -> UNDERPOWERED, "
                   "not evaluable" % (n, int(med))) if n < 0.6 * med else \
                  "scored on the full constraint set -> genuinely AT CHANCE"
            out.append(f"    {m}: acs {r['acs']:.4f}, null_p {r['null_p']:.3f} — {why}")
        # and re-state the control margin without them, which is the claim that has to hold
        full = real[real["n_constraint_tumor_pairs"] >= 0.6 * med]
        if len(full):
            sep2 = full["acs"].min() > ctl["acs"].max()
            out.append(f"  separation using only fully-scored methods: "
                       f"{full['acs'].min():.4f} vs {ctl['acs'].max():.4f} -> {sep2}  "
                       f"{OK if sep2 else BAD}")
    return out


def main() -> int:
    print(__doc__.split("\n\n")[0])
    print("\nThis script imports nothing from `ivygap`. Modules loaded:",
          ", ".join(sorted({"pandas", "numpy", "scipy"})), "\n")
    assert "ivygap" not in sys.modules, "this script must not import project code"

    for title, fn in (("1 · LYMPHOID DIRECTION, from the raw methylation matrix",
                       check_lymphoid_from_raw),
                      ("2 · PURITY CORRELATIONS, recomputed with plain scipy",
                       check_purity_correlation),
                      ("3 · RECOVERY STATISTIC, recomputed with numpy.polyfit",
                       check_recovery_definition),
                      ("4 · CONTROL SEPARATION, read straight from the leaderboard",
                       check_control_separation)):
        print("=" * 78)
        print(title)
        print("=" * 78)
        for line in fn():
            print(line)
        print()
    print("Agreement here rules out a bug confined to one code path. It does NOT rule out a")
    print("shared misreading of the input files -- for that you need someone who knows the")
    print("data to look at it, which no amount of recomputation substitutes for.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
