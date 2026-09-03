"""
agreement.py — THE experiment.

The leaderboard is not the result. A method could satisfy every anatomic constraint and
still be badly calibrated, so ACS is only useful if it *agrees with the truth where the
truth is known*. This module measures that agreement and reports it against a bar fixed
in advance.

    primary result = Spearman rho between
                       (methods ranked by ACS, computed with no ground truth)
                   and (methods ranked by accuracy against real ground truth)

    pre-registered bar: rho >= 0.60, with a bootstrap CI excluding zero

THE THREE YARDSTICKS
--------------------
Ground truth is not one thing, and a rank correlation against a single yardstick could be
an artefact of that yardstick. The protocol names three, each with a different failure
mode, and agreement is reported against each separately as well as pooled:

  synthetic_mixtures  500 pseudobulk mixtures with known composition. Strongest signal,
                      but partly circular: the mixtures are built from the same reference
                      the methods solve against. Its weakness is the reason for the
                      other two.
  absolute_purity     ABSOLUTE DNA-based tumour purity on TCGA-GBM. Orthogonal — it never
                      touches RNA — but it constrains exactly one cell type.
  sc_pseudobulk       Pseudobulk from paired single-cell GBM data, where composition is
                      counted rather than estimated. Independent of the reference if the
                      donors are held out.

BOTH DIRECTIONS ARE RESULTS, AND THE CODE TREATS THEM THAT WAY
--------------------------------------------------------------
`interpret()` returns a verdict for each of the three pre-specified outcomes, including
the null. It never returns "inconclusive, try again" and it never rewrites the bar. A rho
below the threshold is the finding that "matches known biology" is not evidence of
correctness — which is the more useful paper, and the code is written so that outcome is
as easy to report as the other one.

WHAT WOULD INVALIDATE THIS
--------------------------
Ranking on fewer than about six methods makes Spearman almost meaningless — with four
methods the coarsest possible rho values are far apart and the CI spans everything. The
minimum is enforced rather than noted, and controls are excluded from the correlation by
default: they are there to test the constraint set, and leaving them in would inflate rho
by adding two points that every yardstick agrees are terrible.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import stats

from ivygap import config
from ivygap.deconv.controls import CONTROL_NAMES

#: Pre-registered before any agreement number was computed. Do not edit to fit a result.
PREREGISTERED_RHO_THRESHOLD = 0.60
#: Below this many methods a rank correlation is not interpretable.
MIN_METHODS_FOR_CORRELATION = 6

#: The protocol's three yardsticks, in the order it lists them.
YARDSTICKS = ("synthetic_mixtures", "absolute_purity", "sc_pseudobulk")

#: A fourth, clearly separate yardstick used when none of the three is available for
#: enough methods to support a rank correlation. Mixtures are generated from
#: donor-perturbed copies of the frozen signature rather than pooled from real single
#: cells, so it is a simulation. It is named differently, reported separately, and
#: flagged in every artefact — a preliminary reading obtained this way must never be
#: presented as the protocol's yardstick 1.
SIMULATED_YARDSTICK = "simulated_donor_mismatch"
ALL_YARDSTICKS = YARDSTICKS + (SIMULATED_YARDSTICK,)


@dataclass
class AgreementResult:
    yardstick: str
    rho: float
    p_value: float
    ci_low: float
    ci_high: float
    n_methods: int
    #: Methods with a DISTINCT (acs, truth) pair. Degenerate methods — MuSiC on a
    #: collapsed reference is arithmetically NNLS; SCDC ENSEMBLE with one reference is
    #: SCDC — produce byte-identical rows. Counting them twice inflates `n_methods`
    #: without adding an independent observation, and Spearman on duplicated points
    #: understates the sampling variability. The correlation is computed on all rows
    #: (they are real methods and were really run) but this number is what a reader
    #: should judge the correlation's weight by.
    n_distinct: int
    duplicate_groups: list
    meets_threshold: bool
    ci_excludes_zero: bool
    verdict: str
    method_ranks: dict

    def to_dict(self) -> dict:
        return asdict(self)


def _bootstrap_rho(acs: np.ndarray, truth: np.ndarray, n_boot: int,
                   seed: int) -> tuple[float, float]:
    """
    Percentile bootstrap over METHODS.

    Resampling methods is the right unit: the question is whether the ACS ordering of
    *this set of methods* tracks their true ordering, and the sampling variability that
    matters is which methods happened to be in the comparison.
    """
    rng = np.random.default_rng(seed)
    n = len(acs)
    draws = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        a, t = acs[idx], truth[idx]
        # A resample that picks one method repeatedly has no spread and no correlation
        # to measure; skipping it is correct, and how many were skipped is not hidden.
        if len(np.unique(a)) < 2 or len(np.unique(t)) < 2:
            continue
        draws.append(stats.spearmanr(a, t)[0])
    draws = np.array([d for d in draws if np.isfinite(d)])
    if draws.size < 100:
        return (float("nan"), float("nan"))
    return (float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)))


def interpret(rho: float, ci_low: float, ci_high: float, n_methods: int) -> str:
    """The three pre-specified outcomes, and nothing else."""
    if n_methods < MIN_METHODS_FOR_CORRELATION:
        return (f"NOT INTERPRETABLE: {n_methods} methods is below the "
                f"{MIN_METHODS_FOR_CORRELATION} needed for a rank correlation to mean "
                f"anything. Reported for completeness only.")
    if not np.isfinite(rho):
        return "NOT COMPUTABLE: no variation in one of the two rankings."

    excludes_zero = np.isfinite(ci_low) and np.isfinite(ci_high) and ci_low > 0
    if rho >= PREREGISTERED_RHO_THRESHOLD and excludes_zero:
        return ("HEADLINE: ACS ranking tracks true accuracy at the pre-registered bar. "
                "Anatomic concordance is usable to choose a deconvolution method in "
                "tissue with an anatomical atlas, using no ground truth and no outcomes.")
    if rho >= PREREGISTERED_RHO_THRESHOLD and not excludes_zero:
        return (f"BELOW BAR: rho = {rho:.3f} meets the threshold but its CI "
                f"[{ci_low:.3f}, {ci_high:.3f}] includes zero. The pre-registered bar "
                f"required both. Reported as not met.")
    return (f"NULL RESULT: rho = {rho:.3f} is below the pre-registered "
            f"{PREREGISTERED_RHO_THRESHOLD:.2f}. Reproducing known biology is NOT "
            f"evidence that a composition estimate is numerically correct — a direct, "
            f"quantified warning about the informal check the field relies on. This is "
            f"a reportable finding, not a failed run.")


def test_agreement(acs_scores: dict[str, float], truth_scores: dict[str, float],
                   yardstick: str, higher_truth_is_better: bool = False,
                   exclude_controls: bool = True, n_boot: int = 5000,
                   seed: int = config.RANDOM_SEED) -> AgreementResult:
    """
    Correlate the ACS ranking with a ground-truth ranking.

    `truth_scores` is whatever the yardstick produces per method. Set
    `higher_truth_is_better=False` for error metrics such as MAE, where the sign must be
    flipped so that both rankings point the same way — getting that backwards would turn
    a perfect result into a perfectly negative one.
    """
    shared = sorted(set(acs_scores) & set(truth_scores))
    if exclude_controls:
        shared = [m for m in shared if m not in CONTROL_NAMES]

    usable = [m for m in shared
              if np.isfinite(acs_scores[m]) and np.isfinite(truth_scores[m])]

    # Group methods whose (acs, truth) pair is identical to another's.
    seen: dict[tuple, list[str]] = {}
    for m in usable:
        seen.setdefault((round(acs_scores[m], 12), round(truth_scores[m], 12)),
                        []).append(m)
    duplicate_groups = [sorted(g) for g in seen.values() if len(g) > 1]
    n_distinct = len(seen)

    if len(usable) < 3:
        return AgreementResult(
            yardstick=yardstick, rho=float("nan"), p_value=float("nan"),
            ci_low=float("nan"), ci_high=float("nan"), n_methods=len(usable),
            n_distinct=n_distinct, duplicate_groups=duplicate_groups,
            meets_threshold=False, ci_excludes_zero=False,
            verdict=f"NOT COMPUTABLE: only {len(usable)} methods have both scores.",
            method_ranks={})

    a = np.array([acs_scores[m] for m in usable], dtype="float64")
    t = np.array([truth_scores[m] for m in usable], dtype="float64")
    if not higher_truth_is_better:
        t = -t

    with warnings.catch_warnings():
        # A constant ranking has no defined correlation. That is handled below (rho is
        # NaN and `interpret` returns NOT COMPUTABLE); the warning would only add noise
        # to a run whose output is already correct.
        warnings.simplefilter("ignore", stats.ConstantInputWarning)
        rho, p = stats.spearmanr(a, t)
    lo, hi = _bootstrap_rho(a, t, n_boot=n_boot, seed=seed)
    excludes_zero = bool(np.isfinite(lo) and lo > 0)

    verdict = interpret(float(rho), lo, hi, n_distinct)
    if duplicate_groups:
        verdict += (
            f" NOTE: {len(usable)} methods but only {n_distinct} distinct score pairs — "
            f"{'; '.join('='.join(g) for g in duplicate_groups)} are degenerate "
            f"duplicates. Judge the correlation on {n_distinct} points, not {len(usable)}."
        )

    return AgreementResult(
        yardstick=yardstick, rho=float(rho), p_value=float(p),
        ci_low=lo, ci_high=hi, n_methods=len(usable),
        n_distinct=n_distinct, duplicate_groups=duplicate_groups,
        meets_threshold=bool(np.isfinite(rho) and rho >= PREREGISTERED_RHO_THRESHOLD),
        ci_excludes_zero=excludes_zero,
        verdict=verdict,
        method_ranks={
            m: {"acs": float(acs_scores[m]),
                "acs_rank": int(r_a),
                "truth": float(truth_scores[m]),
                "truth_rank": int(r_t)}
            for m, r_a, r_t in zip(
                usable,
                stats.rankdata(-a),          # rank 1 = best ACS
                stats.rankdata(-t),          # rank 1 = best truth, after sign alignment
            )
        },
    )


def run_all_yardsticks(acs_scores: dict[str, float],
                       yardstick_scores: dict[str, dict[str, float]],
                       higher_is_better: dict[str, bool] | None = None,
                       **kwargs) -> pd.DataFrame:
    """
    Every available yardstick, one row each.

    A yardstick with no data is listed as unavailable rather than omitted: "we could not
    obtain ABSOLUTE purity" and "ABSOLUTE purity disagreed" are very different findings
    and must not look the same in the output table.
    """
    higher_is_better = higher_is_better or {}
    rows = []
    for name in ALL_YARDSTICKS:
        scores = yardstick_scores.get(name)
        if not scores:
            rows.append({"yardstick": name, "rho": np.nan, "n_methods": 0,
                         "verdict": "UNAVAILABLE: this yardstick produced no scores "
                                    "in this run."})
            continue
        res = test_agreement(acs_scores, scores, yardstick=name,
                             higher_truth_is_better=higher_is_better.get(name, False),
                             **kwargs)
        row = {k: v for k, v in res.to_dict().items() if k != "method_ranks"}
        if name == SIMULATED_YARDSTICK and np.isfinite(res.rho):
            # Prepend the caveat to the verdict itself, not to a footnote. A verdict
            # that reads "HEADLINE" without it would be quoted without it.
            row["verdict"] = ("[SIMULATED YARDSTICK - NOT the protocol's yardstick 1] "
                              + row["verdict"])
        rows.append(row)
    return pd.DataFrame(rows).set_index("yardstick")


def write_report(table: pd.DataFrame, path) -> dict:
    """Persist the primary result with the bar it was judged against."""
    real = [y for y in YARDSTICKS
            if y in table.index and np.isfinite(table.loc[y, "rho"])]
    report = {
        "primary_result": "Spearman rho between the ACS ranking and the ground-truth "
                          "ranking of the same methods",
        "primary_result_status": (
            "AVAILABLE" if real else
            "NOT YET COMPUTABLE: none of the protocol's three ground-truth yardsticks "
            "covers enough methods for a rank correlation. Any rho below is from the "
            "SIMULATED yardstick and is a machinery check plus a preliminary reading, "
            "not the study's primary result."),
        "yardsticks_with_real_ground_truth": real,
        "preregistered_threshold": PREREGISTERED_RHO_THRESHOLD,
        "preregistered_requirement": "rho >= threshold AND bootstrap CI excluding zero",
        "min_methods_for_correlation": MIN_METHODS_FOR_CORRELATION,
        "controls_excluded_from_correlation": True,
        "by_yardstick": json.loads(table.reset_index().to_json(orient="records")),
    }
    from pathlib import Path
    Path(path).write_text(json.dumps(report, indent=2))
    return report
