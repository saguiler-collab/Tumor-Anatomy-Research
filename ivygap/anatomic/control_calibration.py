"""
control_calibration.py — the negative control as a distribution, not one draw.

THE PROBLEM THIS FIXES
----------------------
`control_shuffled_signature` is specified as a single permutation of the signature's
gene rows, drawn from `config.RANDOM_SEED`. One draw is one sample from a random
variable, and on real Ivy GAP data that variable turns out to have a standard deviation
of about 0.16 — against a total spread of 0.29 across every real method on the
leaderboard.

That is not a rounding concern. Measured on this cohort, the pre-specified draw landed
at the 91st percentile of its own distribution, scoring 0.585 where the mean is 0.361.
Read as "the control", 0.585 ties a real method and beats its own permutation null, and
invites the conclusion that the constraint set is too permissive. Read as one draw from
a distribution centred near chance, it says the opposite. A superseded run of the same
control on a slightly different gene set recorded 0.215 — also a legitimate draw.

So a single-draw control cannot support "the controls behave" in EITHER direction, and
whichever way it happens to land, the reading is an artefact of the seed.

WHAT THIS DOES AND DOES NOT CHANGE
----------------------------------
It does not change the constraint file. It does not change the control's definition, and
it does not change the leaderboard row: the pre-specified single draw is still reported,
exactly as it ran, because that is what was registered.

What it adds is the control's own sampling distribution, so the leaderboard row can be
placed in it, and so the comparison a reader actually wants — "is this method
distinguishable from a meaningless signature?" — can be answered against a percentile
instead of against one arbitrary draw.

Note the direction of travel: replacing a point estimate with a distribution makes the
control HARDER to beat and the test more conservative, never less. That is the only
direction it is safe to move an instrument after seeing what it scored.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from ivygap import config
from ivygap.anatomic import acs, constraints as K

#: Enough draws that the 5th and 95th percentiles are stable to about 0.01, which is
#: below the 1/65 = 0.0154 granularity of ACS on this cohort. More draws cannot make the
#: comparison finer than the metric itself.
DEFAULT_DRAWS = 200


def _draw_shuffled(S: np.ndarray, B: np.ndarray, n_types: int, seed: int) -> np.ndarray:
    """Identical arithmetic to `ShuffledSignatureControl`, with the seed as an argument."""
    rng = np.random.default_rng(seed)
    return np.vstack([nnls(S[rng.permutation(S.shape[0]), :], B[:, j])[0]
                      for j in range(B.shape[1])])


def _draw_random(S: np.ndarray, B: np.ndarray, n_types: int, seed: int) -> np.ndarray:
    """Identical arithmetic to `RandomFractionControl`, with the seed as an argument."""
    from ivygap.deconv.controls import RandomFractionsControl
    rng = np.random.default_rng(seed)
    conc = RandomFractionsControl().concentration
    return rng.dirichlet(np.full(n_types, conc), size=B.shape[1])


#: control name -> the draw it makes. Both controls are single draws as specified, and
#: both therefore need a distribution rather than a point. `control_random` never came
#: close to a real method on this cohort, but "it did not matter this time" is not a
#: reason to leave one of two controls uncalibrated while arguing that single-draw
#: controls are the problem.
DRAWERS = {
    "control_shuffled_signature": _draw_shuffled,
    "control_random": _draw_random,
}


def _to_fractions(raw: np.ndarray, reference, index) -> pd.DataFrame:
    """
    Raw weights -> cell fractions through the pipeline's OWN conversion.

    Imported rather than reimplemented: a calibration that normalised even slightly
    differently from `fit_predict` would produce numbers that cannot be compared against
    the leaderboard row they exist to calibrate.
    """
    from ivygap.deconv.base import project_to_simplex, to_cell_fractions

    types = list(reference.profile.columns)
    cell_size = reference.cell_size.reindex(types).to_numpy()
    rows = []
    for i in range(raw.shape[0]):
        w = project_to_simplex(raw[i])
        if np.isfinite(w).all():
            w = to_cell_fractions(w, cell_size)
        rows.append(w)
    return pd.DataFrame(np.vstack(rows), index=index, columns=types)


def calibrate_one(bulk: pd.DataFrame, manifest: pd.DataFrame, reference,
                  control: str = "control_shuffled_signature",
                  n_draws: int = DEFAULT_DRAWS,
                  seed: int = config.RANDOM_SEED) -> tuple[np.ndarray, pd.Series]:
    """
    Score one control over `n_draws` independent draws.

    Returns (scores, per-constraint satisfied rate). `bulk` and `manifest` must already
    be restricted to the ACS scoring set — both controls treat samples independently, so
    samples that are never scored cannot affect the result and only cost time.
    """
    draw = DRAWERS[control]
    S = reference.profile.to_numpy(dtype="float64")
    B = bulk.to_numpy(dtype="float64")
    n_types = len(reference.profile.columns)

    scores = np.empty(n_draws, dtype="float64")
    per_constraint: list[dict] = []
    for i in range(n_draws):
        raw = draw(S, B, n_types, seed + i)
        est = _to_fractions(raw, reference, bulk.columns)
        table = acs._score_table(acs.tumor_structure_means(est, manifest))
        scores[i] = acs._acs_from_table(table)
        for c in K.CONSTRAINTS:
            sub = table[table["constraint"] == c.id]
            per_constraint.append({
                "constraint": c.id,
                "n_evaluable": len(sub),
                "n_satisfied": int(sub["satisfied"].sum()) if len(sub) else 0,
            })
    pc = pd.DataFrame(per_constraint)
    rates = (pc.groupby("constraint")
               .apply(lambda d: d["n_satisfied"].sum() / max(d["n_evaluable"].sum(), 1),
                      include_groups=False))
    return scores, rates


def calibrate(bulk: pd.DataFrame, manifest: pd.DataFrame, reference,
              observed: dict[str, float] | None = None,
              n_draws: int = DEFAULT_DRAWS,
              seed: int = config.RANDOM_SEED) -> dict:
    """
    Calibrate every negative control, and place each real method against the
    HARDEST of them.

    `observed` maps method name -> ACS. Controls in it are skipped when placing methods,
    but their as-run draw is located within their own distribution.
    """
    per_control = {}
    for name in DRAWERS:
        sc, rates = calibrate_one(bulk, manifest, reference, control=name,
                                  n_draws=n_draws, seed=seed)
        f = sc[np.isfinite(sc)]
        per_control[name] = {
            "n_draws": int(n_draws),
            "acs": {
                "mean": float(f.mean()) if f.size else float("nan"),
                "sd": float(f.std(ddof=1)) if f.size > 1 else float("nan"),
                "min": float(f.min()) if f.size else float("nan"),
                "p05": float(np.quantile(f, 0.05)) if f.size else float("nan"),
                "median": float(np.median(f)) if f.size else float("nan"),
                "p95": float(np.quantile(f, 0.95)) if f.size else float("nan"),
                "max": float(f.max()) if f.size else float("nan"),
            },
            "per_constraint_satisfied_rate": {k: float(v) for k, v in rates.items()},
            "_finite": f,
        }
        if observed and name in observed and np.isfinite(observed[name]) and f.size:
            per_control[name]["as_run_draw"] = {
                "acs": float(observed[name]),
                "percentile_among_draws": float(np.mean(f <= observed[name])),
            }

    # A method must clear the hardest control, not a convenient one.
    #
    # ACS is quantised to multiples of 1/D (1/65 on this cohort), so two controls whose
    # distributions are genuinely different can still share a 95th percentile exactly —
    # it happens here. `max` would then resolve the tie by dict insertion order, which
    # is not a decision anyone made. Break ties explicitly: higher p95, then higher
    # mean, then name, so the choice is reproducible and inspectable.
    hardest = max(per_control,
                  key=lambda n: (per_control[n]["acs"]["p95"],
                                 per_control[n]["acs"]["mean"], n))
    finite = per_control[hardest]["_finite"]
    p95 = float(np.quantile(finite, 0.95)) if finite.size else float("nan")
    report = {
        "what_this_is": (
            "sampling distributions of the negative controls. Each leaderboard control "
            "row is ONE draw from these; nothing here changes a row, a control's "
            "definition, or the constraint file."),
        "n_draws": int(n_draws),
        "n_samples": int(bulk.shape[1]),
        "n_tumors": int(manifest["patient_id"].nunique()),
        "controls": per_control,
        "hardest_control": hardest,
        # kept at the top level so existing readers of `acs` keep working; it is the
        # hardest control's distribution, which is the one a method must clear.
        "acs": per_control[hardest]["acs"],
        "per_constraint_satisfied_rate":
            per_control[hardest]["per_constraint_satisfied_rate"],
    }
    if "as_run_draw" in per_control[hardest]:
        report["as_run_draw"] = per_control[hardest]["as_run_draw"]

    from ivygap.deconv.controls import CONTROL_NAMES
    if observed:
        real = {m: v for m, v in observed.items()
                if m not in CONTROL_NAMES and np.isfinite(v)}
        placed = {}
        for m, v in real.items():
            placed[m] = {
                "acs": float(v),
                "percentile_among_control_draws": (
                    float(np.mean(finite <= v)) if finite.size else float("nan")),
                "distinguishable_from_a_meaningless_signature": bool(v > p95),
            }
        report["methods_vs_control_distribution"] = placed
        report["not_distinguishable_from_noise"] = sorted(
            m for m, d in placed.items()
            if not d["distinguishable_from_a_meaningless_signature"])

    for entry in per_control.values():
        entry.pop("_finite", None)          # numpy arrays are not JSON-serialisable

    nd = report.get("not_distinguishable_from_noise", [])
    report["verdict"] = (
        f"Against the hardest control ({hardest}, 95th percentile {p95:.3f}): "
        + (f"{len(nd)} real method(s) are NOT distinguishable from it — "
           f"{', '.join(nd)}."
           if nd else
           "every real method scores above it.")
    )
    return report
