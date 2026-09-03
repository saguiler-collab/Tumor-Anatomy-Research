"""
acs.py — the Anatomic Concordance Score.

WHAT ACS IS
-----------
For every constraint and every tumour that contributed the structures the constraint
needs, take the paired within-tumour comparison and record whether it points the way the
constraint predicted. ACS is the weighted proportion of satisfied constraint-tumour
pairs.

    ACS = sum over satisfied pairs of w_c  /  sum over evaluable pairs of w_c

Everything is paired within a tumour. A patient is never compared against a different
patient, so between-patient differences in treatment, subtype, immune background and RNA
quality cancel — and Ivy GAP's uneven block counts cannot let a heavily dissected tumour
dominate, because each tumour contributes exactly one verdict per constraint.

WHY THE PERMUTATION NULL IS THE ONE THAT MATTERS
------------------------------------------------
The tempting null is "50%, by coin flip". It is wrong, and it flatters every method.
Cell fractions within a tumour are not exchangeable coin flips: they are compositional
(they sum to one), correlated across types, and unevenly spread across structures. A
method could clear 50% on several constraints through those artefacts alone.

Permuting the STRUCTURE LABELS within each tumour preserves every one of those
properties — the same compositions, the same tumour, the same number of blocks — and
destroys only the correspondence between composition and anatomy. That is the thing the
constraint set claims to detect, so that is the only thing the null should remove.

The permutation is within-tumour for the same reason the scoring is: a global shuffle
would also destroy between-tumour structure, and the resulting null would be too easy to
beat.

WHAT IS DELIBERATELY *NOT* HERE
-------------------------------
Nothing in this module selects a method, and nothing in it can be tuned by looking at a
score. The constraint set is frozen in `constraints.py` and hashed; if a constraint is
edited the hash changes and previously computed scores are invalidated rather than
silently compared against new ones. A constraint that no method satisfies is reported as
such — it is evidence about the constraint, and possibly about the reference, and it is
not permitted to quietly disappear.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ivygap import config
from . import constraints as K

#: A tumour must contribute at least this many samples of a structure for that
#: structure's value to be used. One block is still usable — Ivy GAP dissected each
#: structure once per tumour in many cases — but zero is not, and imputing it would
#: invent the signal being tested.
MIN_SAMPLES_PER_STRUCTURE = 1


@dataclass
class ACSResult:
    """One method's score, with everything needed to argue about it."""
    method: str
    acs: float
    ci_low: float
    ci_high: float
    n_pairs: int
    n_tumors: int
    weight_evaluated: float
    weight_satisfied: float
    null_mean: float
    null_p: float
    per_constraint: pd.DataFrame
    #: One row per tumour: which constraints it could evaluate and how many it satisfied.
    #: The protocol asks for per-tumour results explicitly ("Report per-tumor results,
    #: not just the pool") because with 9 evaluable tumours a pooled ACS can be carried
    #: by one or two of them, and a reader cannot tell from the pooled number alone.
    per_tumor: pd.DataFrame = field(default_factory=pd.DataFrame)
    coverage: dict = field(default_factory=dict)

    def summary(self) -> dict:
        return {
            "method": self.method,
            "acs": self.acs,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "n_constraint_tumor_pairs": self.n_pairs,
            "n_tumors": self.n_tumors,
            "null_mean": self.null_mean,
            "null_p": self.null_p,
            "beats_null": bool(self.null_p < 0.05 and self.acs > self.null_mean),
        }


def tumor_structure_means(estimate: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse to one composition per (tumour, structure).

    A tumour that contributed three cellular-tumour blocks gets one CT row, not three.
    This is where the uneven block counts are neutralised, once, for everything
    downstream.
    """
    df = estimate.copy()
    meta = manifest.reindex(df.index)
    df["tumor_id"] = meta["patient_id"].astype(str)
    df["structure"] = meta["structure"].astype(str)
    return df.groupby(["tumor_id", "structure"]).mean(numeric_only=True)


def _value(ts: pd.DataFrame, tumor: str, structure: str, cell_type: str) -> float:
    try:
        return float(ts.loc[(tumor, structure), cell_type])
    except KeyError:
        return np.nan


def evaluate_constraint(ts: pd.DataFrame, tumor: str, c: K.Constraint) -> float | None:
    """
    Score one constraint for one tumour.

    Returns 1.0 (satisfied), 0.0 (violated), or None when this tumour did not contribute
    the structures the constraint needs. None is not a failure: an unevaluable pair is
    excluded from both numerator and denominator, so a tumour missing MVP does not
    penalise a method for something it was never shown.

    Exact ties score 0. A method that reports identical values in two structures has not
    reproduced the ordering, and giving half credit would let a method that collapses
    every structure to the same composition score 0.5 everywhere.
    """
    if c.kind == "pairwise":
        hi, lo = (_value(ts, tumor, s, c.cell_type) for s in c.structures)
        if not (np.isfinite(hi) and np.isfinite(lo)):
            return None
        return float(hi > lo)

    if c.kind == "maximum":
        target = c.structures[0]
        vals = {s: _value(ts, tumor, s, c.cell_type) for s in c.among}
        if not np.isfinite(vals[target]):
            return None
        others = [v for s, v in vals.items() if s != target and np.isfinite(v)]
        # Requiring "the maximum" of a single comparison would restate C3. Demand a real
        # comparison set before letting this constraint contribute.
        if len(others) < 2:
            return None
        return float(vals[target] > max(others))

    if c.kind == "monotone":
        vals = [_value(ts, tumor, s, c.cell_type) for s in c.structures]
        if not all(np.isfinite(v) for v in vals):
            return None
        return float(all(a < b for a, b in zip(vals, vals[1:])))

    raise ValueError(f"unknown constraint kind {c.kind!r}")


def _score_table(ts: pd.DataFrame) -> pd.DataFrame:
    """Long table: one row per (constraint, tumour) that could be evaluated."""
    tumors = sorted({t for t, _ in ts.index})
    rows = []
    for c in K.CONSTRAINTS:
        for tumor in tumors:
            got = evaluate_constraint(ts, tumor, c)
            if got is None:
                continue
            rows.append({"constraint": c.id, "tumor_id": tumor,
                         "satisfied": got, "weight": c.weight})
    return pd.DataFrame(rows, columns=["constraint", "tumor_id", "satisfied", "weight"])


def _acs_from_table(table: pd.DataFrame) -> float:
    if table.empty:
        return float("nan")
    denom = float(table["weight"].sum())
    if denom <= 0:
        return float("nan")
    return float((table["satisfied"] * table["weight"]).sum() / denom)


def permutation_null_reference(estimate: pd.DataFrame, manifest: pd.DataFrame,
                               n_permutations: int = 1000,
                               seed: int = config.RANDOM_SEED) -> np.ndarray:
    """
    Readable reference implementation of the null, built out of the same pandas pieces
    the real scoring uses.

    This is the definition. `permutation_null` below is a vectorised rewrite that must
    agree with it exactly — `test_fast_permutation_null_matches_the_reference` asserts
    that. Keeping the slow version is what makes the fast one trustworthy: an optimised
    statistic with no independent implementation to check against is a place where a
    subtle indexing error lives forever.
    """
    rng = np.random.default_rng(seed)
    meta = manifest.reindex(estimate.index)
    tumors = meta["patient_id"].astype(str).to_numpy()
    structures = meta["structure"].astype(str).to_numpy()

    by_tumor = {t: np.where(tumors == t)[0] for t in np.unique(tumors)}
    out = np.empty(n_permutations, dtype="float64")

    for i in range(n_permutations):
        shuffled = structures.copy()
        for idx in by_tumor.values():
            shuffled[idx] = rng.permutation(structures[idx])
        fake = pd.DataFrame({"patient_id": tumors, "structure": shuffled},
                            index=estimate.index)
        out[i] = _acs_from_table(_score_table(tumor_structure_means(estimate, fake)))
    return out


def _compile_constraints(cell_types: list[str], structures: list[str]):
    """
    Translate the constraint file into integer index arrays once, so the permutation
    loop does no name lookups.
    """
    ci = {c: i for i, c in enumerate(cell_types)}
    si = {s: i for i, s in enumerate(structures)}
    compiled = []
    for c in K.CONSTRAINTS:
        if c.cell_type not in ci:
            continue
        if c.kind == "pairwise":
            if not all(s in si for s in c.structures):
                continue
            compiled.append(("pairwise", ci[c.cell_type],
                             (si[c.structures[0]], si[c.structures[1]]), c.weight))
        elif c.kind == "maximum":
            if c.structures[0] not in si:
                continue
            among = tuple(si[s] for s in c.among if s in si)
            compiled.append(("maximum", ci[c.cell_type],
                             (si[c.structures[0]], among), c.weight))
        elif c.kind == "monotone":
            if not all(s in si for s in c.structures):
                continue
            compiled.append(("monotone", ci[c.cell_type],
                             tuple(si[s] for s in c.structures), c.weight))
    return compiled


def _acs_from_arrays(means: np.ndarray, has_obs: np.ndarray, compiled) -> float:
    """
    ACS from a (tumour, structure, cell_type) mean array plus a matching observation
    mask. Identical rules to `evaluate_constraint`: strict inequality, and a pair that
    cannot be evaluated is excluded from BOTH numerator and denominator rather than
    scored as a violation.

    `has_obs` is indexed per cell type, not per structure: a method can fail on one cell
    type of a sample while the rest are fine, and treating the whole structure as
    missing would discard usable evidence.
    """
    num = den = 0.0
    n_tumors = means.shape[0]
    for kind, ct, spec, w in compiled:
        if kind == "pairwise":
            hi, lo = spec
            ok = has_obs[:, hi, ct] & has_obs[:, lo, ct]
            if not ok.any():
                continue
            num += w * float((means[ok, hi, ct] > means[ok, lo, ct]).sum())
            den += w * float(ok.sum())
        elif kind == "maximum":
            target, among = spec
            for t in range(n_tumors):
                if not has_obs[t, target, ct]:
                    continue
                others = [s for s in among if s != target and has_obs[t, s, ct]]
                if len(others) < 2:
                    continue
                den += w
                if means[t, target, ct] > max(means[t, s, ct] for s in others):
                    num += w
        else:                                          # monotone
            ok = np.ones(n_tumors, dtype=bool)
            for s in spec:
                ok &= has_obs[:, s, ct]
            if not ok.any():
                continue
            chain = np.ones(int(ok.sum()), dtype=bool)
            for a, b in zip(spec, spec[1:]):
                chain &= means[ok, a, ct] < means[ok, b, ct]
            num += w * float(chain.sum())
            den += w * float(ok.sum())
    return float(num / den) if den > 0 else float("nan")


def permutation_null(estimate: pd.DataFrame, manifest: pd.DataFrame,
                     n_permutations: int = 10_000,
                     seed: int = config.RANDOM_SEED) -> np.ndarray:
    """
    Null distribution of ACS under within-tumour permutation of structure labels.

    Each draw reassigns the structure labels among that tumour's own samples. The
    tumour's compositions, its block count and the number of blocks per structure are all
    preserved; only the correspondence between composition and anatomy is broken.

    Vectorised so the protocol's 10,000 permutations are affordable for every method:
    group means are accumulated with `np.add.at` on integer codes rather than rebuilt by
    a pandas groupby each draw, which is roughly two orders of magnitude faster. The
    readable definition is `permutation_null_reference`, and a test asserts the two
    produce identical values.
    """
    rng = np.random.default_rng(seed)
    meta = manifest.reindex(estimate.index)

    tumor_labels = meta["patient_id"].astype(str).to_numpy()
    struct_labels = meta["structure"].astype(str).to_numpy()

    tumors = sorted(set(tumor_labels))
    structures = sorted(set(struct_labels))
    t_code = np.array([tumors.index(t) for t in tumor_labels])
    s_code = np.array([structures.index(s) for s in struct_labels])

    cell_types = list(estimate.columns)
    values = estimate.to_numpy(dtype="float64")
    n_t, n_s, n_c = len(tumors), len(structures), len(cell_types)

    compiled = _compile_constraints(cell_types, structures)
    by_tumor = [np.where(t_code == i)[0] for i in range(n_t)]

    # A method that fails on a sample emits NaN (see `project_to_simplex`: an all-zero
    # solution has no composition, and returning a uniform vector would be a
    # fabrication). pandas' groupby.mean skips those; a plain np.add.at would propagate
    # them and would additionally count a NaN sample as "structure present". Both
    # divergences are silent and both move the null. So accumulate only finite values,
    # and define presence as "at least one FINITE observation".
    finite = np.isfinite(values)
    safe = np.where(finite, values, 0.0)

    out = np.empty(n_permutations, dtype="float64")
    for i in range(n_permutations):
        shuffled = s_code.copy()
        for idx in by_tumor:
            shuffled[idx] = rng.permutation(s_code[idx])

        sums = np.zeros((n_t, n_s, n_c))
        counts = np.zeros((n_t, n_s, n_c))
        np.add.at(sums, (t_code, shuffled), safe)
        np.add.at(counts, (t_code, shuffled), finite.astype("float64"))

        # Per cell type, because one type can be NaN while others are not.
        has_obs = counts > 0
        means = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=has_obs)
        out[i] = _acs_from_arrays(means, has_obs, compiled)
    return out


def _group_means_array(estimate: pd.DataFrame, manifest: pd.DataFrame):
    """
    The observed (tumour, structure, cell type) means, in the same array form the
    permutation loop uses.

    Exists so the observed statistic and its null are produced by ONE estimator. Two
    implementations that agree on clean data and diverge on NaN would make every
    p-value quietly wrong the first time a method failed on a sample.
    """
    meta = manifest.reindex(estimate.index)
    tumor_labels = meta["patient_id"].astype(str).to_numpy()
    struct_labels = meta["structure"].astype(str).to_numpy()

    tumors = sorted(set(tumor_labels))
    structures = sorted(set(struct_labels))
    t_code = np.array([tumors.index(t) for t in tumor_labels])
    s_code = np.array([structures.index(s) for s in struct_labels])

    cell_types = list(estimate.columns)
    values = estimate.to_numpy(dtype="float64")
    finite = np.isfinite(values)
    safe = np.where(finite, values, 0.0)

    sums = np.zeros((len(tumors), len(structures), len(cell_types)))
    counts = np.zeros_like(sums)
    np.add.at(sums, (t_code, s_code), safe)
    np.add.at(counts, (t_code, s_code), finite.astype("float64"))

    has_obs = counts > 0
    means = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=has_obs)
    return means, has_obs, _compile_constraints(cell_types, structures)


def per_tumor_table(table: pd.DataFrame) -> pd.DataFrame:
    """
    Break the pooled ACS down to one row per tumour.

    ACS is a weighted proportion over (constraint, tumour) pairs, so a tumour that
    happens to supply all five structures contributes more pairs than one that supplied
    two. Reporting only the pool hides both that imbalance and the case where a single
    tumour is responsible for a method's whole margin over another.
    """
    if table.empty:
        return pd.DataFrame(columns=["tumor_id", "n_constraints_evaluable",
                                     "n_satisfied", "weight_evaluated",
                                     "weight_satisfied", "acs", "constraints_satisfied",
                                     "constraints_violated"])
    rows = []
    for tumor, block in table.groupby("tumor_id"):
        w_eval = float(block["weight"].sum())
        w_sat = float((block["satisfied"] * block["weight"]).sum())
        rows.append({
            "tumor_id": tumor,
            "n_constraints_evaluable": int(len(block)),
            "n_satisfied": int(block["satisfied"].sum()),
            "weight_evaluated": w_eval,
            "weight_satisfied": w_sat,
            "acs": w_sat / w_eval if w_eval > 0 else float("nan"),
            "constraints_satisfied": ",".join(
                sorted(block.loc[block["satisfied"] == 1, "constraint"])),
            "constraints_violated": ",".join(
                sorted(block.loc[block["satisfied"] == 0, "constraint"])),
        })
    return pd.DataFrame(rows).sort_values("tumor_id").reset_index(drop=True)


def bootstrap_ci(table: pd.DataFrame, n_boot: int = 2000,
                 seed: int = config.RANDOM_SEED) -> tuple[float, float]:
    """
    Percentile bootstrap over TUMOURS, not over constraint-tumour pairs.

    Resampling pairs would treat one tumour's seven verdicts as seven independent
    observations. They are not — they share a patient, a dissection and a reference
    mismatch — and the interval would come out far too narrow.
    """
    if table.empty:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    tumors = table["tumor_id"].unique()
    if len(tumors) < 2:
        return (float("nan"), float("nan"))

    groups = {t: table[table["tumor_id"] == t] for t in tumors}
    draws = np.empty(n_boot, dtype="float64")
    for i in range(n_boot):
        pick = rng.choice(tumors, size=len(tumors), replace=True)
        draws[i] = _acs_from_table(pd.concat([groups[t] for t in pick], ignore_index=True))
    finite = draws[np.isfinite(draws)]
    if finite.size == 0:
        return (float("nan"), float("nan"))
    return (float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975)))


def score(estimate: pd.DataFrame, manifest: pd.DataFrame, method: str = "method",
          n_permutations: int = 10_000, n_boot: int = 2000,
          seed: int = config.RANDOM_SEED) -> ACSResult:
    """Full ACS for one method: score, bootstrap CI, permutation null, per-constraint."""
    ts = tumor_structure_means(estimate, manifest)
    table = _score_table(ts)

    # The observed statistic MUST come from the same estimator as its null, or the
    # p-value compares two different quantities. The pandas table is kept for the
    # per-constraint breakdown and the bootstrap, but the headline ACS is computed
    # through the identical array path the permutation loop uses, and the two are
    # cross-checked at runtime. A disagreement here is a bug, not a rounding artefact,
    # so it raises rather than being tolerated into a published number.
    means, has_obs, compiled = _group_means_array(estimate, manifest)
    acs = _acs_from_arrays(means, has_obs, compiled)

    acs_pandas = _acs_from_table(table)
    both_nan = not np.isfinite(acs) and not np.isfinite(acs_pandas)
    if not both_nan and not np.isclose(acs, acs_pandas, rtol=0, atol=1e-9):
        raise AssertionError(
            f"the array and pandas ACS estimators disagree "
            f"({acs!r} vs {acs_pandas!r}). They must be identical: the permutation null "
            f"is built with the array path, so any divergence makes the p-value "
            f"meaningless. This is a defect in acs.py, not a data problem."
        )

    lo, hi = bootstrap_ci(table, n_boot=n_boot, seed=seed)

    null = permutation_null(estimate, manifest, n_permutations=n_permutations, seed=seed)
    null = null[np.isfinite(null)]
    # One-sided permutation p with the +1 correction, so a p of exactly zero is never
    # reported from a finite number of permutations.
    null_p = (float((null >= acs).sum() + 1) / (null.size + 1)
              if null.size else float("nan"))

    per_constraint = []
    for c in K.CONSTRAINTS:
        sub = table[table["constraint"] == c.id]
        per_constraint.append({
            "constraint": c.id,
            "description": c.describe(),
            "cell_type": c.cell_type,
            "weight": c.weight,
            "n_tumors_evaluable": int(len(sub)),
            "n_satisfied": int(sub["satisfied"].sum()) if len(sub) else 0,
            "fraction_satisfied": (float(sub["satisfied"].mean()) if len(sub)
                                   else float("nan")),
        })

    return ACSResult(
        method=method,
        per_tumor=per_tumor_table(table),
        acs=acs, ci_low=lo, ci_high=hi,
        n_pairs=int(len(table)),
        n_tumors=int(table["tumor_id"].nunique()) if len(table) else 0,
        weight_evaluated=float(table["weight"].sum()) if len(table) else 0.0,
        weight_satisfied=(float((table["satisfied"] * table["weight"]).sum())
                          if len(table) else 0.0),
        null_mean=float(null.mean()) if null.size else float("nan"),
        null_p=null_p,
        per_constraint=pd.DataFrame(per_constraint),
        coverage=K.coverage_report(),
    )
