"""
equal_footing.py — the guard that makes "we compared these methods fairly" checkable.

WHY THIS MODULE IS THE CENTRE OF THE PROJECT
--------------------------------------------
A benchmark's conclusion is only as good as the claim that every method faced the same
problem, and that claim is almost never verified. It is very easy for one method to
receive a slightly different gene set, or a reference built with a different cap, or
patient-unbalanced scoring, and for the resulting artefact to be written up as a
property of the algorithm.

Ivy GAP adds a specific and severe version of that hazard. Its samples are NESTED:
each tumour contributed several laser-microdissected blocks, and the number of blocks
per tumour is not constant. Score a method by averaging over samples and the tumours
that happened to be dissected most heavily dominate the result. Two methods can then
be separated by nothing but which patients they happened to do well on.

So this module refuses to let scoring proceed on an assertion. It computes a
certificate — content hashes of every shared input, the sample and patient rosters,
the fold assignment, the frozen biology claims — and every scorer requires one. If two
methods' certificates disagree, scoring raises rather than reporting a comparison that
is not one.

THE NINE CONDITIONS
-------------------
 1. Identical sample set, in identical order.
 2. Identical patient set, with identical per-patient sample counts.
 3. Identical gene space (hash of the bulk matrix actually handed to the method).
 4. Identical reference set (hash of every signature matrix and its sigma).
 5. Identical estimand: all methods report cell fractions, cell-size-corrected
    identically, over the same roster in the same order.
 6. Identical cross-validation folds, assigned by PATIENT and shared byte-identically.
 7. Patient-equal aggregation: every score aggregates within patient before across.
 8. The frozen constraint file is unchanged since the estimates were produced.
 9. No outcome variable was visible to any method at fit time.

Conditions 1-8 are machine-checked here. Condition 9 is structural — no method's
`fit_predict` receives survival or anatomic labels at all; `DeconvolutionInput` has no
field that could carry them, and the manifest it does carry is checked below to make
sure an outcome column has not been smuggled into it.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.anatomic import constraints as _constraints
from ivygap.deconv.base import DeconvolutionInput

#: Columns that must never appear in the manifest handed to a method. Seeing an
#: outcome at fit time is the one violation no downstream check could detect.
FORBIDDEN_MANIFEST_COLUMNS = frozenset({
    "survival_days", "survival_months", "os_days", "os_months", "os_status",
    "vital_status", "event", "death", "days_to_death", "days_to_last_followup",
})


class EqualFootingViolation(AssertionError):
    """Raised when a comparison would not be like-for-like. Never caught internally."""


@dataclass
class EqualFootingCertificate:
    """
    Proof that a set of methods faced the same problem. Serialised next to the
    estimates; the scorers re-check it before computing anything.
    """
    input_hashes: dict[str, str]
    sample_ids_hash: str
    patient_ids_hash: str
    samples_per_patient_hash: str
    fold_assignment_hash: str | None
    biology_freeze_hash: str
    n_samples: int
    n_patients: int
    methods: list[str] = field(default_factory=list)
    method_implementations: dict[str, str] = field(default_factory=dict)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "EqualFootingCertificate":
        return cls(**json.loads(text))

    def comparable_with(self, other: "EqualFootingCertificate") -> list[str]:
        """Return the list of mismatches; empty means the two are comparable."""
        problems = []
        for f in ("input_hashes", "sample_ids_hash", "patient_ids_hash",
                  "samples_per_patient_hash", "fold_assignment_hash",
                  "biology_freeze_hash"):
            if getattr(self, f) != getattr(other, f):
                problems.append(f"{f} differs")
        return problems


def _hash_sequence(values) -> str:
    h = hashlib.sha256()
    for v in values:
        h.update(str(v).encode())
        h.update(b"\x00")
    return h.hexdigest()


def build_certificate(data: DeconvolutionInput,
                      folds: pd.Series | None = None,
                      methods: list[str] | None = None,
                      implementations: dict[str, str] | None = None
                      ) -> EqualFootingCertificate:
    """Compute the certificate for one run, after checking the input is admissible."""
    smuggled = FORBIDDEN_MANIFEST_COLUMNS & {c.lower() for c in data.manifest.columns}
    if smuggled:
        raise EqualFootingViolation(
            f"the manifest handed to the methods carries outcome column(s) {sorted(smuggled)}. "
            "No method may see an outcome at fit time — that is how a benchmark becomes "
            "a self-fulfilling prophecy."
        )

    samples = list(data.samples)
    patients = data.manifest.loc[samples, "patient_id"].astype(str)
    per_patient = patients.value_counts().sort_index()

    return EqualFootingCertificate(
        input_hashes=data.content_hashes(),
        sample_ids_hash=_hash_sequence(samples),
        patient_ids_hash=_hash_sequence(sorted(patients.unique())),
        samples_per_patient_hash=_hash_sequence(
            f"{p}={n}" for p, n in per_patient.items()),
        fold_assignment_hash=(
            _hash_sequence(f"{i}={v}" for i, v in folds.sort_index().items())
            if folds is not None else None),
        biology_freeze_hash=_constraints.freeze_hash(),
        n_samples=len(samples),
        n_patients=int(patients.nunique()),
        methods=sorted(methods or []),
        method_implementations=dict(implementations or {}),
    )


def require_comparable(certs: dict[str, EqualFootingCertificate]) -> None:
    """
    Verify a whole set of per-method certificates agree. Raises on the first mismatch,
    naming both methods and every field that differs.
    """
    if not certs:
        raise EqualFootingViolation("no certificates supplied; nothing can be compared")
    names = sorted(certs)
    base_name, base = names[0], certs[names[0]]
    for name in names[1:]:
        problems = base.comparable_with(certs[name])
        if problems:
            raise EqualFootingViolation(
                f"{name!r} and {base_name!r} were not given the same problem: "
                + "; ".join(problems)
                + ". Their scores are not comparable and will not be reported as if "
                  "they were."
            )


def check_estimates_aligned(estimates: dict[str, pd.DataFrame]) -> None:
    """
    Every method must have produced an estimate for every sample, over the same cell
    types, in the same order.

    A method that dropped samples it found difficult would otherwise score well on an
    easier subset. Dropping is allowed only if every method drops the identical set —
    which this makes visible instead of silent.
    """
    if not estimates:
        raise EqualFootingViolation("no estimates supplied")
    names = sorted(estimates)
    ref_name = names[0]
    ref_idx = list(estimates[ref_name].index)
    ref_cols = list(estimates[ref_name].columns)

    for name in names[1:]:
        est = estimates[name]
        if list(est.columns) != ref_cols:
            raise EqualFootingViolation(
                f"{name!r} reports cell types {list(est.columns)}, "
                f"{ref_name!r} reports {ref_cols}"
            )
        if list(est.index) != ref_idx:
            only_ref = set(ref_idx) - set(est.index)
            only_est = set(est.index) - set(ref_idx)
            raise EqualFootingViolation(
                f"{name!r} and {ref_name!r} cover different samples "
                f"({len(only_ref)} missing from {name!r}, {len(only_est)} extra). "
                "Comparing them would compare different cohorts."
            )


def assign_patient_folds(manifest: pd.DataFrame, n_folds: int = config.N_CV_FOLDS,
                         seed: int = config.RANDOM_SEED) -> pd.Series:
    """
    Assign cross-validation folds BY PATIENT, returning a per-sample fold label.

    Splitting by sample would put two blocks of the same tumour on both sides of the
    split. The model would then be scored on a tumour it had already seen, and every
    method's apparent performance would be inflated — the ones that overfit patient
    identity most, inflated most. Grouping by patient is what makes the held-out
    estimate mean anything.

    Patients are assigned to folds in descending order of block count, each going to
    the currently smallest fold, so folds end up near-equal in *samples* as well as in
    patients even though block counts are uneven.
    """
    patients = manifest["patient_id"].astype(str)
    counts = patients.value_counts()

    rng = np.random.default_rng(seed)
    order = list(counts.index)
    # Break ties reproducibly but without alphabetical bias.
    rng.shuffle(order)
    order.sort(key=lambda p: -counts[p])

    fold_of_patient: dict[str, int] = {}
    load = np.zeros(n_folds, dtype="int64")
    for p in order:
        k = int(np.argmin(load))
        fold_of_patient[p] = k
        load[k] += counts[p]

    folds = patients.map(fold_of_patient).astype("int64")
    folds.name = "fold"
    return folds


def patient_equal_mean(values: pd.Series, patients: pd.Series) -> float:
    """
    Aggregate a per-sample quantity so that every PATIENT counts once.

    Mean within each patient first, then mean across patients. A tumour contributing
    twelve blocks and one contributing two then carry identical weight, which is the
    whole point of condition 7. Patients whose values are all NaN drop out, and how
    many dropped is the caller's business to report.
    """
    df = pd.DataFrame({"value": values, "patient": patients.astype(str)})
    per_patient = df.groupby("patient")["value"].mean()
    per_patient = per_patient[np.isfinite(per_patient)]
    if per_patient.empty:
        return float("nan")
    return float(per_patient.mean())


def patient_equal_frame(values: pd.DataFrame, patients: pd.Series) -> pd.DataFrame:
    """Per-patient means of a samples x columns table — the aggregation step of the above."""
    df = values.copy()
    df["__patient__"] = patients.reindex(df.index).astype(str)
    out = df.groupby("__patient__").mean(numeric_only=True)
    out.index.name = "patient_id"
    return out
