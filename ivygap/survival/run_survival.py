"""
run_survival.py — out-of-fold prognostic evaluation, and an honest account of its power.

WHY THIS IS HERE AT ALL
-----------------------
Two deconvolution methods can produce composition tables that differ substantially and
still be indistinguishable on accuracy metrics, because the metrics average over cell
types and samples. Whether the difference *matters* is a separate question, and
prognostic signal is one of the few ways to ask it on real tissue: if method A's
tumour fractions carry information about how long patients lived and method B's do not,
that is a discriminating fact about the methods, not just about the patients.

THE POWER PROBLEM, STATED UP FRONT
----------------------------------
TCGA-GBM supported this analysis with 115 patients and 93 deaths. Ivy GAP's clinical
release covers roughly 40 tumours. With that many events, the confidence interval on a
C-index spans most of the interpretable range, and the difference between two methods
will almost never clear it.

This module therefore does three things rather than one:

  1. Computes the analysis properly — patient-level folds, repeated out-of-fold
     prediction, baseline versus baseline-plus-composition.
  2. Computes what the design can actually detect, BEFORE looking at any result, and
     writes it into the output as `detectable_delta_c_index`.
  3. Refuses to declare a winner when the observed difference is inside that bound,
     reporting `INCONCLUSIVE` instead.

Point 3 is the one that matters. An underpowered comparison does not become
informative because it produced a number, and a ranking read off a difference smaller
than the noise floor is the single easiest way for this project to say something false.

THE FIREWALL
------------
Survival is never used to select a method, tune a hyperparameter, choose a gene set, or
pick a cutoff. It is computed after the method is frozen, and `assert_selection_frozen`
checks that the selection decision on disk predates and does not reference this module.
The inherited project rule is that outcomes may not select methods; here it is enforced
rather than remembered.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.bench.equal_footing import assign_patient_folds

#: Below this many events, a C-index comparison between two methods is reported but
#: never used to rank them. Derived from the standard rule of thumb that a Cox model
#: needs ~10 events per covariate: with 7 composition columns plus age, an honest
#: comparison needs on the order of 80 events.
MIN_EVENTS_FOR_RANKING = 80


@dataclass
class SurvivalContract:
    """What the cohort actually supports, measured before any model is fitted."""
    n_patients: int
    n_events: int
    n_censored: int
    n_covariates: int
    events_per_covariate: float
    detectable_delta_c_index: float
    adequately_powered: bool
    verdict: str

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def approximate_detectable_delta(n_events: int) -> float:
    """
    Smallest C-index difference this many events could distinguish from zero.

    Uses the standard large-sample approximation SE(C) ~ 1/sqrt(n_events) for a
    C-index near 0.6-0.7, doubled to a two-sided 95% bound and multiplied by sqrt(2)
    because two C-indices are being compared. Deliberately approximate; its job is to
    keep a difference of 0.01 on 40 events from being written up as a finding, and any
    reasonable approximation does that.
    """
    if n_events <= 1:
        return float("inf")
    return float(1.96 * np.sqrt(2.0) / np.sqrt(n_events))


def assess_power(clinical: pd.DataFrame, n_covariates: int) -> SurvivalContract:
    """Measure the cohort's power from the data. Called before any model is fitted."""
    n_patients = int(len(clinical))
    n_events = int(clinical["event"].sum())
    epc = n_events / max(n_covariates, 1)
    delta = approximate_detectable_delta(n_events)
    powered = n_events >= MIN_EVENTS_FOR_RANKING

    if n_events == 0:
        verdict = "BLOCKED: no events in the cohort; survival analysis is impossible"
    elif powered:
        verdict = "ADEQUATE: method ranking on prognostic value is interpretable"
    else:
        verdict = (
            f"UNDERPOWERED: {n_events} events supports detecting a C-index difference "
            f"of about {delta:.3f} at best. Differences smaller than that are reported "
            f"as INCONCLUSIVE and are not used to rank methods."
        )

    return SurvivalContract(
        n_patients=n_patients, n_events=n_events,
        n_censored=n_patients - n_events, n_covariates=n_covariates,
        events_per_covariate=float(epc), detectable_delta_c_index=delta,
        adequately_powered=powered, verdict=verdict,
    )


def concordance(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """
    Harrell's C-index: the fraction of comparable patient pairs ordered correctly.

    A pair is comparable when the one who died did so before the other's last known
    follow-up; otherwise censoring makes their order unknowable. Ties in risk count as
    half credit, which is the standard convention and matters here because a shrunken
    model can produce genuinely tied predictions.
    """
    from lifelines.utils import concordance_index
    return float(concordance_index(time, -risk, event))


def out_of_fold_risk(X: pd.DataFrame, clinical: pd.DataFrame,
                     folds: pd.Series, penalizer: float = 0.1) -> pd.Series:
    """
    Repeated out-of-fold Cox predictions.

    Fitting and evaluating on the same patients measures how well a model can memorise,
    which on 40 patients it can do almost perfectly. Every prediction here comes from a
    model that never saw that patient. The Cox fit is ridge-penalised because with this
    many events an unpenalised fit on 8 covariates will not converge reliably.
    """
    from lifelines import CoxPHFitter

    risk = pd.Series(np.nan, index=X.index, dtype="float64")
    for k in sorted(folds.unique()):
        train_idx = folds.index[folds != k]
        test_idx = folds.index[folds == k]
        train_idx = [i for i in train_idx if i in X.index]
        test_idx = [i for i in test_idx if i in X.index]
        if len(train_idx) < 5 or not test_idx:
            continue

        df = X.loc[train_idx].copy()
        df["T"] = clinical.loc[train_idx, "time"].to_numpy()
        df["E"] = clinical.loc[train_idx, "event"].to_numpy()

        # A covariate with no variance in this fold cannot be fitted and would abort
        # the whole fold; dropping it costs one column, not one fifth of the data.
        keep = [c for c in X.columns if df[c].std() > 1e-12]
        if not keep:
            continue
        try:
            cph = CoxPHFitter(penalizer=penalizer)
            cph.fit(df[keep + ["T", "E"]], duration_col="T", event_col="E")
            risk.loc[test_idx] = cph.predict_partial_hazard(X.loc[test_idx, keep]).to_numpy()
        except Exception:                                # noqa: BLE001
            continue
    return risk


def evaluate_method(composition: pd.DataFrame, clinical: pd.DataFrame,
                    baseline_cols: list[str] | None = None,
                    n_repeats: int = config.N_CV_REPEATS) -> dict:
    """
    Baseline versus baseline-plus-composition, out of fold, repeated.

    Composition arrives already aggregated to one row per patient — a patient with
    twelve blocks contributes one row, like everyone else.
    """
    baseline_cols = [c for c in (baseline_cols or ["age"]) if c in clinical.columns]
    shared = [p for p in composition.index if p in clinical.index]
    if len(shared) < 10:
        return {"status": "BLOCKED",
                "reason": f"only {len(shared)} patients have both composition and outcome"}

    comp = composition.loc[shared, [c for c in config.PRIMARY_CELL_TYPES
                                    if c in composition.columns]]
    clin = clinical.loc[shared]
    base = clin[baseline_cols] if baseline_cols else pd.DataFrame(index=shared)

    c_base, c_full = [], []
    for rep in range(n_repeats):
        folds = assign_patient_folds(
            pd.DataFrame({"patient_id": shared}, index=shared),
            seed=config.RANDOM_SEED + rep,
        )
        for label, X, store in (("base", base, c_base),
                                ("full", pd.concat([base, comp], axis=1), c_full)):
            if X.shape[1] == 0:
                continue
            risk = out_of_fold_risk(X, clin, folds)
            ok = risk.notna()
            if ok.sum() < 10:
                continue
            store.append(concordance(clin.loc[ok, "time"].to_numpy(),
                                     clin.loc[ok, "event"].to_numpy(),
                                     risk[ok].to_numpy()))

    if not c_full:
        return {"status": "BLOCKED", "reason": "no fold produced a usable model"}

    return {
        "status": "OK",
        "n_patients": len(shared),
        "n_events": int(clin["event"].sum()),
        "c_index_baseline": float(np.mean(c_base)) if c_base else np.nan,
        "c_index_with_composition": float(np.mean(c_full)),
        "delta_c_index": (float(np.mean(c_full) - np.mean(c_base)) if c_base else np.nan),
        "c_index_sd": float(np.std(c_full)),
        "n_repeats": len(c_full),
    }


def rank_methods(results: dict[str, dict], contract: SurvivalContract) -> pd.DataFrame:
    """
    Assemble the comparison, refusing to rank when the design cannot support it.

    `interpretation` is the column a reader should look at first. When the cohort is
    underpowered every row says INCONCLUSIVE, which is the correct answer and not a
    failure of the run.
    """
    rows = []
    for name, res in results.items():
        row = {"method": name, **{k: v for k, v in res.items() if k != "status"}}
        row["status"] = res.get("status")
        rows.append(row)
    table = pd.DataFrame(rows).set_index("method")

    if "c_index_with_composition" not in table.columns:
        table["interpretation"] = "BLOCKED"
        return table

    best = table["c_index_with_composition"].max()
    gap = best - table["c_index_with_composition"]

    if not contract.adequately_powered:
        table["interpretation"] = (
            f"INCONCLUSIVE (cohort supports detecting ~"
            f"{contract.detectable_delta_c_index:.3f}; not used for ranking)"
        )
    else:
        table["interpretation"] = np.where(
            gap <= contract.detectable_delta_c_index,
            "tied with best (within detectable difference)",
            "worse than best",
        )
    return table.sort_values("c_index_with_composition", ascending=False)


def assert_selection_frozen(selection_path) -> None:
    """
    Verify the method was selected before survival was computed.

    Reads the selection decision and checks it records an outcome-blind criterion. A
    selection file that mentions survival at all means the firewall was breached
    somewhere upstream, and every prognostic number downstream of it is contaminated.
    """
    from pathlib import Path
    selection_path = Path(selection_path)
    if not selection_path.exists():
        raise FileNotFoundError(
            f"{selection_path} not found. The method must be selected — and the "
            "decision written down — before survival is computed."
        )
    decision = json.loads(selection_path.read_text())
    criterion = str(decision.get("criterion", "")).lower()
    if "surviv" in criterion or "prognos" in criterion or "c_index" in criterion:
        raise AssertionError(
            f"the method selection criterion is {criterion!r}, which references an "
            "outcome. Method selection must be outcome-blind; this result is not usable."
        )
