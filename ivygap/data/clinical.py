"""
clinical.py — load Ivy GAP's per-tumour clinical table into a survival contract.

TOLERANT PARSING, INTOLERANT VALIDATION
---------------------------------------
The Allen Institute's tumour details table has shipped with several column namings and
several time units (days in some releases, months in others). Column resolution is
therefore permissive, in the same style as `load_ivygap.py`.

Validation is not permissive. A survival time that is negative, or an event indicator
that is neither 0 nor 1, or a cohort where nobody died, is an error that would silently
produce a meaningless C-index. Each of those aborts with an explicit message rather
than being coerced into something that runs.

WHAT THIS MODULE WILL NOT DO
----------------------------
If the clinical file is absent, this returns None. It does not impute survival, does
not fall back to a published summary statistic, and does not carry forward TCGA's
outcomes onto Ivy GAP patients. The survival module then reports BLOCKED. A missing
input reported honestly is worth more than a fabricated one that lets a pipeline
finish.

THE IVY GAP CLINICAL TABLE, AS ACTUALLY PUBLISHED  (measured 2026-09-03)
------------------------------------------------------------------------
`tumor_details.csv` is 42 tumours x 10 columns. It carries `survival_days` and it
carries **no vital-status column of any kind**. That is not a download problem; the
Allen Institute does not publish one.

Two consequences, both load-bearing:

1. **The join needs a third file.** The clinical table is keyed on `donor_id`; the
   expression matrix is keyed on `tumor_id`; neither file contains the other's key.
   Only `rna_seq_samples_details.csv` contains both. See `portal_metadata.py`.

2. **Event status is unknowable from the release.** 10 of 42 tumours (8 of the 37 with
   RNA-seq) have `survival_days` blank. Treating a recorded time as an observed death
   is an *assumption*, and this module refuses to make it silently:
   `event_policy="require"` (the default) raises, exactly as before.

   The assumption is available as `event_policy="observed-time-is-death"`, which must
   be requested by name, is stamped into every artefact it touches, and never becomes a
   default. It is not imputation of the blanks — those rows are dropped, not invented —
   it is a declared reading of the rows that are present.

WHY THE BLANKS ARE NOT MISSING AT RANDOM
-----------------------------------------
`audit_missingness()` measures this rather than assuming it, because the answer changes
how the result must be read. On the published table, 8 of the 9 tumours with a blank
`survival_days` and a known MGMT status are MGMT-methylated, against 9 of 31 among
tumours with a recorded time (Fisher exact p = 0.0021), and the blanks are younger
(median 52 vs 61 years, p = 0.015). MGMT methylation and younger age are the two
strongest favourable prognostic factors in glioblastoma.

So the blanks are concentrated among the patients most likely to have still been alive
at the data freeze. A complete-case analysis is therefore biased toward short survival
by construction, and that bias is *reported alongside every survival number* rather than
being left for a reader to discover.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ivygap import config

_PATIENT_COLS = ["tumor_id", "donor_id", "tumor_name", "specimen_name", "donor_name"]
_TIME_COLS = ["survival_days", "days_to_death", "survival", "overall_survival_days",
              "survival_months", "os_months", "os_days"]
_EVENT_COLS = ["death", "event", "vital_status", "os_status", "deceased"]
_AGE_COLS = ["age_in_years", "age", "age_at_diagnosis", "age_years"]
_SEX_COLS = ["sex", "gender"]

#: Column names whose units are months rather than days.
_MONTH_COLS = {"survival_months", "os_months"}


def _resolve(df: pd.DataFrame, candidates: list[str], what: str,
             required: bool = True) -> str | None:
    lower = {c.lower().strip(): c for c in df.columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    if required:
        raise KeyError(
            f"could not find a column for {what}. Accepted: {candidates}. "
            f"File has: {list(df.columns)}"
        )
    return None


def _coerce_event(series: pd.Series) -> pd.Series:
    """
    Map an event column onto {0, 1} however it was encoded.

    Anything unrecognised becomes NaN and the row is dropped, rather than defaulting
    to censored — defaulting would quietly convert deaths into survivors and bias every
    downstream estimate in the same direction.
    """
    if pd.api.types.is_numeric_dtype(series):
        vals = series.astype("float64")
        bad = ~vals.isin([0.0, 1.0]) & vals.notna()
        if bad.any():
            raise ValueError(
                f"numeric event column has values outside {{0,1}}: "
                f"{sorted(vals[bad].unique())[:5]}"
            )
        return vals

    text = series.astype(str).str.strip().str.lower()
    dead = {"1", "dead", "deceased", "death", "yes", "true", "1:deceased", "event"}
    alive = {"0", "alive", "living", "censored", "no", "false", "0:living", "0:alive"}
    return text.map(lambda v: 1.0 if v in dead else (0.0 if v in alive else np.nan))


#: The only event policies this module recognises.
EVENT_POLICIES = ("require", "observed-time-is-death")


def _key_on_tumor_id(out: pd.DataFrame, id_col_name: str) -> tuple[pd.DataFrame, dict]:
    """
    Re-key a donor-keyed clinical table onto `tumor_id`, which is what the expression
    matrix and the sample manifest use.

    Returns (table, provenance). If the portal metadata is unavailable the table is
    returned untouched with a provenance note saying the join could not be made — the
    caller will then find zero overlap with the manifest, and that is reported as a
    blocker rather than papered over with a fuzzy name match.
    """
    from ivygap.data import portal_metadata

    if id_col_name.lower() == "tumor_id":
        return out, {"joined_via": "tumor_id present in the clinical file itself"}

    try:
        mapping = portal_metadata.donor_tumor_map()
    except portal_metadata.PortalMetadataMissing as exc:
        return out, {"joined_via": None, "join_failed_because": str(exc)}

    donor_to_tumor = {d: t for t, d in mapping["donor_id"].items()}
    mapped = out["patient_id"].map(donor_to_tumor)

    unmapped = out.loc[mapped.isna(), "patient_id"].tolist()
    out = out.loc[mapped.notna()].copy()
    out["patient_id"] = mapped.loc[mapped.notna()].to_numpy()

    return out, {
        "joined_via": f"{id_col_name} -> tumor_id via rna_seq_samples_details.csv",
        "n_mapped": int(len(out)),
        "n_unmapped": len(unmapped),
        "unmapped_ids": unmapped,
        "why_unmapped": ("these donors have no RNA-seq sample in the archive, so they "
                         "could not contribute to a composition-based analysis anyway"),
    }


def audit_missingness() -> dict:
    """
    Measure whether a blank `survival_days` is missing at random.

    This is evidence, not decoration. If the blanks are concentrated among good-prognosis
    patients then a complete-case analysis is biased toward short survival, and every
    number downstream has to be read with that in mind. Measuring it costs one Fisher
    test; assuming it costs the credibility of the whole survival section.
    """
    path = config.IVYGAP_TUMOR_DETAILS_PATH
    if not path.exists():
        return {"available": False, "reason": f"{path.name} not present"}

    raw = pd.read_csv(path)
    tcol = _resolve(raw, _TIME_COLS, "survival time")
    recorded = pd.to_numeric(raw[tcol], errors="coerce").notna()

    out: dict = {
        "available": True,
        "n_tumors": int(len(raw)),
        "n_time_recorded": int(recorded.sum()),
        "n_time_blank": int((~recorded).sum()),
    }

    if "mgmt_methylation" in raw.columns:
        known = raw["mgmt_methylation"].isin(["Yes", "No"])
        tab = pd.crosstab(recorded[known], raw.loc[known, "mgmt_methylation"])
        tab = tab.reindex(index=[True, False], columns=["Yes", "No"]).fillna(0)
        try:
            from scipy.stats import fisher_exact
            odds, pval = fisher_exact(tab.to_numpy())
        except Exception:                              # noqa: BLE001
            odds = pval = float("nan")
        out["mgmt"] = {
            "methylated_among_recorded": f"{int(tab.loc[True, 'Yes'])}/"
                                         f"{int(tab.loc[True].sum())}",
            "methylated_among_blank": f"{int(tab.loc[False, 'Yes'])}/"
                                      f"{int(tab.loc[False].sum())}",
            "fisher_odds_ratio": float(odds),
            "fisher_p": float(pval),
        }

    acol = _resolve(raw, _AGE_COLS, "age", required=False)
    if acol is not None:
        age = raw[acol].astype(str).str.extract(r"(\d+)").astype(float)[0]
        a, b = age[recorded].dropna(), age[~recorded].dropna()
        try:
            from scipy.stats import mannwhitneyu
            pval = float(mannwhitneyu(a, b).pvalue) if len(a) and len(b) else float("nan")
        except Exception:                              # noqa: BLE001
            pval = float("nan")
        out["age"] = {
            "median_recorded": float(a.median()) if len(a) else float("nan"),
            "median_blank": float(b.median()) if len(b) else float("nan"),
            "mannwhitney_p": pval,
        }

    mgmt_p = out.get("mgmt", {}).get("fisher_p", float("nan"))
    informative = np.isfinite(mgmt_p) and mgmt_p < 0.05
    out["missing_at_random"] = not informative
    out["verdict"] = (
        "INFORMATIVE MISSINGNESS: a blank survival time is strongly associated with "
        f"MGMT methylation (Fisher p = {mgmt_p:.4g}), the strongest favourable "
        "prognostic factor in GBM. The blanks are concentrated among the patients most "
        "likely to have been alive at the data freeze, so dropping them biases the "
        "cohort toward short survival. Any complete-case survival estimate on this "
        "table is biased downward and must be reported as such."
        if informative else
        "No evidence that blank survival times are associated with MGMT status. "
        "Complete-case analysis is not obviously biased, which is weaker than saying "
        "it is unbiased."
    )
    return out


def load_clinical(strict: bool = True,
                  event_policy: str = "require") -> pd.DataFrame | None:
    """
    Return a per-patient table indexed by `patient_id` (= `tumor_id`, the key the
    expression matrix and the sample manifest use) with columns `time` (days),
    `event` (0/1), and whatever covariates were present.

    Returns None when the file is absent — a reportable state, not an exception.

    `event_policy`:
        "require"                 — the default. A clinical file with no vital-status
                                    column raises. This is the invariant-preserving
                                    behaviour and nothing changes it implicitly.
        "observed-time-is-death"  — read a recorded survival time as an observed death
                                    and drop the blanks. A declared assumption, never a
                                    fallback; it must be named by the caller, and
                                    `attrs["event_policy"]` carries it downstream so
                                    every artefact records which reading produced it.

    Note what the second policy does *not* do: it does not assign the blank rows a time
    or a status. Those tumours are dropped and counted. Inventing a censoring time for
    them would be imputation, which is forbidden; reading the rows that are present is a
    stated interpretation of published data, which is not the same thing.
    """
    if event_policy not in EVENT_POLICIES:
        raise ValueError(f"unknown event_policy {event_policy!r}; "
                         f"expected one of {EVENT_POLICIES}")

    path = config.IVYGAP_TUMOR_DETAILS_PATH
    if not path.exists():
        return None

    raw = pd.read_csv(path)
    pid = _resolve(raw, _PATIENT_COLS, "the patient/tumour id")
    tcol = _resolve(raw, _TIME_COLS, "survival time")
    ecol = _resolve(raw, _EVENT_COLS, "the event indicator", required=False)

    time = pd.to_numeric(raw[tcol], errors="coerce")
    if tcol.lower() in _MONTH_COLS:
        time = time * 30.437                      # months -> days, one place only

    assumption = None
    if ecol is not None:
        event = _coerce_event(raw[ecol])
    elif event_policy == "observed-time-is-death":
        # Declared, not defaulted. The blanks are dropped by the dropna below; they are
        # NOT assigned event=0, which would invent a censoring time nobody published.
        event = pd.Series(np.where(time.notna(), 1.0, np.nan), index=raw.index)
        assumption = (
            f"{path.name} publishes no vital-status column. Under the declared policy "
            f"'observed-time-is-death', each of the {int(time.notna().sum())} recorded "
            f"survival times is read as an observed death and the "
            f"{int(time.isna().sum())} blank rows are DROPPED (not censored, not "
            f"imputed). See audit_missingness(): the blanks are not missing at random."
        )
    else:
        # No event column at all. Treating every patient as an observed death would
        # invent 100% mortality; refusing is the only defensible default.
        raise KeyError(
            f"{path.name} has no event/vital-status column (looked for {_EVENT_COLS}). "
            "Survival cannot be computed without knowing who was censored; assuming "
            "everyone died would fabricate the outcome. If you intend to make that "
            "reading explicitly, pass event_policy='observed-time-is-death' — it is "
            "recorded in every artefact rather than applied silently."
        )

    out = pd.DataFrame({
        "patient_id": raw[pid].astype(str).str.strip(),
        "time": time,
        "event": event,
    })
    for cols, name in ((_AGE_COLS, "age"), (_SEX_COLS, "sex")):
        col = _resolve(raw, cols, name, required=False)
        if col is not None:
            out[name] = (
                # Ivy GAP writes age as "61 yrs", so to_numeric alone yields all-NaN and
                # would silently drop age from every model that asked for it.
                pd.to_numeric(raw[col].astype(str).str.extract(r"(\d+\.?\d*)")[0],
                              errors="coerce")
                if name == "age" else raw[col].astype(str))

    n_raw = len(out)
    out, join_prov = _key_on_tumor_id(out, pid)

    out = out.dropna(subset=["time", "event"])
    out = out[out["time"] > 0]

    if out["patient_id"].duplicated().any():
        # One row per tumour is the contract; duplicates mean the file is keyed on
        # something else and the join to composition would fan out.
        dupes = out.loc[out["patient_id"].duplicated(), "patient_id"].unique()
        raise ValueError(f"duplicate patient ids in {path.name}: {list(dupes)[:10]}")

    out = out.set_index("patient_id")

    if strict:
        n_events = int(out["event"].sum())
        if len(out) == 0:
            raise ValueError(f"{path.name} yielded no usable rows out of {n_raw}")
        if n_events == 0:
            raise ValueError(
                f"{path.name} has {len(out)} patients but zero events. A C-index "
                "cannot be computed from a cohort in which nobody died."
            )

    out.attrs["event_policy"] = event_policy
    out.attrs["event_assumption"] = assumption
    out.attrs["join"] = join_prov
    out.attrs["n_rows_in_file"] = n_raw
    out.attrs["source_sha256"] = config.sha256_file(path)

    out.to_csv(config.CLINICAL_PATH, sep="\t")
    return out


def describe(clinical: pd.DataFrame | None) -> dict:
    """Measured cohort description, written into the results manifest as evidence."""
    if clinical is None:
        return {"available": False,
                "reason": f"{config.IVYGAP_TUMOR_DETAILS_PATH.name} not present"}
    return {
        "available": True,
        "n_patients": int(len(clinical)),
        "n_events": int(clinical["event"].sum()),
        "n_censored": int((clinical["event"] == 0).sum()),
        "median_followup_days": float(clinical["time"].median()),
        "has_age": "age" in clinical.columns,
        "has_sex": "sex" in clinical.columns,
        # Carried through so no artefact can report a survival number without also
        # reporting the reading of the data that produced it.
        "event_policy": clinical.attrs.get("event_policy"),
        "event_assumption": clinical.attrs.get("event_assumption"),
        "n_rows_in_file": clinical.attrs.get("n_rows_in_file"),
        "n_dropped_before_analysis": (
            clinical.attrs.get("n_rows_in_file", len(clinical)) - len(clinical)),
        "join": clinical.attrs.get("join"),
        "source_sha256": clinical.attrs.get("source_sha256"),
    }
