"""
Which methods may enter a quantitative comparison, and which merely ran.

WHY THIS IS A MODULE AND NOT A SENTENCE IN A WRITE-UP
-----------------------------------------------------
There are two entirely different reasons a method can produce a poor number:

  1. **It was evaluated properly and did badly.** That is a result about the method.
  2. **Its intended input was not available, so what ran was not the published method.**
     That is a result about the EXPERIMENT, and reporting it as (1) is a false claim about
     someone else's work.

This project's invariant — *never report a degenerate method under its own name* — already
covers (2). What was missing was an enforcement point: the invariant was honoured in prose and
then quietly violated whenever a ranking was computed over "all methods".

So comparability is declared here, once, with the input requirement that fails. Every ranking
imports it. A method excluded here is still REPORTED — with its numbers — under availability,
never deleted and never silently dropped.

THE RULE
--------
A method is comparable when the inputs its published algorithm requires were present. Not when
it happened to return finite numbers.
"""
from __future__ import annotations

#: method -> (the input it requires, why that makes it non-comparable when absent)
#:
#: Keyed by what the method NEEDS, so the same table answers "is this fixable here?".
INPUT_REQUIREMENTS: dict[str, tuple[str, str]] = {
    "bisque": (
        "subjects assayed as BOTH bulk and single cells",
        "Bisque's method IS the assay transform learned from paired subjects. Without overlap "
        "it estimates that transform from marginal distributions, which is its documented "
        "degraded mode. TCGA has no single-cell data, so this is not fixable from any file on "
        "disk -- it needs a genuinely paired cohort.",
    ),
    "music": (
        "cross-donor variance (sigma) in the reference",
        "MuSiC's entire contribution is weighting genes by cross-donor consistency. With a "
        "zero sigma it is residual-weighted NNLS, and returns numbers identical to NNLS.",
    ),
    "scdc_ensemble": (
        "two or more references to weight across",
        "With one reference there is nothing to ensemble; it reduces exactly to SCDC.",
    ),
    "epic": (
        "refProfiles.var (per-gene variance) in the reference",
        "EPIC's published contribution is variance weighting. Without it the package warns and "
        "runs constrained least squares with uniform weights.",
    ),
    "quantiseq": (
        "its own TIL10 signature and a matching gene space",
        "quanTIseq ignores the supplied reference. It is reported separately wherever it runs, "
        "and it returned no finite estimates on either TCGA cohort.",
    ),
    "cibersortx_smode": (
        "the single cells the signature was built from",
        "S-mode's batch correction is defined against the source cells; without them the "
        "method cannot run at all.",
    ),
}


def non_comparable(reference: str, degenerate: set[str] | None = None,
                   failed: set[str] | None = None) -> dict[str, str]:
    """
    Methods excluded from a quantitative comparison on this reference, with the reason.

    `reference` is "frozen" or "h5ad". The frozen signature carries a zero sigma, so it
    disables MuSiC and EPIC and prevents S-mode running. The h5ad-built reference supplies
    sigma, donor profiles and registered cells, which restores those three.

    Bisque, quanTIseq and SCDC ENSEMBLE stay excluded under BOTH references, because their
    missing inputs are not properties of the reference build:

    - Bisque needs subjects assayed as both bulk and single cells. That is a property of the
      COHORT, and TCGA does not have it.
    - quanTIseq needs its own TIL10 signature and a matching gene space.
    - SCDC ENSEMBLE needs two or more references to weight across. **This pipeline supplies
      exactly one** (`DeconvolutionInput(references=(ref,))`) no matter how that one is built,
      so rebuilding from .h5ad does not help. Corrected 2026-09-18: ENSEMBLE was previously
      listed as frozen-only, which would have admitted it to the h5ad ranking under its own
      name while the run itself reported it [DEGENERATE] with numbers identical to SCDC.

    `degenerate` and `failed` are what the run actually observed. They are unioned in rather
    than trusted alone, so a method that degenerates for an unanticipated reason is still
    excluded instead of silently entering the ranking.
    """
    out: dict[str, str] = {}
    always = ["bisque", "quantiseq", "scdc_ensemble"]
    frozen_only = ["music", "epic", "cibersortx_smode"]
    for m in always + (frozen_only if reference == "frozen" else []):
        need, why = INPUT_REQUIREMENTS[m]
        out[m] = f"requires {need}; {why}"
    for m in sorted(degenerate or set()):
        out.setdefault(m, "observed degenerate in this run")
    for m in sorted(failed or set()):
        out.setdefault(m, "failed in this run")
    return out


def comparable(methods, reference: str, degenerate=None, failed=None) -> list[str]:
    """The subset of `methods` that may enter a quantitative comparison, order preserved."""
    bad = non_comparable(reference, degenerate, failed)
    return [m for m in methods if m not in bad]
