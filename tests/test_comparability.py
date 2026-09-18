"""
A method that could not be evaluated must not enter a ranking, and must not vanish either.

The failure this guards against is specific and has happened: the invariant "never report a
degenerate method under its own name" was honoured in prose, then violated whenever a ranking
was computed over every method that returned finite numbers.
"""
from __future__ import annotations

import pytest

from ivygap.deconv import comparability as cmp


def test_bisque_is_excluded_under_both_references():
    """Bisque's requirement is a property of the COHORT, so no reference choice fixes it."""
    for ref in ("frozen", "h5ad"):
        assert "bisque" in cmp.non_comparable(ref), ref


def test_the_h5ad_reference_restores_the_variance_weighted_methods():
    """Only the methods whose missing input WAS the reference are restored by rebuilding it."""
    frozen = cmp.non_comparable("frozen")
    h5ad = cmp.non_comparable("h5ad")
    for m in ("music", "epic", "cibersortx_smode"):
        assert m in frozen, f"{m} should be excluded on the zero-sigma frozen signature"
        assert m not in h5ad, f"{m} should be comparable once sigma is present"
    # and the cohort-level exclusion is untouched by that
    assert "bisque" in h5ad


def test_scdc_ensemble_is_NOT_restored_by_rebuilding_the_reference():
    """CORRECTED 2026-09-18. This test previously asserted the opposite, and was wrong.

    ENSEMBLE needs two or more references to weight across. Rebuilding the one reference from
    .h5ad gives it a BETTER single reference, not a second one -- `DeconvolutionInput` is
    constructed with `references=(ref,)` on both paths. The error was caught empirically: the
    h5ad GBM run reported `scdc_ensemble [DEGENERATE]` with rho identical to `scdc` to four
    decimals, while `non_comparable("h5ad")` would have admitted it to the ranking under its
    own name.

    This is not a test relaxed to accommodate a failure; it is a test that encoded a false
    premise, corrected against what the run measured.
    """
    for ref in ("frozen", "h5ad"):
        why = cmp.non_comparable(ref).get("scdc_ensemble")
        assert why is not None, f"ENSEMBLE must stay excluded on {ref}"
        assert "two or more references" in why


def test_observed_degeneracy_is_unioned_in_on_the_h5ad_reference_too():
    """The static table says what a method REQUIRES; the run says what HAPPENED.

    A method that degenerates for an unanticipated reason must still be excluded, on the
    permissive reference as well as the restrictive one -- otherwise the h5ad ranking, which is
    the defensible one, is the one with no runtime safety net.
    """
    out = cmp.non_comparable("h5ad", degenerate={"dwls"}, failed={"bayesprism"})
    assert "dwls" in out and "bayesprism" in out
    # and it does not quietly displace the requirement-based reason for a method in both
    both = cmp.non_comparable("h5ad", degenerate={"bisque"})
    assert "requires" in both["bisque"], "the specific requirement must survive the union"


def test_every_exclusion_names_the_input_it_needs():
    """A bare 'excluded' is the thing this module exists to prevent."""
    for ref in ("frozen", "h5ad"):
        for m, why in cmp.non_comparable(ref).items():
            assert len(why) > 40, f"{m} on {ref}: reason too thin to act on"
            if m in cmp.INPUT_REQUIREMENTS:
                assert "requires" in why or "observed" in why or "failed" in why


def test_observed_degeneracy_is_unioned_in_not_trusted_away():
    """A method degenerating for an UNANTICIPATED reason must still be excluded."""
    bad = cmp.non_comparable("h5ad", degenerate={"dwls"}, failed={"svr"})
    assert "dwls" in bad and "svr" in bad
    assert cmp.comparable(["nnls", "dwls", "svr"], "h5ad",
                          degenerate={"dwls"}, failed={"svr"}) == ["nnls"]


def test_comparable_preserves_order_and_drops_only_the_excluded():
    ms = ["nnls", "bisque", "svr", "music", "cibersortx"]
    assert cmp.comparable(ms, "frozen") == ["nnls", "svr", "cibersortx"]
    assert cmp.comparable(ms, "h5ad") == ["nnls", "svr", "music", "cibersortx"]
