"""The signature matrix's cell-type columns must contain that cell type's markers.

WHY THIS EXISTS. The study's strongest claim is that 0 of 12 methods reproduce the
methylation-measured T > B ordering. There is exactly one way that result could be an
artefact rather than a finding: if the harness had the **T and B columns transposed**, every
method would appear to put B above T no matter what it actually estimated.

Correlating estimates against methylation truth cannot settle it — estimated-B and
methylation-T both rise with overall lymphocyte abundance, so a positive association appears
with or without a swap. The decisive check is at the reference itself: canonical lineage
markers must score highest in their own column.

Measured 2026-09-23: 9 of 9 B-cell markers land in `B_cell`, by margins of two to four
orders of magnitude (CD79A: B 3038.3 vs T 0.79). No T-cell marker is highest in `B_cell`.
The claim is not a labelling artefact.

The T/NK boundary is deliberately NOT asserted to be clean. CD3D, CD3E, CD3G, CD2 and LCK
score highest in `NK_cell` in this atlas, which is expected rather than wrong: NK and T share
the lymphoid lineage and brain-tumour NK annotations routinely include CD3-positive NKT and
cytotoxic populations. That entanglement is a real property of the reference and is recorded
here so nobody later mistakes it for a discovery -- but it does not touch the T-versus-B
claim, which is what this file guards.
"""
from __future__ import annotations

import pandas as pd
import pytest

from ivygap import config

#: Canonical lineage markers. Chosen for specificity, not for how well they behave here.
B_MARKERS = ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "BLK", "PAX5", "TNFRSF13C"]
T_MARKERS = ["CD28", "IL7R", "CD8A", "THEMIS", "CD3D", "CD3E", "CD3G", "CD2", "LCK"]


@pytest.fixture(scope="module")
def sig():
    p = config.FROZEN_SIGNATURE_PATH
    if not p.exists():
        pytest.skip("frozen signature not present in this checkout")
    return pd.read_csv(p, sep="\t", index_col=0)


def test_the_signature_has_the_lymphoid_columns_the_claim_needs(sig):
    for c in ("T_cell", "B_cell", "NK_cell"):
        assert c in sig.columns, f"{c} missing; the lymphoid claim cannot be posed"


def test_every_b_cell_marker_is_highest_in_the_b_cell_column(sig):
    present = [g for g in B_MARKERS if g in sig.index]
    assert len(present) >= 6, f"only {len(present)} B markers in the signature; too few to judge"
    wrong = [(g, sig.loc[g].idxmax()) for g in present if sig.loc[g].idxmax() != "B_cell"]
    assert not wrong, (
        "a canonical B-cell marker is not highest in B_cell, which would mean the column "
        f"labels cannot be trusted: {wrong}")


def test_no_t_cell_marker_is_highest_in_the_b_cell_column(sig):
    """The specific failure that would fake the headline result."""
    present = [g for g in T_MARKERS if g in sig.index]
    assert len(present) >= 6, f"only {len(present)} T markers present; too few to judge"
    wrong = [(g, float(sig.loc[g, "B_cell"])) for g in present
             if sig.loc[g].idxmax() == "B_cell"]
    assert not wrong, (
        "a T-cell marker scores highest in B_cell — T and B may be transposed, which would "
        f"make the 0-of-12 lymphoid result an artefact: {wrong}")


def test_the_b_column_dominates_by_a_wide_margin_not_a_narrow_one(sig):
    """A marginal win would be consistent with noise; these should not be marginal."""
    present = [g for g in B_MARKERS if g in sig.index]
    thin = []
    for g in present:
        b = float(sig.loc[g, "B_cell"])
        t = float(sig.loc[g, "T_cell"])
        if b <= 0 or b < 10 * max(t, 1e-9):
            thin.append((g, b, t))
    assert not thin, (
        "a B marker's B_cell value is not at least 10x its T_cell value; the separation is "
        f"too thin to license the T-versus-B claim: {thin}")


def test_the_t_nk_entanglement_is_recorded_rather_than_assumed_absent(sig):
    """Documents a known property so a later reader does not rediscover it as a defect.

    This test cannot fail on the entanglement itself -- it exists to fail if someone later
    'cleans' the reference and the recorded description stops matching what is on disk.
    """
    shared = [g for g in ("CD3D", "CD3E", "CD3G", "CD2", "LCK") if g in sig.index]
    if not shared:
        pytest.skip("pan-T markers absent from this signature")
    in_nk = [g for g in shared if sig.loc[g].idxmax() == "NK_cell"]
    assert len(in_nk) >= 1, (
        "no pan-T marker scores highest in NK_cell any more. That may be an improvement, but "
        "tests/test_signature_label_integrity.py and docs/STATISTICAL_VERIFICATION.md both "
        "describe the T/NK entanglement as present — update them to match the reference.")
