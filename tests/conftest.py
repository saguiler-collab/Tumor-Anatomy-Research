"""Shared fixtures. Deliberately small so the suite runs in seconds, not minutes."""

from __future__ import annotations

import pandas as pd
import pytest

from ivygap.data.reference import build_reference, build_synthetic, select_signature_genes
from ivygap.bench import pseudobulk as pb
from ivygap.deconv.base import DeconvolutionInput


@pytest.fixture(scope="session")
def synthetic():
    """A small synthetic single-cell reference with known generative structure."""
    bundle, expression, meta = build_synthetic(
        n_genes=400, n_donors=6, cells_per_type_per_donor=12)
    return bundle, expression, meta


@pytest.fixture(scope="session")
def donor_split(synthetic):
    _, expression, meta = synthetic
    test_set, train_donors, test_donors = pb.build_train_test(expression, meta, n_test=24)
    return test_set, train_donors, test_donors


@pytest.fixture(scope="session")
def problem(synthetic, donor_split):
    """A DeconvolutionInput built the way the real pipeline builds one."""
    _, expression, meta = synthetic
    test_set, train_donors, _ = donor_split

    train_cells = meta.index[meta["donor"].astype(str).isin(train_donors)]
    ref = build_reference(expression[train_cells], meta.loc[train_cells], name="test_ref")
    genes = [g for g in select_signature_genes(ref, n_per_type=25)
             if g in test_set.expression.index]

    manifest = pd.DataFrame(
        {"patient_id": test_set.donors.astype(str),
         "structure": test_set.niche.astype(str)},
        index=test_set.expression.columns,
    )
    data = DeconvolutionInput(bulk=test_set.expression.loc[genes],
                              references=(ref.subset_genes(genes),), manifest=manifest)
    return data, test_set
