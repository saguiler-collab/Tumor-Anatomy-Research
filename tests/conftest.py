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

# =============================================================================
# No test may write into the real results tree
# -----------------------------------------------------------------------------
# Two tests added for the donor-leakage guard called run_benchmark.run() directly.
# That function writes benchmark_summary.csv, method_selection_decision.json and the
# pseudobulk estimates to config.BENCH_DIR / config.ESTIMATES_DIR — the real ones. Run
# alongside a pipeline run, they overwrote its stage-3 artefacts with 10-mixture fixture
# output, and stage 4 reads benchmark_summary.csv as the ground-truth yardstick for the
# study's primary result.
#
# It was caught because the recorded donors were D00-D05 where the run had 110 real
# ones. It would not have been caught by a failing test — every test passed.
#
# The run lock cannot help here: it guards one pipeline run against another, and a
# pytest process is neither. So this is autouse and session-scoped: every test in the
# suite gets its own tree, and a test that writes to a config path writes to a temp
# directory whether or not its author remembered to redirect it.
# =============================================================================

@pytest.fixture(scope="session")
def _test_output_root(tmp_path_factory):
    return tmp_path_factory.mktemp("ivygap_test_outputs")


@pytest.fixture(autouse=True)
def _isolate_output_tree(_test_output_root):
    """
    Function-scoped, not session-scoped, and deliberately so.

    `test_config_separates_synthetic_from_real_paths` calls importlib.reload(config) to
    inspect the module as the file defines it. A reload resets every attribute, wiping a
    session-scoped redirection and leaving the REAL paths in place for every test that
    runs afterwards. Re-applying per test makes the isolation robust against any test
    that mutates, reloads or rebinds the config module — which is the only assumption
    worth making about a module whose whole job is to be mutated.
    """
    from ivygap import config

    root = _test_output_root
    originals = {}
    redirected = {
        "RESULTS_DIR": root / "results",
        "ESTIMATES_DIR": root / "results" / "estimates",
        "BENCH_DIR": root / "results" / "benchmark",
        "ANATOMIC_DIR": root / "results" / "anatomic",
        "SURVIVAL_DIR": root / "results" / "survival",
        "FIGURES_DIR": root / "results" / "figures",
        "RELEASE_DIR": root / "release",
        "PSEUDOBULK_DIR": root / "pseudobulk",
    }
    for name, value in redirected.items():
        originals[name] = getattr(config, name)
        setattr(config, name, value)
        value.mkdir(parents=True, exist_ok=True)

    originals["ALL_OUTPUT_DIRS"] = config.ALL_OUTPUT_DIRS
    config.ALL_OUTPUT_DIRS = list(redirected.values())

    yield root

    for name, value in originals.items():
        setattr(config, name, value)
