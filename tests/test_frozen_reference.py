"""
Tests for the vendored frozen reference (protocol step 1) and for the degeneracy
reporting it makes necessary.

The theme: a collapsed signature matrix cannot support every method in the roster, and
the failure mode that matters is not a crash — it is a method running to completion
while no longer being itself.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.data.reference import load_frozen_reference


@pytest.fixture(scope="module")
def frozen():
    if not config.FROZEN_SIGNATURE_PATH.exists():
        pytest.skip("reference_frozen/signature_matrix.tsv is not present")
    return load_frozen_reference()


def test_frozen_reference_matches_the_project_roster(frozen):
    assert frozen.cell_types == list(config.CELL_TYPES)
    assert frozen.profile.shape[0] > 10_000
    assert (frozen.profile.to_numpy() >= 0).all(), "the mixing model is additive"


def test_cell_sizes_are_positive_and_complete(frozen):
    assert list(frozen.cell_size.index) == list(config.CELL_TYPES)
    assert (frozen.cell_size > 0).all()
    assert frozen.cell_size.notna().all()


def test_frozen_reference_declares_it_has_no_cross_donor_variance(frozen):
    """
    The honest half of loading a collapsed matrix. Zeroing sigma silently would turn
    MuSiC into NNLS while still printing 'MuSiC' on a leaderboard.
    """
    assert frozen.has_cross_donor_variance is False
    assert (frozen.sigma.to_numpy() == 0).all()


def test_the_flag_survives_gene_subsetting(frozen):
    """A subset that quietly re-enabled the flag would defeat the whole mechanism."""
    genes = list(frozen.profile.index[:500])
    assert frozen.subset_genes(genes).has_cross_donor_variance is False


def test_music_reports_itself_degenerate_on_a_collapsed_reference(frozen):
    """The point of the flag: the method must say so, not just behave differently."""
    from ivygap.deconv.base import DeconvolutionInput
    from ivygap.deconv.reference_based import MuSiCDeconvolution

    genes = list(frozen.profile.index[:400])
    ref = frozen.subset_genes(genes)
    rng = np.random.default_rng(0)
    bulk = pd.DataFrame(
        rng.gamma(2.0, 50.0, size=(len(genes), 6)),
        index=genes, columns=[f"s{i}" for i in range(6)])
    mf = pd.DataFrame({"patient_id": ["P0"] * 3 + ["P1"] * 3,
                       "structure": ["CT", "LE", "MVP"] * 2}, index=bulk.columns)

    m = MuSiCDeconvolution()
    m.fit_predict(DeconvolutionInput(bulk=bulk, references=(ref,), manifest=mf))
    assert m.degenerate_ is True
    assert "cross-donor" in m.degeneracy_reason_


def test_music_is_not_degenerate_with_a_real_reference(problem):
    """The flag must discriminate; always-true would be as useless as always-false."""
    from ivygap.deconv.reference_based import MuSiCDeconvolution
    data, _ = problem
    m = MuSiCDeconvolution()
    m.fit_predict(data)
    assert m.degenerate_ is False
    assert m.degeneracy_reason_ is None


def test_provenance_records_hashes_and_limitations():
    """
    The vendored inputs must be traceable and their limits stated. A signature matrix
    committed without provenance is an unsourced number in every downstream result.
    """
    path = Path(config.FROZEN_REFERENCE_DIR) / "PROVENANCE.json"
    if not path.exists():
        pytest.skip("PROVENANCE.json not present")
    prov = json.loads(path.read_text())

    assert prov["files"], "no files recorded"
    for name, info in prov["files"].items():
        assert len(info["sha256"]) == 64, f"{name} has no usable hash"
        f = Path(config.FROZEN_REFERENCE_DIR) / name
        if f.exists():
            assert config.sha256_file(f) == info["sha256"], (
                f"{name} does not match its recorded hash — the vendored input changed"
            )
    assert any("cross-donor" in lim for lim in prov["known_limitations"])


# =============================================================================
# the gene space the frozen path actually hands to the methods
# -----------------------------------------------------------------------------
# `config.USE_SIGNATURE_GENE_SUBSET` was declared, documented as load-bearing, and
# applied only inside run_benchmark — which does not run without a cell-level reference.
# On the configuration this project ships in, the flag was therefore never read and every
# method received the full 16,007-gene intersection. These tests make the flag's effect
# checkable instead of implied.
# =============================================================================

def test_frozen_gene_space_honours_the_subset_switch():
    from ivygap import config as cfg
    from ivygap.data.reference import frozen_gene_space, load_frozen_reference

    frozen = load_frozen_reference()
    bulk_index = frozen.profile.index                    # every gene is "shared"

    genes, prov = frozen_gene_space(frozen, bulk_index)
    assert prov["n_shared"] == len(bulk_index)
    assert prov["n_selected"] == len(genes)
    assert len(genes) < prov["n_shared"], \
        "the subset switch is on but the full intersection came back"
    # top-N per type plus markers, so an order of magnitude smaller, not a trim
    assert len(genes) <= cfg.SIGNATURE_GENES_PER_TYPE * len(cfg.CELL_TYPES) + 64
    assert len(genes) >= cfg.MIN_GENES_SHARED


def test_frozen_gene_space_returns_everything_when_the_switch_is_off(monkeypatch):
    """The negative control: with the switch off the full intersection must come back,
    or the test above would pass for a reason that has nothing to do with the flag."""
    from ivygap import config as cfg
    from ivygap.data.reference import frozen_gene_space, load_frozen_reference

    monkeypatch.setattr(cfg, "USE_SIGNATURE_GENE_SUBSET", False)
    frozen = load_frozen_reference()
    genes, prov = frozen_gene_space(frozen, frozen.profile.index)
    assert len(genes) == len(frozen.profile.index)
    assert prov["strategy"] == "full intersection"


def test_frozen_gene_space_ranks_within_the_shared_space():
    """
    Ranking over the whole frozen profile and intersecting afterwards silently drops
    top-scoring genes the bulk does not carry, so the subset comes out smaller than
    requested for a reason unrelated to how informative the genes are. Rank first over
    the shared space and every selected gene survives.
    """
    from ivygap.data.reference import frozen_gene_space, load_frozen_reference

    frozen = load_frozen_reference()
    half = list(frozen.profile.index[::2])               # pretend the bulk carries half
    genes, prov = frozen_gene_space(frozen, half)

    assert prov["n_shared"] == len(half)
    assert set(genes) <= set(half), "a gene absent from the bulk was selected"
    # nothing was selected and then thrown away
    assert prov["n_selected"] == len(genes)


def test_every_method_sees_the_same_gene_space():
    """Equal footing: the gene set is chosen once, for everyone, and recorded."""
    from ivygap.data.reference import frozen_gene_space, load_frozen_reference

    frozen = load_frozen_reference()
    a, _ = frozen_gene_space(frozen, frozen.profile.index)
    b, _ = frozen_gene_space(frozen, frozen.profile.index)
    assert a == b, "gene selection is not deterministic"
