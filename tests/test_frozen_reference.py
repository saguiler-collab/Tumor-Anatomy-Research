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


# =============================================================================
# reading a real single-cell atlas without densifying it
# -----------------------------------------------------------------------------
# `build_from_h5ad` used to call X.toarray() unconditionally. Core GBmap is 338,564
# cells x ~30,000 genes; dense float32 is about 40 TB, against 9 GB of RAM on the
# machine this was written for. That is not a slow path, it is an impossible one, and
# it made the atlas — which unblocks all four published R tools, method selection and
# two of the three yardsticks — unusable.
# =============================================================================

def _tiny_h5ad(tmp_path, n_cells=800, n_genes=120, n_donors=10):
    import anndata
    from scipy import sparse

    rng = np.random.default_rng(0)
    types = ["AC-like", "MES-like", "TAM-BDM", "TAM-MG", "CD4/CD8", "NK", "B cell",
             "Endothelial", "Oligodendrocyte", "Astrocyte", "NOT-IN-THE-ROSTER"]
    obs = pd.DataFrame(
        {"annotation_level_3": rng.choice(types, n_cells),
         "donor_id": rng.choice([f"D{i}" for i in range(n_donors)], n_cells)},
        index=[f"c{i}" for i in range(n_cells)])
    X = sparse.random(n_cells, n_genes, density=0.15, random_state=0, format="csr") * 100
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])
    path = tmp_path / "atlas.h5ad"
    anndata.AnnData(X=X, obs=obs, var=var).write_h5ad(path)
    return path


def _isolated_ref_dir(tmp_path, monkeypatch):
    from ivygap import config as cfg
    monkeypatch.setattr(cfg, "REFERENCE_DIR", tmp_path)
    monkeypatch.setattr(cfg, "ALL_OUTPUT_DIRS", [tmp_path])


def test_h5ad_is_read_without_densifying_the_whole_matrix(tmp_path, monkeypatch):
    import json as _json
    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)
    path = _tiny_h5ad(tmp_path)

    ref, expr, meta = build_from_h5ad(
        path, name="t", max_cells_per_donor_type=5, max_total_cells=200,
        restrict_to_genes=[f"G{i}" for i in range(80)], export=False)

    samp = _json.loads((tmp_path / "reference_sampling_t.json").read_text())
    assert samp["n_cells_in_file"] == 800
    assert samp["n_cells_after_roster_mapping"] < 800, "unmapped types were not dropped"
    assert samp["n_cells_kept"] <= samp["n_cells_after_roster_mapping"]
    assert samp["n_genes_kept"] == 80 < samp["n_genes_in_file"]
    assert expr.shape == (80, samp["n_cells_kept"])
    # a label outside the roster must never reach the reference
    assert "NOT-IN-THE-ROSTER" not in set(meta["cell_type"])
    assert set(meta["cell_type"]) <= set(config.CELL_TYPES)


def test_subsample_is_donor_balanced_and_seeded(tmp_path, monkeypatch):
    import json as _json
    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)
    path = _tiny_h5ad(tmp_path)

    cap = 4
    _, _, meta = build_from_h5ad(path, name="t", max_cells_per_donor_type=cap,
                                 max_total_cells=10_000, export=False)
    counts = meta.groupby(["donor", "cell_type"], observed=True).size()
    assert counts.max() <= cap, "the per-(donor, cell type) cap was not applied"

    samp = _json.loads((tmp_path / "reference_sampling_t.json").read_text())
    assert samp["strategy"].startswith("donor-balanced")
    assert samp["n_donors_kept"] > 1

    # deterministic: the same seed must select the same cells
    _, _, meta2 = build_from_h5ad(path, name="t2", max_cells_per_donor_type=cap,
                                  max_total_cells=10_000, export=False)
    assert list(meta.index) == list(meta2.index)


def test_atlas_reference_carries_cross_donor_variance(tmp_path, monkeypatch):
    """
    The point of the atlas. A collapsed signature has no cross-donor variance, which is
    why MuSiC runs degenerate and Bisque falls back to spread across cell types. A
    cell-level reference must restore both, or nothing downstream improves.
    """
    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)
    path = _tiny_h5ad(tmp_path)

    ref, _, _ = build_from_h5ad(path, name="t", max_cells_per_donor_type=6,
                                max_total_cells=10_000, export=False)
    assert ref.has_cross_donor_variance, "MuSiC would still be degenerate"
    assert ref.donor_profiles and len(ref.donor_profiles) >= 3, \
        "Bisque would still fall back to spread across cell types"


def test_restrict_to_genes_rejects_a_disjoint_gene_space(tmp_path, monkeypatch):
    """
    The negative control. Handing in Ensembl ids when the atlas uses symbols would
    otherwise select zero genes and fail much later with something unrecognisable.
    """
    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)
    path = _tiny_h5ad(tmp_path)

    with pytest.raises(ValueError, match="shares no gene"):
        build_from_h5ad(path, name="t", restrict_to_genes=["ENSG00000000001"],
                        export=False)



def test_atlas_reference_puts_marker_genes_in_the_right_column(tmp_path, monkeypatch):
    """
    The invariant that catches a scrambled reference.

    `sparse_dataset[...]` returns rows in ascending file order, but the balanced sampler
    returns cells grouped by (donor, cell type). Labelling the columns with the sampler's
    order attached every cell's expression to a DIFFERENT cell's metadata — which
    scrambles every per-type mean while leaving the column ORDER, every shape check and
    every existing test intact.

    It surfaced only as biology: 23 of 24 markers in the wrong column, Oligodendrocyte
    peaking in microvascular proliferation, and a leaderboard where every real method
    scored below a random control. A marker that does not top its own cell type is the
    cheapest possible detector, so it is asserted directly.
    """
    import anndata
    from scipy import sparse

    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)

    # one unmistakable marker per roster type, planted at high expression
    markers = {"AC-like": "EGFR", "TAM-BDM": "CD68", "CD4/CD8": "CD3D",
               "NK": "NKG7", "B cell": "CD79A", "Endothelial": "PECAM1",
               "Oligodendrocyte": "PLP1", "Astrocyte": "AQP4"}
    genes = list(markers.values()) + [f"BG{i}" for i in range(40)]
    labels = list(markers)

    rng = np.random.default_rng(0)
    n_per = 40
    obs_rows, blocks = [], []
    for li, lab in enumerate(labels):
        for k in range(n_per):
            obs_rows.append({"annotation_level_3": lab,
                             "donor_id": f"D{k % 5}"})
        block = rng.random((n_per, len(genes))) * 2.0
        block[:, genes.index(markers[lab])] += 500.0      # the marker dominates
        blocks.append(block)

    X = sparse.csr_matrix(np.vstack(blocks))
    obs = pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(len(obs_rows))])
    var = pd.DataFrame(index=genes)
    path = tmp_path / "planted.h5ad"
    anndata.AnnData(X=X, obs=obs, var=var).write_h5ad(path)

    ref, _, meta = build_from_h5ad(path, name="planted", max_cells_per_donor_type=6,
                                   max_total_cells=10_000, export=False)

    wrong = []
    for label, gene in markers.items():
        expected = config.GBMAP_CELL_TYPE_MAP[label]
        if gene in ref.profile.index:
            top = ref.profile.loc[gene].idxmax()
            if top != expected:
                wrong.append(f"{gene}: expected {expected}, highest in {top}")
    assert not wrong, (
        "the reference is scrambled — cell expression is attached to the wrong "
        "metadata:\n  " + "\n  ".join(wrong))


def test_atlas_reference_cells_keep_their_own_metadata(tmp_path, monkeypatch):
    """
    The same bug stated directly: a cell's expression and its label must travel
    together. Planting one type at a distinctive level makes the pairing checkable
    without going through the signature.
    """
    import anndata
    from scipy import sparse

    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)

    genes = [f"G{i}" for i in range(30)]
    rows, mat = [], []
    rng = np.random.default_rng(1)
    for i in range(160):
        # Endothelial cells are 100x brighter than everything else
        lab = "Endothelial" if i % 8 == 0 else "AC-like"
        vec = rng.random(len(genes)) * (100.0 if lab == "Endothelial" else 1.0)
        rows.append({"annotation_level_3": lab, "donor_id": f"D{i % 4}"})
        mat.append(vec)

    obs = pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))])
    path = tmp_path / "paired.h5ad"
    anndata.AnnData(X=sparse.csr_matrix(np.vstack(mat)), obs=obs,
                    var=pd.DataFrame(index=genes)).write_h5ad(path)

    _, expr, meta = build_from_h5ad(path, name="paired", max_cells_per_donor_type=50,
                                    max_total_cells=10_000, export=False)

    endo = meta.index[meta["cell_type"] == "Endothelial"]
    other = meta.index[meta["cell_type"] == "Tumor"]
    assert len(endo) and len(other)
    # after per-cell CPM every column sums to 1e6, so compare the SHAPE that survives:
    # the planted type's cells must still be the ones labelled Endothelial
    assert expr[endo].mean(axis=1).corr(expr[other].mean(axis=1)) < 0.999, \
        "the two planted populations are indistinguishable; labels are not tracking cells"


def test_cell_size_is_measured_before_normalisation(tmp_path, monkeypatch):
    """
    Normalising every cell to a common library size makes each type's mean total
    identical by construction, so the cell-size factors come out all equal and
    `to_cell_fractions` silently becomes the identity — RNA proportions get reported as
    cell proportions. Observed on the real atlas: every factor exactly 1,000,000.

    A tumour cell and a lymphocyte do not carry comparable amounts of mRNA, and
    config.py names this correction as a requirement, so the factors must differ.
    """
    import anndata
    from scipy import sparse

    from ivygap.data.reference import build_from_h5ad

    pytest.importorskip("anndata")
    _isolated_ref_dir(tmp_path, monkeypatch)

    genes = [f"G{i}" for i in range(24)]
    rng = np.random.default_rng(0)
    rows, mat = [], []
    for i in range(120):
        # Endothelial cells carry 20x the library size of the rest
        lab = "Endothelial" if i % 4 == 0 else "AC-like"
        scale = 20.0 if lab == "Endothelial" else 1.0
        rows.append({"annotation_level_3": lab, "donor_id": f"D{i % 4}"})
        mat.append(rng.random(len(genes)) * scale)

    path = tmp_path / "sizes.h5ad"
    anndata.AnnData(X=sparse.csr_matrix(np.vstack(mat)),
                    obs=pd.DataFrame(rows, index=[f"c{i}" for i in range(len(rows))]),
                    var=pd.DataFrame(index=genes)).write_h5ad(path)

    ref, _, _ = build_from_h5ad(path, name="sizes", max_cells_per_donor_type=50,
                                max_total_cells=10_000, export=False)

    sizes = ref.cell_size.dropna()
    assert sizes.nunique() > 1, (
        f"every cell-size factor is identical ({sizes.iloc[0]:.0f}); the correction is "
        f"a no-op and RNA proportions are being reported as cell proportions")
    assert sizes["Endothelial"] > sizes["Tumor"] * 5, (
        "the planted 20x library-size difference did not survive into the factors")
