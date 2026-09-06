"""
Regression tests for the partial-coverage failure found in run 2026-09-05T2154.

WHAT WENT WRONG THERE
---------------------
quanTIseq models four of the eight roster types and returns NaN for the other four,
which is correct and documented. Those four NaNs then reached `project_to_simplex`,
whose sum went non-finite, so it returned an all-NaN row — destroying the four columns
quanTIseq HAD estimated. Every one of the 122 samples came back entirely NaN, and the
run still recorded `implementation: R:quantiseqr` with `fallback_reason: null`.

A total loss reported as a success is the exact failure mode the project's evidence
rules exist to catch, so it gets tests rather than a comment.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap.deconv.base import (DeconvolutionInput, DeconvolutionMethod,
                                project_to_simplex)


# --- project_to_simplex, the function that did the damage ---------------------

def test_unmodelled_types_do_not_annihilate_the_modelled_ones():
    """The bug, stated directly: 4 structural NaNs must not erase 4 real estimates."""
    w = np.array([np.nan, 0.2, 0.05, np.nan])
    covered = np.array([False, True, True, False])

    out = project_to_simplex(w, covered)

    assert np.isfinite(out[[1, 2]]).all(), "covered entries were destroyed by NaN neighbours"
    assert np.isnan(out[[0, 3]]).all(), "uncovered entries must stay NaN, never zero"


def test_covered_values_keep_their_own_scale():
    """
    Not renormalised to sum 1. quanTIseq's fractions already include an unreported
    "Other" compartment; rescaling four immune types to sum 1 would claim the tumour
    is 100% immune — a composition the method never estimated.
    """
    w = np.array([np.nan, 0.2, 0.05, np.nan])
    covered = np.array([False, True, True, False])

    out = project_to_simplex(w, covered)

    assert out[1] == pytest.approx(0.2)
    assert out[2] == pytest.approx(0.05)
    assert np.nansum(out) == pytest.approx(0.25), "covered mass was rescaled"


def test_a_failed_solve_among_covered_types_still_condemns_the_row():
    """Inside the covered set, NaN keeps its original meaning: that solve failed."""
    w = np.array([np.nan, np.nan, 0.05, np.nan])
    covered = np.array([False, True, True, False])
    assert np.isnan(project_to_simplex(w, covered)).all()


def test_full_coverage_behaviour_is_unchanged():
    """Every other method must be bit-identical to before the change."""
    w = np.array([1.0, 2.0, 1.0])
    assert project_to_simplex(w) == pytest.approx([0.25, 0.5, 0.25])
    assert np.isnan(project_to_simplex(np.array([0.0, 0.0, 0.0]))).all()
    assert np.isnan(project_to_simplex(np.array([np.nan, 1.0, 1.0]))).all()


# --- the declaration reaching the base class ----------------------------------

def _toy_input(n_samples=4, n_genes=6):
    from ivygap.data.reference import ReferenceBundle
    types = ["Tumor", "Macrophage_Microglia", "T_cell", "Endothelial"]
    genes = [f"G{i}" for i in range(n_genes)]
    rng = np.random.default_rng(0)
    profile = pd.DataFrame(rng.random((n_genes, len(types))) + 0.1,
                           index=genes, columns=types)
    samples = [f"S{i}" for i in range(n_samples)]
    bulk = pd.DataFrame(rng.random((n_genes, n_samples)) + 0.1,
                        index=genes, columns=samples)
    ref = ReferenceBundle(name="toy", profile=profile,
                          sigma=pd.DataFrame(0.1, index=genes, columns=types),
                          cell_size=pd.Series(1.0, index=types), n_donors=3)
    manifest = pd.DataFrame({"patient_id": ["P0", "P0", "P1", "P1"][:n_samples],
                             "structure": ["CT", "PAN", "CT", "PAN"][:n_samples]},
                            index=samples)
    return DeconvolutionInput(bulk=bulk, references=(ref,), manifest=manifest,
                              cell_types=tuple(types))


class _PartialMethod(DeconvolutionMethod):
    name = "partial_toy"
    family = "test"
    models_cell_types = frozenset({"Macrophage_Microglia", "T_cell"})

    #: What the solver actually fills, held separately from the DECLARATION so a test
    #: can vary one without the other — that difference is the whole safety property.
    FILLS = ("Macrophage_Microglia", "T_cell")

    def _solve_all(self, data):
        n_s, n_t = len(data.samples), len(data.cell_types)
        out = np.full((n_s, n_t), np.nan)
        idx = [list(data.cell_types).index(c) for c in self.FILLS]
        out[:, idx] = 0.1
        return out


def test_a_declaring_method_survives_fit_predict():
    """End to end: the declaration must reach the projection, not just exist."""
    data = _toy_input()
    est = _PartialMethod().fit_predict(data)

    assert est[["Macrophage_Microglia", "T_cell"]].notna().all().all()
    assert est[["Tumor", "Endothelial"]].isna().all().all()
    assert not est.isna().all(axis=1).any(), "some row came back entirely NaN"


def test_undeclared_partial_output_is_still_condemned():
    """
    The safety property. A method that emits structural NaN WITHOUT declaring coverage
    must keep the old all-NaN behaviour — otherwise this change would start silently
    rescuing genuinely failed solves.
    """
    class _Undeclared(_PartialMethod):
        name = "undeclared_toy"
        models_cell_types = None

    est = _Undeclared().fit_predict(_toy_input())
    assert est.isna().all().all()


def test_a_typo_in_the_declaration_fails_loudly():
    class _Typo(_PartialMethod):
        name = "typo_toy"
        models_cell_types = frozenset({"Macrophage_Microglia", "T_cel"})

    with pytest.raises(ValueError, match="not in the roster"):
        _Typo().fit_predict(_toy_input())


def test_quantiseq_declares_the_four_types_its_signature_covers():
    from ivygap.deconv.reference_based import QuanTIseqDeconvolution
    m = QuanTIseqDeconvolution()
    assert m.models_cell_types == frozenset(
        {"T_cell", "NK_cell", "B_cell", "Macrophage_Microglia"})


def test_the_r_wrapper_carries_coverage_over_from_its_fallback():
    """
    Coverage is a property of the method, not of which implementation ran: quanTIseq
    has no tumour population whether quantiseqr or the reimplementation solves it.
    """
    from ivygap.deconv.r_bridge import RMethod
    from ivygap.deconv.reference_based import QuanTIseqDeconvolution
    wrapped = RMethod("quantiseq", QuanTIseqDeconvolution())
    assert wrapped.models_cell_types == QuanTIseqDeconvolution.models_cell_types


# --- the gene space a self-signature method receives --------------------------

def test_bulk_full_must_be_a_wider_gene_space_not_a_different_sample_set():
    data = _toy_input()
    with pytest.raises(ValueError, match="same samples"):
        DeconvolutionInput(bulk=data.bulk, references=data.references,
                           manifest=data.manifest, cell_types=data.cell_types,
                           bulk_full=data.bulk.iloc[:, :2])


def test_bulk_full_must_contain_every_gene_in_bulk():
    data = _toy_input()
    with pytest.raises(ValueError, match="superset"):
        DeconvolutionInput(bulk=data.bulk, references=data.references,
                           manifest=data.manifest, cell_types=data.cell_types,
                           bulk_full=data.bulk.iloc[2:])


def test_quantiseq_is_the_only_method_that_ignores_the_supplied_signature():
    """
    EPIC reads args$signature — this project's reference — so the shared marker subset
    is right for it. If EPIC ever joins OWN_SIGNATURE_METHODS, equal footing changed
    and the certificate needs to say so.
    """
    from ivygap.deconv.r_bridge import OWN_SIGNATURE_METHODS, SIGNATURE_ONLY_METHODS
    assert OWN_SIGNATURE_METHODS == {"quantiseq"}
    assert OWN_SIGNATURE_METHODS < SIGNATURE_ONLY_METHODS


def test_the_gene_space_decision_is_declared():
    from ivygap.deconv.configs import DECLARED_DEVIATIONS
    params = [d["parameter"] for d in DECLARED_DEVIATIONS["quantiseq"]]
    assert any("gene space" in p for p in params)
    assert any("cell-size" in p or "mRNA" in p for p in params)


# --- the same conflation, one layer up: benchmark metrics ---------------------
#
# `n_failed = any(~isfinite(row))` counted a structurally-unmodelled type as a failed
# SAMPLE, so run 2026-09-05T2154 reported quanTIseq as 200 of 200 benchmark samples
# failed and mae_primary NaN, when it had solved every one of them.

def _truth_and_estimate():
    types = ["Tumor", "Macrophage_Microglia", "T_cell", "Endothelial"]
    idx = [f"S{i}" for i in range(8)]
    rng = np.random.default_rng(1)
    truth = pd.DataFrame(rng.dirichlet(np.ones(4), size=8), index=idx, columns=types)
    est = truth + rng.normal(0, 0.01, truth.shape)
    return truth, est.clip(lower=0)


def test_unmodelled_types_are_not_counted_as_failed_samples():
    from ivygap.bench import metrics
    truth, est = _truth_and_estimate()
    partial = est.copy()
    partial[["Tumor", "Endothelial"]] = np.nan

    naive = metrics.overall_metrics(partial, truth)
    aware = metrics.overall_metrics(partial, truth,
                                    covered=["Macrophage_Microglia", "T_cell"])

    assert naive["n_failed_samples"] == len(truth), "precondition: the old count"
    assert aware["n_failed_samples"] == 0, "a solved sample was still counted as failed"
    assert np.isfinite(aware["mae_primary"])
    # The MAE itself must be computed over the modelled types only, not over a row
    # where pandas quietly skipped the NaNs and averaged a different set per sample.
    assert aware["n_types_scored"] == 2


def test_the_archived_all_nan_failure_is_reproduced_and_fixed():
    """
    The exact shape of the 2026-09-05T2154 artefact: every column NaN in every sample,
    which is what `project_to_simplex` produced before the fix. Both the sample count
    and mae_primary go bad, and both come back.
    """
    from ivygap.bench import metrics
    truth, est = _truth_and_estimate()
    annihilated = est.copy()
    annihilated[:] = np.nan

    naive = metrics.overall_metrics(annihilated, truth)
    assert naive["n_failed_samples"] == len(truth) == 8
    assert np.isnan(naive["mae_primary"]), "the archived NaN mae_primary"

    # With the two upstream fixes the method never emits that table in the first place:
    # only its unmodelled columns are NaN, and those are excluded from the metrics.
    partial = est.copy()
    partial[["Tumor", "Endothelial"]] = np.nan
    aware = metrics.overall_metrics(partial, truth,
                                    covered=["Macrophage_Microglia", "T_cell"])
    assert aware["n_failed_samples"] == 0
    assert np.isfinite(aware["mae_primary"])


def test_partial_metrics_are_flagged_and_carry_their_type_set():
    from ivygap.bench import metrics
    truth, est = _truth_and_estimate()
    partial = est.copy()
    partial[["Tumor", "Endothelial"]] = np.nan

    m = metrics.overall_metrics(partial, truth,
                               covered=["Macrophage_Microglia", "T_cell"])
    assert m["comparable"] is False
    assert m["n_types_scored"] == 2
    assert set(m["types_scored"].split(",")) == {"Macrophage_Microglia", "T_cell"}
    assert metrics.overall_metrics(est, truth)["comparable"] is True


def test_a_partial_method_can_never_be_the_selected_method():
    """
    The selection rule, tested directly. A method scored over two types must not be
    selectable on an MAE that is not the same quantity as everyone else's.
    """
    from ivygap.bench import metrics
    truth, est = _truth_and_estimate()
    partial = est.copy()
    partial[["Tumor", "Endothelial"]] = np.nan
    # Make the partial method look best on its own reduced set.
    partial[["Macrophage_Microglia", "T_cell"]] = truth[["Macrophage_Microglia", "T_cell"]]

    summary = metrics.summarise_methods(
        {"full_method": est, "partial_method": partial}, truth,
        covered={"partial_method": ["Macrophage_Microglia", "T_cell"]})

    assert summary.loc["partial_method", "mae_primary"] < summary.loc["full_method", "mae_primary"]
    assert str(summary.index[0]) == "full_method", "a partial method reached the top row"
    selectable = summary[summary["comparable"]]
    assert "partial_method" not in set(selectable.index)


def test_covered_that_excludes_everything_fails_loudly():
    from ivygap.bench import metrics
    truth, est = _truth_and_estimate()
    with pytest.raises(ValueError, match="excludes every column"):
        metrics.overall_metrics(est, truth, covered=["NotAType"])


def test_covered_with_no_primary_type_fails_loudly():
    """mae_primary over an empty set would be NaN and read as a failure, not a refusal."""
    from ivygap.bench import metrics
    from ivygap import config
    truth, est = _truth_and_estimate()
    non_primary = [c for c in truth.columns if c not in config.PRIMARY_CELL_TYPES]
    if not non_primary:
        pytest.skip("every type in this fixture is a primary type")
    with pytest.raises(ValueError, match="primary"):
        metrics.overall_metrics(est, truth, covered=non_primary)


# --- the R bridge must not report an empty table as a successful run ----------

def test_the_two_coverage_declarations_cannot_drift():
    """
    `r_bridge.OWN_COVERAGE` exists so the all-NaN guard can tell "does not model tumour
    cells" from "returned nothing". If it disagrees with the method's own
    `models_cell_types`, the guard either fires on a healthy method or stays silent on a
    dead one — so the two are pinned to each other here.
    """
    from ivygap.deconv.r_bridge import OWN_COVERAGE
    from ivygap.deconv.registry import build_methods

    for method in build_methods(prefer_r=False):
        declared = getattr(method, "models_cell_types", None)
        r_name = getattr(method, "r_method", method.name)
        if declared:
            assert r_name in OWN_COVERAGE, (
                f"{method.name} declares partial coverage but r_bridge.OWN_COVERAGE "
                f"has no entry for {r_name!r}; the all-NaN guard would reject it")
            assert OWN_COVERAGE[r_name] == declared, (
                f"{r_name}: OWN_COVERAGE {sorted(OWN_COVERAGE[r_name])} != "
                f"models_cell_types {sorted(declared)}")
        else:
            assert r_name not in OWN_COVERAGE, (
                f"r_bridge.OWN_COVERAGE claims {r_name} models only "
                f"{sorted(OWN_COVERAGE.get(r_name, []))}, but the method declares none")


def test_an_all_nan_r_result_is_a_failure_not_a_success():
    """
    The exact artefact from run 2026-09-05T2154: exit 0, a well-formed 122x8 CSV with
    the right index and column names, every cell NaN — recorded as a successful
    `R:quantiseqr` run. It must raise instead, so the fallback machinery engages.
    """
    import ivygap.deconv.r_bridge as rb

    data = _toy_input()
    dead = pd.DataFrame(np.nan, index=data.samples, columns=list(data.cell_types))

    with pytest.raises(rb.RBridgeError, match="no estimate"):
        rb._reject_unusable("music", dead, data, "run_music.R")


def test_a_partial_methods_unmodelled_columns_do_not_trip_the_guard():
    import ivygap.deconv.r_bridge as rb

    data = _toy_input()
    partial = pd.DataFrame(np.nan, index=data.samples, columns=list(data.cell_types))
    partial["Macrophage_Microglia"] = 0.2
    partial["T_cell"] = 0.05

    rb._reject_unusable("quantiseq", partial, data, "run_quantiseq.R")   # must not raise


def test_a_partial_method_that_returns_nothing_still_trips_the_guard():
    import ivygap.deconv.r_bridge as rb

    data = _toy_input()
    dead = pd.DataFrame(np.nan, index=data.samples, columns=list(data.cell_types))
    with pytest.raises(rb.RBridgeError, match="no estimate"):
        rb._reject_unusable("quantiseq", dead, data, "run_quantiseq.R")
