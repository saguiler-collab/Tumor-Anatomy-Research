"""D23: the lymphoid result's two load-bearing measurements, pinned to the shipped artefacts."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

RESULTS = Path(__file__).resolve().parent.parent / "results"


def _load(name):
    p = RESULTS / name
    if not p.exists():
        pytest.skip(f"{name} not built")
    return json.loads(p.read_text())


def test_bisque_returns_the_reference_prior_on_planted_truth():
    d = _load("bisque_anchoring_synthetic.json")
    far, matched = d["far"], d["matched"]
    # the test: a planted mean far from the prior comes back as the prior
    assert far["L1_estimate_to_prior"] < 0.01
    assert far["L1_estimate_to_plant"] > 0.5
    # the control: when plant = prior, Bisque is accurate -- the mechanism, not a broken run
    assert matched["L1_estimate_to_plant"] < 0.15
    # and the ranking survives the level replacement
    assert min(far["per_sample_spearman"].values()) > 0.7


def test_bisque_cohort_means_sit_on_the_reference_prior_in_both_cohorts():
    d = _load("bisque_anchoring.json")
    for c in ("gbm", "lgg"):
        ms = d["cohorts"][c]["methods"]
        b = ms["bisque"]["L1_to_reference_prior_all_types"]
        others = [v["L1_to_reference_prior_all_types"] for k, v in ms.items() if k != "bisque"]
        assert b < 0.15, f"{c}: bisque is {b} from the prior"
        assert b < min(others) / 4, f"{c}: bisque not clearly closer to the prior than others"


def test_the_marker_panels_were_not_revised_after_the_result():
    # Fixed before the first run (2026-09-30). Changing them until the control passes would be
    # selection on the answer; this test fails loudly if anyone tries.
    import sys
    sys.path.insert(0, str(RESULTS.parent / "scripts"))
    from lymphoid_tracking import B_MARKERS, T_MARKERS
    assert T_MARKERS == ["CD3D", "CD3E", "CD3G", "CD5", "CD6", "TRAT1"]
    assert B_MARKERS == ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "BLK", "PAX5", "FCRLA"]


def test_the_failed_positive_control_is_recorded_as_failed():
    d = _load("lymphoid_tracking.json")
    assert d["positive_control_verdict"].startswith("FAILED")
    lgg = d["positive_control_marker_index"]["lgg"]
    assert lgg["rho_vs_truth_TB"] < 0 and lgg["p"] < 0.05


def test_nk_reannotation_baseline_reproduces_the_registered_fit():
    # The variants are only interpretable if the unchanged roster reproduces the registered arm.
    d = _load("nk_reannotation.json")
    for cohort, reg_file in (("gbm", "lymphoid_ordering.json"), ("lgg", "lymphoid_ordering_lgg.json")):
        reg = _load(reg_file)["methods"]
        base = d["cohorts"][cohort]["variants"]["as_registered"]
        for m in ("nnls", "svr"):
            assert base[m]["n_with_lymphoid_signal"] == reg[m]["n_scored"]
            assert abs(base[m]["frac_B_over_T"] - reg[m]["frac_samples_est_B_over_T"]) < 1e-3


def test_removing_or_merging_nk_does_not_restore_t_over_b():
    d = _load("nk_reannotation.json")
    for cohort in ("gbm", "lgg"):
        for kind in ("nk_removed", "t_nk_merged"):
            assert d["cohorts"][cohort]["variants"][kind]["svr"]["frac_B_over_T"] > 0.9


def test_b_is_the_most_tissue_like_lymphoid_profile_in_both_cohorts():
    d = _load("b_profile_tissue_likeness.json")
    for c in ("gbm", "lgg"):
        r = d["cohorts"][c]["r_with_mean_bulk"]
        for key in ("all_genes", "without_rp"):        # holds without ribosomal genes too
            assert r["B_cell"][key] > max(r["T_cell"][key], r["NK_cell"][key]), (c, key)


def test_atlas_b_cells_are_ribosome_rich_but_not_myelin_contaminated():
    d = _load("atlas_cell_fractions.json")["by_type"]
    assert d["B_cell"]["median_rp_fraction"] > 1.4 * d["T_cell"]["median_rp_fraction"]
    # the ambient-myelin explanation was rejected: B carries LESS myelin signal than T
    assert d["B_cell"]["mean_myelin_fraction"] < d["T_cell"]["mean_myelin_fraction"]
    # and B cells are not low-quality: comparable UMIs to T
    assert d["B_cell"]["median_umis"] > 0.8 * d["T_cell"]["median_umis"]


def test_ribosomal_gene_removal_does_not_restore_t_over_b():
    d = _load("ribosomal_test.json")["cohorts"]
    for c in d:
        for m in ("nnls", "svr"):
            base, rm = d[c].get(f"{m}_baseline"), d[c].get(f"{m}_rp_removed")
            if not base or not rm or base["frac_B_over_T"] is None:
                continue
            assert rm["frac_B_over_T"] >= base["frac_B_over_T"] - 0.05, (c, m)
