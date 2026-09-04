"""
Tests for the Anatomic Concordance Score.

Step 4 of the protocol: "Build the scorer, prove it on synthetic anatomy. Composition
that satisfies every constraint must score near 1; permuted labels near chance."

That gate is what these tests implement. The negative cases matter more than the
positive one — a scorer that returns high ACS for everything would pass a
planted-signal test and be worthless.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.anatomic import acs, constraints as K


def _cohort(n_tumors=12, seed=0, satisfy=True, blocks=None, drop=None, flat=False):
    """
    Build composition that satisfies (or violates) every constraint by construction.

    `blocks` gives one tumour extra samples per structure, exercising the uneven-block
    hazard. `drop` removes a (tumour, structure) pair so unevaluable constraints are
    exercised. `flat` makes every structure identical, which must score 0 — ties are not
    half credit.
    """
    rng = np.random.default_rng(seed)
    sign = 1.0 if satisfy else -1.0
    rows, meta, n = {}, {}, 0

    for i in range(n_tumors):
        tid = f"T{i:02d}"
        for s in config.PRIMARY_STRUCTURES:
            if drop and (tid, s) in drop:
                continue
            reps = (blocks or {}).get(tid, 1)
            for _ in range(reps):
                c = dict.fromkeys(config.CELL_TYPES, 0.04)
                c["Tumor"], c["Macrophage_Microglia"] = 0.50, 0.18
                c["Oligodendrocyte"], c["Endothelial"] = 0.08, 0.06
                if not flat:
                    if s == "LE":
                        c["Tumor"] -= sign * 0.35; c["Oligodendrocyte"] += sign * 0.30
                        c["Macrophage_Microglia"] -= sign * 0.08
                    elif s == "IT":
                        c["Tumor"] -= sign * 0.18; c["Oligodendrocyte"] += sign * 0.12
                    elif s == "MVP":
                        c["Endothelial"] += sign * 0.22
                        c["Macrophage_Microglia"] += sign * 0.06
                    elif s == "PAN":
                        c["Macrophage_Microglia"] += sign * 0.14
                v = np.array([c[t] for t in config.CELL_TYPES])
                if not flat:
                    v = np.clip(v + rng.normal(0, 0.004, len(v)), 1e-4, None)
                sid = f"S{n:04d}"
                rows[sid] = v / v.sum()
                meta[sid] = {"patient_id": tid, "structure": s}
                n += 1

    return (pd.DataFrame(rows, index=config.CELL_TYPES).T,
            pd.DataFrame.from_dict(meta, orient="index"))


# --- the constraint file ----------------------------------------------------

def test_constraint_file_is_structurally_valid():
    K.validate()
    assert len(K.CONSTRAINTS) == 7
    assert {c.id for c in K.CONSTRAINTS} == {"C1", "C2", "C3", "C4", "C5", "C6", "C7"}


def test_t_cell_is_excluded_and_no_constraint_scores_it():
    assert "T_cell" in {e["cell_type"] for e in K.EXCLUSIONS}
    assert all(c.cell_type != "T_cell" for c in K.CONSTRAINTS)


def test_no_constraint_rests_on_the_unidentifiable_sidecar():
    for c in K.CONSTRAINTS:
        assert c.cell_type not in config.SIDECAR_CELL_TYPES


def test_freeze_hash_changes_when_a_constraint_changes(monkeypatch):
    before = K.freeze_hash()
    monkeypatch.setattr(K, "CONSTRAINTS", K.CONSTRAINTS[:-1])
    assert K.freeze_hash() != before


def test_coverage_report_names_what_cannot_be_tested():
    cov = K.coverage_report()
    assert "T_cell" in cov["explicitly_excluded"]
    assert cov["structurally_untestable"], "silence about neurons must be explicit"


# --- the scorer: step 4's gate ----------------------------------------------

def test_satisfying_composition_scores_near_one():
    est, mf = _cohort(satisfy=True)
    r = acs.score(est, mf, n_permutations=200, n_boot=200)
    assert r.acs > 0.95
    assert r.n_tumors == 12


def test_violating_composition_scores_near_zero():
    """The negative control for the scorer itself."""
    est, mf = _cohort(satisfy=False)
    r = acs.score(est, mf, n_permutations=200, n_boot=200)
    assert r.acs < 0.10


def test_permutation_null_is_centred_near_chance():
    """
    Step 4's second gate. The null must not sit near 1 (constraints trivially satisfied)
    or near 0 (constraints impossible); it must sit where shuffled anatomy lands.
    """
    est, mf = _cohort(satisfy=True)
    null = acs.permutation_null(est, mf, n_permutations=400)
    assert np.isfinite(null).all()
    assert 0.05 < null.mean() < 0.60, f"null centred at {null.mean():.3f}"


def test_real_signal_beats_its_own_permutation_null():
    est, mf = _cohort(satisfy=True, n_tumors=14)
    r = acs.score(est, mf, n_permutations=500, n_boot=300)
    assert r.acs > r.null_mean
    assert r.null_p < 0.01


def test_flat_composition_gets_no_credit_for_ties():
    """A method that reports every structure identically has reproduced no ordering."""
    est, mf = _cohort(flat=True)
    assert acs.score(est, mf, n_permutations=100, n_boot=100).acs == 0.0


def test_uneven_block_counts_do_not_move_the_score():
    """
    The Ivy GAP hazard: one tumour dissected far more heavily than the others.

    The comparison duplicates T00's EXISTING rows rather than generating new ones, so
    the per-(tumour, structure) means are identical by construction and the only thing
    that changes is how many blocks T00 contributed. Regenerating with a different
    sample count would also change every noise draw, and the test would be measuring the
    RNG rather than the aggregation.
    """
    est, mf = _cohort(n_tumors=10, seed=3)
    balanced = acs.score(est, mf, n_permutations=100, n_boot=100).acs

    dup_rows = mf.index[mf["patient_id"] == "T00"]
    extra_est = pd.concat(
        [est] + [est.loc[dup_rows].rename(index=lambda s, k=k: f"{s}_dup{k}")
                 for k in range(11)])
    extra_mf = pd.concat(
        [mf] + [mf.loc[dup_rows].rename(index=lambda s, k=k: f"{s}_dup{k}")
                for k in range(11)])
    skewed = acs.score(extra_est, extra_mf, n_permutations=100, n_boot=100).acs

    assert skewed == pytest.approx(balanced, abs=1e-9), (
        "one heavily dissected tumour changed the score; the collapse to "
        "(tumour, structure) means is not neutralising uneven block counts"
    )


def test_missing_structure_is_excluded_not_counted_as_failure():
    """
    A tumour that never contributed MVP must not penalise a method on C3/C4/C6. It
    should reduce the number of evaluable pairs, not the score.
    """
    full_est, full_mf = _cohort(n_tumors=10, seed=5)
    drop = {("T00", "MVP"), ("T01", "MVP")}
    part_est, part_mf = _cohort(n_tumors=10, seed=5, drop=drop)

    full = acs.score(full_est, full_mf, n_permutations=100, n_boot=100)
    part = acs.score(part_est, part_mf, n_permutations=100, n_boot=100)

    assert part.n_pairs < full.n_pairs
    assert part.acs == pytest.approx(full.acs, abs=0.02)


def test_c4_requires_a_real_comparison_set():
    """'MVP is the maximum' among one other structure would just restate C3."""
    est, mf = _cohort(n_tumors=6)
    keep = mf.index[mf["structure"].isin(["MVP", "CT"])]
    ts = acs.tumor_structure_means(est.loc[keep], mf.loc[keep])
    assert acs.evaluate_constraint(ts, "T00", K.by_id("C4")) is None


def test_c7_monotone_needs_all_three_and_is_all_or_nothing():
    c7 = K.by_id("C7")
    est, mf = _cohort(n_tumors=6)
    ts = acs.tumor_structure_means(est, mf)
    assert acs.evaluate_constraint(ts, "T00", c7) == 1.0

    # Break only the middle step: the chain must fail as a unit.
    broken = ts.copy()
    broken.loc[("T00", "IT"), "Tumor"] = broken.loc[("T00", "CT"), "Tumor"] + 0.1
    assert acs.evaluate_constraint(broken, "T00", c7) == 0.0


def test_bootstrap_ci_brackets_the_estimate():
    est, mf = _cohort(n_tumors=14, seed=9)
    r = acs.score(est, mf, n_permutations=200, n_boot=600)
    assert r.ci_low <= r.acs <= r.ci_high
    assert r.ci_high - r.ci_low > 0


def test_fast_permutation_null_matches_the_reference():
    """
    The optimised null must be identical to the readable one, not merely similar.

    Both consume the same RNG stream in the same order, so agreement should be exact.
    An optimised statistic with no independent implementation to check against is where
    a subtle indexing error lives forever.
    """
    est, mf = _cohort(n_tumors=8, seed=11)
    fast = acs.permutation_null(est, mf, n_permutations=120, seed=5)
    slow = acs.permutation_null_reference(est, mf, n_permutations=120, seed=5)
    np.testing.assert_allclose(fast, slow, rtol=0, atol=1e-12)


def test_fast_null_handles_a_tumor_missing_a_structure():
    """The equivalence must hold on ragged data, which is what the real archive is."""
    est, mf = _cohort(n_tumors=8, seed=13,
                      drop={("T00", "MVP"), ("T01", "PAN"), ("T02", "LE")})
    fast = acs.permutation_null(est, mf, n_permutations=100, seed=7)
    slow = acs.permutation_null_reference(est, mf, n_permutations=100, seed=7)
    np.testing.assert_allclose(fast, slow, rtol=0, atol=1e-12)


# --- regression tests for two defects found by adversarial probing -----------

def test_fast_null_matches_reference_when_a_method_emits_nan():
    """
    Regression. `project_to_simplex` returns NaN by design for an all-zero solution, so
    a real method CAN fail on samples. The original fast path summed with np.add.at
    (which propagates NaN) while the reference used pandas .mean() (which skips it), and
    additionally counted a NaN sample as 'structure present'. Both divergences were
    silent and both moved the null — measured at 0.32 vs 0.39 on this fixture.
    """
    est, mf = _cohort(n_tumors=8, seed=21)
    est = est.copy()
    est.iloc[[0, 5, 11]] = np.nan                    # three wholly failed samples

    fast = acs.permutation_null(est, mf, n_permutations=80, seed=3)
    slow = acs.permutation_null_reference(est, mf, n_permutations=80, seed=3)
    np.testing.assert_allclose(fast, slow, rtol=0, atol=1e-12)


def test_fast_null_matches_reference_when_one_cell_type_is_nan():
    """
    A method can fail on one cell type while the others are fine. Treating the whole
    (tumour, structure) as missing would discard usable evidence, so the observation
    mask is per cell type — and the two paths must still agree.
    """
    est, mf = _cohort(n_tumors=8, seed=23)
    est = est.copy()
    est.iloc[[2, 7], 3] = np.nan

    fast = acs.permutation_null(est, mf, n_permutations=80, seed=3)
    slow = acs.permutation_null_reference(est, mf, n_permutations=80, seed=3)
    np.testing.assert_allclose(fast, slow, rtol=0, atol=1e-12)


def test_observed_statistic_and_null_use_the_same_estimator(monkeypatch):
    """
    Regression for the subtler defect: the observed ACS came from the pandas path while
    the null came from the array path. They agreed on clean data, so nothing failed —
    but every p-value would have been comparing two different statistics as soon as a
    method emitted NaN.

    Breaking one estimator must now raise rather than produce a plausible-looking
    p-value.
    """
    est, mf = _cohort(n_tumors=6, seed=31)

    real = acs._acs_from_arrays

    def sabotaged(means, has_obs, compiled):
        return real(means, has_obs, compiled) * 0.5     # a wrong but finite answer

    monkeypatch.setattr(acs, "_acs_from_arrays", sabotaged)
    with pytest.raises(AssertionError, match="estimators disagree"):
        acs.score(est, mf, n_permutations=20, n_boot=20)


def test_failed_samples_are_excluded_not_scored_as_violations():
    """
    A tumour whose every block failed for a structure cannot inform a constraint using
    that structure. It must drop out of the denominator, exactly like a tumour that
    never contributed the structure — not be counted as a violation.
    """
    est, mf = _cohort(n_tumors=8, seed=41)
    mvp_rows = mf.index[(mf["patient_id"] == "T00") & (mf["structure"] == "MVP")]

    # NaN out those rows, versus physically removing them from the SAME cohort.
    # Regenerating with `drop=` would change every noise draw and compare two different
    # datasets rather than two ways of expressing one absence.
    broken = est.copy()
    broken.loc[mvp_rows] = np.nan
    dropped_est = est.drop(index=mvp_rows)
    dropped_mf = mf.drop(index=mvp_rows)

    a = acs.score(broken, mf, n_permutations=40, n_boot=40)
    b = acs.score(dropped_est, dropped_mf, n_permutations=40, n_boot=40)
    assert a.acs == pytest.approx(b.acs, abs=1e-12)
    assert a.n_pairs == b.n_pairs


# =============================================================================
# the risks Anatomy_Test.md names, and whether anything actually checks them
# =============================================================================

def test_pre_registration_is_not_claimed_without_a_receipt():
    """
    Risk 1, circularity. The protocol's gate is "registration timestamp precedes every
    result file". A hash proves the constraints have not changed since some moment; it
    does not prove the moment came before the results. The honest default is
    UNREGISTERED, and the project must not describe itself otherwise.
    """
    from ivygap.anatomic import registration

    st = registration.status()
    assert st.state in {"REGISTERED", "UNREGISTERED", "HASH_MISMATCH", "MALFORMED"}
    if not registration.REGISTRATION_PATH.exists():
        assert st.state == "UNREGISTERED"
        assert not st.is_registered
        assert "must not describe itself as pre-registered" in st.verdict


def test_registration_detects_a_constraint_file_edited_after_registration(tmp_path,
                                                                         monkeypatch):
    """The fatal case, and the only one this can actually prove: the hash on record no
    longer matches the live constraint file."""
    import json

    from ivygap.anatomic import registration

    receipt = tmp_path / "REGISTRATION.json"
    receipt.write_text(json.dumps({
        "registry": "OSF", "url": "https://osf.io/xxxxx",
        "registered_utc": "2020-01-01T00:00:00Z",
        "constraint_freeze_hash": "0" * 64,          # not the live hash
    }))
    monkeypatch.setattr(registration, "REGISTRATION_PATH", receipt)

    st = registration.status(check_results=False)
    assert st.state == "HASH_MISMATCH"
    assert not st.is_registered
    assert "edited since registration" in st.verdict


def test_registration_accepts_a_matching_receipt(tmp_path, monkeypatch):
    """The positive control — without it the mismatch test could pass by always failing."""
    import json

    from ivygap.anatomic import constraints as K
    from ivygap.anatomic import registration

    receipt = tmp_path / "REGISTRATION.json"
    receipt.write_text(json.dumps({
        "registry": "OSF", "url": "https://osf.io/xxxxx",
        "registered_utc": "2020-01-01T00:00:00Z",
        "constraint_freeze_hash": K.freeze_hash(),
    }))
    monkeypatch.setattr(registration, "REGISTRATION_PATH", receipt)

    st = registration.status(check_results=False)
    assert st.state == "REGISTERED" and st.is_registered


def test_every_method_config_is_recorded_with_its_deviations():
    """
    Risk 5, method setup fairness. "One shared signature, published configs for every
    method, defaults from each tool's own documentation — no per-method tuning." The
    shared signature was enforced; the configs were recorded nowhere, so "we used the
    defaults" was an assertion a reader had to take on faith.
    """
    from ivygap.deconv import configs
    from ivygap.deconv.controls import build_controls
    from ivygap.deconv.registry import build_methods

    methods = build_methods(prefer_r=False) + build_controls()
    recs = configs.describe(methods)

    assert len(recs) == len(methods)
    assert {r["method"] for r in recs} == {m.name for m in methods}
    for r in recs:
        assert "parameters" in r and "differs_from_class_defaults" in r
        # everything must be JSON-serialisable, or the artefact silently fails to write
        import json as _json
        _json.dumps(r)

    # DWLS's departures from its PUBLISHED defaults must be declared, not implicit
    dwls = next(r for r in recs if r["method"] == "dwls")
    devs = dwls["declared_deviations_from_published_defaults"]
    assert devs, "DWLS runs with opened cutoffs and a cell cap; both must be declared"
    for d in devs:
        assert d["published_default"] and d["used"] and d["why"]
        assert "before any DWLS score existed" in d["decided"]


def test_coverage_reports_what_the_roster_drops_from_the_measurement():
    """
    Risk 3, reference coverage. The claim must come from counting the atlas's labels,
    not from an assumption about which populations are abundant — an earlier draft of
    this module asserted neurons dominate the dropped set, and they do not.
    """
    from ivygap.anatomic.coverage import reference_coverage

    labels = {"AC-like": 100, "MES-like": 100, "Mono": 500, "Neuron": 5, "OPC": 10}
    cov = reference_coverage(atlas_labels=labels)

    assert set(cov["atlas_labels_dropped"]) == {"Mono", "Neuron", "OPC"}
    assert cov["n_cells_dropped"] == 515
    assert cov["fraction_of_atlas_dropped"] == pytest.approx(515 / 715)
    # the myeloid consequence is named because Mono is the largest dropped population
    assert "Macrophage_Microglia" in cov["what_this_means"]
    # and the neuron correction is stated from the count, not asserted
    assert "only 5 neurons" in cov["what_this_means"]
    assert "frozen and hashed" in cov["correction_to_the_frozen_constraint_file"]
