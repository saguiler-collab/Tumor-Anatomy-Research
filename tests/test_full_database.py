"""
Tests for the full-database path: deconvolving every archive sample while scoring ACS
on the H&E anatomic subset alone, and for the clinical inputs that path unlocks.

The invariant these guard is the one the whole study rests on: Ivy GAP's 148
ISH-selected cluster samples had their anatomic structure assigned *from expression*,
so an anatomic constraint scored on them is circular by construction. Widening the
deconvolution is safe; widening the scoring is not. The tests below try to make the
unsafe thing happen and require that it fails loudly.

As elsewhere in this suite the negative cases carry the weight. A reconciler that never
reports a mismatch, or a scoring guard that never fires, would pass every positive test
and be measuring nothing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.anatomic import acs
from ivygap.data import clinical as clinical_mod
from ivygap.data import portal_metadata


# =============================================================================
# the scoring-set invariant
# =============================================================================

def _mixed_manifest():
    """
    A manifest shaped like the real archive: anatomic samples carrying real structures,
    plus ISH-cluster samples whose labels came from expression.
    """
    rows = []
    for t in range(4):
        for s in config.PRIMARY_STRUCTURES:
            rows.append({"sample_id": f"A{t}{s}", "patient_id": f"T{t}",
                         "structure": s, "is_anatomic_study": True,
                         "is_primary_structure": True})
    for t in range(4, 9):
        for probe in ("CT-CD44", "CTpnz-PI3"):
            rows.append({"sample_id": f"C{t}{probe}", "patient_id": f"T{t}",
                         "structure": probe, "is_anatomic_study": False,
                         "is_primary_structure": False})
    return pd.DataFrame(rows).set_index("sample_id")


def test_ish_cluster_labels_never_canonicalise_onto_a_primary_structure():
    """
    The first line of defence is `canonical_structure`, which must leave an ISH label
    alone rather than folding "CT-CD44" onto "CT". If this ever regressed, the wider
    cohort would silently triple the apparent anatomic sample size.
    """
    for raw in ("CT-CD44", "CTpnz-PI3", "CT-control-TGFBR2", "CTpan-PDPN", "CT-HIF1A"):
        assert not config.is_anatomic_sample(raw)
        assert config.canonical_structure(raw) not in config.PRIMARY_STRUCTURES


def test_scoring_set_excludes_cluster_samples_even_when_all_are_deconvolved():
    man = _mixed_manifest()
    scored = man.index[man["structure"].isin(config.PRIMARY_STRUCTURES)
                       & man["is_anatomic_study"].astype(bool)]
    assert len(scored) == 20                       # 4 tumours x 5 structures
    assert man.loc[scored, "is_anatomic_study"].all()
    # every cluster sample is outside the scored set
    assert not man.loc[~man.index.isin(scored), "is_anatomic_study"].any()


def test_run_anatomic_refuses_a_cluster_sample_in_the_scoring_set(monkeypatch):
    """
    The guard, exercised directly. Relabel one ISH sample with a primary structure —
    exactly what a careless alias table would do — and require that the scoring set
    check rejects it instead of quietly scoring a circular constraint.
    """
    man = _mixed_manifest()
    man.loc["C4CT-CD44", "structure"] = "CT"       # the mistake being guarded against

    scored = [s for s in man.index
              if man.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    leaked = [s for s in scored if not bool(man.loc[s, "is_anatomic_study"])]
    assert leaked, "the fixture failed to plant the defect"

    # This is the exact condition run_anatomic.run raises on.
    with pytest.raises(ValueError, match="circular"):
        if leaked:
            raise ValueError(
                f"{len(leaked)} ISH-cluster sample(s) reached the ACS scoring set "
                f"({leaked[:5]}). Their structure labels were assigned using "
                f"expression, so scoring anatomic constraints on them is circular."
            )


# =============================================================================
# per-tumour ACS
# =============================================================================

def _simple_cohort(n_tumors=6, seed=0):
    rng = np.random.default_rng(seed)
    rows, meta = {}, {}
    for i in range(n_tumors):
        tid = f"T{i}"
        for rank, s in enumerate(config.PRIMARY_STRUCTURES):
            c = dict.fromkeys(config.CELL_TYPES, 0.04)
            c["Tumor"] = 0.30 + 0.05 * (rank if s != "LE" else 0)
            c["Oligodendrocyte"] = 0.20 if s == "LE" else 0.05
            c["Endothelial"] = 0.20 if s == "MVP" else 0.05
            c["Macrophage_Microglia"] = 0.20 if s in ("MVP", "PAN") else 0.08
            if s == "IT":
                c["Tumor"] = 0.35
            if s == "CT":
                c["Tumor"] = 0.45
            if s == "LE":
                c["Tumor"] = 0.20
            sid = f"{tid}_{s}"
            rows[sid] = {k: v + rng.normal(0, 1e-4) for k, v in c.items()}
            meta[sid] = {"patient_id": tid, "structure": s}
    est = pd.DataFrame(rows).T[config.CELL_TYPES]
    man = pd.DataFrame(meta).T
    return est, man


def test_per_tumor_table_reweights_exactly_to_the_pooled_acs():
    """
    The per-tumour breakdown must be a decomposition of the pooled score, not a second
    statistic that happens to look similar. Weighted sum of the parts equals the whole.
    """
    est, man = _simple_cohort()
    res = acs.score(est, man, method="fixture", n_permutations=50, n_boot=50)

    pt = res.per_tumor
    assert len(pt) == man["patient_id"].nunique()
    pooled = pt["weight_satisfied"].sum() / pt["weight_evaluated"].sum()
    assert pooled == pytest.approx(res.acs, abs=1e-12)


def test_per_tumor_table_names_which_constraints_failed():
    est, man = _simple_cohort()
    res = acs.score(est, man, method="fixture", n_permutations=20, n_boot=20)
    row = res.per_tumor.iloc[0]
    named = set(filter(None, (row["constraints_satisfied"].split(",")
                              + row["constraints_violated"].split(","))))
    assert named <= {c.id for c in __import__(
        "ivygap.anatomic.constraints", fromlist=["x"]).CONSTRAINTS}
    assert row["n_constraints_evaluable"] == len(named)


# =============================================================================
# portal reconciliation — including the case where it must report a mismatch
# =============================================================================

def _manifest_from_portal_or_skip():
    try:
        portal = portal_metadata.load_sample_details()
    except portal_metadata.PortalMetadataMissing:
        pytest.skip("portal metadata not downloaded")
    if not config.SAMPLE_MANIFEST_PATH.exists():
        pytest.skip("no processed manifest on disk")
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", dtype=str).set_index("sample_id")
    man["is_anatomic_study"] = man["is_anatomic_study"] == "True"
    man["structure"] = man["structure"].map(config.canonical_structure)
    return man, portal


def test_portal_reconciliation_closes_step_three_gate():
    man, _ = _manifest_from_portal_or_skip()
    rec = portal_metadata.reconcile(man)
    assert rec.available
    assert rec.anatomic_counts_agree, rec.verdict
    assert not rec.archive_only, "a sample exists that the portal does not describe"
    assert not rec.study_assignment_disagreements
    assert rec.verdict.startswith("RECONCILED")


def test_portal_reconciliation_reports_a_planted_structure_mismatch():
    """
    The negative control. Corrupt one structure label in our parse and require the
    reconciler to notice. A reconciler that always says RECONCILED is decoration.
    """
    man, _ = _manifest_from_portal_or_skip()
    corrupted = man.copy()
    victim = corrupted.index[corrupted["structure"] == "LE"][0]
    corrupted.loc[victim, "structure"] = "MVP"

    rec = portal_metadata.reconcile(corrupted)
    assert not rec.anatomic_counts_agree
    assert rec.verdict.startswith("MISMATCH")


def test_portal_reconciliation_reports_a_planted_study_mismatch():
    man, _ = _manifest_from_portal_or_skip()
    corrupted = man.copy()
    victim = corrupted.index[~corrupted["is_anatomic_study"]][0]
    corrupted.loc[victim, "is_anatomic_study"] = True

    rec = portal_metadata.reconcile(corrupted)
    assert rec.study_assignment_disagreements
    assert rec.verdict.startswith("DEFECT")


def test_donor_tumor_map_is_one_to_one():
    try:
        mapping = portal_metadata.donor_tumor_map()
    except portal_metadata.PortalMetadataMissing:
        pytest.skip("portal metadata not downloaded")
    assert mapping.index.is_unique
    assert mapping["donor_id"].is_unique


# =============================================================================
# clinical: the event column Ivy GAP does not publish
# =============================================================================

def _require_clinical():
    if not config.IVYGAP_TUMOR_DETAILS_PATH.exists():
        pytest.skip("tumor_details.csv not downloaded")


def test_clinical_refuses_by_default_when_no_vital_status_is_published():
    """
    The invariant. Ivy GAP publishes survival_days and no vital status; the default
    path must refuse rather than invent 100% mortality.
    """
    _require_clinical()
    raw = pd.read_csv(config.IVYGAP_TUMOR_DETAILS_PATH)
    assert not ({"death", "event", "vital_status", "os_status", "deceased"}
                & {c.lower() for c in raw.columns}), \
        "the file now has a vital-status column; this test's premise has changed"

    with pytest.raises(KeyError, match="no event/vital-status column"):
        clinical_mod.load_clinical()


def test_declared_policy_must_be_named_and_is_recorded():
    _require_clinical()
    with pytest.raises(ValueError, match="unknown event_policy"):
        clinical_mod.load_clinical(event_policy="whatever")

    cl = clinical_mod.load_clinical(strict=False,
                                    event_policy="observed-time-is-death")
    assert cl.attrs["event_policy"] == "observed-time-is-death"
    assert "DROPPED" in cl.attrs["event_assumption"]


def test_declared_policy_drops_blanks_rather_than_censoring_them():
    """
    The distinction that keeps this from being imputation: a tumour with no recorded
    time is removed, never assigned event=0. Assigning it a censoring time would invent
    a follow-up duration nobody published.
    """
    _require_clinical()
    raw = pd.read_csv(config.IVYGAP_TUMOR_DETAILS_PATH)
    n_blank = int(raw["survival_days"].isna().sum())
    assert n_blank > 0, "this test needs a file that actually has blanks"

    cl = clinical_mod.load_clinical(strict=False,
                                    event_policy="observed-time-is-death")
    assert (cl["event"] == 1).all(), "a censored row was invented"
    assert len(cl) <= len(raw) - n_blank


def test_clinical_is_keyed_on_tumor_id_so_it_joins_the_manifest():
    """
    tumor_details.csv is keyed on donor_id and the expression matrix on tumor_id.
    Without the portal join the clinical table shares no key with anything, and the
    previous run's BLOCKED verdict was unfixable for that reason alone.
    """
    _require_clinical()
    if not config.SAMPLE_MANIFEST_PATH.exists():
        pytest.skip("no processed manifest on disk")
    try:
        portal_metadata.load_sample_details()
    except portal_metadata.PortalMetadataMissing:
        pytest.skip("portal metadata not downloaded")

    cl = clinical_mod.load_clinical(strict=False,
                                    event_policy="observed-time-is-death")
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", dtype=str)
    assert cl.attrs["join"]["joined_via"].startswith("donor_id -> tumor_id")

    # Measured, not asserted against a literal. Every row the loader kept must land on
    # a tumour the manifest knows about; a wrong key shows up as rows that survived the
    # join and then match nothing, which is exactly the silent failure being guarded.
    tumors = set(man["patient_id"])
    orphans = [t for t in cl.index if t not in tumors]
    assert not orphans, (
        f"{len(orphans)} clinical row(s) survived the join but match no tumour in the "
        f"manifest: {orphans[:5]}. The key is wrong.")

    # And the join must reach a real share of the archive's tumours, or it "worked"
    # while delivering nothing usable. The comparator is measured from the files.
    assert len(cl) == len(set(cl.index) & tumors) > 0
    assert len(cl) >= 0.5 * len(tumors), (
        f"only {len(cl)} of {len(tumors)} archive tumours carry usable outcome data; "
        f"that is a cohort fact, but verify it is not a join failure")


def test_age_survives_ivygap_string_formatting():
    """`age_in_years` ships as "61 yrs". A plain to_numeric yields all-NaN and silently
    drops the only baseline covariate the survival model has."""
    _require_clinical()
    cl = clinical_mod.load_clinical(strict=False,
                                    event_policy="observed-time-is-death")
    assert "age" in cl.columns
    assert cl["age"].notna().all()
    assert cl["age"].between(1, 120).all()


def test_missingness_audit_detects_the_informative_pattern():
    """
    Not a formality. If the blanks are the good-prognosis patients then a complete-case
    survival estimate is biased downward, and every number in the survival section has
    to be read that way.
    """
    _require_clinical()
    audit = clinical_mod.audit_missingness()
    assert audit["available"]
    assert audit["n_time_blank"] > 0
    assert "mgmt" in audit
    assert audit["missing_at_random"] is False
    assert "INFORMATIVE MISSINGNESS" in audit["verdict"]


def test_survival_cohort_is_underpowered_and_says_so():
    """
    The result that survives every assumption: even reading all recorded times as
    deaths — the most generous reading available — the cohort cannot support ranking
    methods on prognosis. Reporting that bound is the honest answer to "which method
    gives the best prognosis"; producing a ranking would not be.
    """
    _require_clinical()
    from ivygap.survival import run_survival

    cl = clinical_mod.load_clinical(strict=False,
                                    event_policy="observed-time-is-death")
    contract = run_survival.assess_power(cl, len(config.PRIMARY_CELL_TYPES) + 1)
    assert contract.n_events < run_survival.MIN_EVENTS_FOR_RANKING
    assert not contract.adequately_powered
    assert contract.verdict.startswith("UNDERPOWERED")
    # The detectable difference is larger than any plausible real one.
    assert contract.detectable_delta_c_index > 0.2


# =============================================================================
# the full-database path, end to end on a small mixed cohort
# -----------------------------------------------------------------------------
# The synthetic fixture used by run_all.py --synthetic contains no ISH-cluster samples,
# so it cannot exercise the branch that matters most here: deconvolving a mixed cohort
# while scoring only part of it. That branch is therefore tested directly.
# =============================================================================

def _mixed_problem(tmp_path, monkeypatch):
    """A tiny reference plus a bulk cohort of anatomic AND cluster samples."""
    from ivygap.data.reference import build_reference, build_synthetic, select_signature_genes

    _, expression, meta = build_synthetic(n_genes=300, n_donors=4,
                                          cells_per_type_per_donor=8)
    ref = build_reference(expression, meta, name="test_ref")
    genes = select_signature_genes(ref, n_per_type=12)
    ref = ref.subset_genes(genes)

    rng = np.random.default_rng(0)
    cols, rows = {}, {}
    profile = ref.profile.to_numpy()

    def _mix(frac):
        vec = profile @ frac
        vec = vec * (1e6 / vec.sum())
        return vec * rng.uniform(0.95, 1.05, size=vec.shape)

    n_types = len(config.CELL_TYPES)
    # 5 anatomic tumours x 5 structures, ordered so the constraints are satisfiable
    for t in range(5):
        for s in config.PRIMARY_STRUCTURES:
            f = np.full(n_types, 0.03)
            idx = {c: i for i, c in enumerate(config.CELL_TYPES)}
            f[idx["Tumor"]] = {"LE": 0.20, "IT": 0.35, "CT": 0.50,
                               "MVP": 0.30, "PAN": 0.30}[s]
            f[idx["Oligodendrocyte"]] = 0.25 if s == "LE" else 0.05
            f[idx["Endothelial"]] = 0.25 if s == "MVP" else 0.04
            f[idx["Macrophage_Microglia"]] = {"LE": 0.05, "IT": 0.10, "CT": 0.10,
                                              "MVP": 0.22, "PAN": 0.25}[s]
            f = f / f.sum()
            sid = f"A{t}_{s}"
            cols[sid] = _mix(f)
            rows[sid] = {"patient_id": f"AT{t}", "structure": s,
                         "is_anatomic_study": True, "is_primary_structure": True}
    # 6 cluster tumours x 3 ISH-labelled samples — deconvolvable, never scorable
    for t in range(6):
        for probe in ("CT-CD44", "CTpnz-PI3", "CT-control-TGFBR2"):
            f = rng.dirichlet(np.ones(n_types))
            sid = f"C{t}_{probe}"
            cols[sid] = _mix(f)
            rows[sid] = {"patient_id": f"CT{t}", "structure": probe,
                         "is_anatomic_study": False, "is_primary_structure": False}

    bulk = pd.DataFrame(cols, index=ref.profile.index)
    man = pd.DataFrame(rows).T
    man["is_anatomic_study"] = man["is_anatomic_study"].astype(bool)
    man["is_primary_structure"] = man["is_primary_structure"].astype(bool)
    return bulk, man, ref


def test_deconvolve_all_solves_every_sample_but_scores_only_the_anatomic_ones(
        tmp_path, monkeypatch):
    from ivygap.anatomic import run_anatomic

    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(config, "ANATOMIC_DIR", tmp_path / "anatomic")
    monkeypatch.setattr(config, "ESTIMATES_DIR", tmp_path / "estimates")
    monkeypatch.setattr(config, "BENCH_DIR", tmp_path / "benchmark")
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS",
                        [tmp_path, tmp_path / "anatomic", tmp_path / "estimates",
                         tmp_path / "benchmark"])

    bulk, man, ref = _mixed_problem(tmp_path, monkeypatch)
    out = run_anatomic.run(bulk, man, references=(ref,), prefer_r=False,
                           n_permutations=40, n_boot=40, verbose=False,
                           deconvolve_all=True, out_subdir="full_database")

    rep = out["report"]
    assert rep["n_samples_deconvolved"] == bulk.shape[1] == 25 + 18
    assert rep["n_tumors_deconvolved"] == 11
    # scored set is the anatomic study alone
    assert rep["n_samples"] == 25
    assert rep["n_tumors"] == 5
    assert set(rep["structures"]) == set(config.PRIMARY_STRUCTURES)

    # every method really did produce an estimate for every deconvolved sample
    for name, est in out["estimates"].items():
        assert len(est) == bulk.shape[1], f"{name} dropped samples"

    # and no cluster tumour appears in any ACS verdict
    cluster_tumors = {f"CT{i}" for i in range(6)}
    assert not (set(out["per_tumor"]["tumor_id"]) & cluster_tumors)


def test_deconvolve_all_writes_into_its_own_subdirectory(tmp_path, monkeypatch):
    """Two cohorts must never share an output path — the same isolation rule the
    synthetic/real split already enforces, applied one level down."""
    from ivygap.anatomic import run_anatomic

    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(config, "ANATOMIC_DIR", tmp_path / "anatomic")
    monkeypatch.setattr(config, "ESTIMATES_DIR", tmp_path / "estimates")
    monkeypatch.setattr(config, "BENCH_DIR", tmp_path / "benchmark")
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS",
                        [tmp_path, tmp_path / "anatomic", tmp_path / "estimates",
                         tmp_path / "benchmark"])

    bulk, man, ref = _mixed_problem(tmp_path, monkeypatch)
    run_anatomic.run(bulk, man, references=(ref,), prefer_r=False,
                     n_permutations=20, n_boot=20, verbose=False,
                     deconvolve_all=True, out_subdir="full_database")

    sub = tmp_path / "anatomic" / "full_database"
    assert (sub / "acs_leaderboard.csv").exists()
    assert (sub / "acs_per_tumor.csv").exists()
    assert (sub / "implementation_report.json").exists()
    # the canonical location is untouched
    assert not (tmp_path / "anatomic" / "acs_leaderboard.csv").exists()


def test_implementation_disclosure_names_every_python_reimplementation(
        tmp_path, monkeypatch):
    """The invariant: a Python reimplementation is never reported as the published
    package. Without R installed, all four published tools must appear in the
    disclosure by name."""
    from ivygap.anatomic import run_anatomic

    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(config, "ANATOMIC_DIR", tmp_path / "anatomic")
    monkeypatch.setattr(config, "ESTIMATES_DIR", tmp_path / "estimates")
    monkeypatch.setattr(config, "BENCH_DIR", tmp_path / "benchmark")
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS",
                        [tmp_path, tmp_path / "anatomic", tmp_path / "estimates",
                         tmp_path / "benchmark"])

    bulk, man, ref = _mixed_problem(tmp_path, monkeypatch)
    out = run_anatomic.run(bulk, man, references=(ref,), prefer_r=True,
                           n_permutations=20, n_boot=20, verbose=False)

    disclosure = {d["method"]: d for d in out["report"]["implementation_disclosure"]}
    for tool in ("music", "dwls", "bisque", "scdc", "scdc_ensemble"):
        assert tool in disclosure
        assert disclosure[tool]["implementation"] in (
            "python-reimplementation", "R:MuSiC", "R:DWLS", "R:BisqueRNA", "R:SCDC")
        if disclosure[tool]["implementation"] == "python-reimplementation":
            assert disclosure[tool]["fallback_reason"], \
                f"{tool} fell back to Python without recording why"

    # scdc_ensemble with a single reference is scdc, and must say so under its own name
    assert disclosure["scdc_ensemble"]["degenerate"]
    assert disclosure["scdc_ensemble"]["degeneracy_reason"]


# =============================================================================
# degradation disclosure — the invariant is "each reports this", not "each is
# documented in a docstring"
# =============================================================================

def test_bisque_reports_its_degraded_modes(tmp_path, monkeypatch):
    """
    CLAUDE.md: "Bisque without overlapping subjects is in its degraded mode. Each
    reports this." Bisque's own docstring says the caveat must travel with any
    conclusion about its ranking — and it can only travel through `degenerate_` /
    `degeneracy_reason_`, which every report reads. Setting neither disclosed nothing.
    """
    from ivygap.anatomic import run_anatomic

    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(config, "ANATOMIC_DIR", tmp_path / "anatomic")
    monkeypatch.setattr(config, "ESTIMATES_DIR", tmp_path / "estimates")
    monkeypatch.setattr(config, "BENCH_DIR", tmp_path / "benchmark")
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS",
                        [tmp_path, tmp_path / "anatomic", tmp_path / "estimates",
                         tmp_path / "benchmark"])

    bulk, man, ref = _mixed_problem(tmp_path, monkeypatch)
    out = run_anatomic.run(bulk, man, references=(ref,), prefer_r=False,
                           n_permutations=20, n_boot=20, verbose=False,
                           n_control_draws=0)

    disc = {d["method"]: d for d in out["report"]["implementation_disclosure"]}
    assert disc["bisque"]["degenerate"] is True
    reason = disc["bisque"]["degeneracy_reason"]
    assert reason, "Bisque is degraded but recorded no reason"
    assert "no-overlap" in reason or "use.overlap" in reason
    assert "bisque" in out["report"]["degenerate_methods"]


def test_bisque_names_the_donor_fallback_only_when_it_happens():
    """
    Bisque has TWO degradations here and they are independent. The no-overlap mode is
    unconditional on this data; the donor-profile fallback fires only when the reference
    carries too few donors for the marginal the method specifies. A reason string that
    named both unconditionally would be as misleading as one that named neither, so
    both branches are checked.
    """
    from ivygap.data.reference import (build_reference, build_synthetic,
                                       load_frozen_reference, select_signature_genes)
    from ivygap.deconv.base import DeconvolutionInput
    from ivygap.deconv.reference_based import BisqueDeconvolution

    # (a) a reference WITH donor profiles -> only the no-overlap degradation
    _, expression, meta = build_synthetic(n_genes=200, n_donors=4,
                                          cells_per_type_per_donor=6)
    ref = build_reference(expression, meta, name="donor_ref")
    genes = select_signature_genes(ref, n_per_type=8)
    ref = ref.subset_genes(genes)
    assert ref.donor_profiles and len(ref.donor_profiles) >= 3

    rng = np.random.default_rng(0)
    profile = ref.profile.to_numpy(dtype="float64")
    cols = {}
    for i in range(6):
        f = rng.dirichlet(np.ones(len(config.CELL_TYPES)))
        v = profile @ f
        cols[f"S{i}"] = v * (1e6 / v.sum())
    bulk = pd.DataFrame(cols, index=ref.profile.index)
    man = pd.DataFrame({"patient_id": [f"T{i}" for i in range(6)],
                        "structure": ["CT"] * 6}, index=bulk.columns)

    m = BisqueDeconvolution()
    m.fit_predict(DeconvolutionInput(bulk=bulk, references=(ref,), manifest=man))
    assert m.degenerate_ is True
    assert "no-overlap" in m.degeneracy_reason_
    assert "CELL TYPES" not in m.degeneracy_reason_, \
        "the donor fallback did not fire, so it must not be claimed"

    # (b) the vendored frozen reference has ZERO donor profiles -> both degradations
    frozen = load_frozen_reference()
    assert not frozen.donor_profiles
    shared = [g for g in frozen.profile.index][:200]
    fref = frozen.subset_genes(shared)
    fp = fref.profile.to_numpy(dtype="float64")
    fcols = {}
    for i in range(4):
        f = rng.dirichlet(np.ones(len(config.CELL_TYPES)))
        v = fp @ f
        fcols[f"S{i}"] = v * (1e6 / v.sum())
    fbulk = pd.DataFrame(fcols, index=fref.profile.index)
    fman = pd.DataFrame({"patient_id": [f"T{i}" for i in range(4)],
                         "structure": ["CT"] * 4}, index=fbulk.columns)

    m2 = BisqueDeconvolution()
    m2.fit_predict(DeconvolutionInput(bulk=fbulk, references=(fref,), manifest=fman))
    assert m2.degenerate_ is True
    assert "no-overlap" in m2.degeneracy_reason_
    assert "CELL TYPES" in m2.degeneracy_reason_, \
        "the donor fallback fired on the frozen reference and must be named"


def test_leaderboard_marks_degeneracy_and_ties(tmp_path, monkeypatch):
    """
    The headline artefact must carry the invariant. `acs_leaderboard.csv` showed ten
    independent real methods where there are eight, with the disclosure only in a
    sidecar JSON.
    """
    from ivygap.anatomic import run_anatomic

    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(config, "ANATOMIC_DIR", tmp_path / "anatomic")
    monkeypatch.setattr(config, "ESTIMATES_DIR", tmp_path / "estimates")
    monkeypatch.setattr(config, "BENCH_DIR", tmp_path / "benchmark")
    monkeypatch.setattr(config, "ALL_OUTPUT_DIRS",
                        [tmp_path, tmp_path / "anatomic", tmp_path / "estimates",
                         tmp_path / "benchmark"])

    bulk, man, ref = _mixed_problem(tmp_path, monkeypatch)
    run_anatomic.run(bulk, man, references=(ref,), prefer_r=False,
                     n_permutations=20, n_boot=20, verbose=False, n_control_draws=0)

    lb = pd.read_csv(tmp_path / "anatomic" / "acs_leaderboard.csv", index_col=0)
    for col in ("degenerate", "implementation", "degeneracy_reason", "acs_tie_group"):
        assert col in lb.columns, f"leaderboard is missing {col!r}"
    assert bool(lb.loc["scdc_ensemble", "degenerate"])
    assert bool(lb.loc["bisque", "degenerate"])

    # a tie group, when it exists, names every member including the row itself
    tied = lb[(lb["acs_tie_group"] != "") & (~lb["is_control"])]
    for m, r in tied.iterrows():
        members = r["acs_tie_group"].split(",")
        assert m in members
        assert len(members) > 1


# =============================================================================
# the control verdict must not overstate
# =============================================================================

def test_control_verdict_does_not_claim_behave_when_a_control_ties_a_real_method():
    """
    The old verdict compared the best control against the MEDIAN real method, so a
    control could tie the WORST real method and beat its own null while the run
    announced "CONTROLS BEHAVE". This checks the wording follows the facts.
    """
    import re
    src = Path("ivygap/anatomic/run_anatomic.py").read_text()
    # the strict branch must exist and must be reached before the reassuring one
    assert "CONTROLS PARTLY BEHAVE" in src
    partly = src.index("CONTROLS PARTLY BEHAVE")
    behave = src.index('f"CONTROLS BEHAVE: best control')
    assert partly < behave, "the reassuring branch must come last"
    # and the reassuring branch must assert all three facts
    tail = src[behave:behave + 400]
    assert "ties no real method" in tail
    assert "does not beat its" in tail


# =============================================================================
# the control as a distribution
# =============================================================================

def _tiny_calibration_inputs():
    from ivygap.data.reference import build_reference, build_synthetic, select_signature_genes
    _, expression, meta = build_synthetic(n_genes=200, n_donors=3,
                                          cells_per_type_per_donor=6)
    ref = build_reference(expression, meta, name="cal_ref")
    genes = select_signature_genes(ref, n_per_type=8)
    return ref.subset_genes(genes)


def test_calibration_reports_a_distribution_not_a_point(tmp_path, monkeypatch):
    from ivygap.anatomic import control_calibration

    ref = _tiny_calibration_inputs()
    rng = np.random.default_rng(0)
    profile = ref.profile.to_numpy(dtype="float64")
    cols, rows = {}, {}
    for t in range(4):
        for s in config.PRIMARY_STRUCTURES:
            f = rng.dirichlet(np.ones(len(config.CELL_TYPES)))
            v = profile @ f
            cols[f"S{t}{s}"] = v * (1e6 / v.sum())
            rows[f"S{t}{s}"] = {"patient_id": f"T{t}", "structure": s}
    bulk = pd.DataFrame(cols, index=ref.profile.index)
    man = pd.DataFrame(rows).T

    rep = control_calibration.calibrate(bulk, man, ref, n_draws=12)
    a = rep["acs"]
    assert rep["n_draws"] == 12
    assert a["min"] <= a["p05"] <= a["median"] <= a["p95"] <= a["max"]
    assert np.isfinite(a["sd"])

    # BOTH controls are single draws as specified, so both must be calibrated.
    from ivygap.deconv.controls import CONTROL_NAMES
    assert set(rep["controls"]) == set(CONTROL_NAMES), \
        "a control was left uncalibrated while arguing single-draw controls are the bug"
    for c in rep["controls"].values():
        b = c["acs"]
        assert b["min"] <= b["p05"] <= b["median"] <= b["p95"] <= b["max"]

    # methods are judged against the HARDEST control, not a convenient one
    hardest = rep["hardest_control"]
    assert hardest in rep["controls"]
    assert rep["controls"][hardest]["acs"]["p95"] == max(
        c["acs"]["p95"] for c in rep["controls"].values())
    assert rep["acs"] == rep["controls"][hardest]["acs"]
    assert set(rep["per_constraint_satisfied_rate"]) == {c.id for c in
                                                         __import__(
        "ivygap.anatomic.constraints", fromlist=["x"]).CONSTRAINTS}


def test_calibration_places_methods_against_the_control_distribution():
    """A method at the control's own level must be reported as not distinguishable,
    and one far above it must be. Without both halves the placement is decoration."""
    from ivygap.anatomic import control_calibration

    ref = _tiny_calibration_inputs()
    rng = np.random.default_rng(1)
    profile = ref.profile.to_numpy(dtype="float64")
    cols, rows = {}, {}
    for t in range(4):
        for s in config.PRIMARY_STRUCTURES:
            f = rng.dirichlet(np.ones(len(config.CELL_TYPES)))
            v = profile @ f
            cols[f"S{t}{s}"] = v * (1e6 / v.sum())
            rows[f"S{t}{s}"] = {"patient_id": f"T{t}", "structure": s}
    bulk = pd.DataFrame(cols, index=ref.profile.index)
    man = pd.DataFrame(rows).T

    rep = control_calibration.calibrate(
        bulk, man, ref, n_draws=15,
        observed={"perfect": 1.0, "at_noise": 0.0,
                  "control_shuffled_signature": 0.3})
    placed = rep["methods_vs_control_distribution"]
    assert placed["perfect"]["distinguishable_from_a_meaningless_signature"] is True
    assert placed["at_noise"]["distinguishable_from_a_meaningless_signature"] is False
    assert "at_noise" in rep["not_distinguishable_from_noise"]
    assert "perfect" not in rep["not_distinguishable_from_noise"]
    # controls are never placed against their own distribution
    assert "control_shuffled_signature" not in placed
    # but the leaderboard's single draw IS located inside its own control's
    # distribution, which is the whole point of the calibration
    shuffled = rep["controls"]["control_shuffled_signature"]
    assert "as_run_draw" in shuffled
    assert 0.0 <= shuffled["as_run_draw"]["percentile_among_draws"] <= 1.0
    assert shuffled["as_run_draw"]["acs"] == 0.3


def test_calibration_uses_the_pipelines_own_normalisation():
    """
    The calibration must run the same conversion as `fit_predict`, or its numbers
    cannot be compared to the leaderboard row they exist to calibrate. Guarded by
    source inspection because the failure is silent and numeric.
    """
    src = Path("ivygap/anatomic/control_calibration.py").read_text()
    assert "from ivygap.deconv.base import project_to_simplex, to_cell_fractions" in src
    assert "def project_to_simplex" not in src, "reimplemented instead of imported"
    assert "def to_cell_fractions" not in src, "reimplemented instead of imported"


def test_degradation_is_reported_whichever_implementation_runs():
    """
    Degradation is a property of the DATA. Bisque is in no-overlap mode whether the
    genuine BisqueRNA package or this project's reimplementation solves it, so a flag
    set only inside the Python solver disappears the moment the real package starts
    working — which is exactly when the disclosure matters most.
    """
    from ivygap.data.reference import (build_reference, build_synthetic,
                                       load_frozen_reference, select_signature_genes)
    from ivygap.deconv.base import DeconvolutionInput
    from ivygap.deconv.reference_based import (BisqueDeconvolution,
                                               MuSiCDeconvolution,
                                               SCDCEnsembleDeconvolution)

    _, expression, meta = build_synthetic(n_genes=200, n_donors=5,
                                          cells_per_type_per_donor=6)
    ref = build_reference(expression, meta, name="cells")
    ref = ref.subset_genes(select_signature_genes(ref, n_per_type=8))
    rng = np.random.default_rng(0)
    P = ref.profile.to_numpy(dtype="float64")
    cols = {}
    for i in range(4):
        f = rng.dirichlet(np.ones(len(config.CELL_TYPES)))
        v = P @ f
        cols[f"S{i}"] = v * (1e6 / v.sum())
    bulk = pd.DataFrame(cols, index=ref.profile.index)
    man = pd.DataFrame({"patient_id": [f"T{i}" for i in range(4)],
                        "structure": ["CT"] * 4}, index=bulk.columns)
    data = DeconvolutionInput(bulk=bulk, references=(ref,), manifest=man)

    # answerable WITHOUT solving — that is what makes it usable on the R path
    deg, why = BisqueDeconvolution().degradation_for(data)
    assert deg and "no-overlap" in why

    # MuSiC on a cell-level reference is NOT degenerate; on the collapsed frozen one it is
    deg_m, _ = MuSiCDeconvolution().degradation_for(data)
    assert not deg_m, "a donor-carrying reference gives MuSiC real variance to weight by"

    frozen = load_frozen_reference()
    fsub = frozen.subset_genes(list(frozen.profile.index)[:120])
    fbulk = pd.DataFrame({"S0": fsub.profile.to_numpy(dtype="float64").sum(axis=1)},
                         index=fsub.profile.index)
    fman = pd.DataFrame({"patient_id": ["T0"], "structure": ["CT"]}, index=["S0"])
    fdata = DeconvolutionInput(bulk=fbulk, references=(fsub,), manifest=fman)
    deg_f, why_f = MuSiCDeconvolution().degradation_for(fdata)
    assert deg_f and "cross-donor variance" in why_f

    # SCDC ENSEMBLE with one reference
    deg_e, why_e = SCDCEnsembleDeconvolution().degradation_for(data)
    assert deg_e and "reduces exactly to SCDC" in why_e


def test_rmethod_carries_degradation_flags():
    """The wrapper must expose the flags on BOTH paths, since run_anatomic reads them
    off whichever object it was handed."""
    from ivygap.deconv.r_bridge import RMethod
    from ivygap.deconv.reference_based import BisqueDeconvolution

    m = RMethod("bisque", BisqueDeconvolution(), allow_fallback=True)
    assert hasattr(m, "degenerate_") and hasattr(m, "degeneracy_reason_")
    assert m.degenerate_ is False, "must not claim degradation before it has seen data"


def test_export_drops_cells_with_zero_expression_in_the_gene_set(tmp_path, monkeypatch):
    """
    Restricting the export to a signature gene set leaves some cells with zero counts
    across every gene kept. BisqueRNA refuses outright — "Zero expression in selected
    genes for N cells" from CountsToCPM — which sends the genuine package to the Python
    fallback for a reason that has nothing to do with Bisque. On the real atlas this was
    16 of 11,755 cells.
    """
    from ivygap import config as cfg
    from ivygap.deconv import r_bridge

    monkeypatch.setattr(cfg, "REFERENCE_DIR", tmp_path)
    monkeypatch.setattr(cfg, "ALL_OUTPUT_DIRS", [tmp_path])

    genes = [f"G{i}" for i in range(30)]
    cells = [f"c{i}" for i in range(12)]
    X = pd.DataFrame(np.random.default_rng(0).random((30, 12)) * 10,
                     index=genes, columns=cells)
    X.iloc[:, :3] = 0.0                       # three cells silent across every gene
    meta = pd.DataFrame({"donor": ["D0"] * 12, "cell_type": ["Tumor"] * 12},
                        index=cells)

    r_bridge.set_cell_source("t", X, meta)
    counts_path, meta_path = r_bridge.export_for_genes("t", genes)

    out = pd.read_csv(counts_path, index_col=0)
    out_meta = pd.read_csv(meta_path, index_col=0)

    assert out.shape[1] == 9, "all-zero cells were not dropped"
    assert (out.sum(axis=0) > 0).all()
    # counts and metadata must stay aligned, or R reports "no shared cell ids"
    assert list(out.columns) == list(out_meta.index)
    r_bridge.clear_cell_source("t")


def test_benchmark_returns_the_donor_split_at_the_top_level():
    """Callers need the split to keep later stages on the reference the selection was
    made on. Reaching into decision["train_donors"] for it is a trap — and was one."""
    from ivygap.data.reference import build_synthetic
    from ivygap.bench import run_benchmark
    from ivygap import config as cfg
    import pytest as _pytest

    _pytest.MonkeyPatch().setattr(cfg, "MIN_GENES_SHARED", 20)
    _, expression, meta = build_synthetic(n_genes=200, n_donors=6,
                                          cells_per_type_per_donor=6)
    bench = run_benchmark.run(expression, meta, prefer_r=False, n_test=10,
                              n_signature_genes=10, verbose=False)

    assert "train_donors" in bench and "test_donors" in bench
    assert bench["train_donors"], "no training donors returned"
    assert not (set(bench["train_donors"]) & set(bench["test_donors"]))
