"""
Loader tests. These target the failure modes that would silently produce a wrong
analysis rather than an error: resolving the wrong column, losing the nested design,
or accepting a clinical table that cannot support survival analysis.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ivygap import config
from ivygap.data.load_ivygap import _normalise_structure, _resolve, build_manifest
from ivygap.data.clinical import _coerce_event


def test_anatomic_labels_normalise_to_the_protocol_vocabulary():
    """
    Verified against the real 2014-11-25 release. The five anatomic structures carry a
    '-reference-histology' suffix; the archive's CTmvp/CTpan become the protocol's
    MVP/PAN.
    """
    assert _normalise_structure("CT-reference-histology") == "CT"
    assert _normalise_structure("IT-reference-histology") == "IT"
    assert _normalise_structure("LE-reference-histology") == "LE"
    assert _normalise_structure("CTmvp-reference-histology") == "MVP"
    assert _normalise_structure("CTpan-reference-histology") == "PAN"


def test_ish_cluster_samples_are_not_mistaken_for_anatomic_ones():
    """
    The bug this guards against would have tripled the anatomic cohort with 148 samples
    whose labels were assigned USING expression — making the whole test circular.
    'CT-CD44' starts with 'CT', so any prefix match silently pulls it in.
    """
    for raw in ("CT-CD44", "CT-control-TGFBR2", "CTpnz-PI3", "CThbv-POSTN",
                "CTmvp-ITGA6"):
        assert not config.is_anatomic_sample(raw)
        assert _normalise_structure(raw) not in config.PRIMARY_STRUCTURES

    # A bare acronym with no suffix is a cluster label too, not an anatomic sample.
    assert not config.is_anatomic_sample("CT")


def test_column_resolution_accepts_alternate_names():
    df = pd.DataFrame(columns=["rna_well_id", "structure_abbreviation", "tumor_id"])
    assert _resolve(df, ["structure_acronym", "structure_abbreviation"], "x") == \
        "structure_abbreviation"


def test_column_resolution_names_what_it_found():
    """A resolution failure must be diagnosable without opening the file."""
    df = pd.DataFrame(columns=["something_else"])
    with pytest.raises(KeyError, match="something_else"):
        _resolve(df, ["structure_acronym"], "the structure")


def test_sample_level_id_mistaken_for_patient_is_rejected():
    """
    The nested design is the premise of the project. Resolving a sample-level column as
    the patient id would silently destroy it and every downstream 'patient-equal'
    aggregation would be a no-op.
    """
    samples = pd.DataFrame({
        "rna_well_id": [f"W{i}" for i in range(6)],
        "structure_acronym": ["CT"] * 6,
        "tumor_id": [f"W{i}" for i in range(6)],       # unique per sample
    })
    with pytest.raises(ValueError, match="sample-level identifier"):
        build_manifest(samples)


def test_manifest_separates_the_two_studies():
    """
    The archive holds two studies in one file: 122 H&E-selected anatomic samples and 148
    ISH-selected cancer-stem-cell cluster samples. Only the first is the instrument.
    """
    samples = pd.DataFrame({
        "rna_well_id": ["W1", "W2", "W3", "W4"],
        "structure_abbreviation": ["CT-reference-histology", "CT-CD44",
                                   "LE-reference-histology", "CTpnz-PI3"],
        "tumor_id": ["T1", "T1", "T2", "T2"],
    })
    m = build_manifest(samples)
    assert m["is_anatomic_study"].tolist() == [True, False, True, False]
    assert m["is_primary_structure"].tolist() == [True, False, True, False]
    assert m["structure"].tolist() == ["CT", "CT-CD44", "LE", "CTpnz-PI3"]


def test_event_coercion_handles_text_encodings():
    s = pd.Series(["Dead", "alive", "1:DECEASED", "0:LIVING"])
    assert _coerce_event(s).tolist() == [1.0, 0.0, 1.0, 0.0]


def test_unrecognised_event_value_becomes_nan_not_censored():
    """
    Defaulting an unknown status to 'censored' would convert deaths into survivors and
    bias every survival estimate in one direction.
    """
    assert pd.isna(_coerce_event(pd.Series(["unknown"]))).all()


def test_numeric_event_outside_zero_one_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        _coerce_event(pd.Series([0.0, 1.0, 2.0]))


def test_config_is_self_consistent():
    config.validate_config()
    assert set(config.MARKER_GENES) == set(config.CELL_TYPES)
    assert config.REFERENCE_STRUCTURE in config.PRIMARY_STRUCTURES


def test_synthetic_runs_cannot_overwrite_real_results(monkeypatch):
    """
    Regression for a defect that cost a completed real run.

    A fixture run and a real run originally wrote to identical paths, so
    `run_all.py --synthetic` silently replaced a 30-minute real result with synthetic
    numbers. The release manifest still said SYNTHETIC FIXTURE — nothing was
    mislabelled — but the real output was gone. Labelling is not isolation.
    """
    import importlib
    from ivygap import config as cfg

    fresh = importlib.reload(cfg)
    real_dirs = {fresh.RESULTS_DIR, fresh.RELEASE_DIR, fresh.ANATOMIC_DIR,
                 fresh.BENCH_DIR, fresh.ESTIMATES_DIR, fresh.SURVIVAL_DIR}

    fresh.use_synthetic_paths()
    synth_dirs = {fresh.RESULTS_DIR, fresh.RELEASE_DIR, fresh.ANATOMIC_DIR,
                  fresh.BENCH_DIR, fresh.ESTIMATES_DIR, fresh.SURVIVAL_DIR}

    assert not (real_dirs & synth_dirs), (
        f"synthetic and real runs share output paths: {real_dirs & synth_dirs}"
    )
    assert fresh.RUN_LABEL == "synthetic"
    assert all("synthetic" in str(d) for d in synth_dirs)

    # Raw inputs are deliberately SHARED — a fixture run should not re-download Ivy GAP.
    assert fresh.RAW_DIR == importlib.reload(cfg).RAW_DIR

    importlib.reload(cfg)          # leave the module as the suite found it


# =============================================================================
# one results tree, one run
# -----------------------------------------------------------------------------
# This project has been bitten twice. First a --synthetic fixture overwrote a completed
# 30-minute real run. Then a background run reported as killed had not actually died —
# only its wrapper had — and was still writing into the tree a fresh run had just
# cleared. Labelling is not isolation, and neither is trusting a "killed" message.
# =============================================================================

def test_run_lock_blocks_a_second_run(tmp_path, monkeypatch):
    from ivygap import config as cfg

    monkeypatch.setattr(cfg, "RESULTS_DIR", tmp_path)
    release = cfg.acquire_run_lock()
    try:
        with pytest.raises(cfg.ResultsTreeBusy, match="held by pid"):
            cfg.acquire_run_lock()
    finally:
        release()


def test_run_lock_is_reclaimed_when_the_holder_is_gone(tmp_path, monkeypatch):
    """A crashed run must not block the next one forever — the lock records a pid, and
    a pid that no longer exists is not a running process."""
    import json

    from ivygap import config as cfg

    monkeypatch.setattr(cfg, "RESULTS_DIR", tmp_path)
    # a pid that cannot be alive
    (tmp_path / ".run.lock").write_text(json.dumps({"pid": 999999, "started": "x"}))

    release = cfg.acquire_run_lock()          # must not raise
    assert (tmp_path / ".run.lock").exists()
    release()
    assert not (tmp_path / ".run.lock").exists()


def test_run_lock_releases_so_a_later_run_can_take_it(tmp_path, monkeypatch):
    from ivygap import config as cfg

    monkeypatch.setattr(cfg, "RESULTS_DIR", tmp_path)
    cfg.acquire_run_lock()()                  # acquire then immediately release
    cfg.acquire_run_lock()()                  # a second run must succeed
