"""
The generated RESULTS.md header must come from the artefacts, not from memory.

This exists because of a specific failure: after the 2026-09-06 rerun, RESULTS.md's body
reported rho = 0.733 on 13 methods while the hand-written header above it still claimed
rho = 0.873 on 10, from a run two generations old. The body was generated; the header was
not. These tests pin the header to the same artefacts as everything under it.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location(
    "summarize_results",
    Path(__file__).resolve().parent.parent / "scripts" / "summarize_results.py")
summarize = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summarize)


@pytest.fixture
def fake_results(tmp_path):
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    (anat / "anatomic_report.json").write_text(json.dumps({
        "constraint_freeze_hash": "deadbeef" * 8,
        "n_samples": 122, "n_tumors": 10, "n_permutations": 10000,
        "structures": ["LE", "IT", "CT", "MVP", "PAN"],
        "deconvolved_all_samples": True,
        "n_samples_deconvolved": 270, "n_tumors_deconvolved": 37,
    }))
    (anat / "equal_footing_certificate.json").write_text(json.dumps({
        "input_hashes": {"n_genes": "657", "references": "gbmap:abc123"}}))
    pd.DataFrame([{
        "yardstick": "synthetic_mixtures", "rho": 0.7329538, "p_value": 0.0043686,
        "ci_low": 0.2709, "ci_high": 0.9437, "n_methods": 13,
        "verdict": "HEADLINE: ACS ranking tracks true accuracy.",
    }]).to_csv(anat / "agreement_test.csv", index=False)
    pd.DataFrame([{"method": "music", "acs": 1.0}]).to_csv(
        anat / "acs_leaderboard.csv", index=False)
    return tmp_path


def test_header_reports_the_rho_the_artefact_holds(fake_results):
    out = summarize.header_section(fake_results)
    assert "0.7330" in out, "the headline rho must come from agreement_test.csv"
    assert "13 methods" in out
    assert "0.271" in out and "0.944" in out, "the bootstrap CI must be reported"


def test_header_cannot_drift_from_the_body(fake_results):
    """The exact regression: header and body disagreeing about the same run."""
    header = summarize.header_section(fake_results)
    agree = pd.read_csv(fake_results / "anatomic" / "agreement_test.csv")
    rho = agree.loc[0, "rho"]
    assert f"{rho:.4f}" in header
    # and a stale value must NOT appear
    assert "0.873" not in header


def test_header_carries_the_freeze_hash_and_cohort(fake_results):
    out = summarize.header_section(fake_results)
    assert "deadbeef" in out
    assert "122" in out and "10 tumours" in out
    assert "270 samples deconvolved" in out


def test_header_survives_a_run_with_no_agreement_test(fake_results):
    """
    A run can legitimately produce no agreement test — too few comparable methods. The
    header must say so rather than raising or inventing a number.
    """
    (fake_results / "anatomic" / "agreement_test.csv").unlink()
    out = summarize.header_section(fake_results)
    assert "No agreement test" in out


def test_the_shipped_results_md_matches_the_generator():
    """
    RESULTS.md in the repo must be what the generator currently produces. If this fails,
    someone edited it by hand or a run happened without regenerating it.
    """
    root = Path(__file__).resolve().parent.parent
    results = root / "results"
    shipped = root / "RESULTS.md"
    if not results.exists() or not shipped.exists():
        pytest.skip("no results tree or RESULTS.md in this checkout")
    regenerated = summarize.render(results)
    assert shipped.read_text().strip() == regenerated.strip(), (
        "RESULTS.md is out of sync with results/. Regenerate it:\n"
        "    python3 scripts/summarize_results.py > RESULTS.md")


# --- the figure blob must come from the artefacts too -------------------------

BFD_SPEC = importlib.util.spec_from_file_location(
    "build_figure_data",
    Path(__file__).resolve().parent.parent / "scripts" / "build_figure_data.py")
bfd = importlib.util.module_from_spec(BFD_SPEC)
BFD_SPEC.loader.exec_module(bfd)


def test_control_band_is_read_from_the_calibration_distribution(tmp_path):
    """
    The bug this pins: `control_calibration.json` keys its controls BY NAME, and reading
    it as a list of dicts yields bare strings, produces an empty band, and the plate
    draws the noise floor as a zero-width line — silently, with every method appearing
    to clear it by miles. The band is the whole point of that figure.
    """
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    (anat / "control_calibration.json").write_text(json.dumps({
        "n_draws": 200,
        "hardest_control": "control_shuffled_signature",
        "controls": {
            "control_shuffled_signature": {
                "n_draws": 200,
                "acs": {"mean": 0.386, "sd": 0.0656, "min": 0.215,
                        "p05": 0.2915, "median": 0.3846, "p95": 0.4923, "max": 0.569},
            },
            "control_random": {"n_draws": 200, "acs": {"mean": 0.40, "p95": 0.55}},
        },
    }))
    pd.DataFrame([
        {"method": "music", "acs": 1.0, "is_control": False, "comparable": True},
        {"method": "weak", "acs": 0.45, "is_control": False, "comparable": True},
        {"method": "control_shuffled_signature", "acs": 0.138,
         "is_control": True, "comparable": True},
    ]).to_csv(anat / "acs_leaderboard.csv", index=False)

    d = bfd.build(tmp_path)
    assert d["control_band"], "the band must not be empty when calibration exists"
    assert d["control_band"]["p95"] == 0.4923
    assert d["control_band"]["n_draws"] == 200, "a band from one draw is not a band"
    assert d["hardest_control"] == "control_shuffled_signature"
    # `weak` sits at 0.45, under the band's p95 — it must be named, not silently passed
    assert d["not_distinguishable"] == ["weak"]


def test_figure_blob_is_json_serialisable_with_no_nan(tmp_path):
    """NaN is valid in pandas and invalid in JSON; an unconverted one breaks the page."""
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    pd.DataFrame([
        {"method": "m", "acs": 0.9, "is_control": False, "comparable": True,
         "acs_tie_group": None, "null_p": float("nan")},
    ]).to_csv(anat / "acs_leaderboard.csv", index=False)
    d = bfd.build(tmp_path)
    text = json.dumps(d)          # must not raise
    assert "NaN" not in text and "Infinity" not in text


def test_blob_counts_comparable_separately_from_real(tmp_path):
    """
    A partial-coverage method is real but not comparable. Collapsing the two counts is
    how quanTIseq ends up quoted as though it were ranked alongside the rest.
    """
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    pd.DataFrame([
        {"method": "music", "acs": 1.0, "is_control": False, "comparable": True},
        {"method": "quantiseq", "acs": 0.6, "is_control": False, "comparable": False},
        {"method": "control_random", "acs": 0.4, "is_control": True, "comparable": True},
    ]).to_csv(anat / "acs_leaderboard.csv", index=False)
    d = bfd.build(tmp_path)
    assert d["n_real_methods"] == 2
    assert d["n_comparable_methods"] == 1
    assert d["n_controls"] == 1


def test_atlas_facts_are_never_the_cohorts_facts(tmp_path, monkeypatch):
    """
    The bug: the rendered caption read "reference 314,700 cells / 10 donors". Both were
    the wrong quantity — 314,700 is atlas cells whose label the roster maps (the build
    subsamples to ~15k), and 10 is Ivy GAP PATIENTS, not atlas donors (there are 110).
    `n_patients` on the anatomic certificate describes the tissue, not the reference.
    """
    from ivygap import config

    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    pd.DataFrame([{"method": "m", "acs": 1.0, "is_control": False, "comparable": True}]
                 ).to_csv(anat / "acs_leaderboard.csv", index=False)
    (anat / "equal_footing_certificate.json").write_text(json.dumps({
        "n_patients": 10, "input_hashes": {"n_genes": "657"}}))
    (anat / "reference_coverage.json").write_text(json.dumps({
        "n_cells_in_atlas": 338564, "n_cells_kept": 314700}))

    refdir = tmp_path / "refdir"
    refdir.mkdir()
    (refdir / "reference_sampling_gbmap.json").write_text(json.dumps({
        "n_cells_kept": 15311, "n_donors_kept": 110, "n_genes_kept": 16758,
        "n_cells_in_file": 338564, "n_cells_after_roster_mapping": 314700}))
    monkeypatch.setattr(config, "REFERENCE_DIR", refdir)

    atlas = bfd.build(tmp_path)["cohort"]["atlas"]
    assert atlas["donors"] == 110, "atlas donors must not come from the cohort's n_patients"
    assert atlas["cells_kept"] == 15311, "must be the built reference, not the mapped atlas"
    assert atlas["cells_after_roster_mapping"] == 314700   # kept, but under its own name


def test_atlas_facts_report_none_rather_than_a_wrong_number(tmp_path, monkeypatch):
    """With no sampling record, report nothing rather than substitute a different quantity."""
    from ivygap import config
    anat = tmp_path / "anatomic"
    anat.mkdir(parents=True)
    pd.DataFrame([{"method": "m", "acs": 1.0, "is_control": False, "comparable": True}]
                 ).to_csv(anat / "acs_leaderboard.csv", index=False)
    (anat / "equal_footing_certificate.json").write_text(json.dumps({"n_patients": 10}))
    empty = tmp_path / "empty"; empty.mkdir()
    monkeypatch.setattr(config, "REFERENCE_DIR", empty)

    atlas = bfd.build(tmp_path)["cohort"]["atlas"]
    assert atlas["donors"] is None and atlas["cells_kept"] is None
    assert "unavailable" in atlas["source"]
