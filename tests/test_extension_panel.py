"""The post-registration extension panel must never leak into the registered analysis.

Methods added on 2026-09-30 (ivygap.deconv.extension) were chosen after the data were seen. They
may be reported; they may not change a statistic registered before the data existed. These tests
pin the separation, and pin the D22 fix (the matrix layer of every re-measurement).
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from ivygap.deconv.extension import (EXTENSION_SPECS, _NoReimplementation,
                                     build_extension_methods, is_extension)
from ivygap.deconv.registry import build_methods

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("prefer_r", [True, False])
def test_no_extension_method_is_in_the_registered_panel(prefer_r):
    registered = {m.name for m in build_methods(prefer_r=prefer_r)}
    leaked = registered & set(EXTENSION_SPECS)
    assert not leaked, f"extension methods leaked into the registered panel: {leaked}"


def test_extension_r_methods_cannot_fall_back():
    for m in build_extension_methods():
        if getattr(m, "requires_r", False):
            assert m.allow_fallback is False, f"{m.name} could fall back to a stand-in"
            assert isinstance(m.fallback, _NoReimplementation)


def test_the_placeholder_refuses_to_produce_numbers():
    with pytest.raises(RuntimeError, match="no reimplementation by design"):
        _NoReimplementation("fardeep")._solve_all(None)


def test_every_extension_spec_carries_a_citation_and_a_literature_reason():
    for name, spec in EXTENSION_SPECS.items():
        assert is_extension(name)
        assert spec.citation and spec.why, f"{name} lacks a citation or a reason"
        assert "ACS" not in spec.why, f"{name}'s reason cites ACS -- ACS may not select methods"


def test_remeasure_passes_the_matrix_layer_D22():
    src = (ROOT / "scripts" / "remeasure_method.py").read_text()
    call = re.search(r"=\s*build_from_h5ad\((.*?)\)", src, re.S)   # the call, not a comment
    assert call and "matrix=matrix" in call.group(1), \
        "remeasure_method.py builds the reference without passing matrix= (D22)"
    assert '"8_matrix_layer"' in src, "the eighth equivalence condition is gone"
    assert "run_provenance.json" in src, "the layer must be read from the run's provenance"


def test_extension_estimates_are_kept_out_of_the_globbed_directory():
    # method_completeness.py and cdseq_matched_comparison.py glob results/estimates/ivygap_*.csv
    src = (ROOT / "scripts" / "remeasure_method.py").read_text()
    assert 'results/extension/ivygap_{M}.csv' in src
    stray = list((ROOT / "results" / "estimates").glob("ivygap_*.csv")) \
        if (ROOT / "results" / "estimates").exists() else []
    bad = [p.name for p in stray if any(n in p.name for n in EXTENSION_SPECS)]
    assert not bad, f"extension estimates in results/estimates/: {bad}"


def test_fardeep_departs_from_defaults_only_in_permn():
    drv = (ROOT / "R" / "run_fardeep.R").read_text()
    call = re.findall(r"fardeep\(X = sig_mat, Y = [^)]*\)", drv)
    assert call, "FARDEEP call not found"
    for c in call:
        args = re.sub(r"X = sig_mat, Y = bulk_mat(\[[^\]]*\])?", "", c)
        assert re.sub(r"[\s,()]|fardeep|permn = 0", "", args) == "", \
            f"FARDEEP called with a non-default argument: {c}"


def test_every_extension_method_carries_a_verified_doi():
    # Three of four citations were first written from memory and were wrong (journal, volume
    # or article number). The DOI is what a reader can resolve; it is required.
    for name, spec in EXTENSION_SPECS.items():
        assert re.fullmatch(r"10\.\d{4,9}/\S+", spec.doi or ""), f"{name}: missing or malformed DOI"


def test_a_gene_subset_fit_never_writes_the_shared_per_sample_file(monkeypatch):
    # D36: Extension E3 called extension_tcga.run with genes_keep and merged its half-gene unmix fits over the
    # registered all-gene rows of results/extension/estimates_full*_extension.csv. A subset fit is not the registered fit.
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import extension_tcga                                                    # noqa: PLC0415
    for v in ("IVYGAP_UNMIX_SHIFT", "IVYGAP_UNMIX_POWER"):
        monkeypatch.delenv(v, raising=False)
    assert extension_tcga.shares_per_sample_file(None)
    assert not extension_tcga.shares_per_sample_file(["GFAP", "OLIG2"])
    assert not extension_tcga.shares_per_sample_file([])
    monkeypatch.setenv("IVYGAP_UNMIX_SHIFT", "1")
    assert not extension_tcga.shares_per_sample_file(None)


def test_a_recide_run_with_its_own_report_keeps_its_estimates_beside_it():
    # D38's rule: only the default run writes the shared estimates file results/extension/ivygap_recide.csv.
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    import run_recide                                                        # noqa: PLC0415
    shared = run_recide.estimates_path(Path(run_recide.DEFAULT_OUT))
    assert shared.name == "ivygap_recide.csv" and shared.parent.name == "extension"
    # the harness passes a RELATIVE --out; the path must still resolve under the project, because the report records
    # it relative to PROJECT_ROOT (that call failed after a 93-minute run on 2026-10-08)
    from ivygap import config                                                # noqa: PLC0415
    rerun = run_recide.estimates_path(Path("results/verification/rerun_outputs/recide_anatomic.json"))
    assert rerun == config.PROJECT_ROOT / "results/verification/rerun_outputs/recide_anatomic.csv"
    assert str(rerun.relative_to(config.PROJECT_ROOT)) == "results/verification/rerun_outputs/recide_anatomic.csv"
    # a relative --out that names the registered report is the registered run (when the suite has not redirected
    # RESULTS_DIR to a temporary tree, the registered report lies under the project)
    default = Path(run_recide.DEFAULT_OUT)
    if default.is_relative_to(config.PROJECT_ROOT):
        assert run_recide.estimates_path(default.relative_to(config.PROJECT_ROOT)) == shared

