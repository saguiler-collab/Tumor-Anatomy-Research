"""The independent (Abdelfattah 2022, GSE182109) reference: integrity of the build, pinned to its provenance.

prespecified/abdelfattah_cluster_mapping.md fixed the cluster-naming rule before any expression was
seen; these tests pin what the build actually did, so a re-run that silently changes it fails.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

REF = Path(__file__).resolve().parent.parent / "data" / "reference" / "abdelfattah_2022"


@pytest.fixture(scope="module")
def prov():
    p = REF / "provenance.json"
    if not p.exists():
        pytest.skip("independent reference not built")
    return json.loads(p.read_text())


def test_log_matrix_was_inverted_to_exact_counts(prov):
    m = prov["matrix"]
    assert m["inverted_to_counts"] is True
    assert m["max_rel_dev_from_integers"] < 1e-9          # whole numbers, to machine precision
    assert m["max_rel_error_totals_vs_library"] < 1e-9    # and they sum to the implied library


def test_reference_is_independent_of_gbmap_and_b_cells_are_not_few_donor(prov):
    assert "not among GBmap's 16 source studies" in prov["independent_of_gbmap"]
    assert prov["patients_per_type"]["B_cell"] >= 10      # GBmap: 78% of B cells from 3 donors


def test_cluster_names_follow_the_prespecified_rule(prov):
    mp = prov["cluster_mapping"]
    assert mp["C3"] == "T_cell" and mp["C11"] == "B_cell"
    assert prov["labelling"].startswith("pre-declared marker rule")
    assert set(prov["types"]) <= {"Tumor", "Macrophage_Microglia", "T_cell", "NK_cell", "B_cell",
                                  "Endothelial", "Oligodendrocyte", "Astrocyte"}


RES = Path(__file__).resolve().parent.parent / "results"
MAIN, IG = RES / "independent_atlas_test.json", RES / "independent_atlas_test_ig_removed.json"


def _orders(path: Path) -> dict:
    d = json.loads(path.read_text())
    return {(c, m): v for c, cc in d["cohorts"].items() for m, v in cc["methods"].items()}


@pytest.mark.skipif(not (MAIN.exists() and IG.exists()), reason="independent-atlas test not run")
def test_independent_atlas_result_is_mixed_and_ig_removal_does_not_change_direction():
    """Pins the reported reading (WHY_B_OVER_T §7i-7j): mixed under the independent atlas --
    NNLS recovers T > B in LGG, SVR does not in either cohort -- and the immunoglobulin-removed
    refit declared before the result agrees cell for cell."""
    main, ig = _orders(MAIN), _orders(IG)
    assert main[("lgg", "nnls")]["T_exceeds_B"] is True
    assert main[("lgg", "svr")]["T_exceeds_B"] is False
    assert main[("gbm", "svr")]["T_exceeds_B"] is False
    for key, v in main.items():
        assert ig[key]["ordering"] == v["ordering"], key
        # GBmap's own answer, carried for comparison, is B above T everywhere
        assert v["gbmap_frozen_for_comparison"]["frac_samples_B_over_T"] > 0.9, key
    assert json.loads(IG.read_text())["ig_removed"] is True


def test_seed_sensitivity_rule():
    """Addendum 2: stable only if every seed reproduces the registered ordering."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ias", Path(__file__).resolve().parent.parent / "scripts" / "independent_atlas_seeds.py")
    ias = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ias)
    assert ias.classify("T>B>NK", ["T>B>NK"] * 3) == "SEED-STABLE"
    assert ias.classify("T>B>NK", ["T>B>NK", "B>T>NK", "T>B>NK"]) == "SEED-DEPENDENT"
    assert ias.classify("T>B>NK", ["T>B>NK", None, "T>B>NK"]) == "BLOCKED"


SEEDS = RES / "independent_atlas_seed_sensitivity.json"


@pytest.mark.skipif(not SEEDS.exists(), reason="seed sensitivity not run")
def test_the_informative_cells_are_seed_stable():
    """Addendum 2 result: LGG NNLS (T > B) and SVR (B > T, both cohorts) hold in every build."""
    d = json.loads(SEEDS.read_text())
    if d.get("missing_runs"):
        pytest.skip("seed runs incomplete")
    cells = d["cells"]
    for k in ("lgg_nnls", "lgg_svr", "gbm_svr"):
        assert cells[k]["reading"] == "SEED-STABLE", k
    assert cells["lgg_nnls"]["registered"].startswith("T>B")
