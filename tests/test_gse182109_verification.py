"""The processed Abdelfattah matrix is GEO GSE182109's raw counts (scripts/verify_gse182109_raw.py)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ART = Path(__file__).resolve().parent.parent / "results" / "gse182109_raw_verification.json"
pytestmark = pytest.mark.skipif(not ART.exists(), reason="run scripts/verify_gse182109_raw.py first")


@pytest.fixture(scope="module")
def v():
    return json.loads(ART.read_text())


def test_every_sample_is_complete_and_named_as_geo_names_it(v):
    assert v["n_gsm_with_complete_triplet"] == v["n_gsm_in_geo_table"] == 44
    assert v["all_names_match_geo_table"] and v["all_dims_agree"]
    assert not v["unexpected_files"]


def test_processed_cells_come_from_their_own_samples(v):
    assert v["processed_cells_found_in_raw"] / v["processed_cells"] > 0.999


def test_every_processed_gene_is_a_raw_feature(v):
    g = v["genes"]
    assert g["processed_genes_in_raw_union"] == g["n_processed_genes"]


def test_recovered_counts_equal_geo_raw_counts(v):
    assert v["exactness_samples"], "no sample was checked for exactness"
    for gsm, r in v["exactness_samples"].items():
        assert r["fraction_entries_exact"] == 1.0, gsm
        assert r["max_abs_count_difference"] < 1e-6, gsm
        assert r["nonzeros_match_count"], gsm           # nothing dropped, nothing added
        assert r["max_rel_library_vs_raw_total_kept_genes"] < 1e-9, gsm
