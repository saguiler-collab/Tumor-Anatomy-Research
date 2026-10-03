"""
Integrity of results/atlas_t_vs_b.json (prespecified/atlas_t_vs_b.md): the artefact must be internally
consistent before any sentence quotes it -- identity check passed, per-patient counts sum to the pooled
counts, the declared-rule T equals strict T plus ties, and the reading follows the rule from the counts.
"""
from __future__ import annotations

import json

import pytest

from ivygap import config

ART = config.PROJECT_ROOT / "results" / "atlas_t_vs_b.json"   # never the conftest-redirected RESULTS_DIR


@pytest.fixture(scope="module")
def art():
    if not ART.exists():
        pytest.skip("results/atlas_t_vs_b.json not produced yet (scripts/atlas_t_vs_b.py)")
    return json.loads(ART.read_text())


def test_identity_check_passed(art):
    assert art["identity_check"]["integral_counts"] is True
    assert art["identity_check"]["library_max_rel_error"] < 1e-6


def test_per_patient_counts_sum_to_pooled(art):
    for k in ("T_strict", "NK", "unresolved_tie", "B", "C12_unmapped", "T_declared_rule", "n_cells"):
        assert sum(r[k] for r in art["per_patient"]) == art["pooled"][k], k
    for r in art["per_patient"]:
        assert r["T_declared_rule"] == r["T_strict"] + r["unresolved_tie"]
    assert art["n_patients"] == len(art["per_patient"])


def test_reading_follows_the_rule_from_the_counts(art):
    """Recomputed here from the per-patient rows, not read from the artefact's own verdict."""
    rows, po = art["per_patient"], art["pooled"]
    n = len(rows)
    t_over_b = sum(r["T_strict"] > r["B"] for r in rows)
    b_over_t = sum(r["B"] >= r["T_strict"] for r in rows)
    if po["T_strict"] > po["B"] and t_over_b > n / 2:
        expected = "SUPPORTS T > B"
    elif po["B"] >= po["T_strict"] or b_over_t > n / 2:
        expected = "CONTRADICTS T > B"
    else:
        expected = "MIXED"
    assert art["reading"] == expected
    assert art["n_patients_T_strict_over_B"] == t_over_b


def test_gbmap_counts_are_consistent():
    """results/gbmap_t_vs_b.json (scripts/gbmap_t_vs_b.py): donor counts partition correctly."""
    p = config.PROJECT_ROOT / "results" / "gbmap_t_vs_b.json"
    if not p.exists():
        pytest.skip("results/gbmap_t_vs_b.json not produced yet")
    g = json.loads(p.read_text())
    assert g["n_donors_T_over_B"] + g["n_donors_B_at_or_over_T"] == g["n_donors_with_any_lymphoid"]
    assert g["n_donors_with_any_lymphoid"] <= g["n_donors"]
    assert set(g["disease"]) == {"glioblastoma"}
