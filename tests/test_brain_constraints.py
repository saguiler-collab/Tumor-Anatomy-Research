"""
The second constraint file — normal human brain.

The GBM set is frozen and registered. This one is not, and these tests exist mainly to
guarantee it CANNOT be used until it is: a constraint file that can be scored before it
is registered is not a pre-registration, and the difference has to be mechanical rather
than remembered.
"""
from __future__ import annotations

import pytest

from ivygap.anatomic import constraints as GBM
from ivygap.anatomic import constraints_brain as B


def test_a_draft_refuses_to_hash():
    """The load-bearing guard. Hashing a draft produces a number that looks like a
    commitment and is not one."""
    assert B.FROZEN is False, "if this file is frozen, the tests below must be revisited"
    with pytest.raises(RuntimeError, match="DRAFT"):
        B.freeze_hash()


def test_the_two_constraint_files_are_independent():
    """
    Editing the brain file must never move the GBM hash. They are separate registrations
    and a shared hash would couple them, so a correction here would void GBM results.
    """
    before = GBM.freeze_hash()
    _ = B.registration_payload()
    assert GBM.freeze_hash() == before


def test_no_gbm_structure_leaks_into_the_brain_set():
    """LE/IT/CT/MVP/PAN are glioblastoma structures. Normal brain has none of them."""
    gbm_structs = {"LE", "IT", "CT", "MVP", "PAN"}
    for c in B.CONSTRAINTS:
        assert not (set(c.regions) & gbm_structs), f"{c.id} names a GBM structure"
        assert not (set(c.among) & gbm_structs), f"{c.id}'s `among` names a GBM structure"


def test_every_region_named_is_declared_and_mapped():
    for c in B.CONSTRAINTS:
        for r in tuple(c.regions) + tuple(c.among):
            assert r in B.REGIONS, f"{c.id} names undeclared region {r!r}"
            assert r in B.REGION_MAP, f"{r!r} has no Allen mapping"
    for region, m in B.REGION_MAP.items():
        assert m.get("verified") is True, f"{region} mapping is unverified"
        assert m.get("allen_node") or m.get("allen_nodes"), f"{region} names no Allen node"


def test_every_cell_type_named_is_on_the_roster():
    for c in B.CONSTRAINTS:
        assert c.cell_type in B.CELL_TYPES, f"{c.id} names {c.cell_type!r}, not on the roster"


def test_monotone_renders_in_the_claimed_direction():
    """
    This rendered backwards on first write — "Neuron: CBC < CTX < WM", the exact
    opposite of B5 — because regions were reversed before joining. The GBM file lists
    lowest-first and joins in order; so does this.
    """
    b5 = next(c for c in B.CONSTRAINTS if c.id == "B5")
    assert b5.kind == "monotone"
    assert b5.describe() == "Neuron: WM < CTX < CBC"
    c7 = next(c for c in GBM.CONSTRAINTS if c.id == "C7")
    assert c7.describe() == "Tumor: LE < IT < CT"


def test_ids_are_unique_and_weights_are_positive():
    ids = [c.id for c in B.CONSTRAINTS]
    assert len(ids) == len(set(ids))
    assert all(c.weight > 0 for c in B.CONSTRAINTS)


def test_the_monotone_chain_carries_more_weight_than_a_single_step():
    """Same rationale as C7: a conjunction is harder to satisfy by chance."""
    mono = [c for c in B.CONSTRAINTS if c.kind == "monotone"]
    pair = [c for c in B.CONSTRAINTS if c.kind == "pairwise"]
    assert mono and pair
    assert min(c.weight for c in mono) > max(c.weight for c in pair)


def test_every_constraint_cites_evidence_that_is_not_a_deconvolution_result():
    """
    The circularity guard. A constraint learned from deconvolution output rigs the test.
    """
    banned = ("deconvolution result", "cibersortx analysis", "our estimate",
              "the leaderboard", "acs")
    for c in B.CONSTRAINTS:
        assert c.rationale.strip() and c.evidence.strip(), f"{c.id} lacks justification"
        low = (c.rationale + " " + c.evidence).lower()
        for b in banned:
            assert b not in low, f"{c.id} cites {b!r} — that would be circular"


def test_the_payload_excludes_prose_so_wording_fixes_do_not_break_a_registration():
    p = B.registration_payload()
    for c in p["constraints"]:
        assert "rationale" not in c and "evidence" not in c
        assert {"id", "kind", "cell_type", "regions", "weight"} <= set(c)


def test_the_payload_excludes_the_region_map():
    """
    A mapping correction is bookkeeping; changing what B1 asserts is not. Only the
    second may move the hash.
    """
    assert "region_map" not in B.registration_payload()


def test_untestable_claims_are_declared_with_reasons():
    assert B.UNTESTABLE
    for u in B.UNTESTABLE:
        assert u.get("fact") and u.get("why_untestable")


def test_summary_names_every_constraint():
    s = B.summary()
    for c in B.CONSTRAINTS:
        assert c.id in s
