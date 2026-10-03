"""
The data inventory (scripts/build_data_inventory.py) -- that it can tell a right
citation from a conflated one, and that what it says about the code is still true.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("inv", ROOT / "scripts" / "build_data_inventory.py")
inv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(inv)

STATUSES = {"used", "derived", "on disk, unused", "duplicate", "considered"}
CACHE = json.loads(inv.CACHE.read_text()) if inv.CACHE.exists() else None


def test_registry_is_well_formed():
    ids = [e["id"] for e in inv.REGISTRY]
    assert len(ids) == len(set(ids))
    for e in inv.REGISTRY:
        assert e["status"] in STATUSES, e["id"]
        assert 0 <= e["group"] < len(inv.GROUPS), e["id"]
        if e["status"] == "used":
            assert e["consumers"], f"{e['id']} is 'used' but names no script that reads it"
        if e["doi"]:
            assert e["first_author"] and e["year"], e["id"]


def test_every_consumer_still_reads_what_it_is_said_to_read():
    """A renamed folder or a refactor that drops a dataset shows up here."""
    missing = [(e["id"], c["file"], c["token"]) for e in inv.REGISTRY
               for c in inv.consumer_check(e["consumers"]) if not c["found"]]
    assert not missing, missing


def test_citation_checker_accepts_a_match_and_rejects_a_conflation():
    rec = {"first_author": "Ravi", "year": 2022, "title": "Spatially resolved ..."}
    assert inv.citation_ok(rec, "Ravi", 2022)[0] is True
    assert inv.citation_ok(rec, "Ruiz-Moreno", 2022)[0] is False          # the D25 error
    assert inv.citation_ok(rec, "Ravi", 2019)[0] is False
    assert inv.citation_ok({"error": "x"}, "Ravi", 2022)[0] is False
    assert inv.citation_ok({"first_author": "Mossi Albiach", "year": 2023},
                           "Mossi Albiach", 2023)[0] is True


@pytest.mark.skipif(CACHE is None, reason="run build_data_inventory.py --verify-dois first")
def test_every_registered_citation_matches_its_crossref_record():
    bad = []
    for e in inv.REGISTRY:
        pairs = ([(e["doi"], e["first_author"], e["year"])] if e["doi"] else []) + \
                list(e.get("extra_dois", []))
        for doi, fa, yr in pairs:
            ok, why = inv.citation_ok(CACHE["crossref"].get(doi), fa, yr)
            if not ok:
                bad.append((e["id"], doi, why))
    assert not bad, bad


@pytest.mark.skipif(CACHE is None, reason="run build_data_inventory.py --verify-dois first")
def test_the_known_conflation_is_rejected_on_real_crossref_data():
    for doi, fa, yr in inv.CONTROL_CONFLATIONS:
        assert inv.citation_ok(CACHE["crossref"].get(doi), fa, yr)[0] is False


def test_gbmap_identifiers_are_the_ones_in_the_file_itself():
    """D25: the inventory must carry the collection the atlas file names, not Siletti's."""
    gbmap = next(e for e in inv.REGISTRY if e["id"] == "D13")
    assert "999f2a15-3d7e-440b-96ae-2c806799c08c" in gbmap["accession"]
    assert "283d65eb" not in gbmap["accession"]
    h5ad = inv.REF / "gbmap_core.h5ad"
    if h5ad.exists():
        import h5py
        with h5py.File(h5ad, "r") as h:
            cit = h["uns"]["citation"][()]
        cit = cit.decode() if isinstance(cit, bytes) else str(cit)
        assert "999f2a15-3d7e-440b-96ae-2c806799c08c" in cit
        assert "861acfd8-25f0-418b-a445-aa96da232827" in gbmap["accession"]
        assert "861acfd8-25f0-418b-a445-aa96da232827" in cit


def test_generated_inventory_reports_no_open_problems():
    md = (ROOT / "docs" / "DATA_INVENTORY.md")
    if not md.exists():
        pytest.skip("inventory not generated")
    assert "**Open problems: none.**" in md.read_text()
