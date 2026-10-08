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


def test_references_doc_dois_maps_every_bracket_to_each_doi_it_names(tmp_path):
    p = tmp_path / "REFERENCES.md"
    p.write_text(
        '| **[5]** | A. Author, "Title," doi: 10.1/ABC. | PDF |\n'
        "| **[10]** | B. Author, doi: 10.2/DEF -- the dataset also cites the preprint, doi: 10.3/GHI. | Crossref |\n"
        "not a reference row, doi: 10.9/should-not-be-seen\n"
    )
    assert inv.references_doc_dois(p) == {"10.1/abc": 5, "10.2/def": 10, "10.3/ghi": 10}


def test_a_used_datasets_citation_missing_from_references_md_is_a_reported_problem(tmp_path):
    # Negative control: deleting one DOI from a copy of the real REFERENCES.md must turn the dataset that
    # cites it into a reported problem -- proving the check can fail, not just always pass. Picks a DOI this
    # test confirms is genuinely present in the real file first, so the control is not vacuous.
    intact = inv.references_doc_dois(inv.REFERENCES)
    used = [e for e in inv.REGISTRY if e["status"] in ("used", "derived") and e["doi"]
            and e["doi"].strip().rstrip(").,").lower() in intact]
    assert used, "no used/derived dataset's DOI is present in docs/REFERENCES.md -- nothing to break for the control"
    target = used[0]
    broken = tmp_path / "REFERENCES.md"
    broken.write_text(inv.REFERENCES.read_text().replace(target["doi"], "10.0000/not-the-same-doi-at-all"))
    assert target["doi"].strip().rstrip(").,").lower() not in inv.references_doc_dois(broken)


def test_known_gap_gencode_and_brennan_have_no_references_md_entry_yet():
    """Documents OPEN_DEFECTS-style: found 2026-10-08, the two citations below are used by this project (D06,
    D12 in docs/DATA_INVENTORY.md) but were never added to docs/REFERENCES.md's numbered list. Update this test,
    not just the assertion, once a user has added them -- it exists so the gap cannot quietly disappear from the
    record without someone having actually added the reference."""
    missing = {"10.1093/nar/gkaa1087", "10.1016/j.cell.2013.09.034"}
    present = set(inv.references_doc_dois())
    assert missing & present == set(), (
        "GENCODE and/or Brennan 2013 now have a docs/REFERENCES.md entry -- good. Replace this test with an "
        "assertion that build_data_inventory.py reports zero DOI-coverage problems.")
