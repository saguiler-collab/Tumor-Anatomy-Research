"""compare_artefacts: reproduced, numeric drift, missing keys, and ignored run metadata."""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("ca", ROOT / "scripts" / "compare_artefacts.py")
ca = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ca)


def test_identical_documents_reproduce():
    d = {"a": 1.0, "b": {"c": [1, 2, "x"]}}
    assert ca.compare(d, d, 1e-9, ca.DEFAULT_IGNORE) == []


def test_a_numeric_change_is_reported():
    """The negative control: a changed result must not pass as reproduced."""
    diffs = ca.compare({"acs": 0.70}, {"acs": 0.71}, 1e-9, ca.DEFAULT_IGNORE)
    assert len(diffs) == 1 and "acs" in diffs[0]


def test_missing_and_extra_keys_are_reported():
    diffs = ca.compare({"a": 1, "b": 2}, {"a": 1, "c": 3}, 1e-9, ca.DEFAULT_IGNORE)
    assert any("missing" in d for d in diffs) and any("extra" in d for d in diffs)


def test_run_metadata_is_ignored():
    assert ca.compare({"elapsed_seconds": 3.0, "out_path": "/a"}, {"elapsed_seconds": 9.0,
                      "out_path": "/b"}, 1e-9, ca.DEFAULT_IGNORE) == []
