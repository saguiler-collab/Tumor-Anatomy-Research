"""The core close-out checker must be able to fail (docs/CORE_CLOSEOUT.md).

A freeze check that passes a changed file, or a pytest check that accepts a failing or stale run, would certify
nothing. These are its negative controls.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import check_core_closeout as cc  # noqa: E402
import freeze_results as fr  # noqa: E402


def _files(root: Path) -> list[Path]:
    return fr.frozen_files(trees=(root / "a",), exclude=())


def test_freeze_passes_an_untouched_tree_and_flags_every_kind_of_change(tmp_path):
    (tmp_path / "a").mkdir()
    for name, text in (("a/x.csv", "1,2\n"), ("a/y.json", "{}"), ("a/.hidden", "not frozen"),
                       ("a/dwls_r_stdout.log", "an R console log, rewritten by every call: not frozen")):
        (tmp_path / name).write_text(text)
    manifest = tmp_path / "manifest.tsv"
    assert [f.name for f in _files(tmp_path)] == ["x.csv", "y.json"]
    fr.write_manifest(manifest, _files(tmp_path), tmp_path)

    clean = cc.compare_manifest(manifest, tmp_path, _files(tmp_path))
    assert (clean["recorded"], clean["changed"], clean["missing"], clean["added"]) == (2, [], [], [])

    (tmp_path / "a/x.csv").write_text("1,3\n")                 # same size, different bytes
    (tmp_path / "a/y.json").unlink()
    (tmp_path / "a/z.txt").write_text("new")
    res = cc.compare_manifest(manifest, tmp_path, _files(tmp_path))
    assert res["changed"] == ["a/x.csv"] and res["missing"] == ["a/y.json"] and res["added"] == ["a/z.txt"]


def test_pytest_evidence_rejects_a_failing_run_a_stale_run_and_an_unmarked_log(tmp_path):
    log, start = tmp_path / "pytest.log", time.time()
    log.write_text(f"START {start}\n....F\n1 failed, 470 passed in 250.0s\n")
    assert not cc.pytest_evidence(log, code_mtime=start - 60)[0]

    log.write_text(f"START {start}\n.....\n471 passed, 1 skipped in 250.0s (0:04:10)\n")
    assert cc.pytest_evidence(log, code_mtime=start - 60)[0]
    assert not cc.pytest_evidence(log, code_mtime=start + 60)[0]          # code changed after the run started

    log.write_text("471 passed in 250.0s\n")                                # no START line: not evidence
    assert not cc.pytest_evidence(log, code_mtime=0)[0]


def test_a_versioned_freeze_refuses_an_unexplained_change():
    previous = {"results/a.json": (10, "h1"), "results/b.csv": (5, "h2"), "results/gone.json": (3, "h3")}
    current = [("results/a.json", 10, "h1"), ("results/b.csv", 6, "h9"), ("results/new.csv", 7, "h4")]
    rows, unexplained = fr.version_changes(previous, current, {"results/b.csv": "recomputed (D39)"})
    assert {(p, k) for p, k, _ in rows} == {("results/b.csv", "changed"), ("results/new.csv", "added"),
                                            ("results/gone.json", "removed")}
    assert sorted(unexplained) == ["results/gone.json", "results/new.csv"]       # no reason given: refused
    rows, unexplained = fr.version_changes(previous, current, {"results/b.csv": "x", "results/new.csv": "y",
                                                                "results/gone.json": "z"})
    assert unexplained == [] and all(w != "NO REASON GIVEN" for _, _, w in rows)

