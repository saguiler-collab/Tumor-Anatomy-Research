"""
The verification harness must be able to FAIL: identical artefacts are REPRODUCED, a planted change
is DIFFERS, sub-1e-6 jitter is NUMERICAL NOISE, and ignored keys (timestamps) never count.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import verify_rerun as vr  # noqa: E402


def _j(p: Path, d: dict) -> Path:
    p.write_text(json.dumps(d)); return p


def test_json_identical_is_reproduced(tmp_path):
    d = {"rho": 0.081, "n": 12, "nested": {"x": [1.0, 2.0]}}
    assert vr.compare_json(_j(tmp_path / "a.json", d), _j(tmp_path / "b.json", d))["status"] == "REPRODUCED"


def test_json_planted_change_is_detected(tmp_path):
    a = _j(tmp_path / "a.json", {"rho": 0.081, "verdict": "INCONCLUSIVE"})
    b = _j(tmp_path / "b.json", {"rho": 0.091, "verdict": "INCONCLUSIVE"})
    r = vr.compare_json(a, b)
    assert r["status"] == "DIFFERS" and r["n_differences"] == 1
    c = _j(tmp_path / "c.json", {"rho": 0.081, "verdict": "SUPPORTED"})
    assert vr.compare_json(a, c)["status"] == "DIFFERS"


def test_json_jitter_is_noise_and_timestamps_are_ignored(tmp_path):
    a = _j(tmp_path / "a.json", {"rho": 0.5, "written_utc": "2026-10-06T00:00:00"})
    b = _j(tmp_path / "b.json", {"rho": 0.5 + 1e-8, "written_utc": "2026-10-07T00:00:00"})
    assert vr.compare_json(a, b)["status"] == "NUMERICAL NOISE"


def test_csv_planted_change_is_detected(tmp_path):
    df = pd.DataFrame({"sample": ["s1", "s2"], "Tumor": [0.9, 0.8], "B_cell": [0.01, 0.02]})
    df.to_csv(tmp_path / "a.csv", index=False)
    df.to_csv(tmp_path / "b.csv", index=False)
    assert vr.compare_csv(tmp_path / "a.csv", tmp_path / "b.csv")["status"] == "REPRODUCED"
    df.loc[1, "B_cell"] = 0.03
    df.to_csv(tmp_path / "c.csv", index=False)
    r = vr.compare_csv(tmp_path / "a.csv", tmp_path / "c.csv")
    assert r["status"] == "DIFFERS" and r["n_cells_over_1e-6"] == 1


def test_every_step_is_well_formed():
    ids = [s["id"] for s in vr.STEPS]
    assert len(ids) == len(set(ids))
    for s in vr.STEPS:
        assert s["tier"] in vr.TIERS and (s["outputs"] or s["verdict"])
        assert Path(s["cmd"][1]).exists(), s["cmd"][1]


def test_added_fields_are_named_but_a_lost_or_changed_value_still_differs(tmp_path):
    a = _j(tmp_path / "a.json", {"rho": 0.081, "n": 12})
    added = _j(tmp_path / "b.json", {"rho": 0.081, "n": 12, "seed": 0})
    r = vr.compare_json(a, added)
    assert r["status"] == "REPRODUCED + NEW FIELDS" and r["n_differences"] == 1
    # an added field never hides a changed value or a field the re-run lost
    changed = _j(tmp_path / "c.json", {"rho": 0.091, "n": 12, "seed": 0})
    lost = _j(tmp_path / "d.json", {"rho": 0.081, "seed": 0})
    assert vr.compare_json(a, changed)["status"] == "DIFFERS"
    assert vr.compare_json(a, lost)["status"] == "DIFFERS"
    # a record from before the category existed is relabelled only if every difference was listed
    full = {"status": "DIFFERS", "n_differences": 1, "examples": ["extra in re-run:   seed = 0"]}
    cut = {"status": "DIFFERS", "n_differences": 30, "examples": ["extra in re-run:   seed = 0"] * 12}
    assert vr._reclassify(full)["status"] == "REPRODUCED + NEW FIELDS"
    assert vr._reclassify(cut)["status"] == "DIFFERS"


def test_a_differing_output_is_moved_aside_and_the_registered_one_restored(tmp_path, monkeypatch):
    snap, res = tmp_path / "snap", tmp_path / "res"
    for d in (snap, res):
        d.mkdir()
    (snap / "x.json").write_text(json.dumps({"rho": 0.081}))
    (res / "x.json").write_text(json.dumps({"rho": 0.5}))
    (snap / "y.json").write_text(json.dumps({"rho": 0.3}))
    (res / "y.json").write_text(json.dumps({"rho": 0.3}))
    monkeypatch.setattr(vr, "RES", res)
    monkeypatch.setattr(vr, "VDIR", tmp_path / "verification")
    st = {"id": "step", "outputs": ["x.json", "y.json"]}
    rec = {"outputs": {o: vr.compare_output(o, snap) for o in st["outputs"]}}
    assert rec["outputs"]["x.json"]["status"] == "DIFFERS" and rec["outputs"]["y.json"]["status"] == "REPRODUCED"
    assert vr.keep_registered(st, rec, snap) == ["x.json"]
    assert json.loads((res / "x.json").read_text()) == {"rho": 0.081}                 # registered back in place
    assert json.loads((tmp_path / "verification" / "rerun_outputs" / "step" / "x.json").read_text()) == {"rho": 0.5}


def test_csv_with_added_rows_is_named_but_a_changed_or_lost_row_still_differs(tmp_path):
    a = pd.DataFrame({"sample_id": ["s1", "s2"], "method": ["nnls", "nnls"], "Tumor": [0.5, 0.7]})
    added = pd.DataFrame({"sample": ["s1", "s2", "s1"], "method": ["nnls", "nnls", "bh"], "Tumor": [0.5, 0.7, 0.9]})
    a.to_csv(tmp_path / "a.csv", index=False)
    added.to_csv(tmp_path / "b.csv", index=False)
    r = vr.compare_csv(tmp_path / "a.csv", tmp_path / "b.csv")
    assert r["status"] == "REPRODUCED + NEW ROWS" and r["new_rows"] == 1 and "bh" in r["new_row_keys"]
    changed = added.copy(); changed.loc[0, "Tumor"] = 0.51
    changed.to_csv(tmp_path / "c.csv", index=False)
    assert vr.compare_csv(tmp_path / "a.csv", tmp_path / "c.csv")["status"] == "DIFFERS"
    lost = added.iloc[1:]
    lost.to_csv(tmp_path / "d.csv", index=False)
    assert vr.compare_csv(tmp_path / "a.csv", tmp_path / "d.csv")["status"] == "DIFFERS"


def test_steps_run_on_the_registered_code_path_and_importing_changes_nothing(tmp_path, monkeypatch):
    # The registered artefacts were made with the repairs off. A step must see IVYGAP_REPAIRED=0 unless the caller set
    # it, and importing the harness must not change the importer's environment (the first version did, and switched
    # the repairs off for every test that ran after it).
    import os
    import subprocess
    assert vr.STEP_ENV == {"IVYGAP_REPAIRED": os.environ.get("IVYGAP_REPAIRED", "0")}
    probe = [sys.executable, "-c", "import os; print('SEEN', os.environ.get('IVYGAP_REPAIRED'))"]
    monkeypatch.setattr(vr, "VDIR", tmp_path / "verification")
    monkeypatch.setattr(vr.config, "PROJECT_ROOT", tmp_path)          # run_step records the log relative to it
    monkeypatch.setattr(vr, "keep_registered", lambda st, rec, snap: [])
    monkeypatch.setattr(vr, "STEP_ENV", {"IVYGAP_REPAIRED": "0"})
    monkeypatch.delenv("IVYGAP_REPAIRED", raising=False)
    rec = vr.run_step({"id": "probe", "cmd": probe, "outputs": [], "verdict": "exit0", "tier": "fast"}, tmp_path)
    assert rec["verdict"] == "PASSED"
    assert "SEEN 0" in (tmp_path / "verification" / "logs" / "probe.log").read_text()
    assert "IVYGAP_REPAIRED" not in os.environ                       # the parent's environment is untouched
    assert subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0, 'scripts'); import verify_rerun, os; "
                           "print(os.environ.get('IVYGAP_REPAIRED'))"], capture_output=True, text=True,
                          cwd=Path(__file__).resolve().parent.parent,
                          env={k: v for k, v in os.environ.items()}).stdout.strip() == "None"


def test_shared_blank_cells_reproduce_but_a_moved_blank_still_differs(tmp_path):
    # A table with the same NaN cells on both sides is identical. It read DIFFERS because the relative scale of a
    # NaN cell was NaN. A blank that moves must still differ.
    a = tmp_path / "a.csv"; b = tmp_path / "b.csv"; c = tmp_path / "c.csv"
    a.write_text("case,x,y\nC1,0.5,\nC2,0.25,1.0\n")
    b.write_text("case,x,y\nC1,0.5,\nC2,0.25,1.0\n")
    c.write_text("case,x,y\nC1,0.5,1.0\nC2,0.25,\n")
    assert vr.compare_csv(a, b)["status"] == "REPRODUCED"
    assert vr.compare_csv(a, c)["status"] == "DIFFERS"

