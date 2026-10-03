"""
The per-method AUC (scripts/anatomic_auc.py) and what keeps it distinct from ACS.

The two statistics share constraints, samples, tumour unit and null, and differ only in
being graded vs thresholded and in tie handling. These tests pin each of those claims,
plus the negative control a scorer must pass before it is believed: an inverted signal
must score below chance, not merely lower.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ivygap import config
from ivygap.anatomic import constraints as K
from ivygap.anatomic.acs import _acs_from_table, _score_table, tumor_structure_means

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("anatomic_auc", ROOT / "scripts" / "anatomic_auc.py")
aa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(aa)

STRUCTS = list(config.PRIMARY_STRUCTURES)
CELLS = list(config.CELL_TYPES)


def _cohort(per_struct: int, signal: float, seed: int = 0, n_tumours: int = 6):
    """Synthetic estimates with the registered orderings planted at strength `signal`
    (negative = inverted). Returns (estimate, manifest)."""
    rng = np.random.default_rng(seed)
    rows, meta = [], []
    rank = {"LE": 0, "IT": 1, "CT": 2, "MVP": 2, "PAN": 2}
    for t in range(n_tumours):
        for s in STRUCTS:
            for k in range(per_struct):
                v = {c: rng.uniform(0.05, 0.15) for c in CELLS}
                v["Tumor"] += signal * rank[s] * 0.1
                v["Oligodendrocyte"] += signal * (0.2 if s == "LE" else 0.0)
                v["Endothelial"] += signal * (0.2 if s == "MVP" else 0.0)
                v["Macrophage_Microglia"] += signal * (0.2 if s in ("PAN", "MVP") else 0.0)
                v = {c: max(x, 1e-6) for c, x in v.items()}
                tot = sum(v.values())
                sid = f"t{t}_{s}_{k}"
                rows.append(pd.Series({c: x / tot for c, x in v.items()}, name=sid))
                meta.append({"sample_id": sid, "patient_id": f"T{t}", "structure": s})
    est = pd.DataFrame(rows)[CELLS]
    man = pd.DataFrame(meta).set_index("sample_id")
    return est, man


def _arrays(est, man):
    return (man.loc[est.index, "patient_id"].astype(str).to_numpy(),
            man.loc[est.index, "structure"].astype(str).to_numpy())


def test_planted_signal_is_recovered():
    est, man = _cohort(per_struct=2, signal=1.0)
    t, s = _arrays(est, man)
    assert aa.weighted(aa.auc_table(est, t, s), "auc") > 0.9


def test_inverted_signal_is_rejected():
    """The negative control that matters: inversion scores BELOW chance."""
    est, man = _cohort(per_struct=2, signal=-1.0)
    t, s = _arrays(est, man)
    assert aa.weighted(aa.auc_table(est, t, s), "auc") < 0.2


def test_no_signal_sits_at_chance_under_the_null():
    est, man = _cohort(per_struct=2, signal=0.0, seed=1)
    t, s = _arrays(est, man)
    null = aa.auc_null(est, t, s, n_permutations=500, seed=0)
    assert abs(null.mean() - 0.5) < 0.02


def test_fast_and_batched_paths_equal_the_readable_definition():
    est, man = _cohort(per_struct=3, signal=0.3, seed=2)
    t, s = _arrays(est, man)
    fast = aa.FastAUC(est, t)
    assert fast.score(s) == pytest.approx(aa.weighted(aa.auc_table(est, t, s), "auc"), abs=1e-12)
    rng = np.random.default_rng(5)
    base = [np.array([fast.code[x] for x in s[idx]]) for idx in fast.idx]
    codes = [rng.permuted(np.tile(c, (25, 1)), axis=1) for c in base]
    batch = fast.score_batch(codes)
    for b in range(25):
        relab = s.copy()
        for k, idx in enumerate(fast.idx):
            relab[idx] = np.array(fast.structures)[codes[k][b]]
        assert batch[b] == pytest.approx(aa.weighted(aa.auc_table(est, t, relab), "auc"),
                                         abs=1e-12)


def _joined(est, man):
    t, s = _arrays(est, man)
    tab = aa.auc_table(est, t, s).set_index(["constraint", "tumor_id"])
    acs = _score_table(tumor_structure_means(est, man)).set_index(["constraint", "tumor_id"])
    return tab.join(acs["satisfied"], how="outer")


def test_both_statistics_evaluate_the_same_constraint_tumour_pairs():
    j = _joined(*_cohort(per_struct=1, signal=0.2, seed=3))
    assert j["auc"].notna().all() and j["satisfied"].notna().all()


def test_with_one_sample_per_structure_pairwise_constraints_are_acs():
    """One sample per (tumour, structure), no ties: each per-tumour AUC of a PAIRWISE
    constraint is the same 0/1 comparison ACS makes."""
    j = _joined(*_cohort(per_struct=1, signal=0.2, seed=3))
    pw = j.loc[[c.id for c in K.CONSTRAINTS if c.kind == "pairwise"]]
    assert pw["tied"].max() == 0
    assert (pw["auc_strict"] == pw["satisfied"]).all()


def test_conjunction_constraints_get_partial_credit_only_from_the_auc():
    """C4 and C7 are where the identity breaks, by design: ACS demands every comparison,
    the AUC averages them. A method placing MVP above three of four structures scores 0
    on ACS's C4 and 0.75 on the AUC's."""
    est, man = _cohort(per_struct=1, signal=0.0, seed=7, n_tumours=1)
    level = {"LE": 0.10, "IT": 0.30, "CT": 0.20, "MVP": 0.25, "PAN": 0.05}
    for sid in est.index:
        est.loc[sid, "Endothelial"] = level[man.loc[sid, "structure"]]
        est.loc[sid, "Tumor"] = {"LE": 0.1, "IT": 0.5, "CT": 0.4}.get(man.loc[sid, "structure"], 0.3)
    j = _joined(est, man)
    assert j.loc[("C4", "T0"), "satisfied"] == 0.0
    assert j.loc[("C4", "T0"), "auc_strict"] == pytest.approx(0.75)
    assert j.loc[("C7", "T0"), "satisfied"] == 0.0            # LE < IT holds, IT < CT does not
    assert j.loc[("C7", "T0"), "auc_strict"] == pytest.approx(0.5)


def test_ties_are_where_the_two_conventions_part():
    """ACS scores an exact tie 0; an AUC scores it 0.5. A method that returns identical
    values in two structures must not get the same credit under both."""
    est, man = _cohort(per_struct=1, signal=0.0, seed=4)
    est[:] = 1.0 / len(CELLS)                       # every structure identical
    t, s = _arrays(est, man)
    tab = aa.auc_table(est, t, s)
    assert aa.weighted(tab, "auc") == pytest.approx(0.5)
    assert aa.weighted(tab, "auc_strict") == pytest.approx(0.0)
    assert _acs_from_table(_score_table(tumor_structure_means(est, man))) == pytest.approx(0.0)


def test_pairs_are_never_formed_across_tumours():
    """Simpson's paradox, built on purpose: each tumour has CT below LE, but the tumour
    sampled mostly at CT sits higher overall, so POOLED samples would show CT above LE.
    Within-tumour pairing must report the inversion; pooling would hide it."""
    rows, meta = [], []
    plan = {"T0": {"LE": (4, 0.30), "CT": (1, 0.20)},     # inverted, low tumour
            "T1": {"LE": (1, 0.90), "CT": (4, 0.80)}}     # inverted, high tumour
    for tum, structs in plan.items():
        for s, (n, level) in structs.items():
            for k in range(n):
                sid = f"{tum}_{s}_{k}"
                rows.append(pd.Series({c: 0.1 for c in CELLS} | {"Tumor": level}, name=sid))
                meta.append({"sample_id": sid, "patient_id": tum, "structure": s})
    est, man = pd.DataFrame(rows)[CELLS], pd.DataFrame(meta).set_index("sample_id")
    t, s = _arrays(est, man)
    c1 = aa.auc_table(est, t, s).query("constraint == 'C1'")
    assert (c1["auc"] == 0.0).all() and len(c1) == 2
    pooled = aa.pair_auc(est.loc[s == "CT", "Tumor"].to_numpy(),
                         est.loc[s == "LE", "Tumor"].to_numpy())[0]
    assert pooled > 0.5          # what pooling across tumours would have reported


def test_harrell_c():
    x = np.arange(10, dtype=float)
    assert aa.harrell_c(x, x)[0] == 1.0
    assert aa.harrell_c(-x, x)[0] == 0.0
    assert aa.harrell_c(np.zeros(10), x)[0] == 0.5


def test_constraint_pairs_cover_the_registered_file():
    got = {c.id: aa.constraint_pairs(c) for c in K.CONSTRAINTS}
    assert got["C1"] == [("CT", "LE")]
    assert got["C7"] == [("IT", "LE"), ("CT", "IT")]
    assert len(got["C4"]) == 4 and all(h == "MVP" for h, _ in got["C4"])


ARTEFACT = config.RESULTS_DIR / "anatomic_auc.json"


@pytest.mark.skipif(not ARTEFACT.exists(), reason="run scripts/anatomic_auc.py first")
def test_artefact_controls_and_the_conclusion_it_supports():
    d = json.loads(ARTEFACT.read_text())
    assert d["controls"]["acs_reproduced_for_every_method"] is True
    assert d["controls"]["registered_agreement_reproduced"]["rho"] == pytest.approx(0.081, abs=1e-3)
    assert d["controls"]["margin_worst_real_minus_best_control"] > 0
    for m, r in d["methods"].items():
        if r["is_control"]:
            assert r["null_p"] > 0.05, f"{m}: a broken-input control beat its null"
    # the conclusion this robustness check supports: a graded score does not rescue
    # anatomy as a predictor of accuracy either
    v = d["agreement_with_auc_in_place_of_acs"]["vs_rho_purity"]
    assert v["meets_threshold"] is False and v["ci_excludes_zero"] is False
    assert d["rank_agreement"]["acs_vs_auc_anat_spearman_comparable_methods"] > 0.8
