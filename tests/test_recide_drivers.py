"""
ReCIDE's drivers: the anatomic run (scripts/run_recide.py) and the TCGA run (scripts/recide_tcga.py)
share one donor draw, one reference build and one R call. These checks pin what must not drift
silently: the implementation label discloses the shim, the draw rule is the one fixed on
2026-10-01, and the TCGA adapter reports itself under the same label. The draw's determinism is a
runtime control in recide_tcga.py (it must reproduce the anatomic run's donors), since it needs the
8 GB atlas.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import recide_tcga  # noqa: E402
import run_recide  # noqa: E402


def test_the_label_discloses_the_vendored_commit_and_the_shim():
    assert "ReCIDE vendored @31bbd6b" in run_recide.IMPLEMENTATION
    assert "shim" in run_recide.IMPLEMENTATION           # never reported as the unmodified package


def test_the_draw_rule_is_the_one_fixed_before_any_result():
    assert "(10,15,20)" in run_recide.DRAW_RULE and ">500" in run_recide.DRAW_RULE
    assert ">20 cells" in run_recide.DRAW_RULE and "project seed" in run_recide.DRAW_RULE


def test_the_tcga_adapter_reports_under_the_same_name_and_label():
    a = recide_tcga.ReCIDEAdapter("gbm", ["d1", "d2"], n_cores=1, budget=None)
    assert a.name == "recide"
    assert a.implementation_ == run_recide.IMPLEMENTATION
    assert recide_tcga.ARM_MATRIX == "raw/X"
