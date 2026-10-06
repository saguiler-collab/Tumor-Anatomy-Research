"""
The Figure 1F measurement (scripts/klemm_figure1f.py; prespecified/klemm_t_vs_b.md).

Positive control: the vector measurement reproduces two numbers the paper prints in its own text
(melanoma-BrM CD8+ T = 33.01% of CD45+; all-BrM lymphocytes = 46.23%). Negative control: the same
geometry read with a wrong calibration (gridlines taken as 50% apart instead of 25%) must FAIL those
controls -- a check that passes under any calibration measures nothing. The recorded reading must
recompute from the recorded per-population values. The paper PDF is gitignored (publisher copyright),
so the geometry tests skip where it is absent; the artefact tests skip where the artefact is absent.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import klemm_figure1f as kf  # noqa: E402

ART = ROOT / "results" / "klemm_t_vs_b.json"            # never the conftest-redirected RESULTS_DIR


@pytest.fixture(scope="module")
def geometry():
    if not kf.config.KLEMM_2020_PDF.exists():
        pytest.skip("Klemm et al. 2020 PDF not on this machine (gitignored)")
    return kf.shapes(kf.config.KLEMM_2020_PDF, kf.PAGE)


def test_measurement_reproduces_the_papers_printed_numbers(geometry):
    c = kf.controls(kf.measure(geometry))
    assert c["melanoma_CD8"]["pass"], c["melanoma_CD8"]
    assert c["brm_lymphocytes_weighted"]["pass"], c["brm_lymphocytes_weighted"]
    assert c["bar_totals_100"]["pass"], c["bar_totals_100"]


def test_a_wrong_calibration_fails_the_printed_number_controls(geometry):
    c = kf.controls(kf.measure(geometry, scale_divisor=50.0))
    assert not c["melanoma_CD8"]["pass"] and not c["bar_totals_100"]["pass"]


@pytest.mark.parametrize("t, b, expected", [(7.28, 0.35, "SUPPORTS T > B"), (0.35, 7.28, "CONTRADICTS T > B"),
                                            (1.000, 1.005, "UNRESOLVED"), (1.005, 1.000, "UNRESOLVED")])
def test_reading_rule(t, b, expected):
    row = {"CD4+ T": t, "Treg": 0.0, "CD8+ T": 0.0, "DNT": 0.0, "B": b, "NK": 0.5}
    assert kf.reading(row)["reading"] == expected


def test_recorded_reading_recomputes_from_recorded_values():
    if not ART.exists():
        pytest.skip("results/klemm_t_vs_b.json not produced yet")
    a = json.loads(ART.read_text())
    assert a["controls"]["all_pass"]
    for key, row in (("IDH_mutant_glioma_n17", "Glioma IDH mut"), ("IDH_wildtype_glioma_n40", "Glioma IDH wt")):
        assert kf.reading(a["rows"][row]) == a[key]
