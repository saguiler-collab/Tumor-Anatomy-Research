"""
registry.py — the roster of methods the benchmark compares.

Ordering is the order results are reported in. Baselines come first so a reader sees
what the newer methods are being compared *against* before seeing their scores.
"""

from __future__ import annotations

from .base import DeconvolutionMethod
from .bayesian import BayesianDeconvolution, HierarchicalBayesianDeconvolution
from .classical import NNLSDeconvolution, SVRDeconvolution
from .elastic_net import ElasticNetDeconvolution
from .r_bridge import RMethod, check as r_check
from .reference_based import (BisqueDeconvolution, DWLSDeconvolution,
                              MuSiCDeconvolution, SCDCDeconvolution,
                              SCDCEnsembleDeconvolution)

#: the four published tools, paired with the Python reimplementation used as fallback
PUBLISHED_TOOLS = ("music", "dwls", "bisque", "scdc", "scdc_ensemble")


def build_methods(prefer_r: bool = True, allow_r_fallback: bool = True
                  ) -> list[DeconvolutionMethod]:
    """
    Construct every method for a run.

    `prefer_r=True` wraps the four published tools so the genuine R package runs when
    it can, falling back to the Python reimplementation otherwise. Set it False to
    force the all-Python path — which is what the test suite and any
    reproduce-without-R workflow should do, so results do not silently depend on
    whether the machine happened to have R.
    """
    methods: list[DeconvolutionMethod] = [
        NNLSDeconvolution(),
        SVRDeconvolution(),
        ElasticNetDeconvolution(),
        BayesianDeconvolution(),
        HierarchicalBayesianDeconvolution(),
    ]

    published = [
        ("music", MuSiCDeconvolution()),
        ("dwls", DWLSDeconvolution()),
        ("bisque", BisqueDeconvolution()),
        ("scdc", SCDCDeconvolution()),
        ("scdc_ensemble", SCDCEnsembleDeconvolution()),
    ]
    for r_name, py_impl in published:
        methods.append(
            RMethod(r_name, py_impl, allow_fallback=allow_r_fallback)
            if prefer_r else py_impl
        )
    return methods


def implementation_report(reference_name: str = "gbmap") -> list[dict]:
    """
    Which of the four published tools will run as the genuine R package, and why not
    when not. Written into the results manifest so the distinction between "MuSiC" and
    "our implementation of MuSiC" is never lost between the run and the write-up.
    """
    return [r_check(m, reference_name).as_dict() for m in PUBLISHED_TOOLS]


def method_names(prefer_r: bool = False) -> list[str]:
    return [m.name for m in build_methods(prefer_r=prefer_r)]
