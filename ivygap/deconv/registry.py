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

#: Published tools whose genuine R package is not run here, with the reason.
#:
#: Empty, and the history is worth keeping. DWLS was disabled here after it ran past a
#: 2,400 s budget that existed to bound it. The budget was not the problem: `subprocess.run`
#: kills the direct child on timeout and then calls `communicate()` again with NO timeout
#: to reap output, which never returns if the child left grandchildren holding the
#: inherited pipe. `_run_bounded` puts the child in its own process group and kills the
#: group, and on the real DWLS call the timeout now fires where it previously did not.
#:
#: A synthetic reproduction of the pipe-holding child timed out correctly under BOTH
#: paths, so the mechanism above is inferred from the real call rather than demonstrated
#: in isolation. What is demonstrated: the same DWLS invocation that ran unbounded under
#: subprocess.run stops at its budget under _run_bounded.
R_PATH_DISABLED: dict[str, str] = {}


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
        if prefer_r and r_name in R_PATH_DISABLED:
            # Deliberately not wrapped: see R_PATH_DISABLED. The reason travels on the
            # object so the disclosure reports it exactly as it reports a fallback.
            py_impl.r_path_disabled_reason_ = R_PATH_DISABLED[r_name]
            methods.append(py_impl)
            continue
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
