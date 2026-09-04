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
#: This is a disclosure, not a preference. DWLS's own signature builder
#: (buildSignatureMatrixMAST) is dominated by a condition-number search over gene
#: counts, and on this hardware it does not terminate predictably. Measured on the real
#: 11,739-cell export:
#:
#:     3,000 cells, cold cache      1,625 s   completed
#:     2,000 cells, per-type cap   >2,100 s   killed, no end in sight
#:     Seurat builder instead      >2,100 s   same — the DE step is not the cost
#:     inside the pipeline         >5,280 s   ran past a 2,400 s subprocess timeout
#:
#: The last line is the deciding one: a wall-clock budget was added precisely to bound
#: this, and it did not fire, so the R path cannot be given a predictable ceiling here.
#: A benchmark that never finishes yields nothing, which is strictly worse than a
#: reimplementation that is labelled as one everywhere it appears.
#:
#: This is not a claim about DWLS. It is a claim about running DWLS's signature build on
#: this machine, and it is recorded so a reader knows the leaderboard's `dwls` row is
#: this project's implementation of the algorithm rather than the published package.
R_PATH_DISABLED: dict[str, str] = {
    "dwls": (
        "the genuine DWLS package is installed and runs, but its signature build "
        "(buildSignatureMatrixMAST) has no predictable ceiling on this hardware — "
        "measured at 1,625 s in the best case and still running past a 2,400 s "
        "subprocess timeout in the worst. Run as this project's Python "
        "reimplementation, which is labelled as such in every artefact. Re-enable by "
        "removing this entry once the R path can be bounded."
    ),
}


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
