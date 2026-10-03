"""
The POST-REGISTRATION EXTENSION PANEL.

These methods were added on 2026-09-30, after the OSF registration (dm2t8, 2026-09-10) and
after the confirmatory run, to widen the panel toward the top-ranked methods of Nguyen et al.
2024 (Briefings in Bioinformatics, 53 methods). They are:

  * run and reported under their own names, in their own artefacts (`results/extension/`);
  * NEVER added to `registry.build_methods`, which defines the registered panel, and never
    counted in any registered statistic (the ACS-vs-accuracy agreement test, the control
    separation, the registered T > B count);
  * genuine published packages only. An extension method has NO Python reimplementation and
    `allow_fallback=False`, so a failure is reported as a failure -- nothing can ever stand
    in under the package's name.

Why a separate panel rather than a longer leaderboard: a method added after the data were
seen can be reported, but it cannot be allowed to change a statistic that was registered
before they were. Labelling the panel is what keeps the confirmatory result confirmatory.
"""
from __future__ import annotations

import contextlib
import functools
import tempfile
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

from ivygap.deconv.base import DeconvolutionInput, DeconvolutionMethod
from ivygap.deconv.r_bridge import RMethod

ADDED = "2026-09-30"


@dataclass(frozen=True)
class ExtensionSpec:
    name: str
    package: str
    route: str                # "R" or "python"
    citation: str
    doi: str                  # verified against Crossref 2026-10-01; the machine-checkable part
    why: str                  # why THIS method, stated from the literature, not from a result


EXTENSION_SPECS: dict[str, ExtensionSpec] = {
    "fardeep": ExtensionSpec(
        name="fardeep", package="FARDEEP 1.0.1 (CRAN)", route="R",
        citation="Hao Y, Yan M, Heath BR, Lei YL, Xie Y. PLoS Comput Biol 2019;15:e1006976",
        doi="10.1371/journal.pcbi.1006976",
        why="Top tier of reference-based methods in Nguyen et al. 2024; robust "
            "least-trimmed-squares regression that drops outlier genes per sample, a "
            "mechanism no method in the registered panel has."),
    "deseq2_unmix": ExtensionSpec(
        name="deseq2_unmix", package="DESeq2 1.52.0 (Bioconductor), function unmix", route="R",
        citation="Love MI, Huber W, Anders S. Genome Biol 2014;15:550",
        doi="10.1186/s13059-014-0550-8",
        why="One of the seven methods Nguyen et al. 2024 rank in the top 10 in all three "
            "scenarios, and one of the two of those seven not yet run here; fits proportions "
            "with an L1 loss in a variance-stabilised space, a loss no other panel method uses. "
            "Configuration fixed before output: prespecified/deseq2_unmix_config.md."),
    "mixture": ExtensionSpec(
        name="mixture", package="MIXTURE 0.0.1 (GitHub elmerfer/MIXTURE, vendored @6332e160)",
        route="R",
        citation="Fernandez EA, Mahmoud YD, Veigas F, Rocha D, et al. Brief Bioinform 2021;22:bbaa317",
        doi="10.1093/bib/bbaa317",
        why="One of the seven methods Nguyen et al. 2024 rank in the top 10 in all three "
            "scenarios, and the last of those seven not yet run here; nu-SVR with recursive "
            "elimination of cell types under a noise bound. Configuration fixed before output: "
            "prespecified/mixture_config.md."),
    "lindeconseq": ExtensionSpec(
        name="lindeconseq", package="LinDeconSeq (GitHub lihuamei/LinDeconSeq, pinned commit)",
        route="R",
        citation="Li H, Sharma A, Ming W, Sun X, Liu H. BMC Genomics 2020;21:652",
        doi="10.1186/s12864-020-06888-1",
        why="Top tier of reference-based methods in Nguyen et al. 2024; weighted robust "
            "regression against a signature matrix."),
    "recide": ExtensionSpec(
        name="recide", package="ReCIDE (GitHub TianLab-Bioinfo/ReCIDE, vendored @31bbd6b)",
        route="vendored",
        citation="Li M, Su Y, Gao Y, Tian W. Brief Bioinform 2024;25:bbae422",
        doi="10.1093/bib/bbae422",
        why="Best performer in Li et al. 2026 (Genome Biology), a real-bulk benchmark; "
            "integrates per-donor signatures. Run by scripts/run_recide.py, not through the "
            "wrapper, because it needs uncapped per-donor cells. Note Li 2026 is the same lab."),
    "aric": ExtensionSpec(
        name="aric", package="ARIC (PyPI)", route="python",
        citation="Zhang W, Xu H, Qiao R, Zhong B, Zhang X, et al. Brief Bioinform "
                 "2022;23:bbab362 (online 2021)",
        doi="10.1093/bib/bbab362",
        why="Top tier of reference-based methods in Nguyen et al. 2024; weighted nu-SVR "
            "with iterative outlier-gene removal and collinearity reduction -- an "
            "outlier-robust variant of the CIBERSORT core the registered panel carries."),
    "rnasieve": ExtensionSpec(
        name="rnasieve", package="rnasieve 0.1.4 (PyPI)", route="python",
        citation="Erdmann-Pham DD, Fischer J, Hong J, Song YS. Genome Res 2021;31:1794-1806",
        doi="10.1101/gr.272344.120",
        why="Top tier of reference-based methods in Nguyen et al. 2024; a likelihood model "
            "that accounts for the single-cell reference's own sampling variance and "
            "estimates the bulk-versus-single-cell profile shift jointly with the "
            "proportions -- the cross-platform gap this cohort has."),
}


class _NoReimplementation(DeconvolutionMethod):
    """Placeholder that refuses to run. An extension method is the genuine package or nothing."""

    family = "reference-based"

    def __init__(self, name: str):
        super().__init__()
        self.name = name

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        raise RuntimeError(f"{self.name} has no reimplementation by design; the genuine "
                           f"package failed and that failure is the result")


def counts_are_integral(counts):
    """True if a matrix recovered as CPM x library_size / 1e6 is on a COUNT scale.

    Relative, not absolute: the cell source is float32, whose representable spacing grows with
    the value (a count near 1e5 is only exact to ~0.01), so an absolute 1e-3 tolerance rejected
    genuine counts (max deviation 0.00614 on Ivy GAP, 2026-10-01). A log or arbitrary scale is off
    by tens of percent, which a 1e-5 relative tolerance still catches.
    """
    import numpy as _np
    c = _np.asarray(counts, dtype="float64")
    nz = c[c > 0]
    if nz.size == 0:
        return True, 0.0
    rel = float((_np.abs(nz - _np.round(nz)) / _np.maximum(nz, 1.0)).max())
    return rel < 1e-5, rel


@contextlib.contextmanager
def _sieve_compat():
    """Scoped shim for rnasieve 0.1.4 on SciPy >= 1.11; see `_PythonPackageMethod._rnasieve`."""
    import scipy.optimize                                           # noqa: PLC0415
    import rnasieve.algo                                            # noqa: PLC0415
    import rnasieve.model                                           # noqa: PLC0415
    orig_min, orig_fm = scipy.optimize.minimize, rnasieve.model.find_mixtures
    scipy.optimize.minimize = lambda fun, x0, *a, **k: orig_min(fun, np.ravel(x0), *a, **k)
    rnasieve.model.find_mixtures = functools.partial(rnasieve.algo.find_mixtures,
                                                     parallelized=False)
    try:
        yield
    finally:
        scipy.optimize.minimize, rnasieve.model.find_mixtures = orig_min, orig_fm


class _PythonPackageMethod(DeconvolutionMethod):
    """A GENUINE published Python package, called through its own public API.

    Not a reimplementation: the package's code runs unmodified. `implementation_` records the
    installed distribution and version so the artefact can never be mistaken for a port.
    """

    family = "reference-based"
    requires_r = False

    def __init__(self, name: str, dist: str):
        super().__init__()
        self.name = name
        self.dist = dist
        self.implementation_ = f"py:{dist} {version(dist)}" + (
            " (+ disclosed SciPy x0-flatten shim; serial path)" if name == "rnasieve" else "")
        self.fallback_reason_ = None
        self.degenerate_, self.degeneracy_reason_ = False, None

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        if self.name == "aric":
            return self._aric(data)
        if self.name == "rnasieve":
            return self._rnasieve(data)
        raise RuntimeError(f"no driver for {self.name}")

    def _rnasieve(self, data: DeconvolutionInput) -> np.ndarray:
        """RNA-Sieve on raw counts recovered exactly from the registered cell source.

        DISCLOSED COMPATIBILITY SHIM. rnasieve 0.1.4 hands `scipy.optimize.minimize` a 2-D x0.
        SciPy < 1.11 flattened it silently; current SciPy raises. No SciPy that old exists for
        this interpreter, so `_sieve_compat` flattens x0 -- exactly what old SciPy did -- for
        the duration of the call and restores the original afterwards. It also routes
        `predict` to the package's own SERIAL path (`parallelized=False`): the default spawns
        10 worker processes, which a 2-core 8 GB machine cannot hold, and spawned workers
        would not inherit the shim. The serial and parallel paths solve the same independent
        per-gene problems, so this changes speed only. Planted-truth gate, 2026-09-30: max
        |est - truth| 0.011, per-type Pearson 0.999-1.000.
        """
        from rnasieve.preprocessing import model_from_raw_counts   # noqa: PLC0415
        from ivygap.deconv import r_bridge                          # noqa: PLC0415
        src = r_bridge.get_cell_source(data.primary.name)
        if src is None:
            raise RuntimeError("rnasieve needs the cell-level reference; none is registered "
                               "(the frozen signature carries no cells). Not run, not imputed.")
        expr, meta = src
        genes = [g for g in data.bulk.index if g in expr.index]
        # Exact counts: the builder stores per-cell CPM over the restricted gene set and the
        # raw library size it divided by.
        counts = expr.loc[genes].to_numpy(dtype="float64") * (
            meta.loc[expr.columns, "library_size"].to_numpy(dtype="float64") / 1e6)
        ok, rel = counts_are_integral(counts)
        if not ok:
            raise RuntimeError(f"recovered counts are not integers (max relative deviation {rel:.3g}); "
                               "the cell source is not on a counts scale")
        counts = np.round(counts)
        ctype = meta.loc[expr.columns, "cell_type"].astype(str).to_numpy()
        cts = list(data.cell_types)
        raw = {c: counts[:, ctype == c] for c in cts}
        empty = [c for c, m in raw.items() if m.shape[1] == 0]
        if empty:
            raise RuntimeError(f"no reference cells for {empty}; rnasieve cannot model them")
        model, psis = model_from_raw_counts(raw, data.bulk.loc[genes].to_numpy(dtype="float64"))
        with _sieve_compat():
            est = model.predict(psis)
        # psis columns are "Bulk 0..n" in input order; alpha_hats columns are sorted labels.
        est.index = list(data.bulk.columns)
        return est[cts].to_numpy(dtype="float64")

    def _aric(self, data: DeconvolutionInput) -> np.ndarray:
        from ARIC import ARIC        # noqa: PLC0415 -- optional dependency, extension only
        ref = data.primary.profile
        genes = [g for g in data.bulk.index if g in ref.index]
        cts = list(data.cell_types)
        with tempfile.TemporaryDirectory(prefix="ivygap_aric_") as tmp:
            tmp = Path(tmp)
            data.bulk.loc[genes].to_csv(tmp / "mix.csv")
            ref.loc[genes, cts].to_csv(tmp / "ref.csv")
            # Package defaults throughout (scale 0.1, delcol_factor 10, iter_num 10,
            # confidence 0.75, w_thresh 10, unknown False). Nothing tuned against ACS.
            ARIC(mix_path=str(tmp / "mix.csv"), ref_path=str(tmp / "ref.csv"),
                 save_path=str(tmp / "out.csv"))
            out = pd.read_csv(tmp / "out.csv", index_col=0)
        # ARIC writes cell types x samples (verified on a planted mixture, 2026-09-30).
        if set(out.index) != set(cts):
            raise RuntimeError(f"ARIC returned rows {list(out.index)[:8]}, expected {cts}")
        return out.T.loc[[str(c) for c in data.bulk.columns], cts].to_numpy(dtype="float64")


def build_extension_methods(names: list[str] | None = None) -> list[DeconvolutionMethod]:
    """Genuine-package wrappers for the extension panel. Never mixed into `build_methods`."""
    out: list[DeconvolutionMethod] = []
    for n, spec in EXTENSION_SPECS.items():
        if names is not None and n not in names:
            continue
        if spec.route == "R":
            out.append(RMethod(n, _NoReimplementation(n), allow_fallback=False))
        elif spec.route == "python":
            out.append(_PythonPackageMethod(n, {"aric": "ARIC", "rnasieve": "rnasieve"}[n]))
    return out


def is_extension(name: str) -> bool:
    return name in EXTENSION_SPECS


def route_of(name: str) -> str:
    return EXTENSION_SPECS[name].route
