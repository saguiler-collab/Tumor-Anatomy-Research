"""
classical.py — the two baselines carried over from the TCGA project: NNLS and nu-SVR.

They are not here to win. They are here because the whole claim of this rebuild is
that regularised and probabilistic methods are *more robust*, and "more robust than
what?" needs an answer computed on identical inputs rather than quoted from a paper.
NNLS is also the solver that failed informatively on TCGA (astrocyte/tumour
collinearity, GFAP dominating the residual), which makes it the right control for
asking whether Elastic Net and the Bayesian prior actually fix that failure.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import nnls
from sklearn.exceptions import ConvergenceWarning
from sklearn.svm import NuSVR

from .base import DeconvolutionInput, DeconvolutionMethod


class NNLSDeconvolution(DeconvolutionMethod):
    """
    Plain non-negative least squares: min ||Sw - b||^2 subject to w >= 0.

    No regularisation at all, so when two reference columns are near-collinear the
    solution is essentially arbitrary along that direction — small changes in the bulk
    sample swing weight between the collinear types. That is precisely the pathology
    the regularised and Bayesian methods below are meant to damp, so keeping an
    unregularised control makes the improvement measurable rather than asserted.
    """

    name = "nnls"
    family = "least-squares"

    def __init__(self, **params):
        super().__init__(**params)
        self.n_hit_iteration_cap_: int = 0

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        self.n_hit_iteration_cap_ = 0
        out = np.zeros((B.shape[1], S.shape[1]))
        for j in range(B.shape[1]):
            out[j], _ = nnls(S, B[:, j])
        return out


class SVRDeconvolution(DeconvolutionMethod):
    """
    CIBERSORT-style nu-support-vector regression.

    The epsilon-insensitive loss ignores residuals inside a tube around the fit, so a
    handful of wildly-off genes cannot drag the whole solution — which is why SVR beat
    NNLS on the TCGA held-out benchmark and became that project's frozen method. nu
    is swept over the standard CIBERSORT grid and the value with the lowest RMSE on
    the *reconstruction* is kept; that criterion uses only the expression data, never
    any outcome, so method selection stays outcome-blind.

    Support-vector coefficients can be negative; they are clipped at zero and
    renormalised in the base class, as CIBERSORT does.
    """

    name = "svr"
    family = "least-squares"

    NU_GRID = (0.25, 0.5, 0.75)

    #: libsvm defaults to unbounded iteration (max_iter=-1). On a well-conditioned
    #: problem that is fine; on a real GBM signature it is not. The reference columns
    #: are strongly collinear — Tumor against Astrocyte above all — and SMO can grind
    #: for hours without converging on exactly the ill-posed geometry this project
    #: exists to study. An unbounded solver inside a pipeline stage is a hazard, not a
    #: nicety: it turns a benchmark into a hang with no diagnostic.
    #:
    #: libsvm returns its best iterate when it hits the cap, so a capped fit is a
    #: slightly-less-converged answer rather than a missing one. Fits that hit the cap
    #: are counted in `n_hit_iteration_cap_` and reported, because "SVR did not converge
    #: on 40% of samples" is a finding about SVR on this reference, not an aside.
    MAX_ITER = 200_000

    def __init__(self, **params):
        super().__init__(**params)
        self.n_hit_iteration_cap_: int = 0

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        self.n_hit_iteration_cap_ = 0

        # SVR is scale-sensitive in a way NNLS is not; CIBERSORT standardises both
        # sides before fitting. Guard against a zero-variance signature column, which
        # would otherwise produce NaN and poison every sample.
        s_mu, s_sd = S.mean(), S.std()
        s_sd = s_sd if s_sd > 0 else 1.0
        S_z = (S - s_mu) / s_sd

        out = np.zeros((B.shape[1], S.shape[1]))
        for j in range(B.shape[1]):
            b = B[:, j]
            b_sd = b.std()
            b_z = (b - b.mean()) / (b_sd if b_sd > 0 else 1.0)

            best, best_rmse, capped = None, np.inf, False
            for nu in self.NU_GRID:
                model = NuSVR(nu=nu, C=1.0, kernel="linear",
                              max_iter=self.MAX_ITER, cache_size=500)
                with warnings.catch_warnings():
                    # ConvergenceWarning here is expected and already counted; letting
                    # it print once per nu per sample would bury the run's real output.
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    model.fit(S_z, b_z)
                if getattr(model, "n_iter_", 0) is not None and \
                        np.max(np.atleast_1d(model.n_iter_)) >= self.MAX_ITER:
                    capped = True
                coef = np.asarray(model.coef_).ravel()
                rmse = float(np.sqrt(np.mean((S_z @ coef - b_z) ** 2)))
                if rmse < best_rmse:
                    best, best_rmse = coef, rmse
            if capped:
                self.n_hit_iteration_cap_ += 1
            out[j] = np.clip(best, 0.0, None)
        return out
