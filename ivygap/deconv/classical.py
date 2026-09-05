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
import pandas as pd
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


def _combat_adjust(x: np.ndarray, batch: np.ndarray) -> np.ndarray:
    """
    Empirical-Bayes batch adjustment, in the shape ComBat uses.

    Genes are rows, samples are columns, `batch` labels each column. For each gene the
    grand mean and pooled variance are estimated, each batch's location and scale shift
    is estimated, those per-gene shifts are shrunk toward the batch's own mean shift
    (the empirical-Bayes step, which is what keeps a single noisy gene from being
    "corrected" by its own noise), and the data are standardised, adjusted, and put back
    on the original scale.

    This is a faithful implementation of the parametric ComBat adjustment, not a call
    into `sva`: neither `sva` (R) nor `pycombat` (Python) is installed here, and adding
    a Bioconductor dependency for one matrix operation is a worse trade than writing the
    fifty lines it actually is. It is labelled as this project's implementation
    everywhere it appears, exactly as the reimplemented solvers are.
    """
    x = np.asarray(x, dtype="float64")
    batch = np.asarray(batch)
    levels = list(dict.fromkeys(batch.tolist()))
    if len(levels) < 2:
        return x

    grand = x.mean(axis=1, keepdims=True)
    var = x.var(axis=1, ddof=1, keepdims=True)
    var[var <= 0] = 1.0
    z = (x - grand) / np.sqrt(var)

    out = z.copy()
    for lv in levels:
        cols = batch == lv
        n = int(cols.sum())
        if n == 0:
            continue
        gamma = z[:, cols].mean(axis=1)                      # location shift per gene
        delta = z[:, cols].var(axis=1, ddof=1) if n > 1 else np.ones(z.shape[0])
        delta[~np.isfinite(delta) | (delta <= 0)] = 1.0

        # empirical-Bayes shrinkage toward the batch's own mean shift
        g_bar, t2 = gamma.mean(), gamma.var(ddof=1) if gamma.size > 1 else 1.0
        t2 = t2 if np.isfinite(t2) and t2 > 0 else 1.0
        gamma_star = (t2 * n * gamma + delta * g_bar) / (t2 * n + delta)

        d_bar, s2 = delta.mean(), delta.var(ddof=1) if delta.size > 1 else 1.0
        s2 = s2 if np.isfinite(s2) and s2 > 0 else 1.0
        a = 2 * s2 / max(d_bar, 1e-12) ** 2 + 2
        b = d_bar * (a - 1)
        delta_star = (0.5 * n * delta + b) / (0.5 * n + a - 1)
        delta_star[~np.isfinite(delta_star) | (delta_star <= 0)] = 1.0

        out[:, cols] = (z[:, cols] - gamma_star[:, None]) / np.sqrt(delta_star)[:, None]

    return out * np.sqrt(var) + grand


class CIBERSORTxDeconvolution(SVRDeconvolution):
    """
    CIBERSORTx, B-mode: nu-SVR deconvolution with bulk-mode batch correction.

    THE BASE ALGORITHM IS ALREADY HERE
    ----------------------------------
    CIBERSORT's core is nu-support-vector regression on standardised inputs, sweeping
    nu over {0.25, 0.5, 0.75} and keeping the fit with the lowest reconstruction error,
    with negative coefficients clipped and the rest renormalised. That is exactly what
    `SVRDeconvolution` does, which is why this class inherits it rather than restating
    it — `svr` in this project's leaderboard IS the published CIBERSORT algorithm, and
    saying so is more useful than shipping a near-identical duplicate row.

    WHAT CIBERSORTx ADDS, AND WHY IT MATTERS HERE SPECIFICALLY
    ----------------------------------------------------------
    Newman et al. (2019), Nat Biotechnol 37:773, add cross-platform batch correction.
    Their B-mode, from the Methods:

        1. deconvolve the mixtures M against the signature S to get fractions F;
        2. build reconstructed mixtures M* = S F in non-log linear space;
        3. log2-adjust M and M*, and apply ComBat to remove technical variation
           between them;
        4. re-estimate the fractions from the adjusted mixtures in linear space.

    The paper is explicit that this is not cosmetic: with signature matrices derived
    from droplet/UMI platforms, "deconvolution may fail ... cell types that are expected
    to be present are observed to 'drop out'". That is this project's exact
    configuration — a 10x Chromium atlas (GBmap) deconvolving 2014 laser-capture bulk —
    so a method that models the platform gap is testing something the others cannot.

    The paper requires at least three mixture samples and recommends ten; with fewer,
    this falls back to the uncorrected fit and says so.

    WHAT IS NOT IMPLEMENTED
    -----------------------
    S-mode, which adjusts the signature rather than the mixtures and is the mode the
    paper actually recommends for droplet-derived signatures. Its algorithm is given in
    Supplementary Note 1, which is not in the PDF available here, and guessing at it
    would produce something that is not S-mode wearing its name. Nor is the hosted
    CIBERSORTx service used: it is web/licence-gated, so no result here comes from
    Stanford's implementation.
    """

    name = "cibersortx"
    family = "least-squares"

    #: The paper's own floor for batch correction.
    MIN_SAMPLES_FOR_BATCH = 3

    def __init__(self, **params):
        super().__init__(**params)
        self.batch_corrected_: bool = False
        self.batch_skip_reason_: str | None = None

    def _fractions(self, raw: np.ndarray) -> np.ndarray:
        w = np.clip(raw, 0.0, None)
        tot = w.sum(axis=1, keepdims=True)
        return np.divide(w, tot, out=np.zeros_like(w), where=tot > 0)

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        raw = super()._solve_all(data)                       # step 1: base CIBERSORT
        S, B = self._as_arrays(data)

        if B.shape[1] < self.MIN_SAMPLES_FOR_BATCH:
            self.batch_corrected_ = False
            self.batch_skip_reason_ = (
                f"{B.shape[1]} mixture sample(s); the published procedure requires at "
                f"least {self.MIN_SAMPLES_FOR_BATCH} and recommends ten. Reported as the "
                f"uncorrected nu-SVR fit.")
            return raw

        F = self._fractions(raw)                             # samples x types
        M_star = S @ F.T                                     # step 2: genes x samples

        # step 3: log2, then ComBat across the two "batches" (observed vs reconstructed)
        lg = lambda m: np.log2(np.clip(m, 0, None) + 1.0)
        stacked = np.hstack([lg(B), lg(M_star)])
        batch = np.array(["M"] * B.shape[1] + ["Mstar"] * M_star.shape[1])
        adj = _combat_adjust(stacked, batch)

        # step 4: back to linear space and re-estimate
        B_adj = np.clip(np.exp2(adj[:, :B.shape[1]]) - 1.0, 0.0, None)

        adjusted = DeconvolutionInput(
            bulk=pd.DataFrame(B_adj, index=data.bulk.index, columns=data.bulk.columns),
            references=data.references, manifest=data.manifest,
            cell_types=data.cell_types,
            apply_cell_size_correction=data.apply_cell_size_correction)

        self.batch_corrected_ = True
        self.batch_skip_reason_ = None
        return super()._solve_all(adjusted)
