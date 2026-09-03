"""
elastic_net.py — regularised deconvolution.

WHAT PROBLEM THIS SOLVES
------------------------
The reference columns of a GBM signature matrix are not independent. Tumour cells in
GBM adopt an astrocyte-like state; TAM-BDM and TAM-MG share most of their myeloid
programme; NK and T cells share cytotoxic machinery. When columns of S are
near-collinear, the least-squares problem is ill-conditioned: many quite different w
vectors reconstruct b almost equally well, so the solution the solver happens to land
on is driven by noise. On TCGA this showed up concretely — Astrocyte was
unidentifiable against Tumor and GFAP alone carried 38.8% of the squared residual.

Elastic Net addresses exactly this. The L2 (ridge) term makes the problem strictly
convex and *shares* weight smoothly between collinear columns instead of letting the
solver pick one arbitrarily. The L1 term encourages genuinely absent cell types to go
to exactly zero rather than to small positive noise, which matters for Ivy GAP's
leading-edge samples where lymphocytes really are near-absent.

    minimise  (1/2n)||Sw - b||^2  +  alpha * [ rho||w||_1 + (1-rho)/2 ||w||^2 ]
    subject to  w >= 0

AN HONEST NOTE ON L1 UNDER A SIMPLEX CONSTRAINT
-----------------------------------------------
Once w is renormalised to sum to 1, ||w||_1 is identically 1, so on the *normalised*
scale the L1 term is a constant and does nothing. It is not useless here, because we
solve on the *unnormalised* scale where w carries magnitude: there L1 does shrink and
does zero out columns, and the pattern of zeros survives renormalisation. But it does
mean rho behaves differently than in ordinary regression, and that a pure-lasso
setting (rho=1) mostly rescales rather than selects. The rho grid below therefore
centres on ridge-leaning values, and the selected (alpha, rho) is recorded per sample
so the choice is auditable rather than buried.

HYPERPARAMETER SELECTION IS OUTCOME-BLIND
-----------------------------------------
alpha and rho are chosen by held-out-*gene* cross-validation on reconstruction error,
per sample. Genes are split into folds, the coefficients are fitted using one subset
of genes, and error is measured on the genes that were held out. This uses only the
expression data — no survival, no anatomic label, no known composition. The project's
inherited rule that outcomes may never be used to select a method applies just as
strongly to selecting a method's hyperparameters.
"""

from __future__ import annotations

import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet

from ivygap import config
from .base import DeconvolutionInput, DeconvolutionMethod


class ElasticNetDeconvolution(DeconvolutionMethod):
    """
    Non-negative Elastic Net deconvolution with per-sample, outcome-blind tuning.

    Parameters
    ----------
    alpha, l1_ratio
        Fix them to skip tuning (useful when the benchmark has already tuned them on
        pseudobulk and wants the same values applied to real samples). Leave both None
        for per-sample held-out-gene selection.
    n_gene_folds
        Folds for the internal gene-wise CV.
    """

    name = "elastic_net"
    family = "regularized"

    #: Number of points on the alpha path. The path itself is computed per sample from
    #: the data (see `_alpha_path`) rather than fixed here: a grid in absolute units
    #: means completely different amounts of shrinkage for a sparse leading-edge sample
    #: and a dense cellular-tumour one, so a fixed grid would silently apply a
    #: different method to different anatomic structures.
    N_ALPHAS = 6
    #: Ratio of the smallest alpha on the path to the largest.
    ALPHA_MIN_RATIO = 1e-3
    #: Ridge-leaning by design — see the module docstring on L1 under a simplex.
    L1_RATIO_GRID = (0.05, 0.2, 0.5)

    def __init__(self, alpha: float | None = None, l1_ratio: float | None = None,
                 n_gene_folds: int = 3, **params):
        super().__init__(alpha=alpha, l1_ratio=l1_ratio,
                         n_gene_folds=n_gene_folds, **params)
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.n_gene_folds = n_gene_folds
        #: (alpha, l1_ratio) actually used for each sample, in sample order
        self.selected_: list[tuple[float, float]] = []

    # -- internals -----------------------------------------------------------

    @staticmethod
    def _fit_one(S: np.ndarray, b: np.ndarray, alpha: float, l1_ratio: float) -> np.ndarray:
        """
        One non-negative Elastic Net fit.

        `fit_intercept=False` is essential: an intercept would absorb a constant
        expression offset that the mixing model attributes to cell content, and the
        model b = Sw has no intercept term to begin with.
        """
        model = ElasticNet(
            alpha=alpha, l1_ratio=l1_ratio, positive=True, fit_intercept=False,
            max_iter=5000, tol=1e-5, selection="cyclic",
            random_state=config.RANDOM_SEED,
        )
        with warnings.catch_warnings():
            # A path point that does not converge is not an error — it is a point on
            # the path that CV should reject. Letting it warn once per fit per sample
            # would bury the run's real output under thousands of lines.
            warnings.simplefilter("ignore", ConvergenceWarning)
            model.fit(S, b)
        return np.clip(np.asarray(model.coef_, dtype="float64"), 0.0, None)

    @classmethod
    def _alpha_path(cls, S: np.ndarray, b: np.ndarray, l1_ratio: float) -> np.ndarray:
        """
        The alpha values worth trying for THIS sample.

        alpha_max is the smallest penalty that drives every coefficient to zero. For
        the elastic-net objective sklearn minimises, that is max|S'b| / (n * l1_ratio).
        Anchoring the path there and descending by a fixed ratio gives a grid that
        spans "everything shrunk away" to "barely penalised" regardless of the sample's
        scale — the standard glmnet construction, and the reason this runs in
        milliseconds where an absolute grid spent most of its time on penalties so
        large that coordinate descent ground against max_iter without moving.
        """
        n = S.shape[0]
        alpha_max = float(np.max(np.abs(S.T @ b))) / max(n * l1_ratio, 1e-12)
        if not np.isfinite(alpha_max) or alpha_max <= 0:
            return np.array([1e-4])
        return np.geomspace(alpha_max * cls.ALPHA_MIN_RATIO, alpha_max, cls.N_ALPHAS)

    def _select(self, S: np.ndarray, b: np.ndarray,
                rng: np.random.Generator) -> tuple[float, float]:
        """Held-out-gene CV over the (alpha, rho) grid. Returns the best pair."""
        n_genes = S.shape[0]
        folds = np.arange(n_genes) % self.n_gene_folds
        rng.shuffle(folds)

        best, best_err = (None, self.L1_RATIO_GRID[0]), np.inf
        for rho in self.L1_RATIO_GRID:
            for alpha in self._alpha_path(S, b, rho):
                errs = []
                for k in range(self.n_gene_folds):
                    tr, te = folds != k, folds == k
                    # A fold that leaves a signature column all-zero cannot inform
                    # that column's coefficient; skip rather than fit a rank-deficient
                    # subproblem and record a meaningless error.
                    if not tr.any() or not te.any() or (S[tr].sum(axis=0) <= 0).any():
                        continue
                    try:
                        w = self._fit_one(S[tr], b[tr], alpha, rho)
                    except Exception:                       # noqa: BLE001
                        errs.append(np.inf)
                        continue
                    errs.append(float(np.mean((S[te] @ w - b[te]) ** 2)))
                if errs:
                    err = float(np.mean(errs))
                    if err < best_err:
                        best, best_err = (float(alpha), rho), err
        if best[0] is None:
            # Every candidate failed; fall back to the least-penalised point rather
            # than returning a sentinel that would silently become alpha=None.
            return float(self._alpha_path(S, b, self.L1_RATIO_GRID[0])[0]), self.L1_RATIO_GRID[0]
        return best

    # -- the interface -------------------------------------------------------

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        rng = np.random.default_rng(config.RANDOM_SEED)

        out = np.zeros((B.shape[1], S.shape[1]))
        self.selected_ = []

        for j in range(B.shape[1]):
            b = B[:, j]

            # Put every sample on a comparable scale so one alpha grid means the same
            # thing for a deeply and a shallowly sequenced sample. Scaling b and S by
            # the same constant leaves the *proportions* untouched, which is all we
            # report, so this is free.
            scale = float(np.linalg.norm(b))
            if scale <= 0:
                out[j] = np.nan
                self.selected_.append((np.nan, np.nan))
                continue
            b_s, S_s = b / scale, S / scale

            if self.alpha is not None and self.l1_ratio is not None:
                alpha, rho = self.alpha, self.l1_ratio
            else:
                alpha, rho = self._select(S_s, b_s, rng)
            self.selected_.append((alpha, rho))

            w = self._fit_one(S_s, b_s, alpha, rho)

            # Heavy shrinkage can drive every coefficient to zero. Rather than emit a
            # meaningless all-zero composition, fall back to the smallest alpha on the
            # grid and record that this happened.
            if w.sum() <= 0:
                weakest = float(self._alpha_path(S_s, b_s, rho)[0])
                w = self._fit_one(S_s, b_s, weakest, rho)
                self.selected_[-1] = (weakest, rho)
            out[j] = w
        return out

    def selection_report(self, samples: list[str]):
        """Per-sample (alpha, rho) actually used — written alongside the estimates."""
        import pandas as pd
        return pd.DataFrame(self.selected_, index=samples,
                            columns=["alpha", "l1_ratio"])
