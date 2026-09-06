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

import dataclasses
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.exceptions import ConvergenceWarning
from sklearn.svm import NuSVR

from .. import config
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

    THIS IS B-MODE, AND B-MODE IS PROBABLY NOT THE RIGHT MODE HERE
    ---------------------------------------------------------------
    Supplementary Table 1d of the same paper records the mode the authors chose for
    every deconvolution they ran, and the pattern is consistent: a signature derived
    from 10x Chromium and applied to a bulk RNA-seq mixture is deconvolved in S-MODE,
    every time. B-mode is what they use for SMART-Seq2- and microarray-derived
    signatures. The split is mechanistic — droplet data carries a strong 3' bias and UMI
    counting, so it sits further from bulk than full-length SMART-Seq2 does.

    GBmap is 87.1% 10x (measured: 214,284 cells 3' v2, 49,262 3' v3, 31,316 5' v1, of
    338,564; Smart-seq2 is 2.7%). So by the paper's own practice this configuration is
    an S-mode configuration. `CIBERSORTxSModeDeconvolution` below implements it; this
    class stays as B-mode and is reported under its own name, so the two modes can be
    compared rather than one silently standing in for the other.

    The hosted CIBERSORTx service is not used: it is web/licence-gated, so no result
    here comes from Stanford's implementation.
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


class CIBERSORTxSModeDeconvolution(SVRDeconvolution):
    """
    CIBERSORTx S-mode: nu-SVR against a batch-ADJUSTED SIGNATURE.

    Where B-mode moves the mixtures toward the signature, S-mode moves the signature
    toward the mixtures. The paper recommends it for signatures derived from
    droplet/UMI platforms, which is this project's configuration — see the note on
    `CIBERSORTxDeconvolution` and `configs.DECLARED_DEVIATIONS["cibersortx"]`.

    THE ALGORITHM, from Supplementary Note 1 (p39 of the Supplementary Information)
    ------------------------------------------------------------------------------
    Given the signature `B` (genes x c types) and the single cells `R` those profiles
    were built from:

      1. mu = each cell type's fractional abundance in R; sigma = 2*mu.
      2. Draw F* (c types x k artificial mixtures) from N(mu, sigma), per type.
      3. Clip negatives in F* to 0; renormalise each mixture's column to sum 1.
      4. Sample single cells per type according to F* and aggregate into k bulk
         profiles in TPM space -> M*.
      5. ComBat on M and M* jointly, in log2 space -> Madj, Madj*.
      6. Back to linear; per gene, NNLS of Madj*[g,:] on F*.T gives that gene's row of
         the adjusted signature -> Badj.
      7. Estimate F from the ORIGINAL mixtures M against Badj.

    Step 7 uses M and not Madj, and that is deliberate in the paper: ComBat moved M*
    toward M's batch, so Badj is already expressed in M's space.

    The paper notes step 6 needs no adaptive noise filtration, because F* is known
    exactly and ComBat is linear in log2 space, so it preserves ordering.

    WHY THIS ONE NEEDS CELLS AND B-MODE DOES NOT
    --------------------------------------------
    B-mode reconstructs its mixtures as S @ F.T, so a collapsed signature is enough.
    S-mode resamples actual single cells, so without a registered cell source there is
    nothing to build M* from. It raises rather than quietly doing something else: an
    S-mode row computed without R would not be S-mode, and this project's rule is that
    a method never reports under a name it did not earn.
    """

    name = "cibersortx_smode"
    family = "least-squares"

    #: Artificial mixtures to build. Step 6 is a per-gene NNLS on a (k x c) design, so
    #: k must exceed c; the paper's fallback for k <= c is to manufacture more
    #: pseudo-mixtures, which at 100 vs 8 types never binds here.
    N_ARTIFICIAL = 100
    #: Cells aggregated per artificial mixture.
    CELLS_PER_MIXTURE = 500

    def __init__(self, **params):
        super().__init__(**params)
        self.signature_adjusted_: bool = False
        self.n_artificial_: int = 0
        self.adjustment_note_: str | None = None

    def _cells(self, data: DeconvolutionInput):
        """
        The registered single cells, aligned to the run's gene space and roster.

        Read through `r_bridge`'s cell-source registry, which is what the benchmark
        stage narrows to TRAINING donors only. Going through it rather than loading the
        atlas directly is what keeps S-mode inside the donor-held-out design instead of
        quietly reading the test donors' cells.
        """
        from ivygap.deconv.r_bridge import get_cell_source

        src = get_cell_source(data.primary.name)
        if src is None:
            raise RuntimeError(
                f"S-mode needs the single cells the signature was built from, and no "
                f"cell source is registered for reference {data.primary.name!r}. "
                f"B-mode reconstructs its mixtures from the signature alone and can run "
                f"without them; S-mode cannot, and will not substitute B-mode under "
                f"S-mode's name."
            )
        expression, meta = src

        shared = [c for c in expression.columns if c in meta.index]
        if not shared:
            raise RuntimeError("no cell ids shared between registered expression and metadata")
        expression, meta = expression[shared], meta.loc[shared]

        genes = list(data.bulk.index)
        missing = [g for g in genes if g not in expression.index]
        if missing:
            raise RuntimeError(
                f"{len(missing)} of the run's {len(genes)} genes are absent from the "
                f"registered cell source; S-mode's artificial mixtures must live on the "
                f"same gene space as the real ones")
        return expression.loc[genes], meta

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        expression, meta = self._cells(data)
        types = list(data.cell_types)
        c = len(types)
        rng = np.random.default_rng(config.RANDOM_SEED)

        labels = meta["cell_type"].astype(str).to_numpy()
        by_type = {t: np.flatnonzero(labels == t) for t in types}
        empty = [t for t in types if by_type[t].size == 0]
        if empty:
            raise RuntimeError(
                f"the registered cell source has no cells for {empty}; S-mode draws its "
                f"artificial mixtures from real cells, so a type with none cannot be "
                f"given a column in the adjusted signature")

        # step 1: mu is each type's fractional abundance in R; sigma = 2*mu
        counts = np.array([by_type[t].size for t in types], dtype="float64")
        mu = counts / counts.sum()
        sigma = 2.0 * mu

        # steps 2-3: F* ~ N(mu, sigma), clipped at zero and renormalised per mixture
        k = self.N_ARTIFICIAL
        F_star = rng.normal(mu[:, None], sigma[:, None], size=(c, k))
        np.clip(F_star, 0.0, None, out=F_star)
        col = F_star.sum(axis=0, keepdims=True)
        # A column can clip to all-zero. Falling back to mu keeps the mixture rather
        # than dropping it, which would silently shrink k below what step 6 was sized for.
        dead = (col <= 0).ravel()
        if dead.any():
            F_star[:, dead] = mu[:, None]
            col = F_star.sum(axis=0, keepdims=True)
        F_star /= col

        # step 4: aggregate real cells into k bulk profiles, in TPM space
        X = expression.to_numpy(dtype="float64")             # genes x cells
        M_star = np.empty((X.shape[0], k), dtype="float64")
        for j in range(k):
            n_per = rng.multinomial(self.CELLS_PER_MIXTURE, F_star[:, j])
            picks = [rng.choice(by_type[t], size=n, replace=True)
                     for t, n in zip(types, n_per) if n > 0]
            idx = np.concatenate(picks)
            M_star[:, j] = X[:, idx].sum(axis=1)
        tot = M_star.sum(axis=0, keepdims=True)
        M_star = np.divide(M_star, tot, out=np.zeros_like(M_star), where=tot > 0) * 1e6

        # step 5: ComBat over the two batches, in log2 space
        _, B = self._as_arrays(data)
        lg = lambda m: np.log2(np.clip(m, 0, None) + 1.0)
        stacked = np.hstack([lg(B), lg(M_star)])
        batch = np.array(["M"] * B.shape[1] + ["Mstar"] * k)
        adj = _combat_adjust(stacked, batch)

        # step 6: linear space, then a per-gene NNLS on the (k x c) design F*.T
        M_star_adj = np.clip(np.exp2(adj[:, B.shape[1]:]) - 1.0, 0.0, None)
        design = F_star.T                                    # k x c
        B_adj = np.empty((M_star_adj.shape[0], c), dtype="float64")
        for g in range(M_star_adj.shape[0]):
            B_adj[g], _ = nnls(design, M_star_adj[g])

        # A gene the adjustment zeroes out everywhere carries no signal and would make
        # the signature rank-deficient; keep the original row for those.
        dead_genes = ~np.isfinite(B_adj).all(axis=1) | (B_adj.sum(axis=1) <= 0)
        S_orig = data.primary.profile.to_numpy(dtype="float64")
        if dead_genes.any():
            B_adj[dead_genes] = S_orig[dead_genes]

        self.signature_adjusted_ = True
        self.n_artificial_ = k
        self.adjustment_note_ = (
            f"signature adjusted by S-mode against {k} artificial mixtures of "
            f"{self.CELLS_PER_MIXTURE} cells each, drawn from {X.shape[1]} registered "
            f"cells; {int(dead_genes.sum())} gene(s) kept their original profile")

        # step 7: deconvolve the ORIGINAL mixtures against the adjusted signature
        ref_adj = dataclasses.replace(
            data.primary,
            profile=pd.DataFrame(B_adj, index=data.bulk.index, columns=types))
        adjusted = DeconvolutionInput(
            bulk=data.bulk, references=(ref_adj,) + tuple(data.references[1:]),
            manifest=data.manifest, cell_types=data.cell_types,
            apply_cell_size_correction=data.apply_cell_size_correction,
            bulk_full=data.bulk_full)
        return super()._solve_all(adjusted)
