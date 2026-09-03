"""
reference_based.py — Python implementations of MuSiC, DWLS, Bisque and SCDC.

WHY THESE FOUR
--------------
They are not four flavours of the same idea; each attacks a different failure mode of
plain NNLS, and Ivy GAP exercises all four:

  MuSiC   weights genes by how *consistent* they are across reference donors, so a
          gene that means different things in different people stops driving the fit.
          Ivy GAP has ~37 tumours and strong between-patient heterogeneity, so
          cross-subject consistency is exactly the right currency.
  DWLS    weights genes by inverse squared fitted value, which stops a handful of
          enormously expressed genes from monopolising the residual and lets rare
          types be seen at all. Ivy GAP's lymphocyte fractions are in the low
          single digits; without this they are numerically invisible.
  Bisque  corrects the platform/assay gap between a single-cell reference and bulk
          RNA-seq with a gene-wise linear transform. Ivy GAP bulk is 2014-era
          laser-capture RNA-seq and the reference is modern droplet scRNA-seq, which
          is about as wide as that gap gets.
  SCDC    runs across several references and ENSEMBLEs them, so no single reference's
          idiosyncrasies decide the answer. This is the only method here that can
          measure its own reference sensitivity.

RELATIONSHIP TO THE PUBLISHED R PACKAGES
----------------------------------------
These are faithful reimplementations of the published algorithms, not bindings. They
exist so the pipeline runs end-to-end with no R toolchain, and so the benchmark can be
reproduced on a machine where an R install has drifted.

When R and the real packages ARE available, `ivygap.deconv.r_bridge` runs them instead
and the registry prefers the genuine implementation. Any run records, per method,
which of the two produced the numbers — see the `implementation` field in the results
manifest. Reporting a Python reimplementation as though it were the published package
would be a straightforward misrepresentation, so the distinction is carried all the
way into the release bundle rather than mentioned in a README.

Known, deliberate simplifications are documented on each class.

REFERENCES
----------
MuSiC   Wang X et al. Nat Commun 10:380 (2019).  doi:10.1038/s41467-018-08023-x
DWLS    Tsoucas D et al. Nat Commun 10:2975 (2019). doi:10.1038/s41467-019-10802-z
Bisque  Jew B et al. Nat Commun 11:1971 (2020).  doi:10.1038/s41467-020-15816-6
SCDC    Dong M et al. Brief Bioinform 22:416 (2021). doi:10.1093/bib/bbz166
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import nnls

from ivygap import config
from .base import (DeconvolutionInput, DeconvolutionMethod, ReferenceBundle,
                   project_to_simplex)

EPS = 1e-12


def _weighted_nnls(S: np.ndarray, b: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    Solve min sum_g w_g (b_g - (Sw)_g)^2 with coefficients >= 0.

    Weighted least squares is ordinary least squares on sqrt(w)-scaled rows, so this
    is one `nnls` call on rescaled inputs rather than a bespoke solver.
    """
    rw = np.sqrt(np.clip(w, 0.0, None))
    try:
        coef, _ = nnls(S * rw[:, None], b * rw)
    except Exception:                                   # noqa: BLE001
        coef, _ = nnls(S, b)                            # fall back to unweighted
    return coef


class MuSiCDeconvolution(DeconvolutionMethod):
    """
    MuSiC: iteratively reweighted NNLS with cross-subject-variance gene weights.

    The weight for gene g at the current estimate p is

        W_g = 1 / ( sum_k p_k^2 * Sigma_gk  +  sigma_eps^2 )

    where Sigma_gk is the variance of gene g's expression in type k *across reference
    donors*, and sigma_eps^2 is the current residual variance. A gene that is stable
    across donors gets a large weight; one that swings gets a small one. The estimate
    and the weights are alternated to convergence.

    Simplification vs. the R package: the recursive tree (`music_prop.cluster`), used
    when cell types are hierarchically related and hard to separate, is not
    implemented. The eight-type GBM roster here is flat by construction, so the tree
    would have nothing to recurse on.
    """

    name = "music"
    family = "reference-based"

    def __init__(self, max_iter: int = 100, tol: float = 1e-5,
                 trim_quantile: float = 0.99, **params):
        super().__init__(max_iter=max_iter, tol=tol, trim_quantile=trim_quantile, **params)
        self.max_iter = max_iter
        self.tol = tol
        self.trim_quantile = trim_quantile
        self.degenerate_: bool = False
        self.degeneracy_reason_: str | None = None

    def _solve_one(self, S: np.ndarray, sigma: np.ndarray, b: np.ndarray) -> np.ndarray:
        p, _ = nnls(S, b)                                # unweighted warm start
        p = project_to_simplex(p)
        if not np.isfinite(p).all():
            return np.zeros(S.shape[1])

        for _ in range(self.max_iter):
            resid = b - S @ p
            sigma_eps = float(np.mean(resid ** 2))
            var_g = sigma @ (p ** 2) + sigma_eps
            w = 1.0 / np.clip(var_g, EPS, None)

            # A near-zero-variance gene can take an astronomically large weight and
            # single-handedly determine the fit. MuSiC trims that tail; without the cap
            # the "weighting" degenerates into fitting one gene.
            cap = np.quantile(w, self.trim_quantile)
            w = np.minimum(w, cap)

            p_new = project_to_simplex(_weighted_nnls(S, b, w))
            if not np.isfinite(p_new).all():
                break
            if np.max(np.abs(p_new - p)) < self.tol:
                p = p_new
                break
            p = p_new
        return p

    def degradation_for(self, data: DeconvolutionInput) -> tuple[bool, str | None]:
        ref = data.primary
        degraded = not (getattr(ref, "has_cross_donor_variance", True)
                        and ref.sigma is not None)
        if not degraded:
            return False, None
        return True, (
            f"reference {ref.name!r} carries no cross-donor variance, so MuSiC's gene "
            f"weighting is constant and the result is arithmetically NNLS. Reported as "
            f"degenerate rather than as a MuSiC result."
        )

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        ref = data.primary
        S = ref.profile.to_numpy(dtype="float64")
        sigma = ref.sigma.to_numpy(dtype="float64")
        B = data.bulk.to_numpy(dtype="float64")

        # With no cross-donor variance the weights collapse to a constant and this is
        # arithmetically NNLS. Recording that is the difference between a benchmark and
        # a leaderboard with a mislabelled row.
        self.degenerate_, self.degeneracy_reason_ = self.degradation_for(data)

        return np.vstack([self._solve_one(S, sigma, B[:, j]) for j in range(B.shape[1])])


class DWLSDeconvolution(DeconvolutionMethod):
    """
    Dampened Weighted Least Squares — built for rare cell types.

    Ordinary least squares minimises *absolute* error, so it is dominated by the most
    highly expressed genes, which belong to the most abundant cell type. A cell type at
    1% contributes almost nothing to the objective and its estimate is essentially
    unconstrained. DWLS reweights by the inverse squared fitted value

        W_g = 1 / (S p)_g^2

    which converts the objective from absolute to *relative* error, putting a rare
    type's genes on equal footing with the tumour's.

    Taken literally this overshoots: a gene with a tiny fitted value gets an enormous
    weight and destabilises everything. Hence "dampened" — weights are capped at
    max(W)/2^j, and j is chosen as the smallest dampening that makes the solution
    stable, measured by re-solving on random gene subsets and taking the j whose
    solution varies least.

    Simplification vs. the R package: the published implementation selects j by a
    multinomial-likelihood criterion over bootstrapped subsets. This uses the same
    subset-resampling idea with a solution-variance criterion, which selects the same
    j on the fixtures in `tests/` but is not guaranteed to agree in every case.
    """

    name = "dwls"
    family = "reference-based"

    def __init__(self, max_iter: int = 100, tol: float = 1e-5,
                 dampen_powers: tuple[int, ...] = (1, 2, 4, 6, 8),
                 n_subsets: int = 5, subset_frac: float = 0.5,
                 select_max_iter: int = 12, **params):
        super().__init__(max_iter=max_iter, tol=tol, dampen_powers=dampen_powers,
                         n_subsets=n_subsets, subset_frac=subset_frac,
                         select_max_iter=select_max_iter, **params)
        self.max_iter = max_iter
        self.tol = tol
        self.dampen_powers = dampen_powers
        self.n_subsets = n_subsets
        self.subset_frac = subset_frac
        # Power selection compares how *stable* solutions are across gene subsets, and
        # that ordering is established long before the fits converge. Running the
        # selection loop to full tolerance multiplied the cost of the method by an
        # order of magnitude and changed the chosen power on none of the fixtures.
        self.select_max_iter = select_max_iter

    @staticmethod
    def _weights(S: np.ndarray, p: np.ndarray, power: int) -> np.ndarray:
        fitted = np.clip(S @ p, EPS, None)
        w = 1.0 / fitted ** 2
        return np.minimum(w, w.max() / (2.0 ** power))

    def _solve_dampened(self, S: np.ndarray, b: np.ndarray, power: int,
                        max_iter: int | None = None) -> np.ndarray:
        p, _ = nnls(S, b)
        p = project_to_simplex(p)
        if not np.isfinite(p).all():
            return np.zeros(S.shape[1])
        for _ in range(max_iter if max_iter is not None else self.max_iter):
            p_new = project_to_simplex(_weighted_nnls(S, b, self._weights(S, p, power)))
            if not np.isfinite(p_new).all():
                break
            if np.max(np.abs(p_new - p)) < self.tol:
                return p_new
            p = p_new
        return p

    def _choose_power(self, S: np.ndarray, b: np.ndarray,
                      rng: np.random.Generator) -> int:
        """Smallest dampening whose solution is stable across random gene subsets."""
        n_genes = S.shape[0]
        k = max(int(n_genes * self.subset_frac), min(50, n_genes))
        subsets = [rng.choice(n_genes, size=k, replace=False)
                   for _ in range(self.n_subsets)]

        best_power, best_var = self.dampen_powers[0], np.inf
        for power in self.dampen_powers:
            sols = []
            for idx in subsets:
                if (S[idx].sum(axis=0) <= 0).any():
                    continue                              # subset cannot see some type
                sols.append(self._solve_dampened(S[idx], b[idx], power,
                                                 max_iter=self.select_max_iter))
            if len(sols) < 2:
                continue
            var = float(np.mean(np.var(np.vstack(sols), axis=0)))
            if var < best_var - 1e-9:
                best_power, best_var = power, var
        return best_power

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        rng = np.random.default_rng(config.RANDOM_SEED)
        out = []
        for j in range(B.shape[1]):
            b = B[:, j]
            out.append(self._solve_dampened(S, b, self._choose_power(S, b, rng)))
        return np.vstack(out)


class BisqueDeconvolution(DeconvolutionMethod):
    """
    Bisque reference-based decomposition: correct the assay gap, then NNLS.

    A single-cell reference and a bulk dataset do not measure genes on the same scale.
    3'-biased droplet chemistry, differing gene-length handling and differing dynamic
    range mean a gene can be systematically inflated in one and deflated in the other.
    Deconvolving without correcting for that pushes the bias straight into the
    composition estimate. Bisque estimates a per-gene linear transform from reference
    space to bulk space and applies it before solving.

    The published method estimates that transform from subjects assayed *both* ways.
    Ivy GAP has no such overlap with any public GBM single-cell atlas, so this runs
    Bisque's documented no-overlap mode (`use.overlap = FALSE`): the transform is
    estimated by matching each gene's mean and standard deviation across the *marginal*
    distributions — across bulk samples on one side, across reference donors on the
    other. This is weaker than the overlap-based version and it is the honest thing to
    report about Bisque on this dataset: it is running in its degraded mode, and any
    conclusion about Bisque's ranking here carries that caveat.
    """

    name = "bisque"
    family = "reference-based"

    def __init__(self, min_donors: int = 3, **params):
        super().__init__(min_donors=min_donors, **params)
        self.min_donors = min_donors
        # Bisque is ALWAYS degraded on this project's data, for two separate reasons,
        # and the docstring above says the caveat has to travel with any conclusion
        # about its ranking. It cannot travel unless it is set here: every report
        # surfaces degradation by reading `degenerate_` / `degeneracy_reason_`, so a
        # class that only explains itself in a docstring reports nothing at all.
        self.degenerate_: bool = False
        self.degeneracy_reason_: str | None = None

    def degradation_for(self, data: DeconvolutionInput) -> tuple[bool, str | None]:
        ref = data.primary
        reasons = ["no subjects are assayed both as bulk and as single cells, so this "
                   "runs Bisque's documented no-overlap mode (use.overlap = FALSE); "
                   "the assay transform is estimated from marginal distributions rather "
                   "than from paired subjects"]
        n_donors = len(ref.donor_profiles or {})
        if n_donors < self.min_donors:
            reasons.append(
                f"the reference carries {n_donors} donor profiles "
                f"(< min_donors={self.min_donors}), so the per-gene reference marginal "
                f"is taken across CELL TYPES instead of across donors — a different "
                f"quantity standing in for the one the method specifies")
        return True, ("; ".join(reasons)
                      + ". Reported as degraded rather than as a Bisque result.")

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        ref = data.primary
        S = ref.profile.to_numpy(dtype="float64")
        B = data.bulk.to_numpy(dtype="float64")

        # Reference-side marginal distribution per gene: across donors if we have real
        # donor-level profiles, otherwise across cell types as a coarse stand-in.
        # Both degradations are declared in degradation_for() so they are reported
        # identically whether this solver or the genuine R package runs.
        self.degenerate_, self.degeneracy_reason_ = self.degradation_for(data)

        have_donors = bool(ref.donor_profiles) and len(ref.donor_profiles) >= self.min_donors
        if have_donors:
            donor_bulk = np.vstack([
                p.to_numpy(dtype="float64").mean(axis=1) for p in ref.donor_profiles.values()
            ]).T                                          # genes x donors
        else:
            # Degradation 2, and a distinct one. The marginal the transform is supposed
            # to match is across reference DONORS. A collapsed signature matrix has
            # none, so the spread across CELL TYPES stands in for it — a different
            # quantity that happens to have the right shape. Silently substituting it
            # would make Bisque look like it ran its no-overlap mode as published.
            donor_bulk = S

        ref_mu = donor_bulk.mean(axis=1)
        ref_sd = donor_bulk.std(axis=1)
        bulk_mu = B.mean(axis=1)
        bulk_sd = B.std(axis=1)

        # Map bulk into reference space gene-by-gene. Where a gene has no variance on
        # either side the transform is undefined; fall back to identity rather than
        # dividing by zero and producing NaN for every sample at once.
        scale = np.where((bulk_sd > EPS) & (ref_sd > EPS), ref_sd / np.maximum(bulk_sd, EPS), 1.0)
        B_t = (B - bulk_mu[:, None]) * scale[:, None] + ref_mu[:, None]
        B_t = np.clip(B_t, 0.0, None)                     # the mixing model is additive

        return np.vstack([nnls(S, B_t[:, j])[0] for j in range(B_t.shape[1])])


class SCDCDeconvolution(DeconvolutionMethod):
    """
    SCDC, single-reference mode: weighted NNLS with iterative outlier down-weighting.

    Close in spirit to MuSiC, but the weights come from the *residual* rather than from
    cross-subject variance, and genes whose residuals are extreme are progressively
    down-weighted rather than trimmed outright. That makes it tolerant of a reference
    that is simply wrong about a handful of genes — which every reference is.
    """

    name = "scdc"
    family = "reference-based"

    def __init__(self, max_iter: int = 50, tol: float = 1e-5,
                 outlier_z: float = 3.0, **params):
        super().__init__(max_iter=max_iter, tol=tol, outlier_z=outlier_z, **params)
        self.max_iter = max_iter
        self.tol = tol
        self.outlier_z = outlier_z

    def _solve_one(self, S: np.ndarray, b: np.ndarray) -> np.ndarray:
        p, _ = nnls(S, b)
        p = project_to_simplex(p)
        if not np.isfinite(p).all():
            return np.zeros(S.shape[1])

        for _ in range(self.max_iter):
            resid = b - S @ p
            mad = np.median(np.abs(resid - np.median(resid)))
            scale = 1.4826 * mad if mad > 0 else (resid.std() or 1.0)
            z = np.abs(resid) / max(scale, EPS)

            # Huber-style: full weight inside the threshold, decaying outside, so an
            # outlier gene is discounted rather than deleted.
            w = np.where(z <= self.outlier_z, 1.0, self.outlier_z / np.maximum(z, EPS))
            w = w / np.clip(np.abs(b) + 1.0, EPS, None)   # relative, not absolute, error

            p_new = project_to_simplex(_weighted_nnls(S, b, w))
            if not np.isfinite(p_new).all():
                break
            if np.max(np.abs(p_new - p)) < self.tol:
                return p_new
            p = p_new
        return p

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        return np.vstack([self._solve_one(S, B[:, j]) for j in range(B.shape[1])])


class SCDCEnsembleDeconvolution(SCDCDeconvolution):
    """
    SCDC ENSEMBLE: run against every reference, then combine.

    This is the only method in the project that can answer "how much does the answer
    depend on which reference I picked?", because it is the only one that sees more
    than one. That question is not academic — reference choice is one of the largest
    sources of disagreement between published deconvolution results.

    Ensemble weights are searched over a grid on the simplex and chosen to minimise
    reconstruction error of the observed bulk, summed over samples. Reconstruction is
    an unsupervised criterion: no known composition, no anatomic label, no outcome.
    The chosen weights are stored in `ensemble_weights_` and written to the results
    manifest, since "which reference won" is itself a finding.

    With a single reference supplied this degenerates to plain SCDC, and says so
    rather than pretending an ensemble happened.
    """

    name = "scdc_ensemble"
    family = "reference-based"
    uses_multiple_references = True

    def __init__(self, grid_step: float = 0.1, **params):
        super().__init__(**params)
        self.params.update(grid_step=grid_step)
        self.grid_step = grid_step
        self.ensemble_weights_: dict[str, float] = {}
        self.degenerate_: bool = False
        self.degeneracy_reason_: str | None = None

    def degradation_for(self, data: DeconvolutionInput) -> tuple[bool, str | None]:
        refs = data.references
        if len(refs) > 1:
            return False, None
        return True, (
            f"one reference supplied ({refs[0].name!r}), so there is nothing to weight "
            f"across: SCDC ENSEMBLE reduces exactly to SCDC. Reported as degenerate "
            f"rather than as an ENSEMBLE result."
        )

    def _simplex_grid(self, n: int) -> list[np.ndarray]:
        """All points on the n-simplex at the configured resolution."""
        steps = int(round(1.0 / self.grid_step))
        pts = []
        for combo in itertools.product(range(steps + 1), repeat=n):
            if sum(combo) == steps:
                pts.append(np.array(combo, dtype="float64") / steps)
        return pts

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        refs: tuple[ReferenceBundle, ...] = data.references
        B = data.bulk.to_numpy(dtype="float64")

        per_ref = []
        for ref in refs:
            S = ref.profile.to_numpy(dtype="float64")
            per_ref.append(np.vstack([self._solve_one(S, B[:, j])
                                      for j in range(B.shape[1])]))

        if len(refs) == 1:
            # Declared in degradation_for() so it is reported identically whether this
            # solver or the genuine R package runs.
            self.degenerate_, self.degeneracy_reason_ = self.degradation_for(data)
            self.ensemble_weights_ = {refs[0].name: 1.0}
            return per_ref[0]

        S_list = [r.profile.to_numpy(dtype="float64") for r in refs]

        best_w, best_err = None, np.inf
        for wts in self._simplex_grid(len(refs)):
            err = 0.0
            for j in range(B.shape[1]):
                b = B[:, j]
                # Reconstruct with each reference's own signature, then blend the
                # reconstructions — blending the *signatures* would invent a chimeric
                # reference that no dataset supports.
                recon = sum(wts[i] * (S_list[i] @ per_ref[i][j]) for i in range(len(refs)))
                denom = np.linalg.norm(b)
                if denom > 0:
                    err += float(np.linalg.norm(b - recon) / denom)
            if err < best_err:
                best_w, best_err = wts, err

        self.ensemble_weights_ = {r.name: float(w) for r, w in zip(refs, best_w)}
        return sum(best_w[i] * per_ref[i] for i in range(len(refs)))
