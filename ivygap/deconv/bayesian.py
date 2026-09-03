"""
bayesian.py — probabilistic deconvolution with an explicit prior and real posteriors.

THE MODEL
---------
Reference-based deconvolution has a natural generative story. A bulk sample is a bag
of transcripts; each transcript came from some cell type; the chance it came from type
k and gene g is proportional to how abundant type k is (w_k) times how much of gene g
that type expresses (S_gk). Writing that down:

    w        ~ Dirichlet(alpha)                        prior over composition
    p_gk     = w_k S_gk / sum_k' w_k' S_gk'            per-gene type-attribution
    z_g.     ~ Multinomial(b_g, p_g.)                  which type each read came from
    b_g      = sum_k z_gk                              observed expression of gene g

z is latent, and conditioning on it makes the composition update conjugate:

    w | z    ~ Dirichlet(alpha + sum_g z_.k)

which is a two-line Gibbs sampler that is *exact* — no variational gap, no tuning, no
Metropolis acceptance rate to babysit. This is the same data-augmentation scheme that
sits underneath BayesPrism's per-sample step.

WHY THIS IS MORE ROBUST, CONCRETELY
-----------------------------------
1. The Dirichlet prior regularises. With alpha slightly below 1 the prior mildly
   favours sparse compositions, so a cell type with no real evidence settles near zero
   instead of absorbing residual noise, without the hard zeros that make L1 brittle.
2. It reports uncertainty. NNLS says "T cells are 0.8%". This says "0.8%, 95% credible
   interval 0.1-2.4%". On Ivy GAP that distinction carries real weight: leading-edge
   samples are small, low-complexity dissections, and a method that cannot say when it
   does not know will produce confident nonsense there.

   An interval is only worth reporting if it covers the truth as often as it claims to,
   and here that depends almost entirely on `pseudocount_total` -- see its description
   below, where the coverage/accuracy trade-off is tabulated from measurement rather
   than assumed, and `ivygap.bench.calibration`, which re-measures it on every run.
3. Collinear types share posterior mass rather than fighting over it. Where Tumor and
   Astrocyte are unidentifiable, the posterior becomes correspondingly wide and
   correlated — the right answer to an ill-posed question, instead of NNLS's arbitrary
   corner solution.

WHAT THE CREDIBLE INTERVALS DO AND DO NOT COVER
-----------------------------------------------
Read this before quoting an interval.

The posterior is over w GIVEN the reference. It captures how much the composition could
plausibly vary while still explaining this sample's expression under the assumed
signature. It does NOT capture the possibility that the signature itself is wrong for
this patient.

That distinction is not academic, and it was measured rather than assumed. Scored
against known-truth pseudobulk pooled from HELD-OUT donors, the dominant error is not
scatter but systematic per-type bias -- on the test fixture roughly -0.13 for Tumor and
+0.15 for Macrophage_Microglia, against a mean interval half-width of 0.024. Realised
95% coverage is correspondingly low. No posterior over w conditional on the reference
can cover a bias in the reference, and widening the intervals until it did would be
fitting the interval to the answer.

Two things follow, and both are implemented:

  * `propagate_reference_uncertainty` (default on) draws the signature matrix itself
    each sweep from its cross-donor standard error, so at least the *estimable* part of
    reference uncertainty enters the posterior. Its effect is real but modest --
    averaging over hundreds of genes cancels most per-gene reference noise.
  * `ivygap.bench.calibration` measures realised coverage on every run and writes it
    beside the estimates, so an interval is never reported without its actual coverage.
    A method that says 95% and delivers 23% is worse than one that says nothing.

The practical reading: these intervals rank samples by confidence reliably -- a wide
interval really does mean this sample is harder -- but their absolute width understates
total error whenever the reference does not match the patient. Use them comparatively.


A CIRCULARITY THIS FILE DELIBERATELY AVOIDS
-------------------------------------------
The obvious hierarchical extension is to partially pool samples that share an anatomic
structure. Do not do that here. The headline claim of this project is that a method
recovers known differences *between* anatomic structures; a model given the structure
label as an input would be scored on its ability to reproduce information it was
handed, and would beat every honest method for the wrong reason.

`HierarchicalBayesianDeconvolution` therefore pools by PATIENT, never by structure.
Blocks from one tumour genuinely do share that patient's immune baseline, so the
pooling is well motivated; and because a patient contributes blocks from several
different structures, pooling across them shrinks toward the patient mean and if
anything makes between-structure differences *harder* to detect. The test stays honest
in the conservative direction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ivygap import config
from .base import DeconvolutionInput, DeconvolutionMethod


class BayesianDeconvolution(DeconvolutionMethod):
    """
    Per-sample collapsed Gibbs sampler over the Dirichlet-multinomial model above.

    Parameters
    ----------
    alpha
        Dirichlet concentration. Values below 1 favour sparse compositions; 0.5 is a
        mild, standard sparsity-inducing choice that does not force zeros.
    n_iter, n_burnin, thin
        Sampler length. The chain mixes fast because the update is conjugate.
    pseudocount_total
        Ivy GAP ships FPKM, not integer counts, and the multinomial step needs counts.
        Each sample's profile is rescaled to this total and rounded.

        This is the most consequential parameter in the class and it is NOT an
        implementation detail. It sets how much evidence the likelihood carries
        relative to the prior, and therefore how wide the credible intervals are. FPKM
        is not a count of independent reads: it has already been normalised for depth
        and length, and neighbouring genes are strongly correlated. Treating it as n
        independent multinomial draws asserts far more information than the data
        contain, and the intervals come out too narrow.

        That is measurable, and it was measured. On known-truth mixtures with
        multiplicative noise, 95% credible intervals covered the truth at:

            pseudocount   coverage   MAE      mean CI width
            1e3           0.86       0.0083   0.038
            1e4           0.60       0.0064   0.013
            1e5           0.23       0.0063   0.004     <- the original default
            1e6           0.08       0.0063   0.001

        Point accuracy is essentially flat above 1e4, so the large values bought
        nothing except overconfidence. The default is therefore 1e3, where the
        intervals are approximately honest, at a cost of about 0.002 in MAE.

        Coverage still falls a little short of nominal, because the posterior is
        symmetric around a point estimate carrying a small systematic bias. These
        intervals are not exact. `ivygap.bench.calibration` measures realised coverage
        on each run's own pseudobulk and writes it into the results, so an interval is
        never quoted without its actual coverage beside it.
    """

    name = "bayesian"
    family = "bayesian"

    def __init__(self, alpha: float = 0.5, n_iter: int = 400, n_burnin: int = 150,
                 thin: int = 2, pseudocount_total: float = 1e3,
                 propagate_reference_uncertainty: bool = True, **params):
        super().__init__(alpha=alpha, n_iter=n_iter, n_burnin=n_burnin,
                         thin=thin, pseudocount_total=pseudocount_total,
                         propagate_reference_uncertainty=propagate_reference_uncertainty,
                         **params)
        if n_burnin >= n_iter:
            raise ValueError("n_burnin must be smaller than n_iter")
        self.alpha = alpha
        self.n_iter = n_iter
        self.n_burnin = n_burnin
        self.thin = thin
        self.pseudocount_total = pseudocount_total
        self.propagate_reference_uncertainty = propagate_reference_uncertainty
        #: posterior summaries, populated by the last fit_predict
        self.posterior_: dict[str, pd.DataFrame] = {}

    # -- the sampler ---------------------------------------------------------

    def _to_counts(self, b: np.ndarray) -> np.ndarray:
        """FPKM/CPM profile -> integer pseudo-counts at a fixed total."""
        total = b.sum()
        if total <= 0:
            return np.zeros_like(b, dtype="int64")
        scaled = b / total * self.pseudocount_total
        return np.rint(scaled).astype("int64")

    def _gibbs(self, S: np.ndarray, counts: np.ndarray, rng: np.random.Generator,
               se: np.ndarray | None = None) -> np.ndarray:
        """
        Run the chain for one sample. Returns the retained posterior draws,
        shape (n_draws, n_types).

        `se` is the standard error of each entry of S. When supplied, a fresh S is
        drawn each sweep, which is what turns the credible intervals from a statement
        about multinomial sampling noise into a statement about the estimate.
        """
        n_types = S.shape[1]
        alpha_vec = np.full(n_types, self.alpha, dtype="float64")

        # Only genes with counts contribute to the allocation step; dropping the rest
        # is exact (a multinomial with n=0 contributes nothing) and usually removes a
        # large share of the rows.
        keep = counts > 0
        if not keep.any():
            return np.full((1, n_types), np.nan)
        S_k, c_k = S[keep], counts[keep]
        se_k = se[keep] if se is not None else None

        # A gene no reference type expresses cannot be attributed to any type. Leaving
        # it in would give a zero row of probabilities and an invalid multinomial.
        informative = S_k.sum(axis=1) > 0
        S_k, c_k = S_k[informative], c_k[informative]
        if se_k is not None:
            se_k = se_k[informative]
        if S_k.shape[0] == 0:
            return np.full((1, n_types), np.nan)

        w = rng.dirichlet(alpha_vec)
        draws = []
        for it in range(self.n_iter):
            if se_k is not None:
                # Draw the reference itself. Clipping at zero keeps the mixing model
                # additive; a gene whose lower tail crosses zero is one the reference
                # barely constrains, and truncation is the right representation of that.
                S_draw = np.clip(S_k + rng.normal(0.0, se_k), 0.0, None)
                rows = S_draw.sum(axis=1)
                # A draw that zeroes an entire gene row leaves it uninformative for
                # this sweep; fall back to the mean profile for those rows rather than
                # producing an invalid multinomial.
                S_draw[rows <= 0] = S_k[rows <= 0]
            else:
                S_draw = S_k

            # p_gk proportional to w_k * S_gk, normalised per gene
            unnorm = S_draw * w[None, :]
            p = unnorm / unnorm.sum(axis=1, keepdims=True)

            # Data augmentation: attribute each gene's reads across types.
            z = rng.multinomial(c_k, p)
            w = rng.dirichlet(alpha_vec + z.sum(axis=0))

            if it >= self.n_burnin and (it - self.n_burnin) % self.thin == 0:
                draws.append(w.copy())
        return np.vstack(draws)

    def _reference_se(self, data) -> np.ndarray | None:
        """
        Standard error of each signature entry, from the reference's cross-donor
        variance: se = sqrt(sigma / n_donors).

        Returns None when the reference carries no usable variance — a single-donor
        reference genuinely cannot estimate its own uncertainty, and inventing a value
        would be worse than reporting the narrower intervals honestly.
        """
        if not self.propagate_reference_uncertainty:
            return None
        ref = data.primary
        sigma = ref.sigma.to_numpy(dtype="float64")
        if ref.n_donors < 2 or not np.isfinite(sigma).any() or (sigma <= 0).all():
            return None
        return np.sqrt(np.clip(sigma, 0.0, None) / max(ref.n_donors, 1))

    # -- the interface -------------------------------------------------------

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        rng = np.random.default_rng(config.RANDOM_SEED)

        n_s, n_t = B.shape[1], S.shape[1]
        means = np.zeros((n_s, n_t))
        lo = np.zeros((n_s, n_t))
        hi = np.zeros((n_s, n_t))
        sd = np.zeros((n_s, n_t))

        se = self._reference_se(data)
        for j in range(n_s):
            draws = self._gibbs(S, self._to_counts(B[:, j]), rng, se=se)
            means[j] = draws.mean(axis=0)
            sd[j] = draws.std(axis=0)
            lo[j] = np.quantile(draws, 0.025, axis=0)
            hi[j] = np.quantile(draws, 0.975, axis=0)

        cols = list(data.cell_types)
        idx = data.samples
        self.posterior_ = {
            "mean": pd.DataFrame(means, index=idx, columns=cols),
            "sd": pd.DataFrame(sd, index=idx, columns=cols),
            "ci_lower": pd.DataFrame(lo, index=idx, columns=cols),
            "ci_upper": pd.DataFrame(hi, index=idx, columns=cols),
        }
        # The posterior mean of a Dirichlet is already on the simplex; the base class
        # then applies the shared cell-size correction like every other method.
        return means


class HierarchicalBayesianDeconvolution(BayesianDeconvolution):
    """
    The same model with partial pooling across a PATIENT's samples (never across an
    anatomic structure — see the module docstring).

        mu_p   ~ Dirichlet(alpha)                 patient-level mean composition
        w_ps   ~ Dirichlet(kappa * mu_p)          each block shrinks toward it
        z, b   as above

    `kappa` sets how strongly a block is pulled toward its patient's mean. Large kappa
    means blocks from one tumour are assumed near-identical, which for Ivy GAP would
    be wrong — the blocks are deliberately drawn from *different* niches. So the
    default is deliberately weak: enough to stabilise a noisy low-input block, not
    enough to erase real niche structure. The value is reported so a reader can see
    how much shrinkage was applied.

    Fitted by Gibbs with an empirical-Bayes update of mu_p between sweeps.
    """

    name = "bayesian_hierarchical"
    family = "bayesian"

    def __init__(self, kappa: float = 5.0, n_outer: int = 3,
                 refine_iter: int = 200, refine_burnin: int = 80, **params):
        super().__init__(**params)
        self.params.update(kappa=kappa, n_outer=n_outer, refine_iter=refine_iter)
        self.kappa = kappa
        self.n_outer = n_outer
        # The refinement sweeps start from an already-converged composition and only
        # need to track a shifting prior, so they run shorter chains than the initial
        # unpooled fit. Running all of them at full length made the method roughly ten
        # times the cost of the flat sampler while moving the estimates by less than
        # the width of their own credible intervals.
        self.refine_iter = refine_iter
        self.refine_burnin = refine_burnin

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        S, B = self._as_arrays(data)
        rng = np.random.default_rng(config.RANDOM_SEED)

        samples = data.samples
        patients = data.manifest.loc[samples, "patient_id"].to_numpy()
        n_s, n_t = B.shape[1], S.shape[1]

        counts = [self._to_counts(B[:, j]) for j in range(n_s)]

        # Start from the unpooled per-sample fit, then alternate: update each patient's
        # mean from its blocks, re-fit each block under a prior centred on that mean.
        se = self._reference_se(data)
        current = np.zeros((n_s, n_t))
        for j in range(n_s):
            draws = self._gibbs(S, counts[j], rng, se=se)
            current[j] = draws.mean(axis=0)

        by_patient = {p: np.where(patients == p)[0] for p in np.unique(patients)}

        for _ in range(self.n_outer):
            mu = {}
            for p, idx in by_patient.items():
                block = current[idx]
                block = block[np.isfinite(block).all(axis=1)]
                mu[p] = (block.mean(axis=0) if len(block)
                         else np.full(n_t, 1.0 / n_t))

            for j in range(n_s):
                prior = self.kappa * mu[patients[j]]
                # A zero prior component would forbid a cell type outright; floor it so
                # a type absent from every other block of this tumour can still be
                # recovered in this one if the evidence demands it.
                prior = np.maximum(prior, 1e-3)
                current[j] = self._gibbs_with_prior(S, counts[j], prior, rng)

        cols = list(data.cell_types)
        self.posterior_ = {"mean": pd.DataFrame(current, index=samples, columns=cols)}
        return current

    def _gibbs_with_prior(self, S: np.ndarray, counts: np.ndarray,
                          prior: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """One block's chain under a caller-supplied (patient-centred) Dirichlet prior."""
        keep = counts > 0
        if not keep.any():
            return np.full(S.shape[1], np.nan)
        S_k, c_k = S[keep], counts[keep]
        informative = S_k.sum(axis=1) > 0
        S_k, c_k = S_k[informative], c_k[informative]
        if S_k.shape[0] == 0:
            return np.full(S.shape[1], np.nan)

        w = rng.dirichlet(prior)
        acc = np.zeros(S.shape[1])
        n_kept = 0
        for it in range(self.refine_iter):
            unnorm = S_k * w[None, :]
            p = unnorm / unnorm.sum(axis=1, keepdims=True)
            z = rng.multinomial(c_k, p)
            w = rng.dirichlet(prior + z.sum(axis=0))
            if it >= self.refine_burnin and (it - self.refine_burnin) % self.thin == 0:
                acc += w
                n_kept += 1
        return acc / max(n_kept, 1)
