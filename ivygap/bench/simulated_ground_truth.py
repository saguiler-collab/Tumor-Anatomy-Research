"""
simulated_ground_truth.py — a ground-truth ranking when no cells are available.

WHY THIS EXISTS, AND WHAT IT IS NOT
-----------------------------------
The agreement test — the study's primary result — correlates the ACS ranking against a
ranking from real ground truth. The protocol names three yardsticks, and on this machine
none of them currently covers enough methods to support a rank correlation:

    synthetic_mixtures   the frozen 500 mixtures exist, but scored only NNLS and SVR.
                         Two methods; Spearman needs at least six.
    absolute_purity      the source project records this as BLOCKED — no ABSOLUTE purity
                         file was ever located.
    sc_pseudobulk        needs paired single-cell data, not present.

So the primary result is genuinely not computable from the real yardsticks, and
`agreement.py` reports exactly that. This module does NOT change that. It is a fourth,
separately named yardstick — `simulated_donor_mismatch` — that makes the machinery
runnable end to end and gives a preliminary reading, and it is labelled as a simulation
everywhere it appears. It must never be reported as the protocol's yardstick 1.

WHY IT IS NOT CIRCULAR
----------------------
The obvious shortcut is to build mixtures directly from the frozen signature and score
methods against it. That is worthless: every method solving S·x = b against the same S
that generated b recovers the answer nearly exactly, all methods tie near zero error, and
the ranking is noise. This is precisely the circularity the protocol names as yardstick
1's weakness, in its most extreme form.

Instead each synthetic donor gets its OWN perturbed signature, and the methods solve
against the unperturbed one. The dominant error is therefore reference mismatch — the
same thing that dominates on real held-out donors — rather than sampling noise. That is
what makes the resulting ranking informative about method robustness rather than about
floating-point arithmetic.

THE GENERATIVE MODEL, AND THE ESTIMAND
--------------------------------------
The frozen signature's columns each sum to about 1e6: they are the expression profile of
an average cell of each type, normalised. `cell_size` (L) is a separate quantity — total
mRNA per cell — so a cell of type k contributes L_k transcripts distributed as
profile_k / colsum_k:

    b  =  sum_k  n_k * L_k * profile_dk                     + Poisson depth noise
    truth_k = n_k / sum_j n_j                               (CELL fractions)

The coefficient on each column is n_k * L_k, because that is exactly what a method
estimates before the shared cell-size correction divides it by L_k. Any other scaling
makes the yardstick measure a bookkeeping mismatch rather than method quality;
`test_zero_perturbation_recovers_truth` pins it by requiring near-perfect recovery when
the donor perturbation is switched off.

Truth is recorded as cell fractions because that is what every method reports after the
shared cell-size correction. Scoring cell-fraction estimates against RNA-fraction truth
manufactures a bias that looks exactly like a method difference.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.bench.pseudobulk import NICHE_TEMPLATES
from ivygap.deconv.base import ReferenceBundle

#: Multiplicative cross-donor spread applied per gene per cell type, on a log scale.
#: 0.35 puts a typical gene within roughly +/-40% of the reference in a given donor,
#: which is the order of magnitude of real between-donor variation in GBM single-cell
#: profiles. The value is a modelling choice and is recorded with the results, because
#: it directly controls how hard the yardstick is.
DONOR_LOG_SD = 0.35

#: Total transcripts per synthetic sample. Sets the sequencing-depth noise floor.
LIBRARY_SIZE = 5e6


@dataclass
class SimulatedSet:
    """Mixtures, their true CELL fractions, and which synthetic donor each came from."""
    expression: pd.DataFrame       # genes x mixtures
    truth: pd.DataFrame            # mixtures x cell types, rows sum to 1
    donors: pd.Series
    niche: pd.Series
    donor_log_sd: float
    n_donors: int

    def __len__(self) -> int:
        return self.expression.shape[1]

    def provenance(self) -> dict:
        return {
            "yardstick": "simulated_donor_mismatch",
            "is_simulation": True,
            "warning": "NOT the protocol's yardstick 1. Mixtures are generated from "
                       "donor-perturbed copies of the frozen signature, not pooled from "
                       "real single cells. Use for machinery validation and a "
                       "preliminary reading only.",
            "n_mixtures": len(self),
            "n_donors": self.n_donors,
            "donor_log_sd": self.donor_log_sd,
            "library_size": LIBRARY_SIZE,
        }


def _perturb(profile: np.ndarray, rng: np.random.Generator,
             log_sd: float) -> np.ndarray:
    """
    One synthetic donor's signature.

    Multiplicative lognormal noise per gene per cell type, then each column renormalised
    to its original total. Renormalising matters: without it a donor whose noise happened
    to run high would also have systematically larger cells, confounding reference
    mismatch with cell size and making the yardstick measure two things at once.
    """
    noise = rng.lognormal(mean=0.0, sigma=log_sd, size=profile.shape)
    out = profile * noise
    original = profile.sum(axis=0)
    new = out.sum(axis=0)
    scale = np.divide(original, new, out=np.ones_like(new), where=new > 0)
    return out * scale


def generate(reference: ReferenceBundle, n_mixtures: int = 500, n_donors: int = 12,
             donor_log_sd: float = DONOR_LOG_SD, niche_fraction: float = 0.5,
             seed: int = config.RANDOM_SEED) -> SimulatedSet:
    """Build the simulated yardstick from a (collapsed) reference bundle."""
    rng = np.random.default_rng(seed)
    types = list(reference.profile.columns)
    genes = list(reference.profile.index)

    profile = reference.profile.to_numpy(dtype="float64")
    L = reference.cell_size.reindex(types).to_numpy(dtype="float64")

    donor_profiles = [_perturb(profile, rng, donor_log_sd) for _ in range(n_donors)]
    niches = list(NICHE_TEMPLATES)

    cols, truths, used_donor, used_niche = {}, [], [], []
    for i in range(n_mixtures):
        d = int(rng.integers(n_donors))
        use_niche = rng.random() < niche_fraction
        niche = str(rng.choice(niches)) if use_niche else None

        if niche is None:
            w = rng.dirichlet(np.full(len(types), config.PSEUDOBULK_DIRICHLET_ALPHA))
            n_zero = int(rng.integers(0, config.PSEUDOBULK_MAX_ZEROED + 1))
            if n_zero:
                w[rng.choice(len(types), size=n_zero, replace=False)] = 0.0
                if w.sum() <= 0:
                    w = rng.dirichlet(np.full(len(types), 1.0))
            w = w / w.sum()
        else:
            template = np.array([NICHE_TEMPLATES[niche].get(t, 0.0) for t in types])
            w = rng.dirichlet(np.clip(template, 1e-4, None) * 60.0)

        total_cells = int(rng.integers(*config.PSEUDOBULK_TOTAL_CELLS))
        n_cells = rng.multinomial(total_cells, w).astype("float64")
        if n_cells.sum() <= 0:
            continue

        # Mixture from THIS donor's signature; the methods solve against the
        # unperturbed one, so reference mismatch is the dominant error.
        #
        # The coefficient on each column must be exactly what the solver estimates, or
        # the yardstick measures a bookkeeping mismatch instead of method quality. A
        # method fits b = S @ x and the base class reports x_k / L_k as the cell
        # fraction, so recovering n_k requires x_k = n_k * L_k, hence
        # b = sum_k n_k * L_k * S_k. An earlier version divided the columns by their
        # totals here; since those totals vary by up to 1.33x across cell types, every
        # method inherited a systematic error of that size and the zero-perturbation
        # control scored 0.057 MAE instead of ~0.
        mix = donor_profiles[d] @ (n_cells * L)

        total = mix.sum()
        if total <= 0:
            continue
        # Poisson thinning at a fixed library size: sequencing depth noise on top of
        # reference mismatch, so a method cannot win purely by being noise-free.
        counts = rng.poisson(mix / total * LIBRARY_SIZE)
        cpm = counts / max(counts.sum(), 1) * 1e6

        name = f"sim{i:05d}"
        cols[name] = cpm
        truths.append(n_cells / n_cells.sum())
        used_donor.append(f"D{d:02d}")
        used_niche.append(niche or "flat")

    if not cols:
        raise ValueError("no mixtures could be generated")

    names = list(cols)
    return SimulatedSet(
        expression=pd.DataFrame(cols, index=genes),
        truth=pd.DataFrame(np.vstack(truths), index=names, columns=types),
        donors=pd.Series(used_donor, index=names, name="donor"),
        niche=pd.Series(used_niche, index=names, name="niche"),
        donor_log_sd=donor_log_sd,
        n_donors=n_donors,
    )
