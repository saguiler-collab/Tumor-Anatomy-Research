"""
base.py — the contract every deconvolution method in this project obeys.

WHY A CONTRACT RATHER THAN NINE SCRIPTS
---------------------------------------
The point of this rebuild is to compare methods, and a method comparison is only
worth reading if the methods were given genuinely identical problems. It is very easy
to accidentally hand one method a slightly different gene set, or a log-scaled matrix,
or a reference built with a different cell cap, and then report the difference as a
finding about the algorithm.

So methods here never load their own data. They receive a frozen `DeconvolutionInput`
that is constructed once, hashed, and passed to all of them. A method that mutates its
input trips an assertion. The hashes travel into the equal-footing certificate that
`ivygap.bench.equal_footing` requires before any score may be computed.

THE ESTIMAND
------------
Every method returns the same thing: **cell fractions**, rows summing to 1.0, over
`config.CELL_TYPES` in that exact order.

That word matters. Raw deconvolution against an expression signature recovers each
type's share of the *RNA*, not its share of the *cells* — a large, transcriptionally
active tumour cell contributes more mRNA than a small lymphocyte. Converting RNA
proportion to cell proportion requires dividing by a per-type cell-size factor. The
conversion is applied once, centrally, in `to_cell_fractions()`, so no method can
quietly report a different estimand from its neighbours.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ivygap import config


# =============================================================================
# INPUTS
# =============================================================================

@dataclass(frozen=True)
class ReferenceBundle:
    """
    Everything a reference-based method may ask of a single-cell dataset.

    Different methods need different summaries of the same reference, and building
    them separately per method is how references drift apart:

        profile    genes x types   mean expression per type (the signature matrix S)
        sigma      genes x types   CROSS-SUBJECT variance of each gene within each
                                   type. MuSiC's central idea: a gene whose expression
                                   in macrophages swings wildly from donor to donor is
                                   a bad basis for estimating macrophage content in a
                                   new donor, and should be downweighted.
        cell_size  types           mean total mRNA per cell of each type, used to turn
                                   RNA proportions into cell proportions.
        donor_profiles  optional {donor -> genes x types}, needed by SCDC's outlier
                                   handling and by Bisque's subject-level transform.
    """
    name: str
    profile: pd.DataFrame
    sigma: pd.DataFrame
    cell_size: pd.Series
    n_donors: int
    donor_profiles: dict[str, pd.DataFrame] | None = None
    #: False when the reference was loaded already collapsed, so `sigma` is a
    #: placeholder rather than a measurement. MuSiC checks this and reports that it is
    #: running degenerate rather than claiming a MuSiC result it cannot produce.
    has_cross_donor_variance: bool = True

    def __post_init__(self) -> None:
        if list(self.profile.columns) != list(self.sigma.columns):
            raise ValueError(f"{self.name}: profile and sigma have different cell types")
        if not self.profile.index.equals(self.sigma.index):
            raise ValueError(f"{self.name}: profile and sigma have different genes")
        if list(self.cell_size.index) != list(self.profile.columns):
            raise ValueError(f"{self.name}: cell_size does not match the cell types")
        if (self.cell_size <= 0).any():
            raise ValueError(f"{self.name}: cell_size must be strictly positive")

    @property
    def cell_types(self) -> list[str]:
        return list(self.profile.columns)

    def subset_genes(self, genes: list[str]) -> "ReferenceBundle":
        return ReferenceBundle(
            name=self.name,
            profile=self.profile.loc[genes],
            sigma=self.sigma.loc[genes],
            cell_size=self.cell_size,
            n_donors=self.n_donors,
            donor_profiles=(
                {d: p.loc[genes] for d, p in self.donor_profiles.items()}
                if self.donor_profiles else None
            ),
            has_cross_donor_variance=self.has_cross_donor_variance,
        )

    def content_hash(self) -> str:
        return config.sha256_frame(
            pd.concat([self.profile, self.sigma], axis=1)
        )


@dataclass(frozen=True)
class DeconvolutionInput:
    """
    The frozen problem statement handed to every method.

    `bulk` is genes x samples on a linear CPM-like scale, restricted to the shared
    gene space. `references` is ordered; `references[0]` is the primary one, and
    single-reference methods use exactly that one so their comparison against the
    multi-reference methods is fair.
    """
    bulk: pd.DataFrame
    references: tuple[ReferenceBundle, ...]
    manifest: pd.DataFrame
    cell_types: tuple[str, ...] = field(default=tuple(config.CELL_TYPES))
    apply_cell_size_correction: bool = True
    #: The same samples on the FULL shared gene space, before marker selection narrowed
    #: it. Optional, and used by exactly one kind of method: one that brings its own
    #: published signature instead of solving against this project's reference.
    #:
    #: The marker subset exists to serve THIS project's signature matrix — top-N genes
    #: per cell type ranked against it, plus its markers. For every method that solves
    #: against that signature, it is the shared gene space and equal footing holds.
    #: quanTIseq does not solve against it: it ships TIL10 and ignores the signature it
    #: is handed. Ranking genes against a different reference kept 34 of TIL10's 138
    #: signature genes, and quanTIseq answered that with macrophages at 100 percent of
    #: cells in the median sample and T cells identically zero in all 122. On the full
    #: space it finds 136 of 138 and returns a GBM-plausible 15.7 percent median
    #: myeloid fraction. Handing it the subset was not equal footing; it was mutilation.
    bulk_full: pd.DataFrame | None = None

    def __post_init__(self) -> None:
        if not self.references:
            raise ValueError("at least one reference is required")
        for ref in self.references:
            if ref.cell_types != list(self.cell_types):
                raise ValueError(
                    f"reference {ref.name!r} has cell types {ref.cell_types}, "
                    f"expected {list(self.cell_types)}"
                )
            if not ref.profile.index.equals(self.bulk.index):
                raise ValueError(
                    f"reference {ref.name!r} gene space differs from the bulk gene space "
                    f"({len(ref.profile.index)} vs {len(self.bulk.index)} genes). "
                    "Align once, up front — never inside a method."
                )
        if (self.bulk.to_numpy() < 0).any():
            raise ValueError("bulk contains negative values; the mixing model is additive")
        missing = set(self.bulk.columns) - set(self.manifest.index)
        if missing:
            raise ValueError(f"{len(missing)} bulk samples are absent from the manifest")
        if self.bulk_full is not None:
            if list(self.bulk_full.columns) != list(self.bulk.columns):
                raise ValueError(
                    "bulk_full must carry the same samples in the same order as bulk "
                    f"({len(self.bulk_full.columns)} vs {len(self.bulk.columns)}). It is "
                    "a wider GENE space, never a different sample set."
                )
            if not set(self.bulk.index).issubset(set(self.bulk_full.index)):
                raise ValueError(
                    "bulk_full must be a superset of bulk's genes; it is the space "
                    "before marker selection, not an unrelated matrix."
                )

    @property
    def primary(self) -> ReferenceBundle:
        return self.references[0]

    @property
    def genes(self) -> list[str]:
        return list(self.bulk.index)

    @property
    def samples(self) -> list[str]:
        return list(self.bulk.columns)

    def content_hashes(self) -> dict[str, str]:
        """Everything the equal-footing certificate needs to prove inputs matched."""
        return {
            "bulk": config.sha256_frame(self.bulk),
            "cell_types": ",".join(self.cell_types),
            "n_genes": str(len(self.genes)),
            "n_samples": str(len(self.samples)),
            "references": ";".join(f"{r.name}:{r.content_hash()[:16]}"
                                   for r in self.references),
            "cell_size_correction": str(self.apply_cell_size_correction),
        }


# =============================================================================
# SHARED NUMERICS
# =============================================================================

def project_to_simplex(w: np.ndarray, covered: np.ndarray | None = None) -> np.ndarray:
    """
    Clip negatives and rescale to sum 1.

    An all-zero solution (every coefficient shrunk away, which heavy L1 can do) has no
    meaningful composition. Returning a uniform vector would be a silent fabrication,
    so we return NaN and let the caller decide — the benchmark counts NaN rows as
    failures rather than scoring them.

    PARTIAL COVERAGE (`covered`)
    ----------------------------
    A NaN normally means "this solve failed", and one NaN condemns the whole row: the
    total is non-finite, so every entry comes back NaN. That is right for a method that
    models all eight roster types and right for the failed samples DWLS reports.

    It is wrong for a method that structurally does not model some of them. quanTIseq
    estimates ten immune populations and puts everything else in an "Other" compartment
    it does not break down, so Tumor, Endothelial, Oligodendrocyte and Astrocyte are NaN
    by construction, in every sample, however well the method worked. Passed through the
    branch above, those four NaNs annihilated the four columns quanTIseq DID estimate,
    and the run recorded `R:quantiseqr` with no fallback reason — a total loss reported
    as a success. `covered` is what keeps the two meanings apart.

    Two things are deliberately NOT done on the covered path:

      * The covered entries are NOT rescaled to sum 1. quanTIseq's fractions are already
        shares of all cells, with the remainder in "Other"; rescaling four immune types
        to sum 1 would assert the tumour is 100 percent immune. Its own scale is the
        honest one, and ordinal constraints — which is all ACS reads — are invariant to
        the monotone rescaling anyway.
      * A NaN among the COVERED entries still condemns the row, because there it carries
        its original meaning: that solve failed.
    """
    w = np.clip(np.asarray(w, dtype="float64"), 0.0, None)
    if covered is None:
        total = w.sum()
        if not np.isfinite(total) or total <= 0:
            return np.full_like(w, np.nan)
        return w / total

    out = np.full_like(w, np.nan)
    sub = w[covered]
    total = sub.sum()
    if not np.isfinite(total) or total <= 0:
        return out
    out[covered] = sub
    return out


def to_cell_fractions(rna_fractions: np.ndarray, cell_size: np.ndarray) -> np.ndarray:
    """
    RNA proportion -> cell proportion.

    A type's share of the mRNA pool is (cells_k * mRNA_per_cell_k) / total, so dividing
    the RNA share by mRNA_per_cell recovers a quantity proportional to cell count.
    Renormalising then gives cell fractions. This is the single place the conversion
    happens for every method in the project.
    """
    out = np.asarray(rna_fractions, dtype="float64") / np.asarray(cell_size, dtype="float64")
    return project_to_simplex(out)


# =============================================================================
# THE METHOD INTERFACE
# =============================================================================

class DeconvolutionMethod(abc.ABC):
    """
    Base class. Subclasses implement `_solve_all`, returning RNA proportions; the
    base class handles simplex projection, cell-size correction, column ordering and
    the input-immutability check.
    """

    #: short, stable identifier used in filenames and result tables
    name: str = "base"
    #: "least-squares" | "regularized" | "bayesian" | "reference-based"
    family: str = "unspecified"
    #: True if the method shells out to R; the registry can then skip it gracefully
    requires_r: bool = False
    #: True if the method genuinely consumes more than one reference
    uses_multiple_references: bool = False
    #: The roster cell types this method is capable of estimating at all, or None when
    #: it models every one of them. Set this ONLY for a method whose published design
    #: omits a type — quanTIseq's TIL10 signature has no tumour, endothelial, glial or
    #: astrocyte population and rolls all of them into one "Other" compartment. It is
    #: not a place to record that a method scored badly on a type, and it is not a
    #: licence to drop a type a method merely estimates poorly.
    models_cell_types: frozenset[str] | None = None

    def __init__(self, **params):
        self.params = params

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r})"

    @abc.abstractmethod
    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        """Return an (n_samples, n_cell_types) array of RNA proportions."""

    def degradation_for(self, data: DeconvolutionInput) -> tuple[bool, str | None]:
        """
        Whether this method is running degraded ON THIS DATA, and why.

        Answerable without solving, because degradation is a property of the inputs, not
        of the implementation. That distinction is load-bearing: Bisque's no-overlap mode
        is degraded whether the genuine R package or this project's reimplementation
        runs it, and a flag set only inside the Python solver silently disappears the
        moment the real package starts working.

        Default: not degraded. Overridden by the methods that can be.
        """
        return False, None

    def fit_predict(self, data: DeconvolutionInput) -> pd.DataFrame:
        """
        Run the method and return a samples x cell_types table of CELL fractions,
        rows summing to 1.0, columns in `config.CELL_TYPES` order.
        """
        # Validated BEFORE solving: a typo in the declaration would otherwise surface
        # only after the method has spent its full runtime.
        covered = self._covered_mask(data.cell_types)

        before = config.sha256_frame(data.bulk)

        raw = np.asarray(self._solve_all(data), dtype="float64")

        after = config.sha256_frame(data.bulk)
        if before != after:
            raise RuntimeError(
                f"{self.name} mutated the shared bulk matrix. Every method receives the "
                "same object; mutating it would corrupt every method that runs after."
            )

        n_s, n_t = len(data.samples), len(data.cell_types)
        if raw.shape != (n_s, n_t):
            raise ValueError(
                f"{self.name} returned shape {raw.shape}, expected {(n_s, n_t)}"
            )

        cell_size = data.primary.cell_size.reindex(list(data.cell_types)).to_numpy()

        # Cell-size correction renormalises, so on a partial-coverage method it would
        # rescale a handful of immune types to sum 1 and assert the tumour is entirely
        # immune. quanTIseq also already applies its own mRNA scaling (scale_mRNA=TRUE,
        # its published default), so ours would be the second correction, not the first.
        apply_cs = data.apply_cell_size_correction and covered is None

        rows = []
        for i in range(n_s):
            w = project_to_simplex(raw[i], covered)
            if apply_cs and np.isfinite(w).all():
                w = to_cell_fractions(w, cell_size)
            rows.append(w)

        out = pd.DataFrame(np.vstack(rows), index=data.samples,
                           columns=list(data.cell_types))
        out.index.name = "sample_id"
        return out

    def _covered_mask(self, cell_types) -> np.ndarray | None:
        """
        Boolean mask over `cell_types`, or None when the method models all of them.

        Raises if the declaration names a type outside the roster: a typo there would
        silently widen the NaN set and quietly shrink what the method is scored on.
        """
        if self.models_cell_types is None:
            return None
        roster = list(cell_types)
        unknown = set(self.models_cell_types) - set(roster)
        if unknown:
            raise ValueError(
                f"{self.name}.models_cell_types names {sorted(unknown)}, which are not "
                f"in the roster {roster}"
            )
        if not self.models_cell_types:
            raise ValueError(f"{self.name} declares it models no cell type at all")
        return np.array([c in self.models_cell_types for c in roster], dtype=bool)

    # -- helpers available to every subclass ---------------------------------

    @staticmethod
    def _as_arrays(data: DeconvolutionInput,
                   ref: ReferenceBundle | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(S, B) as float64 arrays: S is genes x types, B is genes x samples."""
        ref = ref or data.primary
        return (ref.profile.to_numpy(dtype="float64"),
                data.bulk.to_numpy(dtype="float64"))
