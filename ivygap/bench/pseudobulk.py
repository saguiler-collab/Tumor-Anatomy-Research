"""
pseudobulk.py — synthetic mixtures with known composition, used to SELECT the method.

WHY SYNTHETIC MIXTURES AT ALL
-----------------------------
On real Ivy GAP tissue nobody knows the true composition, so "which method is most
accurate" cannot be answered there. Pooling reference cells in proportions we choose
gives mixtures whose truth we do know, and accuracy against that truth is the only
outcome-blind, quantitative basis for choosing a method. The anatomic analysis then
*validates* the chosen method against known biology; it does not select it.

THE SPLIT IS BY DONOR, NEVER BY CELL
------------------------------------
If cells from one donor appear in both the reference and the test mixtures, every
method is being asked to recognise profiles it has already been handed, and the scores
measure memorisation rather than generalisation. Donors are split first; the reference
is rebuilt from training donors only; test mixtures are pooled from held-out donors
only. This is the single most important detail in the file.

ANATOMICALLY SHAPED MIXTURES
----------------------------
Half the mixtures are drawn from a flat Dirichlet, which covers the space evenly. The
other half are drawn around the compositions the five Ivy GAP anatomic structures are
expected to have — leading-edge-like mixtures that are mostly oligodendrocyte with
little tumour, MVP-like mixtures rich in endothelium, and so on. A method can score
well on flat mixtures and still fail in the corner of composition space this project
actually operates in, and that failure is worth surfacing before the real run rather
than after.

Those target compositions come from the same frozen biology used for the anatomic
claims, so nothing new is being invented here; but they are *shapes* for generating
test data, and they are never used to score the real Ivy GAP samples.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ivygap import config

#: Rough composition templates for each anatomic niche, on the primary+sidecar roster.
#: Deliberately coarse: they exist to steer synthetic mixtures into the right region of
#: composition space, not to assert a quantitative truth about Ivy GAP tissue.
NICHE_TEMPLATES: dict[str, dict[str, float]] = {
    "LE":    {"Tumor": 0.15, "Macrophage_Microglia": 0.12, "T_cell": 0.01,
              "NK_cell": 0.005, "B_cell": 0.005, "Endothelial": 0.05,
              "Oligodendrocyte": 0.40, "Astrocyte": 0.26},
    "IT":    {"Tumor": 0.40, "Macrophage_Microglia": 0.18, "T_cell": 0.015,
              "NK_cell": 0.005, "B_cell": 0.005, "Endothelial": 0.06,
              "Oligodendrocyte": 0.22, "Astrocyte": 0.115},
    "CT":    {"Tumor": 0.62, "Macrophage_Microglia": 0.22, "T_cell": 0.02,
              "NK_cell": 0.008, "B_cell": 0.007, "Endothelial": 0.06,
              "Oligodendrocyte": 0.04, "Astrocyte": 0.025},
    "MVP":   {"Tumor": 0.45, "Macrophage_Microglia": 0.20, "T_cell": 0.02,
              "NK_cell": 0.008, "B_cell": 0.007, "Endothelial": 0.26,
              "Oligodendrocyte": 0.03, "Astrocyte": 0.025},
    "PAN":   {"Tumor": 0.52, "Macrophage_Microglia": 0.36, "T_cell": 0.02,
              "NK_cell": 0.008, "B_cell": 0.007, "Endothelial": 0.03,
              "Oligodendrocyte": 0.03, "Astrocyte": 0.025},
}


@dataclass
class PseudobulkSet:
    """Mixtures plus the composition each was built from."""
    expression: pd.DataFrame       # genes x mixtures, linear CPM-like
    truth: pd.DataFrame            # mixtures x cell types, rows sum to 1 (CELL fractions)
    donors: pd.Series              # which held-out donor each mixture drew from
    niche: pd.Series               # the template used, or "flat"

    def __len__(self) -> int:
        return self.expression.shape[1]


def split_donors(meta: pd.DataFrame, test_frac: float = config.PSEUDOBULK_TEST_DONOR_FRAC,
                 seed: int = config.RANDOM_SEED) -> tuple[list[str], list[str]]:
    """Split reference donors into train/test. At least one donor must land in each."""
    donors = sorted(meta["donor"].astype(str).unique())
    if len(donors) < 2:
        raise ValueError(
            f"donor-held-out benchmarking needs at least 2 reference donors, found "
            f"{len(donors)}. With one donor there is no way to distinguish "
            f"generalisation from memorisation."
        )
    rng = np.random.default_rng(seed)
    shuffled = list(rng.permutation(donors))
    n_test = max(1, int(round(len(donors) * test_frac)))
    n_test = min(n_test, len(donors) - 1)          # always leave a training donor
    return sorted(map(str, shuffled[n_test:])), sorted(map(str, shuffled[:n_test]))


def _sample_composition(rng: np.random.Generator, types: list[str],
                        niche: str | None) -> np.ndarray:
    """Draw one target composition, either flat or concentrated around a niche."""
    n = len(types)
    if niche is None:
        w = rng.dirichlet(np.full(n, config.PSEUDOBULK_DIRICHLET_ALPHA))
        # Zero out a few types so structural zeros and limit-of-detection are
        # measurable. Without genuinely absent types, "does this method invent cells
        # that are not there" cannot be asked.
        n_zero = rng.integers(0, config.PSEUDOBULK_MAX_ZEROED + 1)
        if n_zero:
            w[rng.choice(n, size=int(n_zero), replace=False)] = 0.0
            if w.sum() <= 0:
                w = rng.dirichlet(np.full(n, config.PSEUDOBULK_DIRICHLET_ALPHA))
        return w / w.sum()

    template = np.array([NICHE_TEMPLATES[niche].get(t, 0.0) for t in types])
    # Concentration 60 keeps draws recognisably in the niche while still spanning
    # realistic between-patient spread.
    return rng.dirichlet(np.clip(template, 1e-4, None) * 60.0)


def generate(expression: pd.DataFrame, meta: pd.DataFrame, donors: list[str],
             n_mixtures: int, seed: int = config.RANDOM_SEED,
             niche_fraction: float = 0.5) -> PseudobulkSet:
    """
    Pool cells from `donors` into `n_mixtures` mixtures of known composition.

    The truth returned is CELL fractions — the proportion of *cells*, which is what the
    methods report after the shared cell-size correction. Because cell types differ in
    transcript content, the RNA proportions of the same mixture are different numbers;
    scoring cell-fraction estimates against RNA-fraction truth is a classic way to
    manufacture bias that looks like a method difference.
    """
    rng = np.random.default_rng(seed)
    types = list(config.CELL_TYPES)

    meta = meta.loc[meta["donor"].astype(str).isin(donors)]
    if meta.empty:
        raise ValueError(f"no cells for donors {donors}")

    pools: dict[tuple[str, str], list[str]] = {}
    for donor in donors:
        for ctype in types:
            sel = meta.index[(meta["donor"].astype(str) == donor)
                             & (meta["cell_type"] == ctype)]
            if len(sel):
                pools[(donor, ctype)] = list(sel)

    niches = list(NICHE_TEMPLATES)
    cols, truths, used_donor, used_niche = {}, [], [], []

    for i in range(n_mixtures):
        donor = str(rng.choice(donors))
        use_niche = rng.random() < niche_fraction
        niche = str(rng.choice(niches)) if use_niche else None
        target = _sample_composition(rng, types, niche)

        total_cells = int(rng.integers(*config.PSEUDOBULK_TOTAL_CELLS))
        counts = rng.multinomial(total_cells, target)

        picked: list[str] = []
        realised = np.zeros(len(types))
        for k, ctype in enumerate(types):
            pool = pools.get((donor, ctype))
            if not pool or counts[k] == 0:
                continue
            take = rng.choice(pool, size=int(counts[k]), replace=True)
            picked.extend(take)
            realised[k] = counts[k]

        if not picked:
            continue

        # Score against the composition actually realised, not the one requested. A
        # donor missing a cell type cannot contribute it, and holding a method to a
        # target the mixture does not contain would penalise every method for the
        # sampler's shortfall.
        realised = realised / realised.sum()

        mix = expression[picked].sum(axis=1)
        total = mix.sum()
        cols[f"mix{i:05d}"] = (mix / total * 1e6) if total > 0 else mix
        truths.append(realised)
        used_donor.append(donor)
        used_niche.append(niche or "flat")

    if not cols:
        raise ValueError("no mixtures could be built from the supplied donors")

    names = list(cols)
    return PseudobulkSet(
        expression=pd.DataFrame(cols, index=expression.index),
        truth=pd.DataFrame(np.vstack(truths), index=names, columns=types),
        donors=pd.Series(used_donor, index=names, name="donor"),
        niche=pd.Series(used_niche, index=names, name="niche"),
    )


def build_train_test(expression: pd.DataFrame, meta: pd.DataFrame,
                     n_test: int = config.N_TEST_PSEUDOBULK,
                     seed: int = config.RANDOM_SEED
                     ) -> tuple[PseudobulkSet, list[str], list[str]]:
    """
    The benchmark's standard setup: split donors, build test mixtures from the held-out
    ones. Returns (test_set, train_donors, test_donors); the caller rebuilds the
    reference from `train_donors` so no test donor is ever inside the reference.
    """
    train_donors, test_donors = split_donors(meta, seed=seed)
    test_set = generate(expression, meta, test_donors, n_test, seed=seed + 1)
    return test_set, train_donors, test_donors
