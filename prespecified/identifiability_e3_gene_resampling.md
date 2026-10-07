# Extension E3 -- gene-resampling stability: binding rules

**Written 2026-10-07 15:16:26 EDT** (the file's birth time, copied from stat), before any gene-subset fit
was run. Exploratory and post-registration; it selects nothing and changes no registered statistic.

## What had been seen at writing

- E2's results (`prespecified/identifiability_diagnostics.md`, `results/identifiability_diagnostics.json`):
  - loss-scale stability against truth agreement across the 6 primary units: rho 0.886 (exact p 0.017);
    0.943 on matched denominators;
  - agreement between methods: 0.657;
  - the per-unit stabilities and truth agreements.
- The repaired-method and evaluation-matrix results of 2026-10-07.
- **No fit on a gene subset had been run.**

## Why

E2's diagnostic rests on a single perturbation (the loss scale) and 6 units.

- If identifiability belongs to the compartment, a different perturbation should order the compartments the
  same way.
- That perturbation here is which genes inform the fit.
- If it does not, E2's result is a property of the loss scale, not of the data. The paper must then say so.

## Inputs (fixed)

- **Fits:** genuine DESeq2 `unmix` on the TCGA frozen arm, with `scripts/extension_tcga.py`'s inputs unchanged
  except for the gene set. GBM and LGG. The pre-declared setting only: `IVYGAP_UNMIX_SHIFT` and `_POWER` unset,
  as E2's "predeclared".
- **Gene subsets:**
  - K = 10 random subsets, each half of the marker genes, drawn without replacement;
  - drawn with `numpy.random.default_rng(config.RANDOM_SEED + 3)` from the GBM marker list in its stored
    order;
  - the same 10 subsets applied to LGG, each intersected with LGG's marker list.
- **Truth agreement:** E2's, read from `results/identifiability_diagnostics.json` (`truth_rho_unmix`: `unmix`
  at the pre-declared setting on all genes). Not refit.

## Statistics

- **D3 stability** (per compartment × cohort, E2's compartments) is the mean of the 45 pairwise Spearman
  correlations of per-sample estimates between the 10 subset fits.
- A fit that returns a constant for a compartment is excluded from that compartment's pairs. The number
  excluded is reported. With fewer than 3 usable fits, the unit is "not estimable".

## Hypotheses and reading rules (fixed now)

- **H1-E3 (primary).** Across E2's 6 primary units (Tumor, Leukocytes and Lymphoid, in GBM and LGG): Spearman
  between D3 stability and truth agreement, with an exact one-sided permutation p over 720 orderings.
  - **SUPPORTED:** rho > 0 and p < 0.05.
  - **NOT SUPPORTED:** rho <= 0.
  - **INCONCLUSIVE:** otherwise.
- **H1s-E3 (secondary).** The same over E2's 12 units (adding T, B and NK), with a Monte Carlo permutation
  (100,000 orderings).
- **C1 (descriptive).** Do the two perturbations agree? Spearman across the 6 primary units between D1
  (loss-scale) and D3 (gene-resampling) stability.

## Controls

- **E2's truth agreement must be read, not recomputed.** Its Tumor values must equal E2's reported `unmix`
  rows: GBM 0.7801, LGG 0.5916.
- **Each subset fit's tumour rho against ABSOLUTE is reported.** If the median subset fit falls below 0.40,
  the subsets have destroyed the signal, and H1-E3 is reported as uninformative.

## Limits, stated now

- Six primary units.
- One package (`unmix`).
- Exploratory; registered after E2's result was known.

## Implementation

- `scripts/identifiability_e3.py` writes `results/identifiability_e3.json` and per-fit estimates under
  `results/identifiability/e3/`.
- `scripts/extension_tcga.py` gains an optional `genes_keep` argument (default: all genes, which is unchanged
  behaviour).
