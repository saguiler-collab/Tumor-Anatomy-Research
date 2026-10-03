# BayesPrism, authors' configuration vs this project's driver — the prediction and its verdict rule

**Written 2026-10-01 15:20 EDT. At this moment no output of either configuration exists**:
`results/extension/` holds only the authors' run's R log, whose first Gibbs stage is still
sampling (ETA 17:37), and the pipeline configuration has not started. Nothing in this repository
has been committed this session, so the evidence of order is this file's timestamp against the
artefacts' — check `ls -l --time-style=full-iso prespecified/bayesprism_authors_prediction.md
results/extension/bayesprism_tcga_*.json`.

## History, stated plainly

The prediction was first written into the docstring of `scripts/bayesprism_tcga.py`, before the
GBM authors' run started (11:14). That file was edited again at 12:04, after the run had started,
for an unrelated memory change, so its modification time cannot prove the order, and its
docstring cited "WHY_B_OVER_T §7i", a section that did not exist yet and now covers something
else. This file fixes the prediction **and** the exact rule for reading it, before any number
exists.

## The prediction (unchanged from the docstring)

If B acts as an overflow for tumour-specific expression the reference cannot fit, the authors'
configuration (`key = "Tumor"`: BayesPrism learns each tumour's own malignant expression; plus
`cleanup.genes` and protein-coding genes only) should **LOWER the B share of the lymphoid
compartment** relative to this project's driver (`key = NULL`, no gene cleanup).

**Premise, qualified before the result (found 2026-10-01 14:10, `results/b_column_diagnostics.json`).**
The overflow observation is estimator-specific. Removing the B column under SVR sends 92% of B's
estimate to Tumor (GBM), but under NNLS it sends the mass to T. The prediction stands as declared.
Its premise rests on SVR alone, and a failed prediction would be read in that light.

## Verdict rule — fixed now

- **Quantity.** `lymphoid_mean["B_cell"]` as `scripts/bayesprism_tcga.py` writes it: the mean,
  over methylation-matched samples with any lymphoid estimate, of B / (T + NK + B).
- **Comparison.** authors vs pipeline, same cohort, same samples, same reference cells (both runs
  use the registered raw/X arm's inputs). GBM first; LGG only if both LGG runs complete.
- **SUPPORTED** if the authors' B share is lower than the pipeline's; **NOT SUPPORTED** otherwise.
  Exact equality reads as NOT SUPPORTED.
- **Reported alongside, not part of the verdict:** the per-sample paired difference in B share
  (Wilcoxon signed-rank) on samples with lymphoid signal in both runs, `frac_samples_B_over_T`,
  `n_zero_lymphoid`, the purity correlation, and whether either configuration orders T above B.
  If fewer than 10 samples carry lymphoid signal in both runs, the per-sample comparison is
  **INCONCLUSIVE** and only the cohort-level number is read.
- **What it cannot show.** Support makes the overflow account more plausible; it does not establish
  it. Neither outcome bears on the registered statistics, which this run does not enter.
