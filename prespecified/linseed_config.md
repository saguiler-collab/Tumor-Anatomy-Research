# Linseed -- configuration fixed before any output

**Written 2026-10-03.** The time is this file's birth time (`stat`); no clock time is typed here. At
writing, Linseed was not installed and had produced nothing. Post-registration, exploratory: the
reference-free arm's second instrument, beside CDSeq. It selects nothing and changes no registered
statistic.

## Why Linseed

- **Nguyen et al. 2024 [3]:** "Linseed and Deblender consistently have the highest accuracy scores in
  their respective category (reference-free and semi-reference-free)."
- The user asked for it earlier in the project, and said "continue" when it was next.
- CDSeq is the only reference-free method run so far (`scripts/cdseq_anatomic.py`). A second, of a
  different kind -- Linseed finds simplex corners of mutually linear genes, CDSeq is a probabilistic
  topic model -- tests whether the reference-free reading depends on the instrument.

## The package

- linseed 0.99.3, genuine and unmodified, vendored from github.com/ctlab/LinSeed @ 2435a8ce
  (2025-08-22; MIT, LICENSE file present). Built from `vendor/linseed/`, whose C++ code is compiled
  through Rcpp/RcppArmadillo.
- Missing imports are installed first: GEOquery (Bioconductor) and combinat (CRAN), with
  `update = FALSE`.
- **Citation, Crossref-verified:** Zaitsev K, Bambouskova M, Swain A, Artyomov MN. Nat Commun
  2019;10:2209; doi:10.1038/s41467-019-09990-5.

## The pipeline: the package README's tutorial, step for step

```
lo <- LinseedObject$new(matrix, topGenes = 10000)
lo$calculatePairwiseLinearity(); lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100); lo$filterDatasetByPval(0.01)
lo$setCellTypeNumber(k); lo$project("filtered")
lo$smartSearchCorners(dataset = "filtered", error = "norm"); lo$deconvolveByEndpoints()
```

- The seed is `config.RANDOM_SEED`.
- **k = 8, the roster size, fixed a priori -- exactly CDSeq's T.** The tutorial chooses k by looking
  at an SVD plot. A visual choice cannot be pre-declared, so the same fixed rule as the existing
  reference-free arm replaces it.
- **topGenes stays at the tutorial's 10,000.** That needs memory this 8 GB machine has only when no
  other heavy job is running, so the real runs wait for ReCIDE TCGA-LGG to finish. The setting is not
  lowered to fit.

## Inputs (anatomic arm)

The 122 Ivy GAP anatomic samples, as CPM from the per-sample RSEM expected counts -- the source data
CDSeq used -- on gene symbols. No reference, no marker list, and no atlas enter the deconvolution.

## Naming and scoring: CDSeq's, reused unchanged

- **Labelling (a):** correlation of each Linseed signature with the GBmap profile. The atlas is used
  for labelling only.
- **Labelling (b):** enrichment of GBM-specific marker sets (GBMdeconvoluteR/Moreno). No atlas.
- Components mapped to one roster type are summed. A type no component claims is absent (NaN), never
  zero.
- The registered ACS scorer, 10,000 permutations, per labelling.

## Gate before real data (`tests/test_linseed_gate.py`)

- **Planted:** 4 roster profiles from the frozen signature (Tumor, Macrophage_Microglia,
  Oligodendrocyte, Endothelial), 40 samples, Dirichlet proportions, lognormal noise sd 0.05. Linseed
  with k = 4 must recover them: each true type matched to an estimated component by maximum total
  correlation (Hungarian), with minimum matched r >= 0.80. The bar is looser than the
  reference-based 0.90 because a reference-free method must also find the profiles.
- **Negative control:** the same matrix with each gene's values permuted independently across
  samples, which destroys the mixing structure, must FAIL the bar. An error or a missing estimate
  counts as failing.
- **Reading:** if the planted set fails, Linseed does not proceed to real data, and that is reported.

## Real runs (after the gate, when ReCIDE TCGA-LGG has finished)

1. **Ivy GAP anatomic:** ACS under both labellings, with nulls and coverage, set beside CDSeq's (0.846
   profile-labelled, 0.683 marker-labelled). Descriptive.
2. **TCGA-GBM and TCGA-LGG, bulk only:** the same pipeline and labellings.
   - The tumour-labelled component against ABSOLUTE purity (Spearman).
   - The lymphoid ordering against methylation, where lymphoid types are claimed.
   - A reference-free method scored against DNA truth. Descriptive.

## Addendum (written before any Linseed output, after reading `deconvolveByEndpoints` in the source)

`deconvolveByEndpoints` returns signatures W on Linseed's row-normalised scale (each gene divided by
its total across samples) and proportions H (k x samples). CDSeq's labelling functions were written
for gene-expression profiles. Linseed's signatures are therefore passed to both labellings on the
gene-expression scale: each gene's row of W multiplied back by that gene's total in the input. This
is the inverse of the package's own normalisation. Proportions are used exactly as returned.
