# What makes a tumour easier for deconvolution methods to reproduce anatomically?

**Written 2026-10-02, before any factor below was computed.** The per-tumour ACS table
(`results/anatomic/acs_per_tumor.csv`) has been seen (it is a registered output); no association
between it and any factor listed here has been computed. Exploratory, post-registration. It selects
nothing and changes no registered statistic.

## What this can and cannot answer

ACS measures agreement with anatomic expectation, not accuracy (`docs/WHAT_ACS_IS.md`). A factor that
raises per-tumour ACS makes the anatomy *easier to reproduce*. It does not make an estimate more
correct, and must never be reported as a driver of deconvolution accuracy. n = 9 evaluable tumours.

## Unit and outcome

- **Unit:** the tumour (the 9 with at least one evaluable constraint). Samples are never pooled across tumours.
- **Outcome: ACS excess.** For each comparable method (registered panel, excluding the two controls,
  partial-coverage quanTIseq and the degenerate duplicate SCDC ENSEMBLE), per-tumour ACS minus that
  tumour's chance ACS for that method. Chance is the mean of 2,000 within-tumour permutations of
  structure labels. Averaged over methods, it gives one number per tumour.
  **Why chance-correct:** tumours differ in which structures were sampled, hence in which
  constraints are evaluable and how hard they are by chance (C4 and C7 are much harder).
  Raw per-tumour ACS is reported alongside.

## Factors (fixed now)

- **F1 · Anatomy strength in the tumour's own RNA (deconvolution-free).** For each evaluable
  constraint, the contrast between its structures in that tumour's bulk FPKM, signed toward the
  prediction, is mean over panel genes of log2(mean_hi + 1) - log2(mean_lo + 1). C4 is MVP against
  the mean of the other sampled structures; C7 is the mean of its two steps. The per-tumour value is
  the mean over evaluable constraints with ACS weights. Panels are the marker panels fixed in
  `prespecified/abdelfattah_cluster_mapping.md` before any data:
  - Tumor: SOX2, PTPRZ1, EGFR, NES, PDGFRA
  - Oligodendrocyte: MBP, PLP1, MOG, MOBP, MAG
  - Endothelial: PECAM1, VWF, CLDN5, FLT1
  - Macrophage_Microglia: CD14, CD68, AIF1, C1QA, CSF1R
- **F1b · Marker ACS (deconvolution-free).** The constraints scored on the panel means instead of
  on any method's estimates.
- **F2 · Sampling design:** number of anatomic samples; number of structures sampled; evaluated
  constraint weight.
- **F3 · Sequencing quality:** mean total reads and mean % reads aligned to mRNA over the tumour's
  anatomic samples (`data/raw/ivygap_counts/provenance.json`).
- **F4 · Clinical and molecular, descriptive only:** age, initial KPS, MGMT methylation, EGFR
  amplification, molecular subtype, extent of resection, initial vs recurrent surgery. Group means
  only, no p-values, because groups hold 1-5 tumours. **Survival is an outcome and is excluded.**

## Statistics and reading rule

- Continuous factors (F1, F1b, the three F2, the two F3, age, KPS: 9 tests): Spearman rho across
  tumours. The p-value is exact, by full enumeration of the 9! = 362,880 orderings, two-sided, and
  Holm-adjusted across the 9.
- **Negative control:** the same correlations with the two broken-input controls' ACS excess. A
  factor that also "predicts" the controls reflects constraint evaluability, not method behaviour.
- **A factor is called associated only if** Holm-adjusted p < 0.05 **and** its control correlation
  is not significant (unadjusted p >= 0.05).
- **With n = 9 power is low.** An absent association is INCONCLUSIVE, never "no effect".
- **Prediction, stated now:** F1 is positively associated (tumours whose RNA shows the anatomy more
  strongly let methods reproduce it more often). F3 shows no association.
