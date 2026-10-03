# ReCIDE on TCGA -- completing its extension-panel row, and adjudicating two truth-free rankings

**Written 2026-10-02, before any ReCIDE estimate on TCGA exists.** (The clock time comes from the file
system: read it with `stat`, not from this line.) Post-registration and exploratory. It selects
nothing and changes no registered statistic.

## Why

- **ReCIDE finished on the Ivy GAP anatomic cohort on 2026-10-02**
  (`results/extension/recide_anatomic.json`):
  - ACS 0.692 [0.548, 0.829], beating its permutation null;
  - **the lowest anatomic ACS of any real method in the study**.
- **Li et al. 2026 rank it at the top** ([4]; quoted: "ReCIDE and BayesPrism stand out as two robust
  deconvolution methods"). Their criteria are truth-free: consistency with scRNA-seq,
  reproducibility across cohorts, reproducibility of prognostic relevance. The paper calls it "our
  recently developed ReCIDE".
- **Two truth-free criteria therefore disagree.** Only an orthogonal truth can say which ranking, if
  either, reflects accuracy. ReCIDE is the one extension method without a TCGA row.

## Inputs (fixed)

- **Bulk, gene space, samples and purity:** the TCGA raw/X ("h5ad") arm exactly as
  `scripts/extension_tcga.py` builds it for every other extension method in that arm.
- **Reference:** ReCIDE's own, exactly as in its anatomic run (`scripts/run_recide.py`).
  - The same pre-fixed donor-draw rule and project seed, which gives the same donors.
  - Uncapped cells from those donors, restricted to the arm's gene space.
  - Counts recovered and checked integral.
- **Code:** ReCIDE vendored @31bbd6b with the disclosed `summarize_each` shim, through
  `R/run_recide.R`, unbudgeted.
- **Order:** GBM first; LGG afterwards if the machine allows.

## Control (before any ReCIDE number is read)

- The refactored TCGA input path must reproduce an existing raw/X-arm row exactly: `deseq2_unmix` in
  `results/extension/tcga_gbm_h5ad.json`, every field.
- The donor draw must reproduce the anatomic run's 20 donors exactly.
- Either failing aborts.

## Reported (as for every extension method)

- Purity rho against ABSOLUTE and the bias.
- The lymphoid ordering against methylation, B above T per sample, and samples with no lymphoid
  estimate.
- L1 distance of the cohort mean to the reference prior.

## The adjudication (descriptive, one method, no test)

Rank ReCIDE's GBM purity rho among every GBM row of the raw/X arm in the study: the extension
methods' raw/X rows and genuine BayesPrism.

- **Top third:** the reproducibility ranking is consistent with DNA in this cohort, and anatomy's
  lowest-ACS verdict is not.
- **Bottom third:** anatomy's low ranking agrees with DNA here, and the reproducibility ranking does
  not.
- **Middle:** neither truth-free ranking is borne out.

Whatever the result, the study's registered conclusion stands: ACS does not rank methods by accuracy.
This is one method, so the reading is a single observation and is never generalised.

**Also, exploratory:** the extended agreement test gains ReCIDE (n = 17), reported and not tested, as
before.
