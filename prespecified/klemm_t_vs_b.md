# T against B in IDH-mutant gliomas, by flow cytometry (Klemm et al. 2020 [53]) -- does direct measurement settle LGG?

**Written 2026-10-06 16:26:00 EDT**, the file's birth time (`stat`). Exploratory, post-registration.

**What had been seen when this was written:** only text, never a figure or a plotted value. The
user placed the paper and its two supplements in `celldecov_reference_papers/Klemm_et_al_2020/`.
- **Table S2** (`mmc2.pdf`, read in full) holds the flow-cytometry definitions:
  - T cells are CD3+ (CD4+, Treg, CD8+ and double-negative T, each defined);
  - B cells are CD45+ CD11B- CD19/CD20+ CD3-;
  - NK cells are CD3- CD56+.
  - It holds no per-patient values.
- **Table S1** (`mmc1.pdf`, read in full) is the clinical summary of the flow-cytometry cohort:
  - 17 IDH-mutant gliomas (3 oligodendroglioma, 3 anaplastic oligodendroglioma, 7 diffuse
    astrocytoma, 3 anaplastic astrocytoma, 1 IDH-mutant GBM);
  - 40 IDH-wildtype gliomas (38 GBM, 2 anaplastic astrocytoma);
  - 37 brain metastases and 6 non-tumour samples.
  - It holds no immune values.
- **The legends of Figure 1 and Figure S1** (text only):
  - Figure 1F is "Mean of immune cell populations in ... gliomas (nIDH mut = 17, nIDH WT = 40) ... as
    percentage of CD45+ cells".
  - Figure 1D is a per-sample heatmap of the same proportions.
- Earlier (2026-10-03), the paper's text: "the lymphocyte compartment was mostly composed of T cells
  with fewer NK cells and B cells". That statement is pooled over all samples.

## Why

TCGA-LGG is mostly IDH-mutant. In LGG the two methylation instruments disagree: EpiDISH gives
T > B and GIMiCC gives B > T (`gimicc_truth_confirmation.md`). Direct counts so far rest on 2 LGG
patients (`atlas_t_vs_b.md`). Flow cytometry of 17 IDH-mutant gliomas, with gates that separate
CD3+ from CD19/CD20+ cells, is a direct cohort-level measurement of exactly the quantity in dispute.

## The quantity, read from Figure 1F

- **T** = CD4+ T + Treg + CD8+ T + double-negative T, each as % of CD45+, summed.
- **B** = B cells as % of CD45+.
- Both are read for the IDH-mutant glioma group (primary) and the IDH-wildtype group (secondary,
  the GBM-like cohort).
- **How values are read:**
  - If the figure prints numbers, those are used.
  - Otherwise each bar or segment is measured against the figure's own axis.
  - The reading resolution (the smallest difference the panel lets one distinguish) is stated with
    every value.
- NK is read the same way and reported, for the full T / NK / B ordering.

## Reading (Q-K), fixed now

For each glioma group:
- **SUPPORTS T > B** if T exceeds B by more than the stated reading resolution.
- **CONTRADICTS T > B** if B exceeds T by more than the resolution, or B >= T where both are
  printed.
- **UNRESOLVED** if the panel cannot separate them (B not shown on its own, values below
  resolution, or no legible axis).

**Consequence for the study, fixed now:**
- **IDH-mutant group SUPPORTS T > B:** GIMiCC's LGG B > T is contradicted by direct measurement in
  17 IDH-mutant gliomas.
  - The LGG truth is then stated as EpiDISH's direction, supported by cytometry, with GIMiCC's
    disagreement reported as an instrument failure in IDH-mutant tissue.
  - The registered GIMiCC reading (INCONCLUSIVE) does not change.
- **IDH-mutant group CONTRADICTS T > B:** direct measurement agrees with GIMiCC. The LGG headline
  cannot be stated, and the manuscript says so first.
- **UNRESOLVED:** LGG stays contested.

**Secondary, descriptive:** Figure 1D's heatmap, per IDH-mutant sample, T against B, only if the
colour scale makes the comparison legible. A colour reading never overrides Figure 1F.

## Limits, stated now

- These are cohort means of percentages (Figure 1F), not per-patient values. The cohort is external
  (the authors' institutions, not TCGA).
- Dissociation for flow cytometry can lose cells, although T and B cells are both lymphocytes.
- 17 IDH-mutant samples of mixed histology, including oligodendrogliomas.
- A value read from a plot carries reading error, stated with each value.
