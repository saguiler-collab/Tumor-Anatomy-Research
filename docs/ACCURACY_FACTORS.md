# Factors that decide whether deconvolution can be trusted without expensive measurement

Written 2026-10-06. The purpose stated by the user: find the factors that make computational
deconvolution accurate enough to use where single-cell RNA-seq, flow cytometry and other direct
measurements are not affordable.

Every number below is read from an artefact the verification re-run checks
(`docs/VERIFICATION_RERUN.md`).
- It is measured in glioma (Ivy GAP, TCGA-GBM, TCGA-LGG) against DNA truths, a second methylation
  method and direct counts.
- **Updated 2026-10-06 22:55** with an independent cohort and a per-sample whole-genome truth: CPTAC
  glioblastoma, 18 tumours (`results/cptac_wgs_purity.json`).
  - That analysis is **secondary**: registered after the primary CPTAC control failed, and before any
    of its truth values were seen (`prespecified/cptac_wgs_purity_secondary.md`).

**Figures:** `docs/figures/Figure_accuracy_factors` (panels A-D), `Figure_truth_instruments` and
`Figure_cptac_wgs` (A-F).

**Scope:** one tumour type, one family of reference atlases, 12-21 methods. These are measured
factors in glioma, not universal laws.

---

## The factors, strongest first

### 1. What is being estimated (the compartment) -- the strongest factor
Accuracy is measured against DNA truths: Spearman, median across the registered methods
(Figure_accuracy_factors A).

| compartment | GBM | LGG | verdict |
|---|---|---|---|
| tumour content (vs ABSOLUTE, DNA copy number) | 0.41 (best method 0.74) | 0.18 (best 0.55) | usable for **ranking** samples, GBM especially |
| all leukocytes (vs methylation leukocyte fraction) | 0.49 | 0.08 | usable in GBM; weak in LGG for most methods |
| lymphoid total; T, B, NK separately | about 0 (between -0.25 and +0.25, against either methylation truth) | about 0 | **not measured** by bulk RNA under these references |

- The lymphoid failure is not only a matter of noise.
  - In T against B, the methods put B at **51-71%** of lymphocytes.
  - Every direct measurement puts B at **1-9%**: flow cytometry of 17 IDH-mutant and 40 IDH-wildtype
    gliomas; two single-cell atlases (Figure_truth_instruments).
- **Replicated per sample in an independent cohort.** In CPTAC, against whole-genome DNA purity, the
  bulk tumour estimate has median Spearman **0.45** (frozen signature, 13 methods) and **0.43**
  (donor-level reference, 14 methods). TCGA-GBM gives 0.41 and 0.58 for the same methods
  (Figure_cptac_wgs D).
- *What to do:* report tumour content and total immune content. Treat lymphocyte subtypes as
  unmeasured.

### 2. The reference build -- large, and method-specific
The same method, the same tissue and two reference builds give a tumour-content accuracy that moves
by up to **0.5** in either direction (Figure_accuracy_factors B-C).

| GBM, Spearman with DNA purity | frozen signature | donor-level reference from raw counts |
|---|---|---|
| Bisque | 0.03 | 0.52 |
| SCDC | 0.29 | 0.66 |
| MuSiC | 0.43 | 0.67 |
| DWLS | 0.69 | 0.16 |
| NNLS | 0.43 | 0.25 |

- **Methods that model single cells and donors gain from a donor-level reference built from raw
  counts:** MuSiC, SCDC, Bisque.
- **Signature-based methods lose with it** (DWLS, NNLS, EPIC), and so does this project's Python
  BayesPrism reimplementation.
- A reference built from log-transformed data read as if it were linear is simply wrong (OPEN_DEFECTS
  D16).
- The lymphoid ordering also moves with the reference: 0 of 12 methods put T above B on the frozen
  signature, against 3 of 14 on the donor-level reference.
- **CPTAC replicates this only in part** (tumour content against whole-genome purity, frozen →
  donor-level):
  - Bisque gains (-0.15 → 0.43) and MuSiC gains (0.25 → 0.39).
  - SCDC is unchanged (0.45 → 0.43). NNLS loses (0.25 → 0.09).
  - DWLS loses little (0.85 → 0.74), where in TCGA it collapsed (0.69 → 0.16) with the same
    implementation.
  - So the size of the reference effect depends on the cohort as well as the method.
- *What to do:* build the reference to match the method's design, from raw counts. Never assume one
  reference suits every method.

### 3. The method itself
- Tumour-content accuracy across methods runs from **0.03 to 0.74** (GBM, registered signature).
- The best are the Bayesian regression, nu-SVR and CIBERSORTx (0.71-0.74 in GBM, 0.52-0.55 in LGG).
  - All three are this project's implementations: CIBERSORTx is its Python implementation of
    CIBERSORTx's B-mode, not the published web tool.
- The anatomic score does not identify them (factor 6).
- **The ranking transfers between cohorts on the frozen signature.**
  - Each method's accuracy in CPTAC against its accuracy in TCGA-GBM: Spearman **0.68** (p 0.016, 12
    methods).
  - TCGA's top four (Bayesian, nu-SVR, CIBERSORTx, and this project's DWLS reimplementation) are
    four of CPTAC's top five, with 0.59-0.85 in CPTAC.
  - On the donor-level reference the ranking does not transfer: 0.13 (14 methods), or 0.41 without the
    one method whose implementation differed between cohorts (Figure_cptac_wgs E).
- **The ranking does not depend on which DNA truth measures it** (`docs/EVALUATION_MATRIX.md`; post hoc,
  descriptive).
  - Each method's accuracy against ABSOLUTE copy number and against GIMiCC's methylation purity orders
    the methods almost identically: Spearman 0.80-0.87 for tumour content.
  - Against the methylation leukocyte fraction and GIMiCC's immune total it is 0.82-0.99, for immune
    content.
  - This holds in all four arms: glioblastoma and lower-grade glioma, each reference build.
  - So "which method is more accurate" is a reproducible property of the methods, not an artefact of
    one yardstick.
- *What to do:* prefer methods that rank well against a DNA truth in a similar tissue, **on the same
  reference build**. Do not choose by biological plausibility.

### 4. The implementation -- as large as the choice of method
The genuine published package and a reimplementation of the same algorithm are not interchangeable.
Measured 2026-10-07 (`docs/METHOD_REPAIRS.md`, `docs/EVALUATION_MATRIX.md`):

| method, donor-level reference | this project's reimplementation | the genuine package |
|---|---|---|
| BayesPrism, GBM tumour vs DNA | 0.16 | **0.72** |
| BayesPrism, LGG tumour vs DNA | -0.03 | **0.59** |
| BayesPrism, GBM leukocytes | 0.34 | **0.71** |

- **A reimplementation can carry a defect that silently ruins the method.**
  - This project's DWLS capped its weights relative to the largest weight instead of the smallest
    (OPEN_DEFECTS D32), and produced 37% endothelium in glioblastoma.
  - Reimplemented faithfully, it reproduces the package to 5e-11, and its GBM tumour accuracy on the
    frozen signature rises from 0.69 to 0.77.
- **A genuine package can fail for a reason that has nothing to do with its method.** DWLS skipped 73 of
  122 samples because its first solver step was not rescaled for counts-per-million inputs (D31).
  Dividing the inputs by one constant changes no answer and removes every failure.
- **A method can vanish because of how it is scored.** quanTIseq models only immune cells. A scoring
  rule that read only the tumour column dropped it from every TCGA ranking (D33).
- *What to do:*
  - Run the published package.
  - If it fails, find out why before falling back. Report any fallback, with its reason.
  - Check a reimplementation against the package on real inputs before trusting it.
  - Score each method on every compartment it models.

### 5. The tissue: estimates are compressed toward the middle
- Methods recover a median of **33% (GBM) and 23% (LGG)** of the true variation in tumour content,
  among comparable methods.
- They overestimate tumour content in low-purity samples and underestimate it in high-purity ones.
  This replicates in all 12 methods and both cohorts.
- In CPTAC the median slope of estimate on whole-genome purity is 0.22 (frozen) and 0.36
  (donor-level): the same compression.
- No other biological factor tested replicates. Ploidy, genome doubling, subclonal fraction and IDH1
  are null in both cohorts. The mesenchymal subtype is associated in GBM but cannot be tested in LGG.
- *What to do:* use tumour-content estimates to **rank** samples, not as absolute percentages.

### 6. Which ground-truth-free checks work
These are the checks a lab can run with no direct measurement (Figure_accuracy_factors D).

| check | what it ranks | Spearman with real accuracy | reading |
|---|---|---|---|
| anatomic concordance (pre-registered tissue-region constraints) | 12 methods | **0.08** (p 0.80); CPTAC replication **0.12** (p 0.72) | **detects** broken methods (margin +0.31 over two broken controls); does **not** pick the accurate one |
| agreement between methods (DECEPTICON's rule) | methods | 0.28 (mean of 6 units, p 0.035) | weak |
| agreement between methods | 6 compartments | 0.66 (p 0.088) | suggestive |
| **loss-scale stability** (re-fit one method under several loss scales) | 6 compartments | **0.89** (p 0.017) | flags the compartments not to trust |

- *What to do:*
  - Use anatomy, where region-labelled samples exist, to screen out broken methods.
  - Use loss-scale stability to flag compartments that move between fits.
  - Do not read agreement between methods as accuracy.
  - None of these certifies an estimate.

---

### 7. The yardstick: what a one-time validation must be measured against
Measured in CPTAC, where four instruments exist for the same 18 tumours (Figure_cptac_wgs A-C):

| per-sample truth for tumour content | Spearman with whole-genome DNA purity |
|---|---|
| single-nucleus RNA composition (GDC's automated clusters, named by marker panels) | **0.13** (n 15): fails |
| methylation (GIMiCC), all 18 samples on the 769 CpGs complete in every sample | 0.15 |
| methylation (GIMiCC), the 13 samples covering >= 90% of the probes (3,775 CpGs) | **0.92** |
| pathologist's percent tumour nuclei (range only 65-90%) | 0.38 (descriptive) |

- **The expensive measurement was not the valid one.** Single-nucleus composition, as processed here,
  does not measure per-sample tumour content.
  - It fails against both DNA instruments (0.13 each).
  - Nuclear capture differs by cell type, and automated clusters merge populations.
  - The authors' own curated annotation was not tested.
- **A DNA-based purity is.** Whole-genome copy number and complete-coverage methylation agree at 0.92.
- **Input completeness decides whether methylation works** (OPEN_DEFECTS D30).
  - Five samples missing a third or more of the probes shrank the shared probe set to a fifth for everyone,
    and broke the instrument (0.15).
- *What to do:*
  - Validate against a DNA-based purity.
  - Before any methylation deconvolution, check per-sample probe coverage. Set low-coverage samples
    aside rather than shrinking the probe set.

---

## The practical recipe for a lab without single-cell or flow measurement

1. Estimate **tumour content and total immune content only**. Treat lymphocyte subtypes as
   unmeasured.
2. Use the estimates to **rank** samples, not as absolute fractions.
3. Build the reference **from raw counts, matched to the method**, and run the genuine package.
4. **Choose the method from a DNA-truth benchmark in a similar tissue, on the same reference build.**
   - On the frozen signature, four methods led in both TCGA and CPTAC: this project's Bayesian
     regression, nu-SVR, CIBERSORTx B-mode implementation and DWLS reimplementation.
   - A lab running the published packages instead should confirm them once against a DNA truth
     (step 7); factor 4 shows the two are not interchangeable.
5. If region-labelled samples exist, **screen methods with a pre-registered anatomic check** to remove
   broken ones. Do not use it to pick the best.
6. **Re-fit under several loss scales**, and distrust any compartment that moves.
7. **Validate once against a DNA-based purity** on a matched dataset: copy number, or a methylation
   array with complete probe coverage. Do not use single-cell composition as the per-sample truth for
   tumour content.

## Tested on a per-sample truth: CPTAC glioblastoma (2026-10-06)

- **Data:** Wang et al., Cancer Cell 2021, doi:10.1016/j.ccell.2021.01.006. GDC project CPTAC-3, open
  access, 18 tumours.
  - Bulk RNA, single-nucleus RNA from the same cryopulverised material, and methylation.
  - Whole-genome AscatNGS purity from the GDC harmonised pipeline.
- **Primary analysis** (`prespecified/cptac_per_sample_truth.md`): **INCONCLUSIVE.**
  - Its DNA control (nuclei against methylation purity) failed, rho 0.01.
  - Two causes were established afterwards: the methylation instrument was broken by probe coverage
    (D30), and the nuclei do not track DNA.
- **Secondary analysis** (`prespecified/cptac_wgs_purity_secondary.md`, registered after that failure
  and before any of its values were seen):
  - bulk methods recover tumour content against whole-genome DNA (median 0.45 / 0.43);
  - on the frozen signature the method ranking transfers from TCGA (0.68);
  - anatomy does not rank methods by accuracy (0.12);
  - nuclei composition is not a per-sample truth (0.13).
- **Still open:**
  - A per-sample lymphoid truth. CPTAC's nuclei contain B cells in only 2 of 15 tumours, so the T
    against B question (OPEN_DEFECTS D23) stays at the cohort level, where flow cytometry and two
    atlases answer it.
  - A test of the nuclei with the authors' own curated cell annotation.
