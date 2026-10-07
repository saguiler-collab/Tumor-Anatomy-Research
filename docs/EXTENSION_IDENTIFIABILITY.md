# Extension E2 -- Can a user tell, without ground truth, which cell types the bulk data determine?

*Sections 1–6 were written on 2026-10-02 before any analysis in them was run; §0, §2b and §7 were added
after, and say so. Binding rules: `prespecified/identifiability_diagnostics.md` (with Addendum 1) and
`prespecified/agreement_selection_test.md`.*

## 0 · Summary (2026-10-02, after every check below)

**Question.** Without ground truth, can a user tell which cell types a bulk deconvolution actually
determines?

**Design.**
- Two truth-free checks, applied to every compartment in TCGA-GBM and TCGA-LGG:
  - re-fit one genuine package (DESeq2 `unmix`) under seven loss scales and ask whether samples keep
    their order (**stability**);
  - ask the same across the registered methods (**agreement**).
- Then compare both with instruments that share nothing with the RNA: ABSOLUTE purity, the
  methylation leukocyte fraction and EpiDISH.

**Result.**
- **Stability ranks compartments the way DNA does.** rho = 0.886, exact
  p = 0.0167 over six units, as registered. With truth and estimate on
  matched denominators: rho = 0.943, p = 0.0083.
- **Tumour and total leukocyte content are stable and accurate.** Stability >= 0.93. `unmix` vs DNA: tumour 0.78
  (GBM) and 0.59 (LGG); leukocytes 0.78 and 0.60. `unmix`'s leukocyte estimate tracks
  methylation better than every registered method in both cohorts.
- **The lymphoid compartment is the least stable unit in both cohorts, and its T/B/NK split is not
  tracked.** The immune compartment is identified as a whole, not its lymphocyte types. That is a
  sharper form of the study's headline than "no method recovers T above B", and it matches what the
  2024–2026 benchmarks found with truth (DREAM; omnideconv: "all methods struggled to robustly estimate
  tumor-infiltrating lymphocytes"), seen here without truth.
- **Robust:**
  - to sampling: the lymphoid units are the two least stable in 100% of
    2,000 bootstraps;
  - to a near-tie;
  - to three design flaws, each found after the run and checked under a rule written first (a duplicated
    fit, degraded methods in the agreement panel, and truths on a different denominator).
- **Agreement between methods is the weaker signal.** DECEPTICON's selection rule, tested here against DNA
  truth:
  - centrality tracks accuracy only weakly (mean rho 0.28, p 0.035);
  - the selected ensemble beats random pairs in 4 of 6 units;
  - its top pair was a pair of algorithmic siblings in 9 of 10 selections.

**Limits.**
- **Stability flags what not to trust; it certifies nothing.** GBM T cells are stable and do not track
  EpiDISH. Against GIMiCC, a second methylation truth, they do track (+0.64; §7.4), so this
  counter-example depends on which methylation instrument is right.
- **The lymphoid claim is narrower than first stated (§7.4).** Against GIMiCC's tissue fractions,
  `unmix`'s lymphoid total tracks (+0.65 GBM, +0.53 LGG), beyond the leukocyte level (partial +0.36 /
  +0.32). The methods' median does not track. H1 holds under that truth, through agreement (D2 rho 1.00,
  p 0.0014).
- Six primary units.
- One package for the stability arm.
- Post-registration and exploratory throughout. It selects nothing and changes no registered statistic.

## 1 · Why this extension exists

The study's two registered findings leave a gap a method developer cares about.

- **Anatomy (ACS) catches broken methods but does not rank accurate ones** (`docs/WHAT_ACS_IS.md`).
- **No method robustly recovers the T-above-B lymphoid ordering**, while tumour content is
  recovered reasonably (purity rho up to 0.78).

On 2026-10-02 a post hoc ablation inside one genuine package (DESeq2 `unmix`) showed why these differ
(`docs/WHY_B_OVER_T.md` §7l):
- changing only the scale of the fitting loss flips the T-versus-B ordering;
- tumour-content tracking barely moves.

That suggests a general and checkable idea: **a cell type whose estimate changes with an arbitrary
modelling choice is not determined by the data, and a truth-free stability check could tell a user
which estimates to trust.**

ACS answers "is this method broken?". This extension asks a different question of every method:
**"which of its cell types can be believed?"**

## 2 · What the literature says (text-searched in `celldecov_reference_papers/`, quoted verbatim)

- **The problem is known.** Sturm et al. 2019: "deconvolution-based approaches are susceptible to
  background-predictions that might be due to the similarity of the signatures of closely related
  cell types (multicollinearity) and/or to a lower cell-specificity of signature/marker genes."
  Nguyen et al. 2024: "Multicollinearity between cell types also leads to unstable results or multiple
  solutions in cellular deconvolution."
- **The field is moving to truth-free evaluation.** Li et al. 2026 (*Genome Biol*): "traditional
  pseudobulk benchmarks may not always be reliable in real-world settings where absolute cell
  proportions are unknown". They evaluate by "consistency with scRNA-seq, reproducibility across
  cohorts, and reproducibility of prognostic relevance". Gaspard-Boulinc et al. 2025 (*Nat Rev
  Genet*): "For users of deconvolution methods, result validation is a major challenge because the
  ground truth for cell-type composition is often not available."
- **What no source does:** test whether a truth-free robustness criterion predicts accuracy, cell
  type by cell type, against measurements that share nothing with the RNA. Li et al. had no
  absolute truth; Sturm et al. used simulated mixtures. This study has DNA truths (ABSOLUTE purity,
  methylation leukocyte fraction, EpiDISH) for several compartments, in two cohorts.
- **A tension that must be addressed.** Avila Cobos et al. 2020: "Maintaining the data in linear
  scale ... consistently showed the best results (lowest RMSE values) whereas the logarithmic ...
  and VST ... scale led to a poorer performance, with two to four-fold higher median RMSE values".
  Their conclusion: "computational deconvolution should be performed on linear scale".
  `unmix`, whose log-scale loss gave the methylation-consistent ordering here, appears to contradict
  this. It does not: Avila Cobos transformed the **data** before a linear deconvolution, which breaks
  the linear mixture model (the log of a sum is not a sum of logs). `unmix` keeps the mixture
  **linear** and transforms only where the **loss** is measured. This extension makes the
  distinction explicit with an arm that transforms the data instead (secondary arm S1).

## 2b · The 2022–2026 literature, read in full (added 2026-10-02, after the E2 run)

*Nothing in this section was used to design E2. It positions the result.* Nine recent papers were
read in full text (Europe PMC open access; copies in the session scratchpad). Every quoted passage
was text-searched in that copy. References [46]–[49] are in `docs/REFERENCES.md` (Crossref-verified);
the rest are cited here by DOI.

| paper | truth used | finding that bears on E2 (verbatim) |
|---|---|---|
| White et al. 2024, *Nat Commun* 15:7362 — the DREAM challenge [47] | in vitro / in silico admixtures of purified cells | "we find that most methods predict coarse-grained populations well"; "Seven methods showed similar mean limits of detection (3–4%) within the coarse-grained sub-Challenge"; naïve and parental CD4+ T cells: "best-case limits of detection of 6%" |
| Dietrich et al. 2026, *Genome Biol* 27:6 — omnideconv [48] | pseudo-bulk; real lung tumours with IHC CD4+/CD8+ T cells | "all methods struggled to robustly estimate tumor-infiltrating lymphocytes, which are transcriptionally similar and rare compared to the overall sample composition"; a signature's condition number read as its "discriminative power" |
| Huuki-Myers et al. 2025, *Genome Biol* 26:88 [49] | RNAScope/IF in human DLPFC, 22 tissue blocks (orthogonal) | "Bisque and hspe were the most accurate methods"; "the performance rankings of different deconvolution methods have mostly been inconsistent" |
| Deng et al. 2025, *Brief Bioinform* 26:bbaf234 — DECEPTICON [46] | simulations; FACS in blood and non-brain tumours; TCGA via prognosis | selects estimates by agreement: "the most effective analysis strategies for deconvolution demonstrate significant concordance" |
| Li et al. 2026, *Genome Biol* [4] (§2) | none: consistency with scRNA-seq, reproducibility across cohorts and of prognostic relevance | truth-free benchmarking on 5,891 real samples |
| Xu et al. 2025, *Brief Bioinform* 26:bbaf264 (doi:10.1093/bib/bbaf264) | pseudo-bulk with simulated truth | robustness and resilience to perturbed references, scored against simulated truth |
| Ivich et al. 2025, *Genome Biol* 26 (doi:10.1186/s13059-025-03506-9); Fan et al. 2026, *Bioinformatics* 42 (DeconX, doi:10.1093/bioinformatics/btag412) | simulations; real bulk qualitatively | missing reference cell types leave signal in the residuals. Of real bulk: "no ground truth was available" |
| Sutton et al. 2022, *Nat Commun* 13 (doi:10.1038/s41467-022-28655-4) | in silico and RNA mixtures; ~2,000 brain samples | "For datasets without a ground truth, such as bulk brain samples, goodness-of-fit was evaluated" |
| Dai et al. 2024, *Sci Adv* 10 (doi:10.1126/sciadv.adh2588) | matched IHC and sc/snRNA-seq in brain | "the necessity of using real ground truth data in benchmarking studies" |

**Factors the benchmarks name, and what this study measured for each.**
- omnideconv [48] tests five, quoted: "1) cell-type-specific mRNA bias ... 2) unknown cellular
  content ... 3) impact of the resolution of cell type annotations (i.e., from coarse to
  fine-grained); 4) transcriptional similarity between closely-related cell types; 5) variance and
  batch effects affecting reference single-cell datasets".
- Huuki-Myers et al. [49] list eight reasons rankings disagree, including "(3) variability in the
  selection of cell type marker genes", "(5) choice of cell type resolution (fine or broad)" and
  "(8) differences in data normalization and processing".
- DREAM [47]: "The limit of detection by deconvolution is likely to vary from cell type to cell type
  dependent on the uniqueness and strength of their transcriptional signal."

| factor | what this study measured |
|---|---|
| mRNA content per cell | the central cell-size conversion on every mRNA-proportion method (`ivygap/deconv/`) |
| unknown cellular content | the additive model explains a minority of real bulk (MANUSCRIPT §4.5); per-sample fit does not flag untrustworthy samples (Appendix D) |
| resolution, coarse to fine | **E2: all-leukocyte content identified; its lymphoid split not** (§7) |
| similarity of close types | `B_cell`'s profile is closest to `Macrophage_Microglia` (r = 0.497); the `NK_cell` column carries pan-T markers (WHY_B_OVER_T §7) -- properties of the reference, not shown to cause the inversion |
| reference variance / batch | an independent atlas (GSE182109) gives a MIXED answer; the headline is reference-specific (D23) |
| normalisation and processing | **loss scale decides the T-vs-B split inside one package (WHY_B_OVER_T §7l); data scale (S1, §7) changes purity ranking** |
| marker selection | registered fixed rule (`select_signature_genes`), not varied -- not measured |

**Four truth-free criteria are in current use. In the papers read here, none had been checked
compartment by compartment against an orthogonal per-sample truth in tumours:**
1. **Agreement between methods** (DECEPTICON). E2's D2 and `prespecified/agreement_selection_test.md`
   test it directly.
2. **Reproducibility across cohorts or of prognostic association** (Li et al. 2026). It cannot see an
   estimate that is reproducibly wrong; E2's GBM T cells are that case.
3. **Goodness of fit and residuals** (Sutton et al. 2022; Ivich et al. 2025; DeconX). In this study,
   per-sample model fit did not identify untrustworthy samples (MANUSCRIPT Appendix D).
4. **Stability under an arbitrary modelling choice** (E2's D1), here validated against DNA truths.

**Where E2 agrees with the 2024–2026 benchmarks, and what it adds.** DREAM and omnideconv found,
*with* truth, that fine immune subsets and tumour-infiltrating lymphocytes are where estimates fail.
omnideconv attributes this to similarity and rarity. DREAM puts typical detection limits for coarse
types at 3–4%, with 6% best-case for CD4+ T-cell subsets. E2 reaches the same boundary in glioma
from the other side: the failure can be **seen without truth**, because the lymphoid compartment is
the one that moves when an arbitrary modelling choice changes. That tumour content and total
leukocyte content hold is what a user of a glioma deconvolution most needs to know. Whether glioma
lymphocytes sit below DREAM's detection limits cannot be settled here: EpiDISH's blood reference
gives shares of the immune compartment, not of the tissue (Addendum 1, check 7).

## 3 · Design

### Two truth-free diagnostics, per compartment and cohort

**D1 · Loss-scale stability** -- one genuine package, one modelling choice varied.
- DESeq2 `unmix` on the registered TCGA frozen arm, seven settings: the pre-declared shift; forced
  shift 1, 10, 100, 1000 and 10000 at power 1; the pre-declared shift at power 2.
- Per-sample estimates of every cell type are saved for each setting.
- Stability = mean pairwise Spearman correlation of the per-sample estimates across the seven
  settings (do samples keep their order whatever the loss?).

**D2 · Method agreement** -- many published methods, each its own modelling choice.
- The registered panel's saved per-sample TCGA estimates (`results/estimates_full*.csv`).
- Agreement = mean pairwise Spearman across methods.
- Degenerate duplicates are dropped (MuSiC = NNLS on the frozen signature; SCDC ENSEMBLE = SCDC), so
  a method is never counted twice.

### Compartments

Tumor; all leukocytes (macrophage/microglia + T + NK + B); myeloid (macrophage/microglia); lymphoid
(T + NK + B); T; B; NK; endothelial; oligodendrocyte. The Astrocyte sidecar is excluded:
unidentifiable against Tumor by construction.

### Truths -- per sample, sharing nothing with the RNA

| compartment | truth | status |
|---|---|---|
| Tumor | ABSOLUTE purity (DNA copy number) | primary |
| all leukocytes | leukocyte fraction from DNA methylation (Thorsson et al. 2018) | primary |
| lymphoid | EpiDISH T + NK + B (DNA methylation) | primary |
| T, B, NK separately | EpiDISH per type | **secondary**: the per-sample split failed a reference-free positive control (D23) |
| myeloid, endothelial, oligodendrocyte | none in TCGA | stability reported, no truth test |

### Truth agreement

Spearman rho of estimate against truth across samples:
- for D1, at `unmix`'s pre-declared setting;
- for D2, the median across methods.

### The test

Across the compartment-by-cohort units with a primary truth (3 compartments x 2 cohorts = 6 units),
does truth-free stability rank the units the way truth agreement does? D1 and D2 are tested
separately. A secondary analysis adds the T, B and NK units (12 units).

### Secondary arms

- **S1 · data scale versus loss scale.** NNLS on log1p-transformed bulk and signature (the data
  transformed, the Avila Cobos case) against the registered linear NNLS and against `unmix` (the
  loss transformed). This project's code, labelled as such. Expected from Avila Cobos: the log-data
  arm tracks purity worse.
- **S2 · the anatomy arm.** `unmix`'s seven settings on Ivy GAP: is per-constraint satisfaction stable
  across settings for the constraints on each cell type? This links D1 to ACS.

## 4 · What has already been seen (disclosed, so nothing below is mistaken for a prediction)

- `unmix` at six of the seven settings on the frozen arm, GBM and LGG: summary statistics only
  (cohort-level T/B ordering, share of samples with B > T, purity rho). **T and B flip with the loss;
  tumour content does not.** This motivated the hypothesis, so it cannot confirm it.
- Every registered method's purity rho and lymphoid ordering, and its immune-total rho against the
  leukocyte fraction (`results/immune_arm*.json`, September 2026).
- **Not yet computed by anyone:**
  - D1 for any compartment other than Tumor and the T/B ordering;
  - D2 for any compartment;
  - the per-sample truth agreement of `unmix`'s leukocyte and lymphoid totals;
  - the relation between either diagnostic and truth agreement -- the confirmatory question.

## 5 · What each outcome would mean

- **Supported** -- stability tracks truth agreement. A user can screen their own deconvolution, with no
  ground truth, by re-running it under different losses (or methods) and believing only the
  compartments that do not move. That is a practical, generalisable tool, and the natural partner to
  ACS: ACS flags broken methods, stability flags unidentifiable cell types.
- **Not supported** -- stability can coexist with error, for example a consistent bias shared by all
  losses or all methods. Then truth-free robustness, the field's new evaluation direction, is not
  evidence of accuracy either. That would be the same lesson as ACS, one level down, and it is just
  as reportable.
- **Inconclusive** -- with six units, likely for anything short of a clean ordering; reported as such.

## 6 · Limits stated in advance

- Six primary units: a rank test across them can only detect a strong association.
- D1 uses one method family; D2 covers many methods but only the registered ones (the extension
  methods' per-sample TCGA estimates were not saved).
- The truths measure different denominators: purity is a DNA fraction; the leukocyte fraction is a
  methylation estimate of all leukocytes; EpiDISH is a blood reference applied to brain. Each is used
  only for rank agreement, never for level.
- Exploratory and post-registration. Selects nothing, and changes no registered statistic.

---

## 7 · Results (run 2026-10-02 08:47; `results/identifiability_diagnostics.json`)

**Controls passed in both cohorts.**
- `unmix`'s Tumor rho at the pre-declared setting reproduces the reported rows (GBM 0.7801, LGG 0.5916).
- The registered NNLS's reproduces from the saved estimates (0.4271, 0.1635).

| unit (GBM / LGG) | D1 loss-scale stability | D2 method agreement | truth rho, `unmix` | truth rho, median of methods | n with truth |
|---|---|---|---|---|---|
| **Tumor** (ABSOLUTE) | 0.960 / 0.946 | 0.466 / 0.380 | **0.780 / 0.592** | 0.413 / 0.183 | 154 / 510 |
| **Leukocytes** (methylation LF) | 0.954 / 0.927 | 0.561 / 0.354 | **0.780 / 0.603** | 0.493 / 0.079 | 129 / 510 |
| **Lymphoid** (EpiDISH) | 0.713 / 0.494 | 0.277 / 0.184 | **-0.200 / -0.183** | -0.141 / 0.188 | 56 / 510 |
| T (secondary) | 0.893 / 0.749 | 0.108 / 0.037 | -0.119 / -0.273 | -0.083 / -0.027 | 56 / 510 |
| B (secondary) | 0.520 / 0.468 | 0.258 / 0.199 | 0.090 / -0.001 | -0.126 / 0.084 | 56 / 510 |
| NK (secondary) | 0.703 / 0.742 | 0.089 / 0.059 | -0.120 / -0.070 | 0.004 / -0.031 | 56 / 510 |
| Myeloid, Endothelial, Oligodendrocyte (no truth) | 0.936-0.990 | 0.466-0.780 | -- | -- | -- |

- **D2 used all 10 methods for every unit except LGG T.** There, NNLS returns T = 0 for every sample
  and was excluded as constant (rule: reported, not hidden).
- **Exact zeros in `unmix`'s estimates**, as a share of samples averaged over the seven settings
  (GBM / LGG): T 36% / 76%; NK 73% / 91%; B 29% / 33%; lymphoid total 13% / 28%. Rank statistics on
  these compartments rest partly on ties.
- **The GBM lymphoid truth has n = 56**, the EpiDISH subset, so its rho carries an interval of
  roughly +/-0.26.

- **H1 (primary): SUPPORTED.**
  - D1: rho = 0.886 across the 6 primary units, exact one-sided p = 0.0167 (12 of the 720 orderings
    are at least as concordant).
  - D2: rho = 0.657, p = 0.0875 (same direction, not significant alone).
- **H1s (12 units, adding T, B, NK):** D1 rho = 0.587 (p = 0.024); D2 rho = 0.518 (p = 0.044).
- **H2 (leukocytes identifiable): PASS in both cohorts.**
  - Stability is 0.954 (GBM) and 0.927 (LGG).
  - `unmix`'s all-leukocyte estimate tracks the methylation leukocyte fraction at rho = 0.780 and
    0.603. That is **above every registered method in both cohorts**: the best are nu-SVR at 0.688
    (GBM) and DWLS, this project's reimplementation, at 0.590 (LGG). The LGG margin, 0.013, was not
    tested.
  - The per-method values reproduce the September immune arm (`results/immune_arm*.json`) exactly, by
    an independent code path. The degraded-mode methods: Bisque 0.004 / -0.275 and EPIC 0.025 / -0.159.
- **S1 (data scale): the expectation is NOT supported in either cohort.** On identical inputs without
  the cell-size conversion, NNLS on log1p data ranked tumour content better than linear NNLS: 0.715
  vs 0.562 (GBM), 0.557 vs 0.275 (LGG). It also beats the registered NNLS row (0.427 / 0.164).
  That row comes from the registered pipeline's own preprocessing; why the two linear rows differ
  was not examined.

### 7.1 Robustness (Addendum 1 of the binding rules, written before any of it was computed; `results/identifiability_robustness.json`, `scripts/identifiability_robustness.py`)

The control passed: the registered E2 values reproduce from the saved matrices to 1e-16. The seeded
bootstrap gave identical results in two runs.

| check | result |
|---|---|
| 1 · contrast alone | P(the two lowest-truth units are also the two least stable) = 1/15 = **0.067** under the null. The contrast carries most of H1; the order of the four stable units supplies the rest. |
| 2 · near-tie swap | GBM Tumor and GBM Leukocytes truths differ by 0.0003. Swapped: rho 0.829, **p = 0.029** -- still significant. |
| 3 · bootstrap, B = 2000 | Lymphoid units are the two least stable in **100%** of resamples (D1) and **98.8%** (D2). 95% interval of the H1 rho: D1 **0.54 to 1.00**; D2 0.60 to 0.89. Both readings: contrast robust; fine ordering robust to sampling. The bootstrap holds the six units fixed, so it says nothing about how few there are. |
| 4 · H2 wording | `unmix` exceeds every registered method in both cohorts (values above). |
| 5 · **duplicate setting (design flaw)** | In GBM the pre-declared shift is 1, so `predeclared` and `s1` are the same fit (byte-identical). Over 6 distinct fits: GBM D1 Tumor 0.957, Leukocytes 0.950, Lymphoid 0.691; LGG unchanged. **H1 D1 unchanged: rho 0.886, p 0.0167; H1s D1 rho 0.601, p 0.021.** |
| 6 · **degraded methods in D2 (design inconsistency)** | D2 kept Bisque and EPIC, which the registered immune arm excludes. On that arm's 8-method panel: **H1 D2 rho 0.771, p 0.0514** (registered 0.657, 0.0875); H1s D2 rho 0.552, p 0.033. The lymphoid units are still the two lowest. D2's reading does not change: still short of significance. The comparable panel is the cleaner one. |
| 7 · **denominator mismatch in the EpiDISH units (design flaw)** | EpiDISH's blood reference returns shares of the immune compartment (columns sum to 1); E2 compared them with tissue fractions. Matched (each estimate's share of its own leukocyte total; `results/identifiability_denominators.json`): `unmix` lymphoid vs EpiDISH **-0.103 / -0.123** (registered -0.200 / -0.183): **the non-agreement STANDS** (rule: withdrawn only if >= 0.30). The lymphoid share is LESS stable than the tissue fraction (0.552 / 0.341 vs 0.713 / 0.494). **H1 strengthens: D1 rho 0.943, p 0.0083; D2 rho 0.829, p 0.029** (H1s: 0.685, p 0.008; 0.825, p 0.0008). Caveat: matched D2 drops samples where any method returns zero leukocytes -- 26 of 154 (GBM) and **247 of 510 (LGG)**, mostly BayesPrism (reimplementation; 201) and NNLS (156) -- so the LGG D2 value rests on the more immune-rich half. D1 drops none. |

**What this shows.**
1. On real tumours, with DNA truths, the compartments whose estimates do not move when the loss
   changes are the ones the estimates get right. Tumor content and total leukocyte content are stable
   and accurate. The lymphoid compartment and its split are less stable. `unmix`'s lymphoid estimate
   does not track EpiDISH (rho near or below zero, on either denominator). Across the registered
   methods the lymphoid share is at best weakly tracked: median 0.24 in LGG on matched
   denominators, about zero in GBM (n = 56), 0.40 for the agreement-selected ensemble in LGG
   (§7.2). The T/B/NK split is not tracked.
   - A user can run this check with no ground truth.
   - The contrast between {Tumor, Leukocytes} and {Lymphoid} is robust to sampling, to the GBM
     duplicate and to the near-tie.
2. **The immune compartment is identified as a whole but not resolved into lymphocyte types.** That
   is a sharper statement of the study's lymphoid result than "no method recovers T > B".
3. **The anatomic constraints sit on the identified compartments.** The seven registered constraints
   score Tumor, Macrophage_Microglia, Endothelial and Oligodendrocyte. All four are stable here
   (0.936-0.990, both cohorts). T cells were excluded from the constraint file before any output was
   examined: "This pipeline's own synthetic benchmark puts T cells at the detection floor". E2 finds
   the lymphoid compartment the least stable, so that exclusion is now confirmed on real data without
   truth. ACS scores no lymphoid contrast, and the compartments it does score are ones the bulk data
   determine.

**The truth denominators.** The EpiDISH truths are shares of the immune compartment; the registered comparison set them against tissue fractions (check 7). On matched denominators the lymphoid estimate still does not track methylation and the association strengthens. The registered values stay as recorded; the matched ones are the better-posed comparison.

**What it does not show.**
- **Six primary units.** The robust content is a two-group contrast, whose null probability alone
  is 0.067. The rank test cannot say how finely stability grades accuracy.
- **Stability is not sufficient.** T cells in GBM are stable (0.893; 0.804 as a share of leukocytes) and
  do not track EpiDISH (-0.119; -0.006 on matched denominators, so "anti-correlated" was the mismatch): a
  reproducibly wrong estimate, the case a reproducibility criterion cannot detect. A consistent bias
  shared by every loss passes the check, so stability flags what not to trust and certifies nothing.
  Part of T's stability is ties: 36% of GBM samples are exactly zero, and constant zeros agree with
  themselves across settings.
- **S2 (the anatomy arm) ran later, under Addendum 2;** see §7.3.
- **D1 is one method family.** D2, across the registered methods, points the same way; on the cleaner
  8-method panel it falls just short of significance (p 0.051).
- **On S1:** Avila Cobos et al. measured absolute error (RMSE) on simulated mixtures; S1 measures how
  well tumour content is ranked against DNA on real tumours. The two need not agree: a log-scale fit
  can rank samples better while misjudging their levels. S1 did not measure level error, so the paper
  reports the contrast with that difference stated, and does not claim to overturn their result.
- **Post-registration.** The ablation that motivated E2 had already shown the T/B split moving with
  the loss while tumour content held. E2's confirmatory content is the generalisation to leukocytes,
  the per-unit truth agreement, and D2, none of which had been computed.

### 7.2 Does agreement between methods pick the accurate estimate? DECEPTICON's criterion against DNA truth

Rules: `prespecified/agreement_selection_test.md`, written before computing. Code:
`scripts/agreement_selection_test.py`, this project's implementation of the published selection
rule, not the DECEPTICON package. Results: `results/agreement_selection_test.json`.
- Panel: the immune arm's 8 comparable methods.
- Units: E2's, with matched denominators.
- 2000 random-pair draws.
- A 20-draw smoke run of the same code was seen first; P2 does not depend on the number of draws and
  came out identical.

| unit (GBM / LGG) | agreement-selected ensemble | random-pair median (its percentile) | median single method | best single method |
|---|---|---|---|---|
| Tumor | 0.680 / 0.500 | 0.655 (0.60) / 0.427 (0.66) | 0.557 / 0.342 | Bayesian 0.742 / nu-SVR 0.549 |
| Leukocytes | 0.699 / 0.172 | 0.611 (0.90) / 0.172 (0.50) | 0.531 / 0.155 | nu-SVR 0.688 / DWLS (reimpl.) 0.590 |
| Lymphoid share | -0.067 / 0.395 | 0.135 (0.03) / 0.344 (0.99) | -0.050 / 0.242 | DWLS (reimpl.) 0.229 / Bayesian 0.373 |

- **P1 (does the selected ensemble beat chance selection?): INCONCLUSIVE.** It beats the random-pair
  median in 4 of 6 primary units (sign test p 0.34); 5 of 6 on the 10-method panel.
- **P2 (do the methods that agree most agree most with DNA?): SUPPORTED, weakly.** The mean per-unit
  Spearman between centrality and truth agreement is 0.28 (permutation p 0.035); 10-method panel
  0.27, p 0.023. Per unit (GBM / LGG):
  - Tumor 0.55 / 0.74;
  - Leukocytes 0.33 / 0.07;
  - Lymphoid -0.71 / 0.71.

**What drives the agreement.** The highest-correlation pair was a pair of **algorithmic siblings**
in **9 of 10** cell-type selections (12 of all 20 selected pairs):
- this project's `cibersortx` subclasses its `svr` (the same nu-SVR core, plus batch correction);
- the BayesPrism reimplementation starts from NNLS and iterates it against a rescaled reference.

These pairs correlate because they are built the same way, so their agreement measures shared
construction, not shared correctness. The clearest cost is in LGG leukocytes: the selected ensemble
reaches 0.172, no better than random pairs. The most accurate single method (DWLS,
reimplementation, 0.590) ranks fifth of eight by agreement: mean correlation 0.38, against 0.58
for the top method, this project's CIBERSORTx-style nu-SVR, whose truth agreement is 0.20.

**Reading.**
- Agreement carries some information about accuracy on average, which is DECEPTICON's premise.
- It does not reliably pick the accurate estimate. It fails exactly where a panel contains siblings,
  a property of most real panels (one algorithm with several references, or one reference with
  several solvers).
- **GBM lymphoid (-0.71) is not interpretable.** Its truth agreements span -0.18 to +0.23 on n = 56
  samples, all within sampling noise of zero, so no criterion can rank methods there.
- **LGG lymphoid share is the one place where agreement found a better estimate** (0.395, 99th
  percentile of random selections, above every single method). Even the lymphoid total is weakly
  recoverable in LGG, where n = 510. The T/B/NK split is not: the selected ensemble's secondary units lie between
  -0.23 and +0.28. The best single value, 0.39 (SCDC, this project's reimplementation, for B in
  GBM), is the best of eight on n = 56.

**Limits.**
- Eight methods per unit.
- The panel includes this project's reimplementations, labelled as such.
- DECEPTICON's own panel (ten strategies, several templates) differs, so this tests its rule, not
  its package's output.

### 7.3 S2 · Does anatomy see the loss-scale non-identifiability? (Addendum 2, declared before computing)

Declared 2026-10-02 12:55 (file time), after E2's TCGA results were seen and before the run (12:56-13:04). Run on `scripts/remeasure_method.py`'s
exact anatomic inputs: all equivalence conditions pass, and the inputs are built once. Code:
`remeasure_method.py --e2-s2`. Results: `results/identifiability/ivygap_unmix_ablation.json`.

- **Control passed:** the pre-declared setting reproduces the reported row exactly (ACS 0.9077;
  estimates within 1e-16).
- **The pre-declared shift is 1 on Ivy GAP too,** so `s1` duplicates it; there are 6 distinct fits.

| setting | Ivy GAP ACS | constraints that change | TCGA-GBM lymphoid (same setting, §7l) |
|---|---|---|---|
| pre-declared (= shift 1) | 0.908 | -- | T > B; B above T in 24% |
| shift 10 | 0.877 | C7 6/8 | T > B; 33% |
| shift 100 | 0.939 | C7 8/8 | **B > T; 91%** |
| shift 1000 | 0.939 | C6 8/9 | **B > T; 100%** |
| shift 10000 | 0.969 | C6 8/9, C7 8/8 | **B > T; 91%** |
| pre-declared at power 2 | 0.877 | C7 6/8 | B > T; 49% |

- C1–C5 are identical in every setting.
- Every setting beats its permutation null (p ~ 1e-4).

**Reading: INCONCLUSIVE.** The ACS range is 0.092, between the pre-declared 0.05 (SUPPORTED) and
0.10 (NOT SUPPORTED).

**Descriptive, not pre-declared (six fits, no test).** ACS does not favour the settings that reproduce
methylation's lymphoid ordering.
- The two methylation-consistent fits score 0.877 and 0.908.
- The three most inverted fits score 0.939–0.969.
- Spearman between ACS and the TCGA B-above-T share is 0.67.

Anatomy, scored on the stable compartments, moves a little with the loss scale. It moves through the
fine contrasts of C6 (macrophage/microglia) and C7 (the tumour gradient), and not toward getting the
lymphocytes right.

**Stability on Ivy GAP, over 6 distinct fits (descriptive, as declared).** The pattern of E2's TCGA
result repeats on the anatomic cohort:
- **The constrained cell types are stable:** Tumor 0.954, Macrophage_Microglia 0.944, Endothelial
  0.983, Oligodendrocyte 0.965.
- **T, B and NK are less stable:** 0.760, 0.546, 0.551. `unmix` returns exact zeros for 88–93% of
  their per-sample estimates, so their stability rests mostly on ties.

### 7.4 Against a second methylation truth: GIMiCC (2026-10-03; `prespecified/gimicc_truth_confirmation.md` Addenda 3 and 6, written before computing; `results/gimicc_secondary.json`, `results/gimicc_diagnostics.json`)

GIMiCC [52] returns tissue fractions, so its lymphoid units share the estimates' denominator by
construction. That removes check 7's correction. D1 and D2 are unchanged; only the truth column is
replaced (GIMiCC level 3).
- **Reproduction control: exact.** The saved predeclared `unmix` estimates give E2's recorded Tumor
  rhos (0.7801 / 0.5916; median 0.4134 / 0.1830).

| unit | `unmix` vs EpiDISH (E2) | `unmix` vs GIMiCC | methods' median vs GIMiCC |
|---|---|---|---|
| Lymphoid, GBM / LGG | -0.200 / -0.183 | **+0.647 / +0.526** | +0.077 / -0.247 |
| T, GBM / LGG | -0.119 / -0.273 | **+0.643 / +0.491** | +0.180 / +0.096 |
| B, GBM / LGG | +0.090 / -0.001 | +0.253 / +0.296 | -0.126 / -0.219 |
| NK, GBM / LGG | -0.120 / -0.070 | +0.258 / +0.284 | -0.036 / +0.025 |

- **The registered Q3 reading: FAILS.** "The lymphoid estimate does not track methylation" does not
  hold against GIMiCC for `unmix`.
  - That tracking survives controlling for the Thorsson leukocyte fraction: partial **+0.364 / +0.316**
    (D-f). It is more than the leukocyte level.
  - The methods' median does not track (partial +0.121 / -0.047).
- **H1 with GIMiCC's lymphoid truth:** D2 rho **1.00** (exact one-sided p 0.0014); D1 0.771 (p 0.051).
  Truth-free agreement ranks the six primary units exactly as their agreement with this truth.
- **What changes:**
  - The headline of E2 survives, and under this truth it is carried by agreement rather than
    stability.
  - "The lymphoid split tracks nothing" becomes: it tracks nothing against EpiDISH, and the panel's
    median tracks nothing against either instrument. One genuine package (`unmix`) tracks GIMiCC's T
    and lymphoid totals.
  - Which methylation instrument is right is itself open. The two barely agree per sample (T share
    rho 0.08 / 0.11), and they disagree on T versus B in LGG (WHY_B_OVER_T §7n).

## E3 · Gene-resampling stability (registered 2026-10-07 15:16:26, before it was run)

Rule: `prespecified/identifiability_e3_gene_resampling.md`. Script: `scripts/identifiability_e3.py`. Artefact:
`results/identifiability_e3.json`.

- **Method.** Genuine DESeq2 `unmix` at E2's pre-declared setting, re-fit on 10 random halves of the
  marker genes (807 of 1615), GBM and LGG. D3 stability is the mean
  pairwise Spearman correlation between the subset fits.
- **Controls pass.**
  - E2's truth agreement is read as reported (Tumor 0.7801 / 0.5916).
  - The subset fits keep the tumour signal: median rho 0.76
    (GBM) and 0.60 (LGG).
- **H1-E3 (primary, 6 units):** rho +0.714, exact one-sided p 0.068. **INCONCLUSIVE**:
  the same direction as E2, not significant.
  - Secondary, 12 units: rho +0.462 (p 0.066).
  - The two perturbations agree only partly across units: rho +0.60.
- **What differs is the lymphoid total.** Gene resampling rates it stable: 0.89
  (GBM) and 0.76 (LGG). Yet it does not track DNA:
  -0.20 and -0.18. The loss scale
  flagged it: 0.71 and 0.49.
- **Reading for the paper.** Truth-free stability depends on the perturbation, not only on the compartment.
  - The loss-scale perturbation is the informative one here, and E2's claim should be stated for it
    specifically.
  - Gene resampling alone would pass a compartment the DNA says is not determined.

