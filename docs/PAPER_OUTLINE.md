# Paper outline — every number sourced, every limitation stated

**This is the writing scaffold, not a manuscript.** Each claim carries the artefact it comes
from so a number can be checked in one step, and each section carries the limitation that must
travel with it. The LGG validation is complete and filled in; no placeholders remain.

---

## Working title

**Bulk RNA deconvolution recovers only a fifth of tumour-content variation in glioma, and
mesenchymal character predicts the shortfall**

## The one-paragraph abstract shape

Cell-type deconvolution of bulk RNA is widely used to estimate tumour and immune content, but
there is no routine way to know when its output can be trusted. We first asked whether a
tumour's own anatomy could serve as a ground-truth-free check, and found that anatomic
concordance **detects** whether a method is working while failing to **rank** methods against
four independent external yardsticks. We therefore used DNA-measured tumour purity as ground
truth in 154 glioblastomas and found that deconvolution recovers a median of **21.5%** of the
true variation in tumour content, ranging from **−1.1% to 65.6%** across twelve methods.
Deconvolution error is partly predictable from tumour biology: **mesenchymal** transcriptional
character predicts a larger under-call of tumour content, while genomic instability — ploidy,
whole-genome doubling and subclonal fraction — does not. The recovery figure replicates in 496
lower-grade gliomas (coefficient −0.134 against −0.136, 12/12 methods, Holm p = 0.0024) while the
mesenchymal association does not, and IDH mutation is shown at full power not to predict error.
In anatomically
annotated tissue, mesenchymal character is concentrated in the hypoxic perinecrotic niche,
locating the failure spatially.

---

## 1 · Introduction

- Deconvolution is used to estimate immune infiltration and tumour purity from bulk RNA, and
  those estimates inform biological and clinical conclusions.
- **The gap:** benchmarks overwhelmingly use simulated mixtures, which reward methods sharing
  the simulator's assumptions (Avila Cobos 2020; Nguyen 2024). There is no per-sample check.
- **Two questions, in order.** Can plausibility substitute for ground truth? And if not, when
  can any method be trusted?

## 2 · Anatomy as a ground-truth-free check — it detects, it does not rank

*Source: `docs/ENDPOINT.md`, `docs/EXTERNAL_VALIDATION.md`, `results/anatomic/`*

**What works — report first, it is real.**
- Constraints verified in two external single-cell datasets never deconvolved: Albiach, Tumor
  CT>LE **13.34-fold, p = 0.00005**; Darmanis, **4 of 4** gates across 4 patients.
- All 14 methods beat a within-tumour permutation null at **p = 0.0001**;
  `control_random` **0.3385** and `control_shuffled_signature` **0.1385** do not — under **every**
  reference atlas tested.

**What fails — four independent lines.**
| line | result | artefact |
|---|---|---|
| published ground-truth benchmarks | Mann-Whitney **p = 0.50**; DWLS top-tier there, 14th of 14 here | `benchmark_concordance.json` |
| imaging mass cytometry, **pre-registered** | **p = 0.0090**, wrong direction | `imc_anchored_test.json` |
| reference atlas | only **+0.3655 to +0.5099** of the ordering survives | `reference_sensitivity_*.json` |
| marker-set size | ACS moves **0.3846 → 0.6923** on nothing else | `imc_anchored_robustness.json` |

**Two mechanism results that explain why.**
- **ACS is a rank statistic, so it is blind to magnitude.** Restoring EPIC's variance weighting
  moves its tumour estimate by **0.0829** on average and **0.3539** at maximum, and moves its ACS
  by **exactly 0.0000**. `epic_close_out.json`
- **It prefers the reference that violates the mixing model.** Correcting a log-transformed
  reference lowers ACS for **all 14** methods, mean **−0.12**, and holds when each reference
  picks its own markers at matched density. `reference_sensitivity_gbmap_linear_per_arm.json`

**The registered outcome, both ways.** rho = **0.7501** against a GBmap-sharing yardstick;
**+0.0698**, CI [−0.510, +0.769], against DNA purity. Reported **INCONCLUSIVE** on twelve
methods — the interval contains 0.7501, so this is a failure to reproduce, not a refutation.
`yardstick_agreement.json`

> **This section motivates the rest.** If the evaluation cannot say *which* method to trust, ask
> *when* any of them can be.

## 3 · Methods

- Cohorts: Ivy GAP (122 H&E-labelled samples, 10 tumours), TCGA-GBM (**154** with a DNA purity
  call), TCGA-LGG (**510**; **496** with IDH1 status).
- Twelve comparable methods, one frozen GBmap-derived signature, 1,615 markers (200/type).
- Ground truth: ABSOLUTE purity from DNA copy number, PanCanAtlas.
- **Two input files named "counts" contained log2 values** and were inverted before use;
  detected by a scale gate, not assumed. `tcga_*_bulk_provenance.json`
- Pre-specification: `prespecified/biological_failure_factors.md`, committed before any
  error-versus-factor relationship was examined.

## 4 · Result 1 — how much of tumour content is recovered at all

*Source: `results/failure_factors_gbm.json`*

Median **21.5%** of true purity variation recovered; range **−1.1%** (Bisque, r = −0.006) to
**65.6%** (CIBERSORTx). **Bisque's tumour column carries no information while returning a
confident-looking number.**

> **CORRECTED 2026-09-18 — the figure to quote is 33.1% / 23.3%, not 21.5% / 17.2%.**
> The medians above pool **all** methods, including four that were never evaluated under their
> intended inputs (Bisque needs paired bulk/single-cell subjects; MuSiC needs sigma; EPIC needs
> `refProfiles.var`; SCDC ENSEMBLE needs a second reference). Restricting to the **8 comparable**
> methods gives **median recovery 33.1% in GBM and 23.3% in LGG**. Including un-evaluable methods
> understated the figure by about 12 points. `docs/EQUAL_FOOTING.md`,
> `ivygap/deconv/comparability.py`.


**State plainly:** regressing (estimate − purity) on purity is mechanically negative, so the
*consistency* across methods is arithmetic. The **magnitude** is the finding.

## 5 · Result 2 — mesenchymal character predicts the shortfall

Verhaak **Mesenchymal** GBMs are under-called more: **11 of 12** methods (8 of 9 excluding
degenerate), **Holm p = 0.0381**, adjusted for purity and five other factors jointly. Direction
and mechanism pre-specified — MES is transcriptionally myeloid-like, so malignant signal leaks
into the macrophage column. Effect **~0.02**, an order of magnitude below purity's ~0.5, and
reported as such.

## 6 · Result 3 — what does NOT predict error, at equal prominence

**Ploidy** and **subclonal fraction** failed in the direction **opposite** to prediction (3/12
each). **Genomic instability does not break deconvolution in GBM**, and the pre-specified
mechanism was wrong. Genome doublings did not survive Holm. IDH1 and G-CIMP nulls in GBM carry
**no information** — 8 positives each, recorded in advance so they could not be reported as
evidence of absence.

## 7 · Result 4 — external validation in lower-grade glioma (n = 496)

*Source: `results/failure_factors_lgg.json`, `results/cross_cohort_validation.json`*

**Replicated.** Median recovery **33.1% → 23.3%** across the 8 comparable methods (21.5% → 17.2%
if un-evaluable methods are pooled in, which understates it); purity coefficient **−0.13595 → −0.13386**,
**12/12** methods both times, Holm **p = 0.0024**. Coefficients agree to two decimal places across
two grades and 3.4× the sample size. The method ranking survives — CIBERSORTx and SVR top, Bisque
recovering nothing (r = −0.000).

**Did not replicate: mesenchymal character.** 6/12, p = 1.0000. Give all three reasons — the
variable is the declared secondary and not the pre-specified one, LGG spans far less MES range
(median −0.0440 against +0.0117), and the GBM effect was ~0.02 to begin with. **Present MES as a
single-cohort finding awaiting replication, not a validated mechanism.**

**A new informative null.** IDH1 had 8 positives in GBM (recorded in advance as uninformative) and
**388 of 496** in LGG. Result 5/12, p = 1.0 — **IDH mutation does not predict deconvolution
error**, at full power, where it is the majority class.

**Genomic instability null in both** (n = 643 combined), with two GBM factors pointing *opposite*
to prediction. A stronger negative than either cohort alone.

**Every method does worse in LGG, and the reference explanation is as plausible as the biological
one** — a glioblastoma atlas applied to a different disease. Do not claim otherwise.

## 7b · Result 4b — the second ground truth, and a retraction

*Source: `results/immune_arm.json`, `docs/IMMUNE_ARM.md`*

Leukocyte fraction from **DNA methylation** (Thorsson 2018) — a third molecule — on 141 TCGA-GBM
samples. **Both pre-specified predictions failed.**

**Methods UNDER-call immune content** (median error −0.1017; 3 of 12 positive, p = 0.146), and the
under-call **survives** the pre-specified mRNA-content correction (median **0.661×**) that was
declared in advance to work against this conclusion.

**The 4.6–7.9× T-cell over-call is a simulator artefact.** On tissue: **0.0002–0.0085** against
pseudobulk 0.202–0.350. Withdrawn as a clinical claim — the sixth withdrawal in this project, and
the first corrected by a molecule the project had never used.

**The two arms' biases are anti-correlated, ρ = −0.8028, p = 0.0017** — largely a compositional
necessity (estimates sum to one) and reported as such, not as a discovery. Its real force is that
the errors are *not independent*, so no method can be judged on one compartment alone.

**The positive result: CIBERSORTx and SVR are best calibrated against BOTH molecules** —
65.6%/61.5% tumour recovery, 1.03×/0.87× immune fold. A method ranking with orthogonal support,
which the anatomy score could not produce. **Bisque** carries no tumour information (r = −0.006)
while over-calling immune 2.27×.

**Unpredicted observation:** methods place **more B cells than T cells** in GBM (CIBERSORTx 0.046
vs 0.003) while the reference holds 5× more T than B. An identifiability failure, reported as
needing its own test since no direction was pre-specified for it.

## 7c · Result 4c — the anomaly was pre-registered, tested, and CONFIRMED

*Source: `results/lymphoid_ordering_lgg.json`, `results/methylation_celltypes_lgg.json`,
`docs/IMMUNE_ARM.md` §6*

The B-over-T observation in §7b had no per-type truth to check it against. One was then obtained,
and **the prediction was registered before it was computed**, with an explicit falsifier:
*"methylation showing B ≥ T [...] That outcome would be reported and the anomaly withdrawn."*

**Per-type truth:** EpiDISH RPC on `centDHSbloodDMC.m` (**255 of its 333 reference HM450 CpGs**;
64 are all-NA in this matrix and drop out of `complete.cases`), **530 TCGA-LGG
samples**. Both sides renormalised within {T, B, NK}, so the leukocyte-subcomposition vs
all-cell-fraction denominator difference cancels by construction — verified in
`tests/test_lymphoid_ordering.py`, not merely asserted.

**The prediction held.** Methylation: **T 0.476 > NK 0.321 > B 0.203**; T ranks first
in **79.2%** of samples, B in **4.0%**.

**Every method disagrees. 12 of 12 place B above T** — 8 of 8 when restricted to
methods comparable on this reference, so it is not an artefact of degraded stand-ins. **0 of 12
reproduce the true ordering.** `elastic_net` puts **97.7%** of lymphoid signal in B; `music` and
`nnls` return **exactly 0.0000** for T in all 510 samples.

This is a **misassignment, not a proportion error**: lymphocyte signal is assigned to the wrong
column. **The obvious explanation — that the signature cannot separate T from B — was tested and
refuted** (`scripts/identifiability_probe.py`): on mixtures built from the reference itself, plain
NNLS recovers T:B = 2.33 exactly, still recovers T > B at 100% multiplicative noise, and still
recovers it with an entire cell type deleted from the reference; the signature's condition number
is 5.3. The cause lies in the gap between the reference's expression space and real bulk tissue,
and **is not identified**. Its being unanimous across NNLS, SVR, Bayesian and probabilistic-model families
makes it a property of the problem as posed rather than a quirk of one solver — and it
**replicates the GBM observation in a second tissue**, which is what turns an anomaly into a
finding.

**A mechanism was proposed and rejected.** The frozen signature gives T_cell only **5.17%** of
profile mass against B_cell's **13.26%**, and B_cell is the immune profile most correlated with
Macrophage_Microglia (**r = +0.552**). That suggested macrophage spillover into B. The prediction
— B tracking Macrophage more than T does — was tested and **failed: 5 of 10, p = 0.623**. The
explanation is **not adopted**, and is reported as rejected rather than dropped.

## 7d · Result 4d — equal footing changes the winner, not the shortfall

*Source: `results/equal_footing_ranking.json`, `docs/EQUAL_FOOTING.md`*

Four methods in the panel were being measured as something other than themselves: the vendored
signature carries an all-zero sigma, so MuSiC is arithmetically NNLS, EPIC is uniform-weighted,
SCDC ENSEMBLE has nothing to ensemble, and CIBERSORTx S-mode cannot run. Rebuilding the reference
from the atlas's own counts (`matrix="raw/X"`, D16) supplies all of it.

**The ranking does not survive it.** Kendall tau between the two rankings is **+0.214** on the 8
methods rankable under both (**+0.143** excluding the single method that also changed
implementation). **MuSiC goes from excluded-as-NNLS (17.5%) to first at 60.2%**; DWLS falls rank
3 → 12. CIBERSORTx and SVR hold ranks 1→2 and 2→3, so the calibration claim in §7b stands.

**The shortfall does survive it.** Median bias moves from −0.028 to −0.472 and the methods
under-calling tumour go from 7 of 12 to **12 of 12**, against a mean true purity of 0.752.
Recovery is bias-invariant by construction (tested), so these are independent findings, not one
finding seen twice.

**Stated as one sentence:** *giving every method the inputs it was designed to consume changes
which method wins and makes the under-call worse.* Why the calibration moves is **not
established** — the log-vs-counts difference is a plausible cause, not a measured one, and
neither reference is endorsed as correctly calibrated.

## 8 · Result 5 — where the failure lives

*Source: `results/mes_by_niche.json`*

MES score across five niches: LE **−0.0394**, IT **−0.0103**, CT **−0.0034**, MVP **+0.0815**,
PAN **+0.1334**. Monotone, PAN highest as predicted.

**The pre-specified test is not significant (p = 0.3312)** — the statistic was "which structure
is the maximum", which is weak at five categories. Post hoc and labelled: PAN − LE positive in
**6 of 6** tumours, p = 0.0142; trend r = **+0.957**. **Suggestive, not established.**

## 9 · Limitations — none of these belong in a footnote

1. ~~Only the Tumor column has DNA ground truth.~~ **RESOLVED** — the immune compartment now has
   its own independent ground truth, methylation-derived leukocyte fraction (§7b). It produced a
   **retraction**: the pseudobulk 4.6–7.9× T-cell over-call does not exist on tissue, where the
   immune compartment is *under*-called. What remains unresolved is **per-type** immune truth;
   leukocyte fraction is one aggregate number.
2. **One reference for every method.** "The biology breaks deconvolution" cannot be separated
   from "the biology is under-represented in GBmap".
3. **Purity confounds with almost everything**; the joint model adjusts only linearly.
4. **ABSOLUTE purity is itself an estimate** from copy number, not a gold standard.
5. **The spatial layer is nine tumours** and its pre-specified test is null.
6. **The reference is log-transformed** in the archived runs (D16); its direction is measured
   but the full linear re-run is not done.

## 10 · What makes this credible, and it should be said

Five claims were published during this work and then **withdrawn after measurement**; the
retractions are kept in place. Three predictions were pre-registered, and **one failed in the
direction named in advance as most damaging**. The constraint file was frozen and hashed before
any output was seen and was never edited — including when one of its own constraints failed an
external check.

**A result produced under those conditions is worth more than a positive one produced without
them.**
