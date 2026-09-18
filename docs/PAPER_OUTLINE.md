# Paper outline — every number sourced, every limitation stated

**This is the writing scaffold, not a manuscript.** Each claim carries the artefact it comes
from so a number can be checked in one step, and each section carries the limitation that must
travel with it. Placeholders marked **[LGG PENDING]** fill in when the validation run completes.

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
whole-genome doubling and subclonal fraction — does not. **[LGG PENDING]** In anatomically
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

## 7 · Result 4 — external validation in lower-grade glioma

**[LGG PENDING]** — n = 496, IDH1-mutant **78.2%** against GBM's 5.4%. Replication judged on
direction first, then consistency, then significance. `cross_cohort_validation.json`

## 8 · Result 5 — where the failure lives

*Source: `results/mes_by_niche.json`*

MES score across five niches: LE **−0.0394**, IT **−0.0103**, CT **−0.0034**, MVP **+0.0815**,
PAN **+0.1334**. Monotone, PAN highest as predicted.

**The pre-specified test is not significant (p = 0.3312)** — the statistic was "which structure
is the maximum", which is weak at five categories. Post hoc and labelled: PAN − LE positive in
**6 of 6** tumours, p = 0.0142; trend r = **+0.957**. **Suggestive, not established.**

## 9 · Limitations — none of these belong in a footnote

1. **Only the Tumor column has DNA ground truth.** Nothing here speaks to immune accuracy, which
   is what deconvolution is most often used for — and where the pseudobulk arm shows T cells
   over-called **4.6–7.9×**.
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
