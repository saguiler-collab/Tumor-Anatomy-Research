# When does bulk RNA deconvolution fail in glioblastoma?

**Discovery cohort run 2026-09-17 against `prespecified/biological_failure_factors.md`, which
was committed before any error-versus-factor relationship was examined.**

## The headline, and it is a measured quantity rather than a null

> **Across 154 TCGA-GBM samples with a DNA-measured tumour purity, bulk deconvolution recovers
> a median of 21.5% of the true variation in tumour content. The range across twelve methods is
> −1.1% to 65.6%.**

"Recovery fraction" is defined as `1 + slope`, where the slope regresses `estimate − purity` on
purity. A method that tracks purity perfectly scores 1.0; a method carrying **no information at
all** scores 0.0.

| method | slope | recovers | r with purity |
|---|---|---|---|
| **Bisque** *(degenerate)* | −1.011 | **−1.1%** | **−0.006** |
| BayesPrism | −0.897 | 10.3% | +0.305 |
| NNLS / MuSiC / elastic net | −0.819 | 18.1% | +0.352 |
| SCDC | −0.785 | 21.5% | +0.346 |
| EPIC | −0.756 | 24.4% | +0.205 |
| Bayesian | −0.552 | 44.8% | +0.618 |
| DWLS | −0.540 | 46.0% | +0.675 |
| SVR | −0.385 | 61.5% | +0.729 |
| **CIBERSORTx** | **−0.344** | **65.6%** | **+0.699** |

**Bisque's estimates carry no information about tumour content whatsoever** (r = −0.006). It is
degenerate on this reference and labelled as such — but a reader handed its output would have no
way to know its tumour column is noise.

**This is why the purity factor is reported as a magnitude and not as a significance test.**
Regressing `estimate − purity` on purity is mechanically negative: a method with no information
gives exactly −1. The 12-of-12 consistency is arithmetic, not biology, and is recorded as
definitional rather than deleted. What is informative is *how far from −1* each method sits.

## The one genuinely biological factor that survived

> **Mesenchymal glioblastomas have their tumour content under-called more than other subtypes —
> consistently, across 11 of 12 methods (8 of 9 excluding degenerate ones), Holm-corrected
> p = 0.0381, adjusted for purity and five other factors in a joint model.**

The direction and the mechanism were both **fixed in advance**: mesenchymal GBM is
transcriptionally myeloid-like, so malignant signal should leak into the macrophage column and
depress the tumour estimate. It did.

Unlike purity, this factor is **not** subtracted from the outcome, so there is no definitional
component — and because the model is joint, the coefficient is already adjusted for purity.

The effect is **small in absolute terms**: roughly 0.02 in tumour fraction between MES and
non-MES, against purity's ~0.5 across its range. Consistent and real, an order of magnitude
smaller than the purity effect, and it should be reported that way.

## What failed, at the same prominence

Five of seven candidate factors did not survive. Two failed in the direction **opposite** to
what was predicted, which the pre-specification defines as a failed prediction rather than a
discovery.

| factor | predicted | observed | verdict |
|---|---|---|---|
| **ploidy** | \|error\| increases | 3/12 in the predicted direction (1/9 without degenerate) | **wrong direction** — aneuploidy does not degrade the estimate; if anything it improves it |
| **subclonal fraction** | \|error\| increases | 3/12 | **wrong direction** — subclonal heterogeneity does not degrade the estimate |
| genome doublings | \|error\| increases | 9/12, unadjusted p = 0.146 | not significant after Holm |
| IDH1 mutation | \|error\| increases | 5/12 | **underpowered** — 8 positives of 148 |
| G-CIMP | \|error\| increases | 5/12 | **underpowered** — 8 positives of 151, and 0.868 correlated with IDH1 |

**The ploidy and subclonal results are informative failures.** The pre-specified mechanism —
that aneuploidy and subclonal heterogeneity should break a model assuming one expression profile
per cell type — is **not supported**, and the point estimates lean the other way. Whatever makes
deconvolution fail in GBM, it is not genomic instability.

**The IDH1 and G-CIMP nulls carry no information**, exactly as recorded before the run: eight
positives cannot support a test. This was written down in advance precisely so it could not be
reported as evidence of absence afterwards.

## What this does and does not establish

**Does:** deconvolution error in GBM is partly predictable from tumour biology — specifically
from transcriptional subtype, and not from genomic instability. And the absolute performance
figure, 21.5% of purity variation recovered, is a usable number for anyone deciding whether to
trust a GBM deconvolution.

**Does not:** this is one cohort, one reference, and only the *Tumor* column has DNA ground
truth. Nothing here speaks to immune-cell accuracy. The reference is GBmap-derived for every
method, so "the biology breaks deconvolution" cannot be separated from "the biology is
under-represented in GBmap".

## Robustness

The pre-specified repeat excluding the three degenerate methods (Bisque, MuSiC, SCDC ENSEMBLE)
agrees with the primary analysis on **all seven factors**. No conclusion depends on a method
that had silently reduced to another.

Artefacts: `results/failure_factors_gbm.json`, `results/absolute_purity_per_sample.csv`.

---

# VALIDATION IN LOWER-GRADE GLIOMA — n = 496

**Run 2026-09-17. The factors, directions, model and statistic are those fixed in
`prespecified/biological_failure_factors.md` before the GBM discovery, and nothing was added,
dropped or re-signed after seeing GBM.**

Replication is judged on **direction first**, then consistency across the twelve methods, then
significance. A factor can be significant in both cohorts and still have failed to replicate if
it points the other way.

## The headline quantity replicates, and tightly

| | GBM (n = 147) | LGG (n = 496) |
|---|---|---|
| **median recovery of true purity variation** | **21.5%** | **17.2%** |
| purity coefficient | **−0.13595** | **−0.13386** |
| methods in the predicted direction | **12 / 12** | **12 / 12** |
| Holm-corrected p | 0.0034 | **0.0024** |

**The coefficients agree to two decimal places across two cohorts, two grades and 3.4× the
sample size.** Bulk deconvolution recovers roughly a fifth of the real variation in tumour
content in glioma, and that is now a replicated measurement rather than a single-cohort number.

**Per method**, the ranking survives while every value drops:

| method | GBM | LGG | change |
|---|---|---|---|
| CIBERSORTx | 65.6% | **43.6%** | −22.0% |
| SVR | 61.5% | **41.2%** | −20.3% |
| DWLS | 46.0% | 31.0% | −15.1% |
| Bayesian | 44.8% | 29.4% | −15.4% |
| EPIC | 24.4% | **35.0%** | **+10.6%** |
| SCDC | 21.5% | 11.0% | −10.5% |
| NNLS / MuSiC | 18.1% | 17.2% | −0.9% |
| elastic net | 18.1% | 14.4% | −3.7% |
| BayesPrism | 10.3% | 7.2% | −3.1% |
| **Bisque** | −1.1% | **−0.0%** | +1.0% |

**Bisque recovers nothing in either cohort** (r = −0.000 in LGG). EPIC is the only method that
improves, and no claim is made about why.

**Every other method performs worse in LGG, and that is interpretable two ways.** It may be a
genuine property of lower-grade glioma — lower malignant fraction, different microenvironment —
or simply that the reference is a *glioblastoma* atlas applied to a different disease. The joint
model adjusts for purity but cannot separate these. **The reference explanation is at least as
plausible as the biological one and is not dismissed.**

## The mesenchymal finding did NOT replicate

| | GBM | LGG |
|---|---|---|
| variable | Verhaak Mesenchymal *(pre-specified)* | continuous MES score *(declared secondary)* |
| result | **11/12, Holm p = 0.0381** | **6/12, p = 1.0000 — chance** |

**This is reported as a failure to replicate, not explained away.** Three things are true at once
and all three belong in the paper:

1. **The variable is not the same one.** LGG has no Verhaak class, so the test used the continuous
   MES score declared in advance as a secondary. A secondary failing is weaker evidence than a
   primary failing.
2. **LGG has far less mesenchymal character to detect.** Median MES contrast is **−0.0440** in
   LGG against **+0.0117** in GBM. A factor cannot show an effect across a range the cohort does
   not span.
3. **The GBM effect was small to begin with** — about 0.02 in tumour fraction, an order of
   magnitude below purity's ~0.5. Small effects are exactly what fails to replicate.

**The honest conclusion: the mesenchymal association is GBM-specific, under-powered in LGG, or
both, and these data cannot distinguish those.** It should be presented as a single-cohort
finding awaiting replication, not as a validated mechanism.

## IDH1 is now an INFORMATIVE null — and that is new

In GBM, IDH1 had **8 positives of 148** and its null was recorded in advance as carrying no
information. In LGG it has **388 of 496 (78.2%)**.

**Result: 5/12 methods, β = −0.00236, p = 1.0 — a clean null at full power.**

> **IDH mutation status does not predict bulk deconvolution error**, in a cohort where IDH-mutant
> is the majority class. That is a real negative result rather than an absence of data, and it
> could only be obtained in LGG. It also settles the IDH/G-CIMP axis that GBM could not test.

## Genomic instability: null in both cohorts

Ploidy, genome doublings and subclonal fraction are null in **both** cohorts, and in GBM two of
them pointed the *opposite* way to prediction. **Genomic instability does not break bulk
deconvolution in glioma.** The pre-specified mechanism — that aneuploidy and subclonal
heterogeneity should violate a one-profile-per-cell-type model — is not supported at n = 643
across two cohorts. That is a stronger negative than either cohort alone.

## Summary of the two-cohort study

| claim | status |
|---|---|
| deconvolution recovers ~a fifth of tumour-content variation | **REPLICATED** — 21.5% / 17.2%, coefficients to 2 dp |
| the method ranking by DNA agreement | **REPLICATED** — CIBERSORTx/SVR top, Bisque nothing, both cohorts |
| mesenchymal character predicts error | **DID NOT REPLICATE** — GBM-specific, under-powered, or both |
| genomic instability predicts error | **NULL in both** — mechanism not supported |
| IDH mutation predicts error | **NULL at full power** in LGG — informative |

Artefacts: `failure_factors_gbm.json`, `failure_factors_lgg.json`,
`cross_cohort_validation.json`, `absolute_purity_per_sample{,_lgg}.csv`.

