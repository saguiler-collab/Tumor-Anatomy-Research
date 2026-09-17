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
