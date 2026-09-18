# The immune arm: a second ground truth, and a retraction of this project's own claim

**Run 2026-09-18 against `prespecified/immune_failure_factors.md`, committed before any immune
error existed. Both pre-specified predictions FAILED.**

Ground truth: **leukocyte fraction from DNA methylation** (Thorsson et al. 2018, *Immunity*
48(4)) — 2,000 differentially methylated loci, mixture model, 141 of 174 TCGA-GBM samples. A
third molecule, independent of the RNA the methods consume and of the copy number the tumour arm
uses.

---

## 1 · P1 FAILED: methods UNDER-call immune content on tissue

Predicted: over-call, median error > 0. **Observed: 3 of 12 methods positive, p = 0.146, median
error −0.1017.**

| method | immune estimate | leukocyte fraction | fold | after mRNA correction |
|---|---|---|---|---|
| Bisque | 0.428 | 0.189 | 2.27× | 2.46× |
| Bayesian | 0.318 | 0.189 | 1.68× | 2.05× |
| **CIBERSORTx** | **0.195** | 0.189 | **1.03×** | 1.25× |
| **SVR** | **0.165** | 0.189 | **0.87×** | 1.06× |
| DWLS | 0.123 | 0.189 | 0.65× | 0.86× |
| SCDC | 0.087 | 0.189 | 0.46× | 0.66× |
| elastic net | 0.086 | 0.189 | 0.46× | 0.54× |
| NNLS / MuSiC | 0.053 | 0.189 | 0.28× | 0.34× |
| EPIC | 0.052 | 0.189 | 0.28× | 0.31× |
| **BayesPrism** | 0.031 | 0.189 | **0.16×** | **0.20×** |

**The under-call survives correction.** The pre-specification named two biases that push the
estimate downward — the roster drops 7.0% of the atlas including some leukocytes, and immune cells
carry ~2.245× less mRNA than tumour cells. Converting mRNA share to cell share with the
**pre-specified** vector (`prespecified/mrna_content_gbmap_10x.csv`) leaves **4 of 12** methods
over-calling and a **median 0.661×**. Methods still under-call by about a third after the
correction that was declared in advance to work against this conclusion.

---

## 2 · The 4.6–7.9× T-cell over-call is a SIMULATOR ARTEFACT

This project's most-quoted clinical claim was that every method over-calls T cells by 4.6–7.9× in
high-purity tumour, with the missing tumour mass becoming T cells. That came from **synthetic
pseudobulk built from GBmap's own cells**. On real tissue:

| method | T_cell on TCGA-GBM tissue | T_cell in pseudobulk |
|---|---|---|
| EPIC | **0.0002** | — |
| MuSiC / NNLS | **0.0006** | 0.202 |
| SVR | 0.0008 | — |
| BayesPrism | 0.0009 | 0.234 |
| elastic net | 0.0010 | — |
| Bayesian / CIBERSORTx | 0.0026 | 0.339 |
| DWLS | 0.0085 | — |
| SCDC | 0.0279 | 0.350 |
| Bisque | 0.0774 | — |

> **On tissue these methods place T cells at 0.0002–0.0085 — two orders of magnitude below their
> pseudobulk estimates of 0.202–0.350.** The over-call does not exist outside the simulator.

**This is withdrawn as a clinical claim.** It was real *in the benchmark* and the benchmark was
measuring its own construction: in pseudobulk, T cells absorbed the mass that methods failed to
assign to tumour. Given a real tumour they do not.

**It is the sixth claim withdrawn in this project after measurement**, and the first one where the
correcting evidence came from a molecule the project had never used before.

---

## 3 · P2 FAILED — but the relationship it looked for exists in a different form

Predicted: signs opposite between arms — tumour under-called, immune over-called. **Observed: 3 of
12.** Not a systematic asymmetry.

**Instead the two biases are strongly anti-correlated across methods: Spearman −0.8028,
p = 0.00167.** A method that over-calls tumour under-calls immune, and vice versa.

**This is mostly a compositional necessity and is reported as one, not as a discovery.** The
estimates sum to one, so mass assigned to tumour is mass unavailable to everything else. What it
does establish is that **the two errors are not independent**, so a method cannot be judged on one
compartment alone — which is an argument for having two ground truths rather than a finding
produced by them.

---

## 4 · The cross-validated ranking — the positive result

**CIBERSORTx and SVR are the best-calibrated methods against BOTH molecules.**

| | tumour recovery (copy number) | immune fold (methylation) |
|---|---|---|
| **CIBERSORTx** | **65.6%** | **1.03×** |
| **SVR** | **61.5%** | **0.87×** |
| DWLS | 46.0% | 0.65× |
| BayesPrism | 10.3% | 0.16× |
| **Bisque** | **−1.1%** | 2.27× |

Two independent DNA-based measurements, two compartments, and the same two methods come out on
top. **That is a method ranking with orthogonal support** — and it is exactly what the ACS
ordering could not produce.

**Bisque remains the clearest warning in the panel:** no tumour information at all (r = −0.006)
while over-calling immune content 2.27×.

---

## 5 · NEW, and not predicted: methods place more B cells than T cells in GBM

Unprompted by any hypothesis, and visible only once per-type estimates were persisted:

| | T_cell | B_cell |
|---|---|---|
| **GBmap reference composition** | **0.165** | **0.032** |
| CIBERSORTx estimate | 0.0026 | **0.0462** |
| Bayesian estimate | 0.0026 | **0.1444** |
| MuSiC / NNLS estimate | 0.0006 | 0.0048 |

**The reference contains 5× more T cells than B cells. Most methods estimate the reverse.** In
glioblastoma, where T cells are the dominant lymphocyte and B cells are rare, this is backwards
against both the biology and the reference the methods were given.

This is an **identifiability failure**, not a proportion error: the solvers are placing lymphocyte
signal in the wrong column. It was not pre-specified and no direction was fixed for it, so it was
reported here as an observation requiring its own test.

**That test was then pre-registered and run. See section 6 — the anomaly is confirmed against
DNA methylation and replicates into a second tissue.**

---

## 6 · The B-over-T anomaly is CONFIRMED against an orthogonal per-type measurement

Section 5 recorded the B-over-T inversion as "an observation requiring its own test, not an
established mechanism". That test was then **pre-registered before it was run**
(`prespecified/immune_failure_factors.md`, per-cell-type addendum), with the prediction fixed in
advance and an explicit falsifier:

> **PREDICTION: DNA methylation will show T cells > B cells.**
> **What would falsify it:** methylation showing **B ≥ T**. That would mean the methods may be
> right and my "anomaly" was an artefact [...] **That outcome would be reported and the anomaly
> withdrawn.**

### The result: the prediction held

DNA methylation (EpiDISH RPC, `centDHSbloodDMC.m`) on **530 TCGA-LGG samples**. EpiDISH ships
333 reference CpGs; **64 are all-NA across every sample in this matrix**, so `epidish()` — which
is given `b[complete.cases(b), ]` — actually ran on **255 of 333**. That is the number that
belongs in the sentence, and it is read from EpiDISH's own stdout rather than asserted:

**T > B in 93.2% of samples** — mean T 0.4769 vs B 0.2040.

### And every method disagrees with it

| method | est. T | est. B | calls higher |
|---|---|---|---|
| bayesian | 0.0237 | 0.8330 | **B** |
| bayesprism | 0.3038 | 0.4265 | **B** |
| bisque | 0.1511 | 0.3412 | **B** |
| cibersortx | 0.0316 | 0.8787 | **B** |
| dwls | 0.4656 | 0.4842 | **B** |
| elastic_net | 0.0108 | 0.9768 | **B** |
| epic | 0.2854 | 0.4208 | **B** |
| music | 0.0000 | 0.9575 | **B** |
| nnls | 0.0000 | 0.9575 | **B** |
| scdc | 0.2793 | 0.6521 | **B** |
| scdc_ensemble | 0.2793 | 0.6521 | **B** |
| svr | 0.0008 | 0.9356 | **B** |
| **methylation (truth)** | **0.4769** | **0.2040** | **T** |

**12 of 12 methods put B above T. The orthogonal measurement puts T above B in
93.2% of samples.** The anomaly is not withdrawn; it is confirmed, and it
replicates from GBM into a second tissue.

### Why this comparison is like-for-like

Both sides are renormalised **within the same three types** (T_cell, B_cell, NK_cell) before
comparison (`scripts/methylation_celltypes.py`). The methylation estimate is a leukocyte
sub-composition and the deconvolution estimate is a fraction of all cells, but after
renormalisation within the lymphoid triple both are *relative composition among lymphocytes*.
**The denominator difference is therefore removed by construction**, which is what makes the
directional comparison legitimate at all.

### What this is, and what it is not

This is a **misassignment**, not a proportion error. The solvers are not merely mis-scaling
lymphocyte abundance; they are assigning lymphocyte signal to the **wrong column**. A method can
be given a reference holding 5× more T than B, run on tissue whose methylation says T outnumbers
B roughly 2.3:1, and still return B as the dominant lymphocyte — unanimously, across twelve
methods spanning NNLS, SVR, Bayesian, and probabilistic-model families.

That it is **unanimous across method families** is the important part. It is not a quirk of one
solver's regularisation.

> **CORRECTED 2026-09-18.** This paragraph previously called the result "a property of the
> problem as posed — a direction in gene space along which T-cell and B-cell signal are not
> separable". **That explanation was tested and is wrong**, and it is withdrawn. See below.

### The "not separable" explanation was tested and REFUTED

"The signature cannot separate T from B" is a claim about the **matrix**, testable with no tissue
at all: build mixtures *from the reference itself* with a known T:B ratio and see whether a solver
recovers it. `scripts/identifiability_probe.py` does exactly that, with the planted truth chosen
to match the tissue (T:B = {d['planted_T_over_B']}, the methylation value) rather than to produce an outcome.

**The signature separates T from B in every condition tested:**

| condition | est T | est B | T:B | verdict |
|---|---|---|---|---|
| multiplicative noise CV = 0.0 | 0.0700 | 0.0300 | 2.33 | recovered |
| multiplicative noise CV = 0.01 | 0.0701 | 0.0301 | 2.33 | recovered |
| multiplicative noise CV = 0.05 | 0.0702 | 0.0302 | 2.33 | recovered |
| multiplicative noise CV = 0.1 | 0.0698 | 0.0302 | 2.31 | recovered |
| multiplicative noise CV = 0.25 | 0.0671 | 0.0299 | 2.24 | recovered |
| multiplicative noise CV = 0.5 | 0.0627 | 0.0386 | 1.63 | recovered |
| multiplicative noise CV = 1.0 | 0.0711 | 0.0537 | 1.32 | recovered |
| reference missing `Tumor` | 0.2272 | 0.0878 | 2.59 | recovered |
| reference missing `Macrophage_Microglia` | 0.1158 | 0.0551 | 2.10 | recovered |
| reference missing `NK_cell` | 0.0966 | 0.0343 | 2.81 | recovered |
| reference missing `Endothelial` | 0.0707 | 0.0315 | 2.24 | recovered |
| reference missing `Oligodendrocyte` | 0.0779 | 0.0279 | 2.79 | recovered |
| reference missing `Astrocyte` | 0.0690 | 0.0288 | 2.40 | recovered |

Plain NNLS recovers **T:B = 2.33 exactly** on noiseless data, still recovers T > B at **100%
multiplicative noise**, and still recovers it when an entire cell type is **deleted from the
reference** — the standard explanation for misassigned signal. The condition number of the
eight-type signature is **5.29**, which is well-conditioned, not marginal.

**So the inversion is not caused by the signature's conditioning, not by noise, and not by a
missing cell type.** Whatever produces it lives in the gap between this reference's expression
space and real bulk tissue — platform and normalisation differences between a single-cell-derived
profile and TCGA bulk, or genuine glioma biology the reference does not represent. **That cause is
not identified, and this document does not claim one.**

This matters practically: the failure is invisible to the checks a careful analyst would actually
run. You cannot catch it by inspecting the condition number of your signature matrix, by adding
noise, or by worrying about a missing cell type. It only appears when you have an orthogonal
per-type measurement to check against — which is the argument for obtaining one.

### CORRECTED 2026-09-18 — "12 of 12" was a statement about MEANS, and the denominator matters

The count above is computed on per-method **means**. For several methods that mean is taken over
a small minority of the cohort, because `renormalise` sends a sample with no lymphoid signal at
all to NaN and `.mean()` skips it. Reporting "12 of 12 put B above T" without that denominator
overstates a clean result and hides a second, more severe failure. Separated properly, there are
**two distinct failure modes**:

#### Mode 1 — ABSENCE. Four methods report no lymphocytes at all in most samples.

These return **exactly 0.0000** for T_cell, B_cell *and* NK_cell — verified as exact zeros, not
small values near zero:

| method | GBM: samples with zero lymphoid | LGG: samples with zero lymphoid |
|---|---|---|
| `bayesprism` | 51 / 56 (91.1%) | 409 / 510 (80.2%) |
| `elastic_net` | 48 / 56 (85.7%) | 286 / 510 (56.1%) |
| `music` | 55 / 56 (98.2%) | 443 / 510 (86.9%) |
| `nnls` | 55 / 56 (98.2%) | 443 / 510 (86.9%) |

For these methods the B-versus-T question is **vacuous**. They are not placing B above T; they
are reporting that the tumour contains no lymphocytes. Against tissue whose methylation says T
cells are the dominant lymphocyte in ~95% of samples, **that is a more severe failure than
getting the ratio backwards**, and it was invisible for as long as the comparison was made on
means.

`epic` is a near-miss in the same direction: no exact zeros, but a median lymphoid sum of
**7.6e-05** in LGG.

#### Mode 2 — MISASSIGNMENT. Among methods that do return lymphocytes, B is placed above T.

Per-sample, paired within each sample rather than compared as two cohort averages:

| method | GBM: B>T | GBM discordant | LGG: B>T | LGG discordant |
|---|---|---|---|---|
| `bayesian` | 100.0% of 56 | 94.6% | 100.0% of 510 | 93.5% |
| `svr` | 90.2% of 51 | 84.3% | 96.5% of 482 | 90.9% |
| `scdc` | 67.9% of 56 | 64.3% | 93.5% of 510 | 87.6% |
| `scdc_ensemble` | 67.9% of 56 | 64.3% | 93.5% of 510 | 87.6% |
| `cibersortx` | 89.6% of 48 | 83.3% | 90.3% of 462 | 85.3% |
| `epic` | 48.2% of 56 | 46.4% | 57.2% of 510 | 54.7% |
| `bisque` | 40.0% of 50 | 40.0% | 55.3% of 468 | 53.2% |
| `dwls` | 54.0% of 50 | 50.0% | 52.0% of 323 | 47.7% |

**8 of 8 in LGG and 6 of 8 in GBM place B above T on a majority of the samples they
score.** "Discordant" is the fraction of samples where the method says B>T *and* methylation says
T>B in that same sample — the honest pairing, since it never credits a method for agreeing on a
sample methylation also got the other way.

#### What survives, stated exactly

- **The pre-registered prediction held.** Methylation puts T above B in **94.8%** (GBM) and
  **93.2%** (LGG) of samples. Not withdrawn.
- **No method reproduces it.** Zero of twelve, in either cohort, on the mean; and among the
  eight that return lymphoid signal at all, the majority place B above T per-sample.
- **The finding replicates across two cohorts** with independent methylation matrices, and the
  true ratio is *larger* in GBM (T:B = 4.76) than in LGG (2.34).
- **What is no longer claimed:** that this is one clean unanimous inversion. It is two failures —
  four methods find no lymphocytes, the rest put them in the wrong column — and conflating them
  made the result look tidier than it is.

### Does equal footing FIX the inversion? Partly — and it makes the other failure worse

The inversion above was measured against the frozen signature, where four methods were degraded
stand-ins. So the obvious question is whether giving every method its intended inputs repairs it.
Measured on GBM, same 155-sample methylation truth, same 56 matched samples, changing only the
reference (`--reference h5ad`, sigma + donor profiles + registered cells):

| | frozen signature | h5ad (sigma) |
|---|---|---|
| methods scored | 12 | 14 |
| **agree with methylation on T > B** (mean) | **0 of 12** | **3 of 14** |
| reproduce the full T>NK>B ordering | 0 of 12 | 2 of 14 |
| **return ZERO lymphoid signal in most samples** | **4 of 12** | **6 of 14** |

**Both things move, in opposite directions.**

**The misassignment partly repairs.** Zero methods agreed with methylation on the frozen
signature; **3 do once the reference carries what they were designed to consume**, and two
reproduce the full three-way ordering that nothing reproduced before. The clearest case is
`bisque`, which goes from placing B above T on 40.0% of samples to **3.6%** — i.e. it gets the
lymphoid order right on 96.4% of samples once given donor structure. (`bisque` remains
non-comparable for the tumour-content ranking, because TCGA cannot supply its paired-subject
requirement; that exclusion is about a different question and does not apply here.)

**The absence gets worse.** Methods returning *exactly zero* T, B and NK across most samples rise
from **4 of 12** to **6 of 14** — now including `scdc`, `scdc_ensemble` and `svr`, which had
returned lymphoid signal on every sample under the frozen signature. `svr` goes from scoring all
56 samples to **13**.

So equal footing does not simply improve matters. It trades one failure for another: **fewer
methods put lymphocytes in the wrong column, and more methods report no lymphocytes at all.**

**This is not a comfortable result and it is not being smoothed.** It means the lymphoid
compartment of a glioma is not reliably recoverable by any method in this panel under *either*
reference — the failure merely changes shape. And it is a caution against the reading that the
h5ad arm is simply "the fixed version": on the question the previous section asked, it is better;
on whether a method reports lymphocytes at all, it is worse.

Artefacts: `results/lymphoid_ordering.json` (frozen), `results/lymphoid_ordering_h5ad.json`.

### Secondary, NOT registered: the full lymphoid ordering

The same measurement yields NK at no extra cost, so the three-way ordering is reported here.
**It was not pre-registered, no direction was fixed for it, and it is therefore exploratory** —
only T > B is confirmatory. It is included because it shows the inversion is not a narrow T/B
quirk.

Methylation truth within {T, B, NK} on 530 samples:
**T 0.4758 > NK 0.3208 > B 0.2034**. Per sample, **T ranks first in
79.2%** and **B ranks first in 4.0%**.

| method | T | NK | B | ordering |
|---|---|---|---|---|
| bayesian | 0.0237 | 0.1433 | 0.8330 | B>NK>T |
| bayesprism | 0.3038 | 0.2697 | 0.4265 | B>T>NK |
| bisque | 0.1511 | 0.5077 | 0.3412 | NK>B>T |
| cibersortx | 0.0316 | 0.0897 | 0.8787 | B>NK>T |
| dwls | 0.4656 | 0.0502 | 0.4842 | B>T>NK |
| elastic_net | 0.0108 | 0.0124 | 0.9768 | B>NK>T |
| epic | 0.2854 | 0.2939 | 0.4208 | B>NK>T |
| music | 0.0000 | 0.0425 | 0.9575 | B>NK>T |
| nnls | 0.0000 | 0.0425 | 0.9575 | B>NK>T |
| scdc | 0.2793 | 0.0687 | 0.6521 | B>T>NK |
| scdc_ensemble | 0.2793 | 0.0687 | 0.6521 | B>T>NK |
| svr | 0.0008 | 0.0636 | 0.9356 | B>NK>T |
| **methylation (truth)** | **0.4758** | **0.3208** | **0.2034** | **T>NK>B** |

**0 of 12 methods reproduce the true ordering. Every method places B first.**

Two of them place it there almost totally: `elastic_net` puts **97.7%** of lymphoid signal in B
and `music`/`nnls` put **95.8%** there while returning **exactly 0.0000** for T. A solver
returning a hard zero for the dominant lymphocyte in the tissue is not making a quantitative
error — the T-cell column is not identifiable for it at all.

Restricted to the 8 methods comparable on this reference (`ivygap/deconv/comparability.py`),
the count is unchanged: **0 agree on T > B**. The inversion is not an artefact of
including degraded stand-ins.

Artefact: `results/lymphoid_ordering_lgg.json`. Verified by
`tests/test_lymphoid_ordering.py`, which includes an inverted-signal control and a test that
adding an arbitrarily large Tumor column does not move the lymphoid ordering — i.e. the
denominator really is removed by construction rather than merely claimed to be.

### A mechanism was proposed, tested, and REJECTED

The inversion invites an explanation, and the frozen signature offers an obvious one. On the
1,615-gene marker space the eight profiles are **not** balanced in magnitude, even though marker
*ownership* is balanced by construction (200 genes each):

| cell type | share of profile mass |
|---|---|
| Oligodendrocyte | 22.17% |
| Astrocyte | 16.93% |
| Endothelial | 15.26% |
| **B_cell** | **13.26%** |
| NK_cell | 12.01% |
| Macrophage_Microglia | 10.65% |
| **T_cell** | **5.17%** |
| Tumor | 4.55% |

**The T_cell column carries 2.6× less mass than B_cell** and is the weakest immune column in the
matrix. And `B_cell` is the immune profile most correlated with `Macrophage_Microglia`
(**r = +0.552**, the largest off-diagonal among immune types) — the dominant immune population in
glioma.

That suggested a mechanism: **B absorbs macrophage/microglia spillover.** It makes a testable
prediction — across samples, the B estimate should track the Macrophage estimate more closely
than the T estimate does — so it was tested rather than asserted.

**It failed.** 5 of 10 methods show r(B,M) > r(T,M); sign-test **p = 0.623**. That is a
coin flip.

**The spillover explanation is therefore not adopted.** The inversion is a measured fact; the
profile-mass asymmetry is a measured fact about the reference; the causal link between them is
**not established by these data**, and is recorded here as a rejected hypothesis rather than
quietly dropped. Whether the inversion is a property of *this reference* or of the deconvolution
problem itself is tested directly by rebuilding the reference from `gbmap_core.h5ad` — that is
the `--reference h5ad` arm, and it is the decisive experiment, not this correlation.

One incidental confirmation did come out of it: `music` and `nnls` return **zero variance** in
T_cell across all 510 samples (the correlation is undefined), which independently confirms the
hard-zero reading above — for those solvers the T-cell column is not merely small, it is never
used.

Artefact: `results/spillover_lgg.json`.

### What "the reference holds 5× more T than B" does and does not mean

Three different quantities get called "the reference composition", and the argument above depends
on which one is meant. All three are measured here so the claim cannot drift:

| quantity | T_cell | B_cell | T:B |
|---|---|---|---|
| **GBmap core atlas**, all 338,564 cells | 0.1724 (54,257) | 0.0040 (1,250) | **43.4×** |
| **Balanced subsample** the reference is built from (15,311 cells, cap 50/donor/type) | **0.1653** (2,531) | **0.0320** (490) | **5.17×** |
| **Profile mass** on the 1,615-gene marker space, within {T, B, NK} | 0.3027 | 0.3350 | **0.90×** |

The **5.17×** quoted above is the middle row — the subsample, which is the right one, because
those are the cells whose mean expression became the reference the methods were handed.

Two things follow, and they pull in opposite directions:

1. **The comparison is conservative.** The atlas itself holds T cells **43×** more abundant than
   B cells; the balanced sampler caps cells per (donor, type) and compresses that to 5.17× by
   design. So the reference *understates* how dominant T cells are in glioma, and the methods
   still invert the order.
2. **Cell abundance is not profile brightness.** On the marker space the T_cell profile carries
   slightly *less* mass than B_cell (0.3027 vs 0.3350 within the lymphoid triple; 5.17% vs 13.26%
   across all eight types). A rare cell type with a bright, distinctive profile is easier to
   detect than an abundant one with a dim profile — which is the most likely shape of an
   explanation here, and is **not** the same claim as the spillover hypothesis that was tested and
   rejected above.

*(The recurrence of "5.17" in two rows of this section is a coincidence: 5.17× is the T:B cell
ratio and 5.17% is T_cell's share of profile mass across all eight types. They are unrelated
quantities and neither is evidence for the other.)*

### Limits, as declared in advance

- EpiDISH's reference is a **blood** reference applied to **brain tumour** tissue. It is used here
  only for the *lymphoid sub-composition*, where the CpGs are lineage markers, not for absolute
  scale.
- `Macrophage_Microglia` is **not measurable by this route** and is not compared: microglia are
  brain-resident and appear in no blood reference, and substituting `Mono` would invent a
  correspondence the reference does not have.
- This is **LGG**. The anomaly was first seen in GBM. The pre-specification required checking
  whether it exists in LGG at all before testing it — it does, in 12 of 12 methods, which is what
  makes this a replication rather than a single observation.

Artefacts: `results/methylation_celltypes_lgg.json`, `results/methylation_celltypes_lgg.csv`.

---

## What this arm cannot do, as declared in advance

Leukocyte fraction is **one aggregate number**. It cannot verify a per-type magnitude. (Methylation
can, for the lymphoid types only — that is section 6, and it is a different measurement from the
leukocyte fraction discussed here.) The T-cell
decomposition above is therefore evidence about **estimates**, not about per-type truth: it shows
the methods' tissue T-cell values are two orders of magnitude below their pseudobulk values, which
is enough to withdraw the over-call claim, but it does not establish what the true T-cell fraction
is.

Artefacts: `results/immune_arm.json`, `results/estimates_full.csv`.
