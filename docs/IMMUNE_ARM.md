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
signal in the wrong column. It is reported as an observation requiring its own test, not as an
established mechanism — it was not pre-specified and no direction was fixed for it.

---

## What this arm cannot do, as declared in advance

Leukocyte fraction is **one aggregate number**. It cannot verify a per-type magnitude. The T-cell
decomposition above is therefore evidence about **estimates**, not about per-type truth: it shows
the methods' tissue T-cell values are two orders of magnitude below their pseudobulk values, which
is enough to withdraw the over-call claim, but it does not establish what the true T-cell fraction
is.

Artefacts: `results/immune_arm.json`, `results/estimates_full.csv`.
