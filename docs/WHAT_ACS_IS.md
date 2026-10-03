# What the Anatomic Concordance Score (ACS) is -- and what it is not

*Written 2026-10-02. Every number below is read from an artefact named beside it.*

## In one sentence

**ACS is the fraction of anatomic predictions -- written down before any method was run -- that a
deconvolution method's estimates agree with.** It asks "do this method's cell-type estimates rise
and fall between the regions of a tumour the way a neuropathologist would expect?"

## What it is NOT

**ACS is not an accuracy measure of tumour-microenvironment composition.** It never compares an
estimate with a true cell fraction, because in this cohort there is no true cell fraction to
compare with. A method can score ACS = 1.0 and still be wrong about how many of each cell type
there are. Testing whether ACS can stand in for accuracy is the question this whole study asks,
and the answer it found is **no**: ACS separates broken methods from working ones, but it does not
pick out the accurate one.

## The ingredients

**1. Tissue with anatomy (Ivy GAP).** 122 laser-microdissected samples from 10 glioblastomas, each
cut from one of five anatomic structures a pathologist outlined on the slide:

| structure | what it is |
|---|---|
| LE -- leading edge | the outermost rim: mostly normal brain with a few infiltrating tumour cells |
| IT -- infiltrating tumour | the transition zone |
| CT -- cellular tumour | the dense malignant core |
| MVP -- microvascular proliferation | regions of florid new blood-vessel growth |
| PAN -- pseudopalisading cells around necrosis | tumour cells lining dead tissue, hypoxic |

**2. Seven predictions (the constraints), fixed and hashed before any method was scored**
(`ivygap/anatomic/constraints.py`, registered at OSF `dm2t8`):

| id | prediction | why a pathologist expects it | weight |
|---|---|---|---|
| C1 | Tumor: CT > LE | the core is dense tumour, the edge is mostly brain | 1 |
| C2 | Oligodendrocyte: LE > CT | the edge keeps normal white matter; the core has displaced it | 1 |
| C3 | Endothelial: MVP > CT | MVP is *defined* by blood-vessel proliferation | 1 |
| C4 | Endothelial: MVP is the highest of all five | stricter form of C3 | 1 |
| C5 | Macrophage/microglia: PAN > LE | the hypoxic niche recruits myeloid cells | 1 |
| C6 | Macrophage/microglia: MVP > CT | the perivascular niche is a myeloid reservoir | 1 |
| C7 | Tumor: LE < IT < CT | a three-step chain, much harder to satisfy by luck | 2 |

Only four of the eight cell types appear. **T, NK and B cells and astrocytes carry no constraint**
(T cells were declared out of scope before any output was examined, because they sit at the
detection floor). So ACS cannot see an error in the lymphoid compartment at all.

**3. The scoring, step by step** (`ivygap/anatomic/acs.py`):

1. Run a deconvolution method on all 122 samples -> one composition per sample.
2. **Collapse to one composition per (tumour, structure)** by averaging that tumour's blocks of
   that structure. A tumour with three CT blocks contributes one CT value, not three.
3. For every constraint and every tumour that sampled the structures it needs: score **1 if the
   predicted ordering holds, 0 if not.** Exact ties score 0, so a method that gives every region
   the same composition cannot collect half credit.
4. **ACS = weighted share of those (constraint, tumour) pairs that score 1.** 57 pairs are evaluable
   across 9 tumours (one of the 10 lacks the structure pairs).
5. **Is that more than luck?** Shuffle the structure labels *within each tumour* 10,000 times and
   rescore. This keeps each tumour's compositions and block counts and breaks only the link to
   anatomy. Chance comes out near **0.37**, not 0.5: ties score 0, and C4 and C7 are hard to hit
   by accident (median null mean over methods, `results/anatomic/acs_leaderboard.csv`).
6. A 95% confidence interval by resampling tumours.

## A worked example from this study

SVR (this project's CIBERSORT implementation) in tumour 703393, which sampled all five structures
three times each (`results/estimates/ivygap_svr.csv`; per-structure means):

| | LE | IT | CT | MVP | PAN |
|---|---|---|---|---|---|
| Tumor | 0.247 | 0.670 | 0.702 | 0.328 | 0.880 |
| Oligodendrocyte | 0.242 | 0.086 | 0.030 | 0.001 | 0.000 |
| Endothelial | 0.021 | 0.023 | 0.044 | 0.471 | 0.002 |
| Macrophage/microglia | 0.022 | 0.020 | 0.130 | 0.161 | 0.057 |

- C1: CT 0.702 > LE 0.247 -- yes
- C2: LE 0.242 > CT 0.030 -- yes
- C3: MVP 0.471 > CT 0.044 -- yes
- C4: MVP 0.471 is the maximum -- yes
- C5: PAN 0.057 > LE 0.022 -- yes
- C6: MVP 0.161 > CT 0.130 -- yes
- C7: LE 0.247 < IT 0.670 < CT 0.702 -- yes

All seven hold, so this tumour's ACS for SVR is 8/8 = **1.0**. The same SVR, run on TCGA
glioblastoma, places B cells above T cells in **90%** of samples, where DNA methylation says T
cells outnumber B cells (`results/lymphoid_ordering.json`). It is perfect on anatomy and wrong on
the lymphoid compartment, which ACS cannot see.

## "A tumour's ACS"

The registered ACS belongs to a **method**, pooled over tumours. Per-tumour ACS -- the same score
inside one tumour -- tells you how often a method reproduced that tumour's anatomy. On 2026-10-02
we asked which tumours methods reproduce best and why (`prespecified/tumour_acs_factors.md`,
`results/tumour_acs_factors.json`, 13 comparable methods):

- **Every tumour sits well above its own chance level.** ACS minus chance runs +0.30 to +0.57.
- **None of the factors explained which tumours score higher** (all Holm-adjusted p = 1.0, n = 9:
  INCONCLUSIVE). That covers how strongly the anatomy shows in the tumour's own RNA (the
  prediction made in advance, which was not supported: rho = -0.23), sequencing depth, alignment
  rate, age and performance score. MGMT, EGFR and subtype group differences were small, with
  groups of 2-5 tumours. Every tumour had complete primary resection.
- **The negative control found a real measurement property.** In tumours with fewer sampled
  structures, even the broken-input controls score above their chance level (control rho with the
  number of structures = -0.81, p = 0.016). Per-tumour ACS is least reliable where a tumour was
  sampled most sparsely.

## Why ACS cannot be an accuracy measure -- five reasons, each measured here

1. **It never sees a true composition.** It scores orderings against expectation, not values
   against truth.
2. **It scores contrasts, not levels.** Bisque's TCGA cohort mean sits 0.087 (L1) from its
   single-cell reference's composition (`results/bisque_anchoring.json`). On planted truth it
   returns the reference prior exactly while ranking samples correctly. A method that gets the
   level entirely wrong but the differences between regions right scores well
   (`OPEN_DEFECTS.md` D23).
3. **It ignores half the roster.** A lymphoid error is invisible to it. FARDEEP has the highest
   ACS of any method in the study (0.9846) and places B above T in 100% of TCGA GBM samples on the
   frozen signature (`results/extension/`).
4. **It cannot tell a fraction from a score.** Any transformation that preserves each cell type's
   order across regions leaves ACS unchanged, so a marker score scores like a true proportion.
5. **Measured directly, it does not track accuracy.** The ACS ranking against the ranking by
   correlation with DNA-measured tumour purity gives Spearman rho = 0.081 (n = 12 methods; the
   registered bar was 0.60; `results/yardstick_agreement.json`). A graded AUC version gives -0.06
   (`results/anatomic_auc.json`). A protein-level test (imaging mass cytometry) saw ACS rank two
   marker sets in the reverse order (`results/imc_anchored_test.json`). With 12 methods the test is
   underpowered -- a true 0.7 would be detected 65% of the time (`results/agreement_power.json`) --
   so the honest wording is "not shown to predict accuracy", not "proven unrelated".

## What ACS IS good for

**Detecting broken methods.** The two broken-input controls (random proportions; a signature with
its gene labels shuffled) score 0.4000 and 0.3692, inside their own chance range. Every real method
scores 0.7077 or more and beats chance at p < 1e-4. A method that cannot reproduce basic anatomy is
broken, and ACS reliably catches it **without any ground truth.** That is half of the paper's
title -- *detects* -- and it is a real, usable result.

## Where this study's accuracy numbers come from instead

| truth | what it measures | file |
|---|---|---|
| ABSOLUTE purity from DNA copy number (TCGA) | true tumour content per sample, sharing nothing with RNA | `results/absolute_purity_yardstick*.json` |
| EpiDISH on DNA methylation (TCGA) | the T / NK / B ordering, at cohort level | `results/methylation_celltypes*.json` |
| synthetic mixtures with planted composition | exact truth, but built from the same atlas the methods use | `reference_frozen/tcga_benchmark/` |

## What this means for replacing other measurements

Bulk-RNA deconvolution estimates **which cells are in a sample and in what proportions**. It
competes with single-cell RNA sequencing, flow cytometry, imaging mass cytometry and histology,
which measure composition directly. It does **not** compete with CT, MRI or ultrasound, which image
where a tumour is, how large it is and what it displaces: no composition estimate carries that
information. On the evidence here it cannot yet replace single-cell measurement for the immune
compartment of glioma (no method robustly recovered the T-above-B ordering). And choosing among
methods still needs a ground truth. Anatomy can screen out broken methods; it cannot certify an
accurate one.
