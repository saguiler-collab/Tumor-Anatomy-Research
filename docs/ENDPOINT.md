# The endpoint: what this project found, and what the paper says

**Status 2026-09-15. Everything below is measured and has an artefact behind it.** This is the
synthesis document — the claim, the evidence for it, the evidence against the original
hypothesis, and what is deliberately left undone.

---

## The one-sentence finding

> **Anatomic concordance is a valid detector of whether a deconvolution method is working, and
> an invalid ranker of which method is working best.**

The tumour's anatomy does carry real, independently verifiable information about cell-type
composition, and methods that recover it are doing something genuine. But the *ordering* it
produces among working methods is not stable enough, and not externally valid enough, to choose
a method with.

**That is a negative result on the registered hypothesis, and it is the more useful finding.**
"Use biological plausibility instead of ground truth" is an increasingly popular idea across
deconvolution, spatial transcriptomics and single-cell annotation. This project tested it
properly — pre-registered, with negative controls, against three independent external yardsticks
— and found where it breaks. A field-wide caution supported by measurement is worth more than
one more benchmark that agrees with the others.

---

## Part 1 — What works. The detector claim, and it is strong.

### The constraints are real biology, verified in tissue that this project did not deconvolve

The seven pre-registered anatomic constraints were checked directly against **measured**
single-cell composition, with no deconvolution method involved:

| source | what it is | result |
|---|---|---|
| **Albiach et al. 2023** | 135,482 annotated cells with anatomic zones | **3 of 4 testable constraints satisfied.** C1 (Tumor CT > LE) **13.34-fold**, p = 0.00005. C2 (Oligodendrocyte LE > CT) **3.53-fold**, p = 0.00075. C7 by ordering. C5 violated and marked WEAK. |
| **Darmanis et al. 2017** | 4 patients, FACS-gated neoplastic fraction, core vs periphery | **4 of 4 gates support C1**, two at p = 0.0001 with 3 of 3 patients. The only cross-patient constraint test in the project. |
| **Ivy GAP ISH** | 18,778 rows, 480 genes, 899 sub-blocks, 42 donors | Markers declared before values were read; C2 explicitly **not checkable** because the panel has no oligodendrocyte-lineage marker, and OLIG2 was refused on the grounds that in glioma it measures tumour. |

**The constraint file was frozen and hashed before any of this** (`2d1fb47c…5807a`) and was not
edited when a constraint failed. C5's violation is reported, not repaired.

### ACS separates real methods from controls decisively

| | ACS | beats permutation null? |
|---|---|---|
| 14 real methods | 0.7385 – 1.0000 | **all yes, p = 0.0001** |
| `control_random` | **0.3385** | no, p = 0.78 |
| `control_shuffled_signature` | **0.1385** | no, p = 0.999 |

This separation **holds under every reference atlas tested** — GBmap, GBmap-linear, Neftel,
Darmanis, and both GBmap assay subsets. It is the most robust result in the project.

### The registered primary outcome passed its pre-registered bar

**Spearman rho = 0.7501** between the ACS ranking and the accuracy ranking, CI [0.328, 0.957],
p = 0.0020. The bar, fixed in advance, was rho ≥ 0.60 with a CI excluding zero.

**It passed. And Part 2 is about what it actually measures.**

---

## Part 2 — What fails. The ranker claim, refuted four independent ways.

### 2a. It failed against every external ground truth tested

**Three independent anchors, none of which this project controls:**

| anchor | what it is | result |
|---|---|---|
| **Avila Cobos 2020** (*Nat Commun* 11:5650) | benchmark against known composition | **Mann-Whitney p = 0.50.** No evidence their top-tier methods score higher on ACS. **DWLS is in their top tier and ranks 14th of 14 here.** |
| **Ajaib 2023** (*Neuro-Oncology* 25:1236) | **imaging mass cytometry**, 33 antibodies, protein, 10 GBM with matched bulk | **Pre-specified prediction FAILED**, p = 0.0090 in the wrong direction. ACS ranked GBmap-derived markers (0.6667) above the GBM-specific markers IMC ranks first (0.5128). |
| **Sturm 2019** (*Bioinformatics* 35:i436) | flow cytometry + immune mixtures | Descriptive only — 3 overlapping methods. EPIC ranks 2nd here, consistent with their recommendation. |

The IMC test is the important one because it was **pre-registered**: the prediction, the arms,
the preprocessing and the constraint subsets were committed before the run, and the outcome that
occurred was named in advance as the most damaging one.

### 2b. The ordering is substantially a property of the reference atlas, not the methods

Rebuilding the same atlas from its own counts to break a confound (D16):

| comparison | atlas | expression space | Spearman |
|---|---|---|---|
| GBmap_10x vs GBmap_SS2 | same | same | +0.8170 |
| GBmap_log vs GBmap_linear | same | different | **+0.9161** |
| GBmap_linear vs Neftel | different | same | **+0.5099** |
| GBmap_linear vs Darmanis | different | same | **+0.3655** |

Changing platform or expression space within one atlas costs little. **Changing the atlas leaves
about half the ordering.** And the registered primary outcome compares two rankings that
**both** use GBmap — the ACS arm deconvolves against it, the accuracy arm builds its mixtures
from its cells (C4). Part of rho = 0.7501 is two GBmap-based rankings agreeing about GBmap.

### 2c. It cannot resolve the top of its own leaderboard

| rank | method | ACS |
|---|---|---|
| 1 | MuSiC | 1.0000 |
| 2= | NNLS, SVR, elastic net, EPIC | **0.9846 — a four-way tie** |
| 6= | CIBERSORTx, CIBERSORTx S-mode, Bisque | **0.9692 — a three-way tie** |

**MuSiC's win is one satisfied pair out of 57.** The question "why did MuSiC win?" is
unanswerable at this resolution, and saying so is more honest than constructing a story around a
single pair.

### 2d. It is movable by choices the method never made

- **One constraint carries the ranking.** Drop C6 and the ordering changes: rho 0.678, and 13 of
  14 methods move (D8).
- **Marker-set size moves ACS by 0.31.** Same atlas, same algorithm, same bulk, same constraints,
  same tumours — only markers-per-type changed, and ACS went **0.3846 → 0.6923**, non-overlapping
  CIs, from indistinguishable-from-null to significant, **non-monotonically**. The whole
  leaderboard spans ~0.60 to 1.00. A nuisance parameter moves ACS across three-quarters of the
  range it is being used to rank fifteen methods across.

### 2e. The reference-free test — and it cuts BOTH ways, which is the useful part

Run 2026-09-15. Every other method in the panel is handed GBmap, so none of them can separate
"the methods recover the anatomy" from "GBmap recovers the anatomy". **CDSeq estimates cell types
de novo from the bulk alone**, which makes it the only instrument that can. It was blocked until
the Ivy GAP read counts were obtained (**C9**, a retraction).

Genes were chosen by **variance in the bulk**, not from the leaderboard's 657 markers — those
were selected on GBmap, and using them would reintroduce the dependence this arm exists to
remove. Every method is scored on the **same 41 constraint–tumour pairs**.

| | ACS | p |
|---|---|---|
| MuSiC, CIBERSORTx S-mode | 1.0000 | 0.0001 |
| elastic net, EPIC, NNLS, SVR | 0.9756 | 0.0001 |
| BayesPrism, SCDC | 0.9268 | 0.0001 |
| Bisque | 0.8780 | 0.0001 |
| **DWLS** | **0.7561** | 0.0002 |
| **CDSeq — atlas for labelling only** | **0.7561** | 0.0002 |
| **CDSeq — no atlas at all** | **0.6829** | 0.0009 |
| `control_random` | 0.5122 | **0.2233 — fails its null** |
| `control_shuffled_signature` | 0.1951 | 0.9977 |

**Reading one: the constraints are recoverable without any atlas.** A method that never saw
GBmap — not for estimation, not even for naming its own output — reaches **0.6829 at
p = 0.0009**, decisively above a random control that does not beat its own null. The anatomic
signal is in the bulk, not manufactured by the reference. This is the strongest support the
**detector** claim has.

**Reading two, and it is a genuine counterweight to this project's own thesis: the reference
earns its place.** CDSeq lands at the *bottom* of the real methods, tying DWLS exactly. Thirteen
of fourteen reference-using methods beat it. So the reference is not merely a source of
arbitrary ordering — it contributes real signal that bulk-only estimation does not recover.

**That matters for honesty about C4.** "The ordering is substantially a property of GBmap" is
established. "GBmap contributes nothing but its own bias" is **not**, and this measurement is
what rules it out. Both sentences belong in the paper.

**The caveat is severe and is the reason the two CDSeq rows differ.** CDSeq returns anonymous
cell types, and mapping them onto the roster changed the answer: correlating with the GBmap
profile called four of eight estimated types **Tumor**; marker enrichment called **none** of them
Tumor and three of them Astrocyte. The de novo types are not cleanly identifiable against this
eight-type roster, and a reference-free method that cannot name its own output is of limited
clinical use whatever it scores.

Artefacts: `results/cdseq_anatomic.json`, `results/cdseq_matched_comparison.json`.

---

## Part 3 — The clinical findings, which stand on their own

These do not depend on the ranking question and are directly reportable.

- **Every method under-calls tumour by 0.33–0.51 at high purity** — exactly where real tissue
  sits. Error grows monotonically with tumour content.
- **The missing mass becomes T cells**, over-called **4.6–7.9×**. This is the largest error mode
  in the panel.
- **ACS is structurally blind to it.** `T_cell` carries no constraint, by pre-registered design,
  because no anatomic fact about T-cell distribution was defensible enough to register.
- **1 of 15 methods reports uncertainty at all**, and its intervals cover **9.25%** against a
  nominal 95%.
- **Split-conformal intervals on tumour content span ±0.49 (Bisque) to ±0.67** — the latter is
  the width of the whole simplex, i.e. no information.
- **Only cross-donor variance weighting buys accuracy** (5 of 6 outcomes). Regularisation, batch
  correction and donor hierarchy buy nothing on this cohort.

**Clinical bottom line:** no method in this panel is fit to report a patient's immune content,
and the evaluation framework most likely to be used to choose between them cannot see the error.

---

## Part 4 — The methodological contribution

Seventeen defects were found, measured and recorded, several of which invalidate work that
looked finished. The ones a reader should see:

| | what it was | why it matters beyond this project |
|---|---|---|
| **D12** | the central mRNA→cell conversion was **the identity**, applied nowhere | proven twice, the second time without reading any code: applying the transform and skipping it gave bit-identical output |
| **D16** | the reference was built from **log-transformed** values treated as linear | the accuracy arm is internally consistent and survives; the **anatomic arm is not**. One defect explains why the benchmark looked healthy while tissue accuracy was poor |
| **D14** | the ordering does not survive a change of atlas | the first attribution ("it is the atlas") was confounded with expression space; corrected after measurement, magnitude withdrawn |
| **D17** | variant builds silently overwrote the primary reference's provenance record | a figure script read that exact path |
| **C9** | claimed Ivy GAP counts were not public | **retracted the same day**: GEO says "raw data not provided", the Allen portal serves per-sample RSEM counts. Wrong because one authoritative source was mistaken for an exhaustive search |

**Three claims in this project were published and then withdrawn after measurement** (D1 three
times, D14's magnitude, C9). Each retraction is kept in place rather than edited away. That
record is part of the result.

---

## What the paper claims, precisely

1. Anatomic constraints in glioblastoma are **real and independently verifiable** — confirmed in
   two external single-cell datasets that were never deconvolved.
2. Deconvolution methods **satisfy them far above chance**, and negative controls do not. ACS is
   a working **detector**.
3. **ACS is not a valid method selector.** It fails against published ground-truth benchmarks
   (p = 0.50), fails a pre-registered test against protein ground truth (p = 0.0090, wrong
   direction), cannot resolve the top of its own leaderboard (a four-way tie), rests largely on
   one constraint, and moves by 0.31 on marker-set size alone.
4. The apparent success of the registered outcome (rho = 0.7501) is **substantially a shared
   dependence on one reference atlas**, and every analysis that breaks that dependence weakens
   it.
5. Separately and robustly: **current deconvolution methods are not clinically usable for immune
   content in GBM**, and anatomy-based evaluation is structurally blind to that failure.

**The honest headline is not "anatomy works" and not "anatomy fails". It is: anatomy tells you
whether a method is working, not which one is best — and the difference between those two
questions is where the field has been assuming rather than measuring.**

---

## What is deliberately NOT done, and why that is fine

- **The full linear re-run (B1).** D16 means the anatomic arm should ideally be recomputed on a
  linear reference. It is declared as a limitation with its direction measured (+0.9161 — the
  expression space barely moves the ordering), so the paper does not depend on it.
- **CDSeq on the anatomic arm.** Unblocked and feasibility-measured (~30 min), not run. It would
  strengthen Part 2b; it is not required for any claim made.
- **A second tissue.** Nine evaluable tumours is the statistical ceiling and is stated as one.
- **Retuning ACS so it performs better.** Forbidden by the project's own invariants, and the
  reason the negative result is credible.
