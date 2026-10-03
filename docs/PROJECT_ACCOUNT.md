# The Anatomy Test — a complete account of the project

What was asked, what was built, what was measured, what went wrong, and what the study
concludes. Written to be read start to finish by someone who has not seen the repository.

Every figure in this document is read from a frozen artefact. The archive is
`results_archive/2026-09-23T1221`, 242 files with a recorded hash each.

---

## 1 · The problem

Reference-based deconvolution estimates what fraction of a bulk RNA sample came from each cell
type. It is used because the direct measurement — single-cell sequencing — is not available in
most settings. The estimate is only as good as the method, and the methods disagree.

The field selects between them by benchmarking against known composition: synthetic mixtures
assembled from a reference, or flow cytometry on the same specimen. Both require exactly the
measurement deconvolution is used in place of. A laboratory that cannot afford single-cell
sequencing also cannot validate the method it is using instead.

That circularity is the problem this study addresses.

---

## 2 · The proposed test, and why glioblastoma

Glioblastoma is not homogeneous. A neuropathologist names its regions on an H&E slide, and
some compositional facts about those regions are not in dispute — microvascular proliferation
is *defined* by proliferating endothelium; the leading edge is infiltrated brain and therefore
retains oligodendrocytes that dense tumour has displaced.

Those facts are independent of any deconvolution algorithm. If a method's estimates reproduce
them, that is evidence about the method — obtained from a slide and a pathologist's label,
both of which a low-resource setting already has.

**The registered question**, verbatim from the protocol:

> Can a tumor's own anatomy stand in for ground truth when choosing a cell-type deconvolution
> method — and does a method that gets the anatomy right also get the biology right?

Two clauses, and they are different questions. The first is about *selection*. The second is
about *implication*. A study can answer one and not the other, and this one does.

---

## 3 · The instrument: the Anatomic Concordance Score

### 3.1 The constraint file

Seven constraints, fixed and hashed before any method was scored
(`ivygap/anatomic/constraints.py`, SHA-256 `2d1fb47c…`):

| | type | cell type | claim | why it is safe to assert |
|---|---|---|---|---|
| **C1** | pairwise | Tumor | CT > LE | Near-definitional. The leading edge is brain infiltrated at low tumour density; cellular tumour is dense tumour. |
| **C2** | pairwise | Oligodendrocyte | LE > CT | The leading edge retains white matter; cellular tumour has displaced it. |
| **C3** | pairwise | Endothelial | MVP > CT | Microvascular proliferation is *defined* by florid endothelial proliferation. |
| **C4** | maximum | Endothelial | MVP exceeds all four others | Stricter form of C3. |
| **C5** | pairwise | Macrophage_Microglia | PAN > LE | The hypoxic peri-necrotic niche recruits myeloid cells. |
| **C6** | pairwise | Macrophage_Microglia | MVP > CT | The perivascular niche is a documented myeloid reservoir. |
| **C7** | monotone | Tumor | LE < IT < CT | A three-step chain, far harder to satisfy by chance, scored as one unit. |

**One exclusion was registered, and it matters to how this story ends.** T cells were
*deliberately left out* of the constraint set, with this reason recorded in advance:

> *"This pipeline's own synthetic benchmark puts T cells at the detection floor, so a T-cell
> constraint would score noise."*

The study knew, before scoring anything, that T-cell estimates were not trustworthy enough to
build a constraint on. It later found — through a completely different instrument — that this
was correct in a way nobody anticipated. See §7.

### 3.2 Scoring

ACS is the fraction of constraint/tumour pairs a method satisfies. Two design decisions carry
the result:

**The null is a within-tumour permutation, not a coin flip.** Cell fractions are compositional
and correlated, so 50% is the wrong null. Structure labels are permuted *within* each tumour
across 10,000 draws, preserving both the tumour's composition and the structure sizes.

**Samples are nested within tumours**, so scoring collapses to (tumour, structure) before
aggregation. Ignoring that would inflate the effective sample size.

### 3.3 The cohort

Ivy GAP holds 270 region-labelled RNA-seq samples. **Only 122 of them are usable here**, and
the reason is the single most important thing the data loader does.

| study | samples | how structure was assigned | used? |
|---|---|---|---|
| Anatomic structures | **122** (10 tumours) | H&E histology, by a neuropathologist, with no reference to RNA | **yes** |
| Cancer-stem-cell clusters | 148 (34 tumours) | in-situ hybridisation — i.e. *using expression* | **no** |

Scoring the second group would be circular: their structure labels were derived from the same
kind of signal the deconvolution reads. They are deconvolved and never scored.

Of the 10 anatomic tumours, **9 contribute at least one evaluable constraint pair**.

---

## 4 · What was tested, and against what

### 4.1 The panel — 15 estimators, and they are not all the same kind of thing

| tier | n | which |
|---|---|---|
| Published package, vendor's own code | 6 | MuSiC, SCDC, SCDC ENSEMBLE, EPIC, Bisque, quanTIseq |
| Published algorithm, reimplemented here | 3 | CIBERSORTx (B- and S-mode), nu-SVR |
| Published package, Python reimplementation used | 2 | BayesPrism, DWLS |
| Classical regression baseline, not a published tool | 4 | NNLS, elastic net, Bayesian, hierarchical Bayesian |

All fifteen decompose **bulk RNA**. None reads another modality — which is the point, because
every ground truth deliberately does.

### 4.2 Two negative controls

One estimator receives **random cell proportions**. The other receives the signature matrix
with its **gene labels shuffled**. If a scorer rates real and broken inputs alike it is
measuring nothing, and that cannot be discovered by looking at real inputs.

### 4.3 Three independent ground truths, sharing nothing with the anatomic arm

| truth | instrument | what it resolves | n |
|---|---|---|---|
| ABSOLUTE purity | DNA copy number | tumour content | 147 GBM / 496 LGG |
| DNA methylation, EpiDISH RPC | HM450 arrays | T / NK / B sub-composition | 155 GBM / 530 LGG |
| Synthetic pseudobulk | mixtures from the atlas | full composition | 14 methods |

The third shares its atlas with the anatomic arm. That turns out to matter enormously.

---

## 5 · Result 1 — anatomy detects

Every real method clears its own permutation null at **p < 1e-4**. Both controls fail theirs.

| | ACS |
|---|---|
| best method (`music`) | 0.9692 |
| lowest fully evaluable method (`bisque`) | **0.7077** |
| best control (`control_random`) | **0.4000** |
| shuffled-signature control | 0.3692 |

**Margin: +0.31 on a 0–1 scale.** The constraint set discriminates. Whatever else follows,
this part works.

*One methodological point that changes the number.* quanTIseq scores 0.600 and looks like the
worst method. It is not: it can only be scored on **15 of 57** constraint/tumour pairs, because
it ignores the supplied reference and models only four immune types. On that reduced set,
chance alone reaches 0.502, and its result is indistinguishable from chance (p = 0.31). It
**could not be evaluated**; it did not perform badly. Reporting it as the floor understated the
study's own margin by a third.

---

## 6 · Result 2 — anatomy does not rank

The registered criterion: Spearman ρ ≥ 0.60 between the ACS ordering and a ground-truth
ordering, **and** a bootstrap confidence interval excluding zero.

| yardstick | shares with the ACS arm | ρ | p | meets bar |
|---|---|---|---|---|
| Synthetic pseudobulk | **the GBmap atlas** | **+0.6372** | 0.0143 | **yes** |
| DNA tumour purity | nothing | **+0.0810** | 0.80 | no |
| …excluding degenerate methods | nothing | −0.1255 | 0.75 | no |

**The registered outcome passed — against the one yardstick built from the same atlas the
anatomic arm deconvolves against.** Every yardstick sharing nothing with ACS fails. That is a
pattern, not one bad run.

**This result is reported as INCONCLUSIVE, not as a refutation.** Twelve comparable methods
cannot separate "no relationship" from "a relationship this test is too small to see," and the
confidence interval [−0.63, +0.69] is wide enough to contain the other arm's value. Calling it
a refutation would overstate it.

---

## 7 · Result 3 — the finding the paper is built on

The study's second clause — *does a method that gets the anatomy right also get the biology
right?* — was answered by a separate pre-registered prediction.

**The prediction, and its timestamp.** On **2026-09-17 at 23:37:27**, before the data existed
in analysable form, this was committed to `prespecified/immune_failure_factors.md`:

> **PREDICTION: DNA methylation will show T cells > B cells.**
>
> *What would falsify it:* methylation showing B ≥ T. That would mean the methods may be right
> and my "anomaly" was an artefact. **That outcome would be reported and the anomaly
> withdrawn.**

The LGG methylation measurement was produced at **2026-09-18 01:12:23** — ninety-five minutes
later. The prediction genuinely precedes the data and names its own falsifier.

**The measurement.** DNA methylation places T above B in **94.8%** of glioblastomas
(Wilcoxon p = 1.1 × 10⁻²⁴) and **93.2%** of lower-grade gliomas (p = 1.1 × 10⁻⁷⁶).

**The result.**

| | GBM *(discovery)* | LGG *(registered replication)* |
|---|---|---|
| methods reproducing T > B | **0 of 12** | **0 of 12** |
| published tools and algorithms only | **0 of 9** | **0 of 9** |
| classical baselines only | 0 of 3 | 0 of 3 |

The two cohorts play different roles. The anomaly was *observed* in glioblastoma, where no
per-cell-type truth existed to check it. The prediction was then registered, and confirmed in
an independent cohort of a different tumour type at 3.4× the sample size.

**Two failure modes, which must not be merged.**

| mode | what happens | who |
|---|---|---|
| **ABSENCE** | returns *exactly* 0.0000 for T, NK and B in most samples | `music`, `nnls`, `bayesprism`, `elastic_net` — 4 of 12 in each cohort |
| **MISASSIGNMENT** | reports lymphocytes, then orders them wrongly | 6 of 8 (GBM), 8 of 8 (LGG) of the rest |

A method reporting no lymphocytes has not inverted the compartment; it has declined to
estimate it. `music` and `nnls` return zero T in **55 of 56** GBM samples.

**And this is where the registered T-cell exclusion comes back.** The constraint file left T
cells out because the pipeline's own synthetic benchmark put them at the detection floor. An
entirely separate instrument — DNA methylation on a different molecule — later showed that
T-cell estimates are not merely noisy but *systematically inverted*. The exclusion was
conservative and correct, and the study found out why by a route it had not planned.

---

## 8 · Why B above T? — five mechanisms excluded, two suspects, no answer

| # | proposed explanation | outcome |
|---|---|---|
| 1 | The signature cannot separate T from B | **rejected** — NNLS recovers a planted T:B = 2.333 exactly, surviving 100% noise and deletion of a whole cell type |
| 2 | B absorbs macrophage/microglia spillover | **rejected** — 5 of 10 methods, sign-test p = 0.62 |
| 3 | High tumour purity explains the absence mode | **rejected** — GBM 2/8, LGG 1/8; several run the other way |
| 4 | Per-sample model fit is a usable trust signal | **rejected** — 2 LGG methods significant in the *wrong* direction |
| 5 | Misfit in other markers propagates to lymphoid coefficients | **rejected** — refit on lymphoid-only markers still gives T > B in 0.0%; mass moves to NK |

**The result is not a numerical failure.** Over the 634 genes actually solved on, the full
signature has condition number **30.3** and the lymphoid block only **3.78**. The reference
holds ~5× more T than B. The signature demonstrably separates them. And every method still
inverts it on real tissue — which localises the problem to the gap between reference and real
bulk, not to the algebra.

**Two structural properties of the reference are suspects, neither proven:**

- The `NK_cell` column carries the *defining* T-cell markers. Five of nine canonical pan-T
  markers are strongest there, not in `T_cell`: CD3D 2.29×, CD3E 2.20×, CD3G 4.30×, CD2 1.42×,
  LCK 5.83×. CD3 is the T-cell marker and NK cells are CD3-negative by definition, so this
  atlas's NK population very likely contains CD3⁺ cells. **This is exactly what mechanism 5
  observed** when the mass moved to NK.
- `B_cell`'s profile is closest to `Macrophage_Microglia` (r = 0.497) — the population that
  dominates the glioma immune compartment.

Neither is sufficient. The CD3 story predicts T→NK confusion, but among the ten methods with a
computable ratio, NK/T exceeds truth in only 5 (GBM) and 7 (LGG), while **B/T exceeds truth in
all ten, in both cohorts**.

**The mechanism is unknown, and the study says so.** The next experiment is named: re-annotate
the NK population with CD3⁺ cells removed and re-run.

---

## 9 · What else was measured

| result | finding |
|---|---|
| Tumour-content recovery | median **33.1%** (GBM) / **23.3%** (LGG) among comparable methods — most of the true variation is not recovered |
| The mixing model's premise | median **63.8%** (GBM) / **76.4%** (LGG) of marker-space variance unexplained, bounded by a shuffled-label floor (−0.008) and a top-8 SVD ceiling (0.980) |
| Purity as a failure factor | the one biological factor that **replicates**: β −0.136 vs −0.134, 12/12 methods, Holm p = 0.0024 |
| Everything else | ploidy, genome doubling, subclonal fraction, IDH1 — null in both cohorts, reported at equal prominence |

**Extension E2 — which cell types can be believed without ground truth (post-registration,
exploratory; `docs/EXTENSION_IDENTIFIABILITY.md`).**
- **Design.** One genuine package (DESeq2 `unmix`) is re-fitted under seven loss scales, and each
  compartment's stability across those fits is set against its agreement with DNA truths.
- **Result.** Stability ranks compartments the way DNA does: rho **0.886** (exact p 0.017), and
  **0.943** with truth and estimate on matched denominators.
  - Tumour and total leukocyte content are stable and accurate.
  - The lymphoid compartment is the least stable, and its T/B/NK split tracks nothing against
    EpiDISH. Against GIMiCC (2026-10-03), one package's (`unmix`) T and lymphoid totals track; the
    panel's median does not.
  - **The immune compartment is identified as a whole, not its lymphocyte types.**
- **Disclosed and checked:** three flaws in the registered design, each under a rule written before
  it was computed.
- **Agreement between methods is weaker.** Tested here, DECEPTICON's agreement rule tracks accuracy only
  weakly (mean rho 0.28). Its top pairs are algorithmic siblings in 9 of 10 selections.
- **Limit.** Stability flags what not to trust; it certifies nothing.

---

## 10 · What went wrong, and how it was caught

Twenty-one defects are logged with dates and evidence. Four are worth knowing because each
changed what the paper says:

**D16 — the reference was log-transformed data treated as linear.** Deconvolution solves a
linear system; the atlas layer being read had already been log-transformed. Rebuilding from
genuine counts moved the top score from 1.000 to 0.969 and one method by 0.26. *No conclusion
changed* — but the published numbers had been indefensible, so they were replaced.

**D21 — the orthogonal yardstick's correlation was published with its sign inverted.** One
line of code declared every yardstick an error metric; one of them was a correlation, so its
truth vector was negated. Found because two artefacts disagreed about the same number, and
settled by recomputing from the per-method values. Conclusion unaffected; a published sign was
wrong.

**A verification script that could not fail.** A checker was written to confirm every statistic
appears in a result file. It harvested 755,490 numbers and asked whether each quoted value was
among them — at which density essentially everything matches. It reported a clean bill of
health on an invented figure *and* on a flipped sign. Rewritten to pin each statistic to the
field that produces it, it immediately found eight stale values.

**A title claim that was false.** A draft read *"the best-scoring method inverts the lymphoid
compartment."* The best-scoring method is `music`, which returns exactly zero T, NK and B in 55
of 56 GBM samples — **absence, not inversion**, and the wrong failure mode for the method
named.

The pattern across all four: **a check that only ever passes is indistinguishable from no
check.** Every verification in this project now ships with a deliberate failure.

---

## 11 · How the numbers were verified

Four layers, each independent of the one below:

1. **335 automated tests** over the analysis code.
2. **Document checkers.** `check_doc_numbers.py` verifies every ACS in every document against
   the *current* run, not merely against some run. `check_doc_statistics.py` pins each headline
   statistic to the artefact field that produces it.
3. **Independent recomputation.** `independent_verification.py` imports no project code and
   recomputes four headline quantities with plain scipy — including reading the raw methylation
   matrix and solving it with NNLS instead of EpiDISH RPC. All four agree; the largest
   disagreement is 0.000049.
4. **Statistical re-derivation.** Every correlation, p-value and confidence interval recomputed
   from the underlying data. p-values rest on an asymptotic approximation questionable at
   n = 9–14, so each was checked against a permutation p — **exact by full enumeration at
   n = 9 (362,880 permutations)**. No difference exceeded 0.004. Both bootstrap intervals
   reproduce bit-for-bit from the recorded seed.

All 51 references are verified: read from a PDF on disk or resolved against Crossref. One,
[10] (the GBmap citation), carried another paper's title and DOI until 2026-10-01 and is
corrected (OPEN_DEFECTS D25). Results are frozen in a hash-verified archive.

---

## 12 · What the study concludes

**Clause 1 — can anatomy substitute for ground truth in method selection?**
Partly. It *detects*: the margin over two negative controls is +0.31 and both controls fail
their own nulls. It does not *rank*: the ordering does not predict accuracy against an
orthogonal instrument, and that result is underpowered at n = 12, so it is reported as
inconclusive.

**Clause 2 — does getting the anatomy right imply getting the biology right?**
No. The method ranked first by anatomic concordance returns no lymphoid compartment at all in
98% of samples, and **0 of 12 methods reproduce a pre-registered ordering that an orthogonal
molecular instrument resolves at p = 10⁻²⁴**.

**Therefore:** concordance with established biology is not evidence that a composition estimate
is quantitatively correct. That informal check is what the field currently relies on, and this
study measures how far it gets you — which is far enough to catch a broken method, and not far
enough to choose a working one.

---

## 13 · What is bounded, and what is open

> ⚠ SUPERSEDED (OPEN_DEFECTS D22, 2026-09-30): the genuine-package re-measurements quoted below ran on GBmap's log `X` layer while the leaderboard ran on `raw/X`, so these figures are cross-matrix and are withdrawn. Current values: MANUSCRIPT §4.7 and `results/*_remeasured.json` (DWLS: no estimate for 73 of 122 samples on raw/X; Bisque and EPIC reproduce their leaderboard rows).

**Bounded and stated:** one tumour type; twelve comparable methods as a power ceiling; ReCIDE
absent, published during this study; no deep-learning methods, excluded because they train on
simulated mixtures and would compound a limitation the study already discloses; two genuine R
packages that fell back to reimplementations on this hardware, with the difference measured
(DWLS scores 0.7846 as the package against 0.7231 as the reimplementation).

**Open:** the mechanism. Five candidates excluded, two suspects named, one experiment
specified.
