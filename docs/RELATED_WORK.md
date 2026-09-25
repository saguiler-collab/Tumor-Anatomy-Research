# Related work — where this study agrees, disagrees, and what is new

Five papers, read 2026-09-10, compared against the confirmatory run
(`results_archive/2026-09-10T2039`). Full citations in `sources/SOURCES.md`.

The short version: **every prior benchmark ranks methods against known ground truth.
This study asks whether methods can be ranked without any.** That is the gap, and one of
the reviews below names it explicitly.

---

## 1 · The papers

| | what it is |
|---|---|
| **Avila Cobos et al. 2020**, *Nat Commun* 11:5650 | Benchmarks deconvolution pipelines on simulated pseudobulk. doi:10.1038/s41467-020-19015-1 |
| **Sturm et al. 2019**, *Bioinformatics* 35:i436 | Benchmarks immune-cell quantification for immuno-oncology; produced `immunedeconv`. doi:10.1093/bioinformatics/btz363 |
| **Nguyen et al. 2024**, *Nucleic Acids Res* 52:4761 | "Fourteen years of cellular deconvolution" — reviews and benchmarks **53 methods** across 283 cell types, 30 tissues, 63 individuals. doi:10.1093/nar/gkae267 |
| **Gaspard-Boulinc et al. 2025**, *Nat Rev Genet* 26:828 | Deconvolution for **spatial** transcriptomics. doi:10.1038/s41576-025-00845-y |
| **Liu, Qian & Ma 2025**, *bioRxiv* | **DNA-methylation-based** deconvolution of the brain-tumour microenvironment. doi:10.1101/2025.01.19.633794 |
| **Li et al. 2026**, *Genome Biology* 27:38 | Benchmarks **5 methods** on **18 real bulk cohorts, 5,891 samples, 9 cancer types**, scoring *reproducibility* of differentially-proportioned cell types instead of accuracy against known proportions. doi:10.1186/s13059-026-03942-1 |

*(The file `Nguyen et al., 2014.pdf` is misnamed — it is the 2024 NAR review.)*

---

## 2 · Where our results agree with theirs

**Least-squares and SVR methods rank near the top — confirmed twice over.**

Avila Cobos: *"the five best bulk deconvolution methods (OLS, nnls, RLR, FARDEEP, and
CIBERSORT)... regression methods (RLR, FARDEEP) and support vector regression
(CIBERSORT) consistently showed the smallest RMSE and highest Pearson correlation."*

Ours: `nnls` and `svr` — and `svr` **is** CIBERSORT's published nu-SVR core — sit in the
top tier at ACS 0.9846, and rank 3rd and 6th on accuracy (MAE 0.0561, 0.0610).

Avila Cobos reached that ranking with **known ground truth on simulated mixtures**. We
reach a compatible placement on **real tissue using no ground truth at all**, from
pre-registered anatomy.

**But do not overstate this — section 4b tests it formally and the overall result is
null (p = 0.5).** What holds is narrower: the three methods whose *identical published
implementation* both studies test — MuSiC, NNLS and CIBERSORT — are our top three. The
agreement is real where the comparison is clean, and the formal test fails because DWLS
is a reimplementation here and because two methods outside their top tier score as high
as ours.

**EPIC ranks highly — agrees with Sturm.**

Sturm: *"due to a robust overall performance, we recommend EPIC and quanTIseq for general
purpose deconvolution."* Ours: EPIC at ACS 0.9846, tied for second.

**Linear scale, not log — agrees with Avila Cobos.**

Avila Cobos: *"the most relevant factors affecting the deconvolution results are: (i) the
data transformation, with linear transformation outperforming the others."* This project
normalises to CPM and never log-transforms for deconvolution, because the mixing model is
additive on the linear scale. Independently arrived at, independently confirmed.

**Reference completeness matters — agrees with Avila Cobos.**

Avila Cobos: *"failure to include cell types in the reference"* degrades results. This
project measures exactly that and reports it: `reference_coverage.json` records that the
eight-type roster drops **7.0% of the atlas**, that the largest dropped population is
myeloid, and that the absorbed signal therefore concentrates on C5 and C6 rather than
spreading evenly.

---

### Li et al. 2026 independently confirm the objection this study raises against its own primary outcome

This is the closest paper to ours in motivation, published contemporaneously, and it reaches
our central methodological conclusion from a completely different direction.

Their opening premise is the one our correction **C4** makes about our own registered outcome:

> *"current benchmarking studies for deconvolution methods invariably lean on pseudobulk data
> or flow cytometry, assuming known absolute cell type proportions. In real bulk RNA-expression
> deconvolution, such precision is a mirage."*

And they measured it. On **GSE176078**, a breast-cancer cohort profiled by *both* scRNA-seq and
bulk RNA-seq on the same samples, they report that pseudobulk yields **systematically higher
sample-level correlations than real bulk**, and — the part that matters for us — that **method
rankings transfer poorly between the two**: ranking consistency is significantly lower across
pseudobulk-versus-real-bulk than within either modality alone (p < 0.001).

That is an independent, 5,891-sample demonstration of the exact pattern our two yardsticks
produce. Our synthetic-mixture arm clears the pre-registered bar at **ρ = 0.6372**; our
orthogonal DNA-purity arm fails it at **ρ = 0.0810**. We attributed the gap to the pseudobulk
yardstick sharing its atlas with the ACS arm. Li et al. show the gap is general, not an artefact
of our cohort or our reference.

**What this changes for us: nothing in the numbers, and a great deal in the framing.** A
reviewer's most natural objection to our study is *"your pseudobulk arm passed — maybe your
orthogonal arm is the broken one."* The answer is no longer only our own internal argument.

### The MuSiC convergence, which is the sharper point

Li et al. find **"MuSiC showed limited performance with real-world data, potentially due to its
sensitivity to noise."**

In our study **MuSiC ranks first of fifteen on anatomic concordance (ACS 0.9692)** — and
recovers only **18.1%** (GBM) and **17.2%** (LGG) of true tumour-content variation, and places B
cells above T cells in **100%** of GBM and **98.5%** of LGG samples.

So the method our ground-truth-free metric likes *best* is the one an independent real-data
benchmark of 5,891 samples singles out as performing *worst*. That is "anatomy detects, it does
not rank" confirmed from outside this project, on different cancers, by different authors, using
a different criterion. It belongs in the discussion.

---

## 3 · Where our results disagree, and why

### BayesPrism: they call it robust; our orthogonal truth puts it near the bottom

Li et al. conclude that **ReCIDE and BayesPrism** are the two robust methods across all three of
their scenarios. Our results do not look like that:

| | Li et al. 2026 | this study |
|---|---|---|
| BayesPrism | one of two most robust across 3 scenarios | ACS **0.8000** (9th of 15); recovery **10.3%** GBM / **7.2%** LGG; places B over T in **80.0%** / **57.4%** of samples |

**The disagreement is real and should not be smoothed over — but it is probably not a
contradiction, because the two studies measure different things.** Li et al. score the
*reproducibility* of which cell types shift between conditions, and say so explicitly: they
focus on *"the presence or absence of DP cell types"* rather than magnitudes, having found the
qualitative criterion more stable than the quantitative one. We score *quantitative* recovery of
a continuous truth and *per-cell-type ordering correctness*.

A method can reproducibly detect **which** populations move while getting **how much** they move,
and their relative ranking, wrong. Those are compatible findings about different properties, and
the practical implication is the interesting part: **BayesPrism may be the right choice for
"which cell types differ between my two groups" and the wrong choice for "what fraction of this
tumour is T cells"** — a distinction no single benchmark score expresses.

The honest caveat on our side: BayesPrism runs here as a Python reimplementation (the genuine R
package scores **0.8154**, still mid-table), and our cohort is glioma only.

---


### DWLS: they rank it best among scRNA-reference methods; we rank it last

Avila Cobos: *"DWLS performed best among the deconvolution methods that use scRNA-seq
data as input."* Ours: DWLS is **last among real methods** — ACS 0.7385, MAE 0.0763.

**RESOLVED 2026-09-10: it is a real disagreement, not an artefact.**

The confirmatory run's DWLS row is the Python reimplementation — it exceeded its 2,400 s
budget in both anatomic stages — so the comparison was initially unresolvable.
`scripts/remeasure_dwls.py` therefore ran the **genuine R package** on the same cohort
with a 4-hour budget. It completed in **3,873 s (1.08 h)**, above the pipeline budget and
well inside the new one.

| | ACS | 95% CI | pairs | null_p | genes | donors |
|---|---|---|---|---|---|---|
| **genuine `R:DWLS`** | **0.7846** | [0.657, 0.906] | 57 | 0.0001 | 657 | 88 train |
| Python reimplementation | 0.7385 | — | 57 | 0.0001 | 657 | 88 train |

**Measured 2026-09-14 on verified-equivalent inputs.** Seven conditions are checked and
recorded in the report, and the run aborts if any fails: gene space by hash, the 88 training
donors by name, the 22 held-out donors by name, no silent sample loss, identical sample IDs
in order, identical cell-type ordering, and the same normalisation and scale.

> **An earlier version of this table reported 0.7077 and is VOID.** That run used the
> 1,591-gene fallback space *and* built its reference from all 110 donors, so the method saw
> the 22 held-out donors' cells. Its conclusion — "the genuine package scores lower than the
> reimplementation, and stays last" — **reversed** once the inputs were made equivalent. See
> `docs/OPEN_DEFECTS.md` D10.

**The genuine package scores HIGHER than the reimplementation**, by 0.046, which is three
weighted constraint–tumour pairs out of 57. At 0.7846 DWLS is no longer last among the real
methods: it passes both Bayesian models (0.7692). So part of this study's disagreement with
Avila Cobos *was* an artefact of substituted software, and the part that remains is smaller
than reported.

So **the substituted-software explanation accounts for part of the disagreement, and the
comparison can now be made.** Avila Cobos et al. rank DWLS best among single-cell-reference
methods; this study ranked it last at 0.7385, and the genuine package measures 0.7846 — above
both Bayesian models, no longer last, and still far from first. The two readings below remain
the candidate explanations for the residual gap. They are not exclusive, and neither is
established:

**DWLS is built for rare cell types.** Its dampening suppresses the dominance of highly
expressed genes belonging to abundant populations — the regime where a rare type would
otherwise be swamped. This roster is the opposite: tumour cells are roughly 50-60% of
every sample, and the constraints that discriminate (C5, C6, C7) concern myeloid and
tumour populations that are not rare. Avila Cobos's simulated mixtures were not
tumour-dominated in this way.

**And the yardsticks differ.** Avila Cobos ranks by RMSE against simulated composition;
this ranks by ordinal agreement with anatomy on real microdissected tissue. A method can
be well calibrated in magnitude and still order structures wrongly, which is precisely
the distinction this study exists to examine.

The leaderboard row still carries the reimplementation, because that is what the
archived confirmatory run produced and the archive stands as recorded. The genuine
measurement is reported **alongside** it in `results/dwls_remeasured.json`, never
substituted in — and that file now carries a `COMPARABILITY: NOT COMPARABLE` header.

**Any statement about DWLS in the paper must name the software AND the gene space.**
Quoting 0.7077 beside 0.7385 without both is the error this section made. The same applies
to BayesPrism, re-measured at 0.8923 on the same 1,591-gene space against a leaderboard row
of 0.8769 on 657.

A second, independent factor: DWLS is built for **rare** cell types, damping the
dominance of highly expressed genes belonging to abundant populations. This roster is
dominated by tumour cells at roughly 50–60% of every sample. That is close to the
opposite of the regime DWLS was designed for.

### quanTIseq: Sturm recommends it; we cannot rank it

Sturm recommends quanTIseq for general-purpose deconvolution. Ours reports it at ACS
0.600 on **15 constraint–tumour pairs** rather than 57, flagged `comparable = False`, and
excluded from the ranking, the tie groups and the primary correlation.

**Different question, not a different answer.** Sturm's benchmark is immuno-oncology —
immune cells in tumours. quanTIseq ships TIL10, which models ten immune populations and
rolls everything else into "Other". Our roster is eight types including Tumor,
Endothelial, Oligodendrocyte and Astrocyte, four of which quanTIseq does not model at
all, in any sample, however well it worked.

So quanTIseq is being asked a question it was not built to answer. Reporting it with its
own denominator and refusing to rank it is the correct handling, and it is what
`comparable = False` exists for.

---

## 4 · What is original here

### What remains original AFTER Li et al. 2026 — read this before writing the introduction

Li et al. is the nearest neighbour this study has, and pretending otherwise would be the
fastest way to lose a reviewer. Both papers reject pseudobulk as the arbiter and both build a
benchmark that needs no known cell proportions. The differences are real and worth naming
precisely rather than asserting novelty in general terms.

| | Li et al. 2026 | this study |
|---|---|---|
| ground-truth-free criterion | **reproducibility** — do the same cell types shift across independent cohorts? | **anatomy** — does the estimate obey constraints a pathologist fixed in advance? |
| what it asks of a method | *consistency* | *correctness against a named prior* |
| minimum data needed | two or more cohorts sharing a disease contrast | **one** cohort with region labels |
| is it checked against orthogonal truth? | **no** — reproducibility is the endpoint | **yes** — DNA methylation per cell type, DNA copy-number purity |
| negative controls | not reported | **two**, and the margin over them licenses every claim |
| scale | 5,891 samples, 9 cancers, 5 methods | 122 anatomic + 664 truth-linked samples, 1 cancer type, 15 estimators |

Three things follow, and they are the introduction's argument:

1. **Reproducibility is not correctness, and this study supplies the missing half.** Li et al.
   establish that ReCIDE and BayesPrism reproduce their DP cell-type calls across cohorts.
   Nothing in that design can detect a method that is *reproducibly wrong* — and we exhibit
   exactly that failure: **0 of 12 methods** reproduce a T > B lymphoid ordering that DNA
   methylation measures at p = 1.1 × 10⁻²⁴ (GBM) and 1.1 × 10⁻⁷⁶ (LGG), and they fail
   *consistently*, in both cohorts. A reproducibility benchmark would score that consistency
   favourably. An orthogonal-truth benchmark scores it as the failure it is.

2. **Their design needs a contrast; ours needs a slide.** DP cell types are defined by
   condition-level change, so the method requires at least two groups and, for the strongest
   scenario, several cohorts. Anatomic concordance is computed within a single tumour set from
   histology labels. For the low-resource setting this project is motivated by — one cohort,
   no matched single-cell data, no second site — that difference is the whole point.

3. **Neither benchmark can select a method, and only one of us can show it.** Li et al. rank
   methods and stop. We rank methods, then check the ranking against independent truth, and
   report that it does not survive (ρ = +0.081, n = 12, p = 0.80 — inconclusive) and that the
   top-ranked method inverts the lymphoid compartment. **The negative result about our own
   metric is the contribution**, and it is a kind of result their design has no way to produce.

### What we should concede, plainly

- **ReCIDE is not in our panel.** Li et al.'s best performer is absent here. Say so; it is a
  stated limitation, not a hidden one.
- **One cancer type against their nine.** Our replication is GBM → LGG, both glioma. Their
  generalisation across nine cancers is broader than anything we claim.
- **Their sample size is 8× ours.** The power ceiling we disclose (12 comparable methods) is
  the one thing more data would genuinely fix, and they have more data.


### The gap, named by the field's own review

Nguyen et al. 2024 enumerate five validation strategies. Four require known ground truth —
simulation, purified samples, flow cytometry, matched scRNA-seq. Of the fifth:

> *"The fifth approach relies on **domain experts to interpret the deconvolution results**
> to indirectly assess the performance of deconvolution methods."*

That is exactly the practice this project argues is unfalsifiable: the expert looks at the
output and judges whether it is sensible, **after** seeing it. A method cannot fail such a
check, because any result can be rationalised once observed.

**This study's contribution is to convert that fifth approach into a test.** The expert
expectation is written down first, frozen, cryptographically hashed
(`2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`) and publicly
registered (<https://osf.io/dm2t8>) before any deconvolution output is inspected. The
same judgement then becomes a prediction a method can miss.

To our knowledge no prior deconvolution benchmark pre-registers its expectations.

### And a criticism of the standard approach that this design answers

The same review says of simulation-based benchmarking — the dominant strategy, and this
project's own accuracy arm:

> *"simulation is subjected to bias because simulated data is generated based on some
> assumptions which are usually identical with the assumptions made in designing the
> approach. **Presumably, any algorithm would be the best, when applied to data that was
> simulated based on the same set of assumptions.**"*

This is a serious charge against the field's default yardstick, and this study should
quote it rather than avoid it, because the design has an answer: **anatomic concordance
does not depend on the simulator's assumptions.** It is measured on real microdissected
tissue against claims derived from histology. The two arms fail in different ways, which
is why correlating them is informative — and it is why our own observation that the
pseudobulk arm cannot reward CIBERSORTx S-mode (see `OPEN_DEFECTS.md` D5) is an instance
of exactly the bias the review describes, found independently in our own data.

### Four things that are unusual even among careful benchmarks

1. **Negative controls run through the identical pipeline, calibrated over 200 draws.**
   Random compositions and a shuffled signature matrix. The protocol commits in advance to
   publishing the constraint file unchanged even if the controls score well. None of the
   five papers above runs a shuffled-signature control this way.
2. **The null is not a coin flip.** Cell fractions are compositional and samples nest
   within tumours, so significance comes from permuting structure labels *within* each
   tumour, 10,000 draws.
3. **A registered failure branch.** rho < 0.60 was pre-committed to mean that reproducing
   known biology is *not* evidence of numerical accuracy — reported with equal prominence.
4. **The yardstick itself is checked.** Ivy GAP's own ISH independently supports C3
   (ESM1, p = 0.005), C4 (ESM1, p = 0.021) and C6 (CD163, p = 0.0085), from data that
   never touched a deconvolution method. Benchmarks validate methods; this validates the
   measuring stick.

---

## 4b · Concordance test, run 2026-09-10 — and it is not a clean win

`scripts/benchmark_concordance.py`. The published claim, the method set and the statistic
were all fixed in the script before it ran, because otherwise this is trivially fakeable.

**Avila Cobos 2020.** Their top tier, in our ACS ranking of 14 comparable methods:

| their top tier | our ACS | our rank |
|---|---|---|
| MuSiC | 1.0000 | **1** |
| nnls | 0.9846 | **2** |
| CIBERSORT (= our `svr`) | 0.9846 | **2** |
| SCDC | 0.9538 | 8 |
| DWLS | 0.7385 | **14** (last) |

| evaluated, not top tier | our ACS |
|---|---|
| EPIC | 0.9846 |
| elastic net | 0.9846 |
| Bisque | 0.9231 |

**One-sided Mann-Whitney p = 0.5. No evidence the top tier ranks higher overall.**

That is the result and it should be reported as such. Two things drive it:

1. **DWLS sits last.** In the confirmatory run it ran as a Python reimplementation after
   exceeding its wall-clock budget, so this row is not the package Avila Cobos evaluated.
   It is kept in the statistic rather than dropped, because removing the one method
   expected to disagree is precisely the move this design exists to prevent. But the
   correct reading is that **this study does not test the published DWLS**, so the
   comparison is weakened at exactly the point where it most needs to be clean.
2. **EPIC and elastic net score 0.9846** — as high as our top tier — while sitting
   outside Avila Cobos's top tier.

**What survives.** Of the five overlapping methods, the three whose *identical published
implementation* both studies test — MuSiC, NNLS, CIBERSORT — are our **top three**. The
agreement is real where the comparison is clean and breaks where it is not. That is a
narrower claim than "ACS agrees with published benchmarks", and it is the one the data
supports.

**Sturm 2019** cannot be tested statistically: only EPIC and CIBERSORT overlap in a way
that permits a group comparison, and quanTIseq is excluded from the statistic in advance
because its ACS is computed over 15 constraint-tumour pairs rather than 57. EPIC, which
Sturm recommends, is at rank 2 here.

### The finding underneath the finding

**The two published benchmarks disagree with each other about EPIC.** Sturm recommends it
for general-purpose deconvolution; Avila Cobos does not place it in the top tier. So there
is no single published consensus ranking for ACS to agree with, and any paper claiming
"our ranking matches the literature" should be asked *which* literature.

That is worth stating plainly, because it also reframes what this study is for. If
established benchmarks using real ground truth disagree with one another about a method,
then the value of a ground-truth-free check is not that it settles the question — it is
that it is a *cheap, independent* axis that can be applied where no ground truth exists at
all, which is almost everywhere.

---

## 5 · What these papers say we should do next

- **Nguyen 2024** benchmarks 53 methods; this study has 15. The review's list is where to
  find methods that differ in *kind* — the 14 comparable methods here produce only 8
  distinct ACS values, so another least-squares variant adds nothing.
- **Gaspard-Boulinc 2025** covers spatial deconvolution. Relevant if the second tissue
  becomes Visium rather than Allen bulk: it would keep the GBM constraint file unchanged,
  at the cost of a 55-micron spot being a handful of cells rather than a tissue block.
- **Liu, Qian & Ma 2025** deconvolve the brain-tumour microenvironment from **DNA
  methylation** — an orthogonal modality that shares no failure mode with RNA-based
  deconvolution. That makes it a candidate third yardstick alongside ABSOLUTE purity, and
  it is specifically about brain tumours.
- **Avila Cobos 2020** and **Sturm 2019** both publish per-method rankings. Correlating
  our ACS ranking against theirs would ask the study's own question with *someone else's*
  ground truth on *someone else's* tissue. It needs no new data. See
  `docs/EXTERNAL_CHECKLIST.md` section F2 — and note it must be pre-specified which
  methods and which published metric, before looking, or it is cherry-picking.
