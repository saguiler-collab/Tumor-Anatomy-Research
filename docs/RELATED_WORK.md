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

## 3 · Where our results disagree, and why

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
