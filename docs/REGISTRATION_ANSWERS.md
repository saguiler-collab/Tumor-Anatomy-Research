# OSF registration — field-by-field answers

Copy-paste text for the OSF Preregistration form. Every number here is measured from
the artefacts, not typed from memory; regenerate any of them with the commands noted.

**Read this first — honesty about status.** Two complete analyses have already been run
(`results_archive/2026-09-05T2154`, `results_archive/2026-09-06T1103`). This registration
is therefore NOT a preregistration of an unrun analysis, and must not be presented as
one. What it timestamps is narrower and still worth having: the constraint file was
written and frozen **before any deconvolution output was inspected**, and that is the
design claim the whole study rests on. Wherever OSF asks about data collection or
analysis status, say that data are archival and already collected, and that pilot
analyses have been completed. Then re-run after registering so the reported results
postdate the registration.

---

## The hash

```
2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a
```

SHA-256 of the canonical JSON payload of the seven pre-registered anatomic constraints —
their ids, cell types, structures, directions, kinds and weights. It is **not** a hash of
the source file, so reformatting or re-commenting `constraints.py` does not change it,
while altering a single claim does. Regenerate:

```bash
python3 -c "import sys;sys.path.insert(0,'.');from ivygap.anatomic import constraints as K;print(K.freeze_hash())"
```

`ivygap/anatomic/registration.py` compares this against `REGISTRATION.json` and reports
`REGISTERED`, `UNREGISTERED` or `HASH_MISMATCH`. `HASH_MISMATCH` is the fatal state: it
means the constraints changed after registration, and every result computed under them
is void.

---

## Page 3 — Research Design

### Study type

**Simulation study.** OSF defines it as "using existing or synthetic data to assess
performance of a model ... including demonstration of methods", which is precisely this
study: it assesses the performance of 15 cell-type deconvolution methods against both a
synthetic yardstick with known truth and real archival tissue.

> If the control accepts more than one, also tick **Non-randomized study** — the
> anatomic arm is an observational secondary analysis of existing human tissue data,
> with no assignment of any kind.

### Intention for causal interpretation

**No causal relationship inferred.**

The primary estimand is a rank correlation between two orderings of the same methods. It
asks whether one ranking predicts another, not whether anything causes anything.

### Blinding of experimental treatments

Tick: **"Subjects will not be aware of the assigned treatment during data collection
(either because the subjects are not human participants or because of blinding
procedures)."**

Subjects are archival tumour tissue specimens from a public release. There are no living
participants, no assignment, and no treatment.

### Additional blinding during research or analysis

> There is no treatment to blind, but the study's validity depends on three masking
> procedures, all enforced in code rather than by convention.
>
> First, the constraint file was written, frozen and SHA-256 hashed
> (2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a) before any
> deconvolution output was inspected. It is never edited in response to a result: where
> a stated rationale was later found to be incomplete, the correction is recorded
> alongside the file rather than applied to it, because editing would change the hash and
> void every result computed under it.
>
> Second, no deconvolution method receives the anatomic structure label. Structure is
> used only in scoring, after every method has produced its estimate. An equal-footing
> check verifies that the manifest handed to the methods carries no outcome column, and
> aborts if it does.
>
> Third, method selection is blind to outcome by construction. The selection criterion is
> recorded before selection and an assertion aborts the survival stage if the recorded
> criterion mentions an outcome. Selection is likewise never made on the Anatomic
> Concordance Score itself, since ACS is the quantity under test and choosing on it would
> assume the conclusion.

### Study design

> **Design.** A two-arm methods-evaluation study with a paired, within-method
> comparison. Fifteen cell-type deconvolution methods plus two negative controls are run
> on identical inputs, and each is scored twice: once against known ground truth on
> synthetic mixtures, and once against pre-registered anatomic expectations on real
> tumour tissue. The primary analysis is the rank correlation between those two
> orderings. Every method appears in both arms, so the comparison is fully within-method
> (paired); there is no between-subject factor and no assignment.
>
> **Arm 1 — accuracy against known truth (the yardstick).** Pseudobulk mixtures are
> synthesised from a single-cell atlas (GBmap, 338,564 cells) with composition known by
> construction. The split is donor-held-out: mixtures for testing are built only from
> donors whose cells were excluded from the reference the methods solve against, so no
> method is scored on a donor it was trained on. Accuracy is mean absolute error on the
> primary cell types, aggregated donor-equally.
>
> **Arm 2 — anatomic concordance (the quantity under test).** The same methods
> deconvolve bulk RNA-seq from anatomically annotated GBM tissue (Ivy GAP), sampled by
> H&E-guided microdissection into five structures: leading edge (LE), infiltrating tumour
> (IT), cellular tumour (CT), microvascular proliferation (MVP), and pseudopalisading
> cells around necrosis (PAN). Each method's estimates are collapsed to (tumour,
> structure) means and scored against seven ordinal constraints fixed in advance:
>
> | id | claim | kind | weight |
> |---|---|---|---|
> | C1 | Tumor higher in CT than LE | pairwise | 1 |
> | C2 | Oligodendrocyte higher in LE than CT | pairwise | 1 |
> | C3 | Endothelial higher in MVP than CT | pairwise | 1 |
> | C4 | Endothelial highest in MVP among all five structures | maximum | 1 |
> | C5 | Macrophage/microglia higher in PAN than LE | pairwise | 1 |
> | C6 | Macrophage/microglia higher in MVP than CT | pairwise | 1 |
> | C7 | Tumor increases monotonically LE < IT < CT | monotone | 2 |
>
> The Anatomic Concordance Score (ACS) is the weighted proportion of satisfied
> (constraint x tumour) pairs. C7 carries double weight because a three-step monotone
> chain is far harder to satisfy by chance than a single pairwise step. Constraints a
> method structurally cannot address are excluded from numerator and denominator alike
> rather than counted as failures.
>
> **Unit of analysis and nesting.** Multiple blocks are sampled per tumour, so samples
> are nested within tumours. Estimates are collapsed to (tumour, structure) before any
> aggregation, and every interval and permutation resamples or permutes **within
> tumours**. The evaluable set is 122 H&E-selected anatomic samples across 9-10 tumours;
> a wider 270-sample cohort spanning 37 patients is deconvolved as a sensitivity check
> but scored on the anatomic subset only.
>
> **A deliberate exclusion.** Ivy GAP also contains 148 samples whose structure label was
> assigned from expression rather than histology. Scoring anatomic constraints on those
> would be circular, so they are excluded from ACS by an assertion that requires the
> H&E-study flag, not by convention.
>
> **Negative controls.** Two controls run through the identical pipeline: random
> compositions, and deconvolution against a shuffled signature matrix. These are the
> load-bearing check. A constraint set that a control can satisfy is measuring nothing,
> so the controls are calibrated over 200 independent draws each and reported whatever
> they score. If controls score highly, that is published as the finding rather than
> answered by retuning the constraints.
>
> **Comparability.** Methods that model only part of the cell-type roster are reported
> with their own denominator and explicitly excluded from the ranking, the tie groups,
> the control-check median and the primary correlation, because a score over two
> constraints is not on the same scale as a score over seven.

### Randomization

> Not a randomized study; no subjects are assigned to conditions. Randomness enters only
> as analysis machinery, seeded from a fixed constant so every artefact is reproducible:
> (a) the within-tumour permutation null, 10,000 draws, which permutes structure labels
> inside each tumour and leaves the nesting and the compositional dependence intact;
> (b) the bootstrap over tumours, 2,000 resamples, for interval estimates;
> (c) donor-balanced subsampling of the single-cell atlas when building the reference;
> (d) the donor-held-out split defining which donors build the reference and which build
> the test mixtures.

---

## Page 4 — Sampling

### Existing data / data collection status

**Data already collected; analysis already begun.** Both datasets are public archival
releases. No new data are collected and no participants are recruited.

### Data source

> **Ivy Glioblastoma Atlas Project (Ivy GAP)**, RNA-Seq release of 2014-11-25 (Allen
> Institute for Brain Science; Puchalski et al., *Science* 2018). Bulk RNA-seq of
> laser-microdissected anatomic structures from human glioblastoma.
>
> **GBmap Core** single-cell atlas, 338,564 cells, used only to build the cell-type
> reference and the synthetic mixtures — never scored.

### Sample size, and how it was arrived at

> Sample size is fixed by the archive; nothing was chosen. The release contains **270
> RNA-seq samples from 37 tumours**. Of these, **122 samples across 10 tumours** carry
> H&E-guided anatomic structure labels and are the pre-registered scoring set. The
> remaining 148 had their structure assigned from expression and are excluded from
> scoring as circular.
>
> After collapsing to (tumour, structure) and requiring both structures of a constraint
> to be present in a tumour, **9 tumours** contribute at least one evaluable pair, giving
> **57 constraint-tumour pairs** and a weighted denominator of **65**.
>
> The synthetic arm uses **500 pseudobulk mixtures** built from **22 held-out donors**,
> with the reference built from the other **88**.

### Stopping rule

> None applicable. Both datasets are closed archival releases of fixed size; there is no
> sequential collection and no interim analysis that could stop it.

### Power / precision, stated honestly

> This is the study's principal limitation and it is reported as a result rather than
> buried. Nine evaluable tumours and a weighted denominator of 65 mean ACS can take only
> 66 distinct values, so several methods necessarily tie, and the bootstrap intervals of
> the top methods overlap heavily. The study is powered to detect a strong rank
> association (the pre-registered bar is rho >= 0.60) and is **not** powered to separate
> adjacent methods. Any claim that one top-ranked method beats another is reported as
> INCONCLUSIVE. Widening this requires a second tissue, which no analysis inside Ivy GAP
> can substitute for.

---

## Page 5 — Variables

### Manipulated variables

> None. Nothing is manipulated; no condition is assigned.

### Measured variables

> **Per method, per sample — the estimate.** A composition over eight cell types
> (Tumor, Macrophage/Microglia, T cell, NK cell, B cell, Endothelial, Oligodendrocyte,
> Astrocyte), non-negative and summing to one. Converted from mRNA share to cell share by
> a per-type mRNA-content factor, applied once and centrally.
>
> **Primary outcome — Spearman rho** between two rankings of the same methods: their ACS
> ranking and their accuracy ranking on held-out synthetic mixtures.
>
> **ACS** — the weighted proportion of satisfied (constraint x tumour) ordinal pairs,
> per method.
>
> **Accuracy** — mean absolute error against known composition, donor-equally
> aggregated, excluding the Astrocyte column.
>
> **Reported alongside, never as the headline:** a within-tumour permutation p-value per
> method, a bootstrap interval over tumours, and the two negative controls' scores.

### Indices, and one deliberate exclusion

> The Astrocyte column is excluded from the primary accuracy metric. It is known to be
> unidentifiable against the Tumor column in this reference, and including it would let a
> method's handling of an acknowledged-unresolvable column move the headline number. It
> is still estimated and reported; it is simply not permitted to carry a definitional
> claim.

---

## Page 6 — Analysis Plan

### Statistical models

> **Primary.** Spearman rank correlation between each method's ACS and its accuracy on
> donor-held-out synthetic mixtures, across all methods scored over the full cell-type
> roster. Pre-registered decision rule, fixed in `agreement.py` before any result:
>
> - rho >= 0.60 **and** the bootstrap CI excludes zero -> anatomic concordance is usable
>   for method selection;
> - rho >= 0.60 with a CI including zero -> underpowered, reported as such;
> - rho < 0.60 -> the finding is that "matches known biology" is **not** evidence of
>   accuracy. That is a real result, not a failure, and is reported with equal prominence.
> - fewer than 6 comparable methods -> no correlation is reported at all.
>
> **Inference on ACS.** The null is **not** a coin flip. Cell fractions are compositional
> and correlated, and samples are nested in tumours, so a 50% null would be wrong in both
> directions. The null is instead a within-tumour permutation of structure labels,
> 10,000 draws, which preserves each tumour's composition and the number of blocks per
> structure while destroying only the structure assignment. Intervals come from a
> bootstrap over tumours, 2,000 resamples.
>
> **Ties.** A paired bootstrap over donors identifies which methods are statistically
> indistinguishable from the best. Methods with identical ACS are marked as a tie group
> so a reader counting points for the rank correlation counts the right number.

### Transformations

> Counts are normalised to CPM, never log-transformed for deconvolution — the mixing
> model is additive on the linear scale. Genes are restricted to a marker subset selected
> from the reference alone, with no reference to any outcome or structure. One method
> (quanTIseq) ships its own fixed signature and receives the full gene space instead; the
> departure and its reason are declared in the artefacts.

### Inference criteria

> Two-sided, alpha = 0.05, with the pre-registered rho >= 0.60 as the substantive bar. A
> permutation p-value below 0.05 is necessary but not sufficient for any claim: a method
> must also exceed both negative controls.

### Data exclusion

> Stated in advance and enforced in code, not by judgement:
> - the 148 expression-labelled samples are excluded from ACS scoring (circular);
> - a (tumour, constraint) pair is excluded where the tumour lacks a required structure;
> - a sample a method fails to solve is NaN and counted, never imputed;
> - no sample is excluded for being an outlier, for any definition of outlier.

### Missing data

> Nothing is imputed. Where a required input is missing the analysis reports BLOCKED
> rather than proceeding. This is load-bearing for the prognostic question: the Ivy GAP
> release publishes survival times with **no vital-status column**, so who was censored
> is unknowable and no C-index can be computed. Of 42 tumours, 32 have a recorded time
> and 10 are blank, and the blanks are **not** missing at random — a blank is strongly
> associated with MGMT methylation (Fisher p = 0.0021), the strongest favourable
> prognostic factor in GBM, so dropping them would bias the cohort toward short survival.
> A survival analysis is therefore reported only under an explicitly declared assumption,
> into a separately labelled directory, and never as the headline.

### Exploratory analyses

> Everything not named above is exploratory and labelled as such, including: the wider
> 270-sample cohort, the per-constraint and per-tumour breakdowns, the comparison of
> CIBERSORTx B-mode against S-mode, and any survival analysis.

---

## Page 7 — Other

> **Prior work by the authors on this dataset.** Two complete analyses have already been
> run and are archived with per-file hashes (`results_archive/2026-09-05T2154`,
> `results_archive/2026-09-06T1103`). They are pilot runs. The constraint file predates
> both and is unchanged, which is what the hash above certifies. Results reported as
> confirmatory will come from a run executed after this registration.
>
> **Deviations already known.** Several methods are run with documented departures from
> their published defaults, each recorded with its reason and the evidence the decision
> was made from, all of it outcome-independent. Two published packages exceeded their
> wall-clock budget in one stage and fell back to reimplementations; the fallback is
> disclosed per method rather than absorbed. A Python reimplementation is never reported
> under a published package's name.
>
> **Code and artefacts.** The full pipeline, the constraint file, the negative controls
> and every artefact required to reproduce these numbers accompany this registration.
