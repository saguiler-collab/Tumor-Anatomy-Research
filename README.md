# The Anatomy Test

**Anatomic Concordance: Neuropathology as a Ground-Truth-Free Benchmark for Cell-Type
Deconvolution in Glioblastoma**

> Can a tumour's own anatomy stand in for ground truth when choosing a cell-type
> deconvolution method — and does a method that gets the anatomy right also get the
> numbers right?

`Anatomy_Test.md` is the protocol and takes precedence over this file. This README is
the orientation: what the project is, why it exists, what it found, and what it cannot
claim.

---

## 1. The problem, in plain terms

### What deconvolution is

When you sequence a piece of tumour, you do not get one cell's worth of RNA. You get the
average of every cell in the sample — cancer cells, immune cells, blood vessels, normal
brain — blended together. **Deconvolution** is the statistical attempt to run that
blending backwards: given the mixed signal and a reference describing what each cell type
looks like on its own, estimate what fraction of the sample each type made up.

This matters because composition is often the biology. Whether a tumour is full of
macrophages or nearly free of them changes what it is and how it might be treated.

### Why nobody knows if it works

Here is the difficulty that motivates this entire project: **for almost any real tissue,
nobody knows the true answer.** You cannot count every cell in a tumour. So when a
deconvolution method reports "23% macrophages", there is nothing to check it against.

Dozens of methods exist and they disagree with each other, sometimes wildly. Faced with
no ground truth, the field falls back on two habits:

1. **Consensus** — run several methods, trust what they agree on. But methods sharing an
   assumption can be wrong together, and agreement is not accuracy.
2. **"It matches known biology"** — the estimates look like what we expect, so the method
   must be working.

The second habit is the one this project is about, because it has a structural flaw.

### The circularity

Suppose you deconvolve tumour regions and find more immune cells near dead tissue. You
report this as biologically sensible, and take it as evidence the method worked.

But you decided the pattern was sensible **after seeing it**. Had the method produced the
opposite, you could as easily have called that "an unexpected finding". The check can
never fail, so it is not a check.

```
HOW IT IS USUALLY DONE
    method ──solves──> fractions per region ──read as──> anatomy      (a CONCLUSION)

THE ANATOMY TEST
    constraints written in advance ──> anatomy ──scores──> method     (a TEST)
```

**Write the expected ordering down first, hash it, and the same pattern becomes a
prediction the method can miss.** That is the entire idea. Everything else is machinery
to make it honest.

---

## 2. Why glioblastoma, and why its anatomy is unusual

Glioblastoma (GBM) is the most aggressive primary brain tumour; median survival is
roughly 15 months. It is also, for this purpose, unusually well suited — because a
pathologist looking at a GBM slide can name distinct anatomic structures **by sight**,
and those structures are *defined by their cell content*.

The Ivy Glioblastoma Atlas Project (Ivy GAP) did something rare: it cut those structures
apart under a microscope with laser microdissection, then sequenced each one separately.
Five structures:

| | structure | what it is |
|---|---|---|
| **LE** | Leading edge | The outer margin. Brain tissue infiltrated by tumour but still largely intact — normal cells, including the myelin-making oligodendrocytes of white matter, are still present. |
| **IT** | Infiltrating tumour | The intermediate zone. More tumour than LE, less than CT. |
| **CT** | Cellular tumour | The dense core. Tumour cells have crowded out the normal brain. |
| **MVP** | Microvascular proliferation | Florid, disordered blood-vessel growth. This is one of the two features that *define* GBM under WHO diagnostic criteria. |
| **PAN** | Pseudopalisading cells around necrosis | Tumour cells crowded in dense rows around a patch of dead tissue. The tissue is hypoxic — starved of oxygen. The other WHO-defining feature. |

### The key move

Because these structures are **defined by** what they contain, some compositional facts
are known without ever measuring them:

- MVP is *diagnosed by* proliferating blood vessels. So endothelial cells — the cells
  lining blood vessels — must be more abundant there than in the tumour core. That is
  not a hypothesis. It is close to a tautology.
- LE is *defined as* preserved brain tissue. So oligodendrocytes must be more abundant
  there than in the core, where the brain has been displaced.
- Hypoxia recruits myeloid cells, so the peri-necrotic zone should be macrophage-rich
  relative to the leading edge.

**Anatomy gives free, partial ground truth** — not exact fractions, but reliable
statements about which way the ordering must go. And unlike real ground truth, it costs
nothing and already exists in any tissue with an anatomical atlas.

---

## 3. The seven constraints

Fixed before any deconvolution output was inspected, and SHA-256 hashed so it can be
proven they were not edited afterwards:

```
2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a
```

| id | claim | kind | weight |
|---|---|---|---|
| C1 | Tumour cells higher in **CT** than **LE** | pairwise | 1 |
| C2 | Oligodendrocytes higher in **LE** than **CT** | pairwise | 1 |
| C3 | Endothelial higher in **MVP** than **CT** | pairwise | 1 |
| C4 | Endothelial **highest in MVP** of all five structures | maximum | 1 |
| C5 | Macrophage/microglia higher in **PAN** than **LE** | pairwise | 1 |
| C6 | Macrophage/microglia higher in **MVP** than **CT** | pairwise | 1 |
| C7 | Tumour cells increase **LE < IT < CT** | monotone | 2 |

Notes that matter:

- **C7 carries double weight** because a three-step monotone chain is much harder to
  satisfy by chance than a single pairwise comparison.
- **C4 is a strictly harder version of C3.** A method can put endothelial cells above CT
  and still not make MVP the maximum. Scoring both separates "roughly right" from "right".
- **T cells were excluded in advance.** This project's own synthetic benchmark puts them
  at the detection floor, so a T-cell constraint would be scoring noise. Excluding it
  before looking is honest; excluding it after would not be.
- **Neuronal content is untestable here** — and the recorded reason was later found
  incomplete. The roster has no neuron column, but the deeper problem is that the
  reference atlas contains **22 neurons**, so a neuron column could not have been
  estimated anyway. The constraint file is frozen, so this correction is recorded
  *alongside* it rather than by editing it. Editing would change the hash and void every
  result computed under it.

**Anatomic Concordance Score (ACS)** = the weighted proportion of (constraint × tumour)
pairs a method satisfies.

---

## 4. What makes this a test and not a demonstration

Four safeguards do the real work. Without them this would be a story, not a study.

### The negative controls are the point

Two fake "methods" run through the identical pipeline:

- **`control_random`** — makes up compositions at random.
- **`control_shuffled_signature`** — deconvolves against a randomly permuted reference.

If a nonsense method can satisfy the constraints, the constraints are measuring nothing,
and the whole leaderboard is decoration. Each control is run as a **distribution over 200
draws**, not a single lucky roll.

The protocol's rule is explicit and uncomfortable: **if the controls score well, publish
that.** Do not retune the constraints and re-report. A benchmark that only reports
flattering results is not a benchmark.

### The null is not a coin flip

It is tempting to say a method satisfying a constraint "by chance" would do so 50% of the
time. That is wrong twice over: cell fractions are **compositional** (they must sum to
one, so they are inherently correlated), and samples are **nested** — several tissue
blocks come from each tumour, so they are not independent observations.

Instead the null shuffles **structure labels within each tumour**, 10,000 times. This
preserves each tumour's composition and its number of blocks per structure, destroying
only the thing under test.

### Nothing is ever selected on an outcome

- Method selection never uses survival. An assertion aborts the analysis if the recorded
  selection criterion so much as mentions an outcome.
- Method selection never uses ACS either — ACS is the quantity being tested, and
  selecting on it would assume the conclusion.

### The circular samples are excluded by code

Ivy GAP contains 148 further samples whose structure label was assigned **from
expression** rather than from histology. Scoring anatomic constraints on those would be
viciously circular. An assertion requires the H&E-study flag and aborts otherwise — the
rule is enforced, not remembered.

---

## 5. How the experiment actually works

Two arms. Every method appears in both, so the comparison is paired.

**Arm 1 — accuracy against real ground truth.** Synthetic mixtures are built from a
single-cell atlas (GBmap, 338,564 cells) by pooling known quantities of known cell types.
Here the true composition *is* known, by construction. Crucially the split is
**donor-held-out**: test mixtures come only from donors whose cells were excluded from
the reference, so no method is graded on people it was trained on.

**Arm 2 — anatomic concordance.** The same methods deconvolve real Ivy GAP tissue and are
scored on the seven constraints.

**The experiment** is the comparison of the two rankings:

> Does the ordering produced by anatomy match the ordering produced by truth?

One number — Spearman ρ — judged against a bar fixed in code before it was computed:

| outcome | meaning |
|---|---|
| ρ ≥ 0.60, CI excludes zero | Anatomy **can** be used to choose a method without ground truth. |
| ρ ≥ 0.60, CI includes zero | Underpowered. Reported as such. |
| ρ < 0.60 | Reproducing known biology is **not** evidence of numerical accuracy. |

**All three are publishable and the code treats them that way.** `agreement.py` returns
the null result as a first-class outcome, not an error. Registering a rule that can
return an unwelcome answer is the point.

---

## 6. What was found

From `results_archive/2026-09-06T1103` (97 files, each hash-verified).

### The controls behave

Best control **0.400** against a median real method of **0.954**. Neither control beats
its own permutation null (p = 0.43 and p = 0.999). The constraint set discriminates —
so the leaderboard means something.

### The leaderboard (122 samples, 9 evaluable tumours)

| method | ACS | 95% CI | implementation |
|---|---|---|---|
| MuSiC | 1.000 | [1.000, 1.000] | R:MuSiC |
| NNLS / SVR / Elastic Net / EPIC | 0.985 | [0.953, 1.000] | mixed |
| CIBERSORTx (B-mode) | 0.969 | [0.930, 1.000] | this project |
| SCDC / SCDC ENSEMBLE | 0.954 | [0.911, 0.986] | R:SCDC |
| Bisque | 0.923 | [0.866, 0.983] | R:BisqueRNA |
| BayesPrism | 0.877 | [0.786, 0.968] | reimplementation |
| Bayesian / Hierarchical | 0.769 | [0.627, 0.921] | this project |
| DWLS | 0.738 | [0.600, 0.867] | reimplementation |
| *control_random* | *0.400* | | |
| *control_shuffled_signature* | *0.138* | | |
| quanTIseq | 0.600 | [0.333, 0.882] | **not ranked** — partial coverage |

### The headline

**ρ = 0.733** (95% CI [0.271, 0.944], p = 0.0044, 13 methods) — above the pre-registered
0.60 bar with the interval excluding zero. On the wider 270-sample cohort, **ρ = 0.750**.
The two ACS rankings agree with each other at **ρ = 0.987**.

> **Anatomic concordance ranks deconvolution methods approximately the way ground truth
> does — using no ground truth and no outcomes.**

### And the honest caveat that goes with it

**ACS cannot resolve the top of its own leaderboard.** MuSiC "wins" by exactly **one
constraint–tumour pair out of 57** — it and NNLS differ on C5 alone (6/6 versus 5/6) and
are identical on the other six constraints. Six methods have confidence intervals that
include MuSiC's score. Nine tumours and a weighted denominator of 65 mean ACS can take
only 66 distinct values.

So the correct claim is *"anatomy separates good methods from bad ones"*, not *"anatomy
identifies the best method"*. The second sentence is not supported and is not made.

---

## 7. Limitations, stated plainly

- **One tissue is a case study, not a method.** Nine evaluable tumours is the ceiling and
  no analysis inside Ivy GAP raises it. A second tissue is the only fix.
- **No prognostic claim is possible.** Ivy GAP publishes survival times with **no
  vital-status column**, so who was censored is unknowable and no C-index can be computed.
  The analysis reports **BLOCKED** rather than guessing. Worse, of 42 tumours the 10 with
  blank times are **not missing at random** — a blank is associated with MGMT methylation
  (Fisher p = 0.0021), the best favourable prognostic factor in GBM, so dropping them
  would bias the cohort toward short survival. Under an explicitly declared assumption the
  result is **UNDERPOWERED** (29 events, detectable difference ≈ 0.515) and reported as
  INCONCLUSIVE.
- **The Astrocyte column is unidentifiable** against Tumor in this reference. It is
  estimated and reported but never permitted to carry a definitional claim.
- **Two published packages fell back** to reimplementations in one run after exceeding a
  wall-clock budget. Disclosed per method; a reimplementation is never reported under a
  published package's name.
- **7.0% of the reference atlas is dropped** by the eight-type roster. Those cells do not
  vanish from the tissue, only from the model, so their signal is absorbed by the nearest
  retained type. The largest dropped population is myeloid, and the nearest retained
  column is Macrophage/Microglia — which carries C5 and C6. The absorption is therefore
  concentrated on two of seven constraints rather than spread evenly.
- **This registration is retrospective.** Analyses were run before registering. What the
  hash certifies is narrower and still the load-bearing claim: the constraint file
  predates any inspected deconvolution output.

---

## 8. Repository map

| path | what it is |
|---|---|
| `Anatomy_Test.md` | The protocol. Takes precedence over everything else. |
| `ivygap/anatomic/constraints.py` | The pre-registered constraint file. Hashed. |
| `ivygap/anatomic/acs.py` | ACS scoring, permutation null, bootstrap. |
| `ivygap/anatomic/agreement.py` | The experiment, and the pre-registered decision rule. |
| `ivygap/deconv/` | The methods: baselines, published packages, R bridge, controls. |
| `ivygap/config.py` | The control panel. No path or cell-type name is hard-coded elsewhere. |
| `docs/METHODS.md` | Every method and every deviation from a published algorithm. |
| `docs/DATA_SOURCES.md` | What the archive contains, verified against the release. |
| `docs/EXTERNAL_ACTIONS.md` | The five things that must be done outside this repo. |
| `docs/REGISTRATION_ANSWERS.md` | Field-by-field text for the OSF registration. |
| `results_archive/` | Read-only, hash-verified snapshots of completed runs. |
| `scripts/run_all.py` | The whole pipeline, six stages. |
| `scripts/fetch_allen.py` | Allen Human Brain Atlas fetcher, for the second tissue. |

```bash
pytest tests/ -q                     # ~10 minutes, 198 tests
python3 scripts/run_all.py --synthetic   # end-to-end on fixtures, proves nothing scientific
python3 scripts/run_all.py --full-database --control-draws 200
python3 scripts/archive_run.py --verify <stamp>
```

---

## 9. What comes next

1. **Register the constraint file** (`docs/EXTERNAL_ACTIONS.md` item 1). Free, and it is
   what turns "we say we wrote this first" into something a reader can verify.
2. **A second tissue.** The Allen Human Brain Atlas is downloadable today, but the data
   is not the blocker: the constraints, the cell-type roster and the platform all have to
   be dealt with first, and a normal-brain run needs its own pre-registered constraint
   file written before any of its output is seen.
3. **A run that includes CIBERSORTx S-mode**, which postdates both archived runs.
