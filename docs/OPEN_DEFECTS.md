# Open defects

Things known to be wrong, written down so they survive a context compaction and so
nobody has to rediscover them. Each entry states what was **verified** versus what is
**suspected**, because the difference decides whether a fix is safe.

Closed defects live in git history and in `docs/METHODS.md`; this file is only what is
still open.

---

## D1 · Cell-size correction is applied twice to some R methods

**Severity: high.** It affects the accuracy numbers that the primary result is computed
from, and one of the methods involved is the top-ranked one.

### The rule this breaks

`docs/METHODS.md` opens with it, and the OSF registration (<https://osf.io/dm2t8>) states
it in *Measured variables*:

> "Converted from mRNA share to cell share by a per-type mRNA-content factor, **applied
> once and centrally**."

So this is not a deviation to declare. It is the code failing to do what was registered,
and fixing it makes the code match the record rather than departing from it.

### Why it happens

`DeconvolutionMethod.fit_predict` converts RNA share to cell share for every method whose
coverage is full:

```python
apply_cs = data.apply_cell_size_correction and covered is None
```

That is correct for a solver that returns a share of the **transcript pool** — NNLS, SVR,
elastic net, the Bayesian pair — because those solve `bulk ≈ S @ w` against an expression
profile and `w` is an mRNA share by construction.

It is wrong for a package that already did the conversion internally. quanTIseq was
caught (`scale_mRNA = TRUE` is the same correction, and the code skips ours for it). The
others were not.

### EPIC — CONFIRMED, verified empirically 2026-09-10

EPIC returns **both** quantities as separate elements, and they differ:

```
EPIC returns: mRNAProportions, cellFractions, fit.gof

sample 1     Bcells   CAFs  CD4_T   CD8_T  Endo   Macrophages  NK   otherCells
mRNA         0.0099  2e-04  0.0292  0.0089 0.0248  0.0411      0     0.8858
cell         0.0102  2e-04  0.0304  0.0093 0.0255  0.0119      0     0.9124
identical?   FALSE
```

Macrophages move by a factor of ~3.5 between the two. `R/run_epic.R:44` takes
`est$cellFractions` — the already-corrected one — and `fit_predict` then applies
`to_cell_fractions()` on top.

**Fix:** take `est$mRNAProportions` so the central correction is the only one. Unambiguous,
because the package hands you the choice explicitly.

### All five now VERIFIED — 2026-09-10

Each package was run on a purpose-built reference where the two cell types differ 3x in
mRNA content, mixed **50/50 by cell count**. That separates the two possible answers
cleanly and puts neither at a boundary:

- if a package returns **0.50**, it already divided out cell size — it reports a CELL share
- if it returns **0.75**, it reports an mRNA share and our correction is the first one

(An earlier version of this test used a 10x ratio, which put the mRNA answer at 0.988 —
both answers crowd the ceiling and the result is unreadable. The ratio matters.)

| method | returned | reports | our correction is |
|---|---|---|---|
| **MuSiC** | **0.500** | CELL share | **a SECOND correction — defect** |
| **Bisque** | **0.500** | CELL share | **a SECOND correction — defect** |
| **EPIC** | *(offers both)* | CELL share, as taken | **a SECOND correction — defect** |
| BayesPrism | 0.750 | mRNA share | correct, the first |
| SCDC | 0.653 | between the two | **partial** — see below |

#### That table is now reproducible — added 2026-09-12

**The five numbers above were measured once by a script that was never saved, and no
artefact recorded them.** They existed only as this table and its copy in
`ROAD_TO_PAPER.md`. That breaks the project's own rule — *measure and record, never
validate against hard-coded literals* — and a blocking design decision rested on numbers
neither a reviewer nor this project could re-check.

`scripts/verify_cell_size_semantics.py` now runs the probe and writes
`results/cell_size_semantics.json`. It is not a benchmark: a package is not better for
using either convention, and the script reports no ranking.

**Status of the re-measurement:**

| row | prose value | re-measured | agrees |
|---|---|---|---|
| BayesPrism | 0.750 | **0.7500** (reimplementation) | yes |
| MuSiC | 0.500 | pending — needs R | — |
| Bisque | 0.500 | pending — needs R | — |
| EPIC | offers both | pending — needs R | — |
| SCDC | 0.653 | pending — needs R | — |

The four R rows are pending only because the cores were busy with the BayesPrism
re-measurement; nothing about them is in doubt yet. **Until they are measured, treat the
prose values as unverified rather than as wrong.**

#### What building the probe clarified about the defect's scope

The probe's first version passed *raw* single cells to `build_reference`, which put each
type's library size into the profile columns. Plain NNLS then returned 0.500 — apparently
"already a cell share". That is an artefact of the probe, not a property of NNLS, and it
matters because production does the opposite: `build_from_h5ad` normalises every cell to
1e6 **before** averaging into the profile and captures the raw library sizes separately as
the cell-size factors. Every profile column therefore carries the same total.

The consequence, which was implicit before and is now measured: because the signature's
columns are equal-total, a least-squares solve against it returns **mRNA share**, and the
central conversion is that family's *first* conversion and is correct. NNLS comes out at
exactly 0.7500 — the theoretical value, which is what certifies the probe is built right.

**So the defect is narrower and better defined than "the correction is applied twice".**
It is specific to packages that perform their *own* internal size conversion, using their
own estimate. It does not touch the least-squares family (NNLS, SVR, both CIBERSORTx
modes, Elastic Net), which is measured above and confirmed. That does not shrink the
defect's *impact* — 15 of 16 methods still move when the correction is applied twice,
because the renormalisation is composition-dependent and a changed column changes every
other column in the sample.

**Three of the top methods are double-corrected, including the top-ranked one.** MuSiC
leads the ACS leaderboard at 1.000 and the accuracy benchmark at MAE 0.0508, and its
estimates have had a per-type cell-size factor applied twice.

**SCDC is the awkward one.** At 0.653 it sits between 0.50 and 0.75, closer to the mRNA
answer. It does *some* internal size handling — `ct.cell.size` defaults to library size
derived from the data — without fully converting. So neither "apply ours" nor "skip ours"
is right for SCDC, and the third design below (pass this project's factors in via
`ct.cell.size`, skip the central step) is the only one that handles it cleanly.

### What this changes — CORRECTED 2026-09-11

**An earlier version of this section claimed the double correction "cannot change ACS,
because a per-type constant cannot reorder a column across structures". That claim was
wrong, was repeated several times, and is retracted here.** It was asserted from the
shape of the formula and never tested. When tested, it fails.

`to_cell_fractions` divides by cell size **and then renormalises each sample to sum 1**:

```python
out = rna_fractions / cell_size
return project_to_simplex(out)          # <- the renormalisation
```

The renormalisation divisor is `sum_j(w_j / c_j)`, which **depends on that sample's own
composition**. So the effective scaling applied to cell type k is not `1/c_k` but
`(1/c_k) / sum_j(w_j/c_j)`, and that denominator differs from sample to sample. Two
samples with the same true ordering of type k can come out ordered differently. The
per-column-constant argument would hold for the division alone; it does not survive the
renormalisation.

**Measured, by applying the correction a second time to each archived estimate table:**

| | |
|---|---|
| methods whose ACS changes | **15 of 16** |
| largest single change | **0.0462** (Bisque, 0.9231 -> 0.8769) |
| methods whose RANK moves | **7 of 14** |
| Spearman(before, after) | 0.968 |

DWLS rises from 14th to 12th; EPIC and SVR rise to joint first with MuSiC; NNLS and
Elastic Net fall from 2nd to 4th.

**So the severity is higher than previously stated.** The defect affects:

- **ACS**, and therefore the leaderboard ordering;
- **the accuracy arm**, where a doubled correction is an error in magnitude;
- **rho**, which correlates the two.

Both arms of the agreement test move. Nothing about the *design* is invalidated — the
controls, the permutation null, the constraint file and the ISH validation are untouched,
because none of them depends on the cell-size conversion. But every number in the
leaderboard should be treated as provisional until the fix lands.

### Why it was not fixed on discovery

Found mid-run, while the confirmatory run (the first to postdate registration) was in
stage 4a. Changing method behaviour then would have produced a run that was half one
definition and half another, which is worse than a known, documented bias applied
uniformly. The run finishes as-is.

### What the fix must include

- Verification per method, not a blanket change.
- **Before/after numbers for every affected method**, published, so the effect on the
  leaderboard is visible rather than silent. The correction changes accuracy values, which
  changes the ranking, which changes rho. That must be shown, not absorbed.
- A note in `RESULTS.md` and in the paper that the pre-correction runs
  (`2026-09-05T2154`, `2026-09-06T1103`, and the confirmatory run) carry the bias.

### What it does NOT invalidate

**This section previously restated the retracted claim and contradicted the correction
above. Rewritten 2026-09-12.**

What the defect does **not** touch, because none of these depends on the cell-size
conversion:

- **The constraint file** and its hash — registered, and independent of any method.
- **The negative controls.** `control_random` is unaffected entirely (delta 0.0000 when
  the correction is applied a second time); the shuffled-signature control moves by
  0.0154. Neither comes near a real method either way.
- **The permutation null**, which shuffles structure labels within a tumour and is
  computed on whatever estimates it is given.
- **The ISH validation**, which involves no deconvolution at all.
- **The design** — pre-registration, equal footing, the nesting rule, outcome-blind
  selection.

What it **does** touch is stated above: ACS, the leaderboard ordering, the accuracy arm,
and rho. Every leaderboard number is provisional until the fix lands.

---

## D2 · SCDC ENSEMBLE runs degenerate for want of a second reference

**Severity: medium.** Disclosed correctly, so it is honest — but a published tool is
running with its entire contribution switched off.

`scdc_ensemble`'s recorded degeneracy reason:

> "one reference supplied ('gbmap'), so there is nothing to weight across: SCDC ENSEMBLE
> reduces exactly to SCDC. Reported as degenerate rather than as an ENSEMBLE result."

**The fix is a second GBM single-cell reference**, passed as `sc.eset.list`. Neftel et al.
2019 (doi:10.1016/j.cell.2019.06.024) is the obvious candidate and is already cited in
`docs/DATA_SOURCES.md`.

**A tempting non-fix, to be avoided.** `SCDC_ENSEMBLE` also accepts `prop.input`, "a list
of SCDC_prop objects ... allows users to ENSEMBLE results calculated from other
deconvolution methods." That would un-degenerate it without a second reference, and it
must not be used here: it would make SCDC ENSEMBLE a consensus of the other leaderboard
entries, destroying its independence as a peer entry, making its score depend on which
methods were fed to it, and instantiating the exact "consensus of methods" practice this
project's opening argument criticises.

---

## D3 · SCDC is not given the project's cell-size factors

**Severity: low, and entangled with D1.**

`R/run_scdc.R` does not pass `ct.cell.size`, so SCDC computes library sizes from the data
rather than using `reference_frozen/cell_size_factors.csv`. Whether that is wrong depends
on D1's resolution: if SCDC's output is taken as an mRNA share and corrected centrally,
not passing `ct.cell.size` is right. Resolve D1 first.

**Verified not a defect:** `truep` is never passed to SCDC. Passing it would leak the
answer into the estimate on the benchmark arm, where true proportions are known.

---

## D4 · The Python fallback is unbounded

**Severity: medium, operational.**

`r_bridge.timeout_for` bounds the R path. The Python reimplementation that replaces it on
timeout is bounded by nothing. In the 2026-09-10 run DWLS took **28,387 s** — 2,400 s of
correctly-killed R plus roughly 26,000 s of unbounded substitute — against 3,471 s for
identical work four days earlier.

The proximate cause that run was memory pressure from concurrent work, not the algorithm.
The structural problem is real regardless: a method stopped for being slow can spend eight
hours in its replacement.

**Partly addressed:** wall-clock per method now reaches `implementation_report.json` with
a note when a method exceeds the budget its R path was held to, so it is at least visible.

**Not addressed:** the fallback is still unbounded. Capping it is not obviously right — a
truncated method has no result at all, which is worse than a slow one — so this needs a
decision rather than a patch.

---

## D5 · The synthetic yardstick cannot reward platform correction

**Severity: interpretive, not a bug. Belongs in the paper.**

The benchmark arm builds pseudobulk from the GBmap atlas and scores methods against it.
The reference is built from that same atlas. So the synthetic mixtures have **no
cross-platform gap** — the thing CIBERSORTx S-mode exists to correct.

Observed in the confirmatory run: S-mode scores `mae_primary` 0.0693 against B-mode's
0.0608, i.e. **worse on the synthetic yardstick**. On real tissue it is the opposite —
S-mode puts myeloid at 17.5% against B-mode's 7.7%, and quanTIseq, independently and from
a different signature, says 15.5%.

Both can be true, and the reason matters: donor-held-out pseudobulk tests generalisation
across **donors** but not across **platforms**, so it systematically undervalues a method
whose contribution is platform correction. Any claim that S-mode is "less accurate" has to
be qualified by what the yardstick can and cannot see.

This is a limitation of the accuracy arm and should be stated as one.

---

## D6 · ACS cannot distinguish methods that disagree about composition by 20 points

**Not a bug. A property of the design, and the most important limitation to state in the
paper.** Found 2026-09-12 while auditing the ties.

Four methods score **exactly** ACS 0.9846 — `nnls`, `svr`, `elastic_net`, `epic`. Their
estimates are not close:

| method | mean Tumor | mean Macrophage | mean Endothelial |
|---|---|---|---|
| svr | **0.590** | 0.078 | 0.142 |
| nnls | 0.574 | 0.057 | 0.147 |
| elastic_net | 0.494 | 0.070 | 0.160 |
| epic | **0.395** | 0.073 | 0.157 |

Mean absolute difference per sample, svr vs epic: **Tumor 0.206, Oligodendrocyte 0.125,
Astrocyte 0.106.** Largest single-sample difference across the four: **0.72**.

So the anatomic test declares four methods equally good while they disagree about what
fraction of the tissue is cancer by twenty percentage points — a relative difference of
about 50%.

### Why this happens, and why it is not fixable by tuning

ACS scores **order**, not magnitude. Every constraint asks "is X higher in structure A
than in structure B", and a method can get every ordering right while being badly
calibrated in absolute terms. Two methods that both put endothelium highest in MVP score
identically whether one says 14% and the other says 30%.

This is inherent to any ordinal criterion, and it is the honest cost of the thing that
makes the design work at all: orderings are knowable from anatomy without ground truth,
and magnitudes are not.

### What it means for the study's claim

The registered question is whether anatomic concordance can be used **to choose a
deconvolution method**. This bounds the answer:

- ACS **can** separate methods that get orderings wrong from methods that get them right —
  the controls at 0.138 and 0.400 against a 0.969 median establish that, and it is a real
  result.
- ACS **cannot** choose between methods that order correctly but differ substantially in
  composition. On this cohort that is six of fourteen methods.

A paper that says "anatomy identifies the best method" would be wrong. One that says
"anatomy separates methods that reproduce known biology from methods that do not, and
cannot rank within that set" is supported.

### This compounds the resolution ceiling, it is not the same as it

The 66-value granularity (weighted denominator 65) explains why *scores* collide. This is
the stronger statement: even where scores collide, the underlying *estimates* differ by
amounts that would change any biological conclusion drawn from them. Reporting the tie
groups is necessary but not sufficient — the composition spread within a tie group should
be reported too.

**Suggested for the paper:** a figure showing mean composition per method within the
0.9846 tie group. It makes the limitation visible in one image and is more honest than a
footnote.

---

## D7 · Degenerate per-sample solves are neither detected nor disclosed

**Severity: medium. A real defect, with a bounded and measured effect.** Found 2026-09-12.

Ten estimates across three methods put **100% of a sample into a single cell type** —
always Tumor:

| method | samples at 100% Tumor |
|---|---|
| DWLS | 6 |
| BayesPrism | 2 |
| SCDC / SCDC ENSEMBLE | 2 |

A glioblastoma block containing no macrophages, no endothelium and no oligodendrocytes is
not a composition; it is a solver that failed and returned a corner of the simplex. The
pipeline currently treats these as valid estimates. Nothing flags them, nothing counts
them, and they enter the (tumour, structure) means that constraints are scored on.

This sits against the project's own standard — *"no failed folds, no silently dropped
samples, no method quietly missing from a ranking"*. These samples are not dropped, but
they are silently **wrong**, which the standard does not currently name.

### Where they land, and why that is the interesting part

**Seven of the eight distinct degenerate samples are PAN** — pseudopalisading cells around
necrosis. That is not random. PAN is necrotic tissue: degraded RNA, low library
complexity, and a composition genuinely dominated by tumour cells. It is the hardest
structure in the cohort to deconvolve, and it is also **the structure C5 depends on**
(Macrophage/microglia: PAN > LE).

### The effect is real but dose-dependent — measured, not assumed

| tumour | method | degenerate PAN blocks | C5 outcome |
|---|---|---|---|
| 705803 | DWLS | **2 of 3** | **violated** |
| 705803 | SCDC | 1 of 3 | satisfied |
| 703393 | SCDC | 1 of 3 | satisfied |
| 705803 | BayesPrism | 0 | violated — a different cause |

So a degenerate block corrupts the structure mean only when it is a **majority** of that
structure's blocks, because the mean over the remaining good blocks survives one bad one.
And C5 also fails for reasons unrelated to degeneracy, so the two must not be conflated.

**What is established:** degenerate solves exist, concentrate in PAN, are undetected, and
can flip a constraint when they dominate a structure's blocks.
**What is not established:** that they changed the leaderboard ordering. On this cohort
the one clear case (DWLS, 705803) affects a method already last.

### The fix

Detect and disclose, do not silently drop. A per-sample estimate whose maximum component
exceeds some threshold — 0.99, or exactly 1.0 — should be counted and reported per method
in `implementation_report.json`, the same way wall-clock and fallbacks are. Whether such
samples should also be *excluded* from the structure mean is a study-design decision and
must not be made quietly: excluding them would change ACS, and the exclusion rule would
have to be declared.

**Do not** fix this by clipping or smoothing the estimates. The degenerate value is the
method's actual output and the honest record of what it produced.


---

## D8 · The ranking below the top rests largely on one constraint, and that constraint sits on the roster's weakest column

**Status: MEASURED and DISCLOSED. Not a code defect — a limit on how finely the
leaderboard can be read. Recorded 2026-09-12.**

`scripts/constraint_sensitivity.py` recomputes ACS with each constraint and each tumour
left out in turn, from the archived (constraint x tumour) satisfaction matrix. The matrix
reproduces all 17 published ACS values to 1e-16 before any variant is reported, so this is
arithmetic on the recorded result, not a re-run that could differ for other reasons.

**The headline is robust.** MuSiC is top under every single-constraint exclusion and every
single-tumour exclusion. No one claim and no one tumour is responsible for it.

**The ordering below the top is not.**

| excluded | rho vs the full ranking | max rank move |
|---|---|---|
| C6 — Macrophage_Microglia: MVP > CT | **0.678** | **6 positions** |
| C7 — Tumor: LE < IT < CT | 0.929 | 5 |
| C5 — Macrophage_Microglia: PAN > LE | 0.973 | 2 |
| C1, C2, C3, C4 | 1.000 | 0 |
| worst single tumour (703393) | 0.926 | 5 |

Dropping C6 alone moves 13 of the 14 ranked methods and takes rank agreement to 0.678. So
the separation among the middle-ranked methods is substantially one constraint's work.

### Why this is worse than it first looks

C6 and C5 are the two Macrophage_Microglia constraints, and that column is the one the
reference is least able to support. `RESULTS.md` §1 records it independently: the roster
drops Mono (14,215 cells), DC (3,961) and RG (2,807) — 7.0% of the atlas — and the
dominant dropped signal is myeloid, whose nearest retained column is
Macrophage_Microglia. The absorbed signal is therefore concentrated on exactly the two
constraints that carry most of the ranking.

These two facts were recorded separately and neither pointed at the other. Together they
say: **the part of the leaderboard that discriminates most is the part resting on the
roster's most heavily loaded column.**

### What follows, and what does not

**Does not follow:** that the constraint set should change. C6 stays. The protocol is
explicit that constraints are never edited in response to scores, and a constraint that
carries the ranking is the last one it would be honest to remove after the fact. C3
likewise stays despite contributing nothing to the ordering — it earns its weight against
the controls, which satisfy it only 56% of the time.

**Does follow:**

1. **Report rank differences among the middle of the leaderboard as weak.** The tie
   structure already says this — ACS takes 66 distinct values on a weighted denominator of
   65, and D6 shows it cannot separate methods differing by 20 composition points. This is
   the third independent reason for the same caution.
2. **The second tissue must not reuse a roster with the same weakness**, or the same
   column carries the ranking twice and the replication is not independent of it. See
   `docs/SECOND_TISSUE.md`.
3. **A reference that retains the myeloid subtypes would test this directly.** If C6's
   discriminating power survives un-collapsing Mono and DC, the concern is answered; if it
   does not, the finding is about the roster and must be stated that way.

### What is established, and what is not

**Established:** the top method is robust to any single exclusion; C6 carries most of the
ordering below it; C6 and C5 sit on the column absorbing 7.0% of the atlas.
**Not established:** that the myeloid absorption *causes* C6's discriminating power. The
two are linked by the column they share, and that is a reason to test it, not a result.
