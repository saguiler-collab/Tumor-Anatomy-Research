# Open defects

Things known to be wrong, written down so they survive a context compaction and so
nobody has to rediscover them. Each entry states what was **verified** versus what is
**suspected**, because the difference decides whether a fix is safe.

Closed defects live in git history and in `docs/METHODS.md`; this file is only what is
still open — plus entries whose *investigation* concluded but whose reasoning is worth
keeping, marked RESOLVED at the top.

| | status |
|---|---|
| **D1** cell-size correction | **SUPERSEDED BY D12.** Package conventions measured for all five; Bisque alone returns cell share. But the conversion itself is the identity on this reference, so nothing was corrected twice — or once |
| **D2** SCDC ENSEMBLE degenerate | open, disclosed in the leaderboard |
| **D3** SCDC `ct.cell.size` | RESOLVED with D1 — not passing it is correct; SCDC returns full-space mRNA share |
| **D4** unbounded Python fallback | open, operational |
| **D5** yardstick cannot reward platform correction | open, design limit |
| **D6** ACS cannot separate 20-point composition differences | open, design limit |
| **D7** degenerate per-sample solves undetected | open |
| **D8** one constraint carries the ranking below the top | measured, disclosed |
| **D9** false pre-registration violation published | FIXED; residual weakness disclosed |
| **D10** re-measurements used a different gene space | mechanism FIXED; the two measurements must be redone |
| **D11** EPIC returns the within-subset mRNA share, not the full-space one | **OPEN** — but see D12: the conversion it is applied to is the identity, so the practical effect on the published numbers is nil |
| **D12** the cell-size correction has never been applied — it is the identity | **OPEN, high** — supersedes most of D1; no number changes, but the registration says otherwise |
| **D13** EPIC runs with `refProfiles.var` unset, so its gene weighting is off; its output is mislabelled as a cell fraction | **OPEN, high** — EPIC is joint-2nd at 0.9846 |

---

## D1 · Cell-size correction — investigated as a DOUBLE correction, measured to be a SUBSTITUTED one

> **STATUS, 2026-09-12 (third revision today — read all of it).** This entry has now been
> wrong in both directions, and the sequence matters more than any one claim in it:
>
> 1. It first said the double correction **cannot change ACS**. Wrong; retracted.
> 2. It then said the double correction **happens and every leaderboard number is
>    provisional**. Asserted from prose with no artefact.
> 3. A measured probe then said the packages return mRNA share, so the correction is
>    applied **once** and D1 is resolved. Withdrawn — the probe exported its full gene
>    space, and production exports a marker subset.
> 4. The probe was rebuilt with an asymmetric subset matching the real atlas's structure and
>    **re-run. The "applied once" conclusion holds for MuSiC, SCDC and the least-squares
>    family**, for a reason the earlier version had not identified: MuSiC's own conversion
>    corrects the *subset bias*, not cell size, and what it leaves is exactly what this
>    project's conversion expects.
>
> **Current state: RESOLVED for MuSiC, SCDC and least squares. NOT a double correction.**
> Unmeasured for Bisque and BayesPrism, which the probe cannot read. **And it turned up a
> separate, real defect in EPIC — see D11.** Every step above is kept because the reasoning
> was wrong twice and the record of how is worth more than the conclusion.

**Severity as first recorded: high.** It affects the accuracy numbers that the primary result is computed
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

### MEASURED 2026-09-12 — the double correction does not happen in production

**The premise of everything above is refuted for three of the five packages.** The probe is
now a script (`scripts/verify_cell_size_semantics.py`) rather than prose, it was run against
the genuine R packages, and it was run twice: once with the cell export production actually
writes, and once with raw per-cell library sizes, which is what each package's published
usage assumes. Only the export differs between the two — the reference profile is built the
production way in both — so a package that moves between them moved because of its own
internal cell-size estimate and nothing else.

| package | production export (1e6 CPM/cell) | raw export | moves? |
|---|---|---|---|
| **MuSiC** | 0.7501 — **mRNA share** | 0.5001 — CELL share | **yes** |
| **SCDC** | 0.7501 — **mRNA share** | 0.5336 — partial | **yes** |
| **EPIC** | 0.7501 — **mRNA share** | 0.7501 — mRNA share | no |
| Bisque | UNREADABLE | UNREADABLE | — |

Artefacts: `results/cell_size_semantics_normalized.json`,
`results/cell_size_semantics_raw.json`, `results/cell_size_semantics_comparison.json`.
The probe's true mRNA fraction is 0.7500 and its true cell fraction 0.5000, so 0.7501 is
the mRNA answer to within the solver's own error.

**The mechanism.** MuSiC and SCDC *do* convert cell size when they can estimate it —
`music_basis` computes `M.S`, the mean library size per cell type, and `SCDC_basis` derives
one the same way. Production hands them cells normalised to 1e6 each, so every type's mean
library size is **identical by construction**, their estimate comes out flat, and their own
conversion becomes the identity. They return mRNA share. `reference.py` states this exact
hazard in a comment about why `cell_totals` must be captured before normalisation; what was
missed is that it applies to the R packages' internal estimates too, not only to ours.

EPIC does not move because it never reads the cells — it consumes the signature matrix, and
its `mRNA_cell` argument is left at the default, which supplies no factor for cell-type
names outside its own built-in roster.

**So the project's conversion is the FIRST and ONLY one for these three, and it lands on
the truth.** Taking MuSiC's 0.7501, dividing by the cell-size factors (3000, 1000) and
renormalising gives **0.5001** against a true cell fraction of **0.5000**.

### WITHDRAWN: the resolution above was premature — measured 2026-09-12, same day

**The probe exported its full gene space. Production exports a marker subset. That
difference decides the answer, and it invalidates the inference above.**

`run_r_method` calls `export_for_genes(ref_name, list(data.bulk.index))`, so the R packages
receive only the genes in the bulk — 657 of roughly 27,625 in the cell source. Cells are
normalised to 1e6 **across all genes**; they are emphatically not equal across a subset of
them, because cell types differ in how much of their transcriptome falls inside a marker
panel.

Measured on the real atlas export that is on disk (`sc_counts_gbmap_c3d46536…`, 1,591 genes
x 15,311 cells, CPM-normalised then subset — the production construction):

| cell type | mean library size over the exported subset | relative |
|---|---|---|
| Oligodendrocyte | 65,003 | 2.17 |
| Astrocyte | 61,263 | 2.04 |
| Macrophage_Microglia | 60,253 | 2.01 |
| NK_cell | 51,646 | 1.72 |
| Endothelial | 51,140 | 1.70 |
| T_cell | 43,705 | 1.46 |
| B_cell | 40,871 | 1.36 |
| **Tumor** | **29,999** | **1.00** |

**Spread 2.17x — not flat.** `music_basis` computes `M.S` from exactly this quantity, so
MuSiC's own cell-size estimate in production is *not* the identity. It is a real conversion,
applied with factors that are an artefact of which genes the marker panel happens to
contain, and it is then followed by this project's central conversion.

So the probe's clean 0.7501 shows only that MuSiC returns mRNA share **when its input gives
it no per-type variation to find**. Production gives it 2.17x of variation. The claim that
the correction is applied once does not follow, and the earlier "the leaderboard is not
provisional for this reason" is withdrawn with it.

**Why this was caught.** The probe was built symmetric — every cell type owning an equal
block of genes with identical concentration — so restricting it to a subset preserved
flatness by construction. Real cell types are not symmetric: a tumour cell puts 30,000 of
its million into the marker space while an oligodendrocyte puts 65,000. A synthetic probe
that mirrors production in the variable under test but not in the structure that modulates
it will produce a confident wrong answer, which is what happened here.

### The experiment was run — MEASURED, biased subset, 2026-09-12

`--subset-mode biased` keeps a third of the BIG type's marker block and all of every other
block, giving per-type library sizes over the export that span **1.92x with the BIG type
lowest** — the structure measured on the real atlas, where Tumor sits at 29,999 against
Oligodendrocyte's 65,003 (2.17x). The probe's truths, computed from the cells themselves:

| quantity | value |
|---|---|
| true CELL fraction of the big type | **0.5000** |
| true mRNA share over the FULL gene space | **0.7500** |
| true mRNA share over the EXPORTED SUBSET | **0.6094** |

Two different right answers are therefore available, and which one a package returns is the
whole question.

| package | measured | which estimand | after the central cell-size conversion | error vs 0.5000 |
|---|---|---|---|---|
| **MuSiC** | 0.7500 | full-space mRNA share | **0.5000** | **+0.0000** |
| NNLS (control) | 0.7506 | full-space mRNA share | 0.5008 | +0.0008 |
| **SCDC** | 0.7469 | full-space mRNA share | 0.4959 | −0.0041 |
| **EPIC** | **0.6104** | **within-SUBSET mRNA share** | **0.3431** | **−0.1569** |

**MuSiC and SCDC come out correct, and the double-correction concern is refuted under the
condition that actually matters.** MuSiC's `M.S` is genuinely non-flat here — that part of
the withdrawal was right — but what it corrects is the **subset bias**, not cell size. It
*cannot* correct cell size: normalising each cell to 1e6 destroyed that information before
MuSiC ever saw the cells. The residual after its own conversion is exactly the
full-transcriptome mRNA share, which is precisely the quantity this project's central
conversion expects. The two compose correctly, and they compose correctly *because* of the
normalisation rather than in spite of it.

So the answer to "is the cell-size correction applied twice?" is **no, for MuSiC, SCDC and
the least-squares family** — established now on the production-like subset, not on the
symmetric full-space probe whose result was rightly withdrawn.

**EPIC is a different matter, and it is a new defect. See D11.**

Bisque and BayesPrism remain unreadable on this probe in every mode (≥31% of their mass
lands on the roster types absent from the mixture), so their convention is still unmeasured
and must not be asserted.

### ALL FIVE MEASURED — the table is complete, 2026-09-13

The probe could not read Bisque or BayesPrism because the two-type mixture left six roster
types absent, and those two put 31-75% of their mass on absent columns. `--mixture all_types`
removes the cause: all eight types appear at equal cell counts, BIG still carrying 3x the
mRNA, so there are no absent types to leak into.

| package | returned (BIG:SMALL) | off-pair mass | convention | this project's conversion is |
|---|---|---|---|---|
| **MuSiC** | 0.7502 | 0.600 | mRNA share | the FIRST — correct |
| **SCDC** | 0.7503 | 0.600 | mRNA share | the FIRST — correct |
| **EPIC** | 0.7502 | 0.600 | mRNA share | the FIRST — correct *(but see D11)* |
| **BayesPrism** | 0.7502 | 0.600 | mRNA share | the FIRST — correct |
| **Bisque** | **0.5000** | **0.750** | **CELL share** | **a SECOND correction — defect** |

**Each verdict is confirmed twice over, by two independent quantities.** The pair ratio says
what convention a method used; the *off-pair* mass says the same thing without reference to
the pair, because a method reporting mRNA share must leave 1-(3+1)/(3+7) = **0.600** outside
the pair, while one reporting cell share must leave 1-2/8 = **0.750**. Every method's two
readings agree. A coding error in the classifier could move the first number; it could not
move both into agreement.

Artefact: `results/cell_size_semantics_normalized_all_types.json`.

**So the original prose table was right about Bisque and wrong about MuSiC, SCDC and EPIC.**
Bisque is the one genuine double correction in the panel, and it is the only one. The
sequence of revisions in this entry was not wasted: the answer is package-specific, and
nothing short of measuring each package on a probe that mirrors production would have got it
right.

### Bisque: what the double correction means, and what must NOT be done about it

Bisque's `ReferenceBasedDecomposition` returns cell proportions, and `to_cell_fractions` then
divides by cell size again. Its leaderboard row (ACS 0.9231) and its benchmark numbers carry
that.

This matters more than its rank suggests, because Bisque is the panel's **most reliable
method on the tumour compartment** by two independent measures: the smallest tumour bias
(-0.0256 against MuSiC's -0.0757) and the narrowest conformal interval (+/-0.489 against
MuSiC's +/-0.557). A correction applied twice to the method that is otherwise best calibrated
on the clinically decisive quantity is worth fixing properly rather than quickly.

**The fix is to skip the central conversion for Bisque**, exactly as it is already skipped for
quanTIseq, and to declare it the same way. That is a mechanism change, not a tuning: the rule
is "apply the conversion once", and for Bisque applying it centrally is applying it twice.

**What must not be done:** do not adopt the fix *because* it may improve Bisque's numbers.
The decision rests on the measurement above, which is independent of any ACS or benchmark
score — the probe is synthetic and its truth is known by construction. Record the before and
after for every affected number, as the D1 fix section already requires.

### What this means for the leaderboard

**The "15 of 16 methods change ACS, 7 of 14 ranks move" measurement above describes a
counterfactual, not a bias in the published numbers.** It was produced by applying the
correction *a second time* to each archived estimate table, which is the right experiment
for the hypothesis that the correction had already been applied once by the package. That
hypothesis is now measured to be false for MuSiC, SCDC and EPIC. Applying the correction
twice does change ACS — that finding stands and the retraction above of the
"cannot change ACS" claim stands with it — but **nothing in the pipeline applies it twice.**

**The leaderboard is therefore not provisional for this reason.** The earlier statement that
"every number in the leaderboard should be treated as provisional until the fix lands" is
withdrawn. There is no fix to land for these three.

### What remains genuinely open

1. **Bisque is unresolved.** The probe cannot read it: Bisque assigns **75%** of its mass to
   the six roster types absent from the mixture, in both scalings, so the ratio between the
   two present types measures nothing about a cell-size convention. Its documented
   behaviour (`ReferenceBasedDecomposition` returning cell proportions) suggests cell share,
   but this project has not measured it and must not assert it. A probe that fits Bisque
   would need a roster with no absent types.
2. **BayesPrism is unresolved under production scaling** for the same reason — 31% leakage.
   Its Python reimplementation measures 0.7500, mRNA share, cleanly.
3. **The earlier prose numbers are unexplained but consistent.** MuSiC 0.500 matches the
   raw-export measurement almost exactly; SCDC 0.653 sits between this probe's raw 0.5336
   and mRNA 0.7501. The likeliest reading is that the unsaved original probe used raw cells
   — the same probe-versus-production mismatch that this script's own first version fell
   into and that its comments now record. Not established, and it does not need to be: the
   production measurement is what governs.

### The real defect this uncovered — UNDECLARED DEVIATION

Not a double correction. A **substituted** one.

Production silently neutralises MuSiC's and SCDC's own cell-size machinery by normalising
every cell to a common library size before export, and then supplies this project's factors
instead. That is arguably the *better* choice — one set of factors, computed one way,
applied identically to every method, which is what equal footing demands, and it is clearly
what the design intends, since `cell_totals` is deliberately captured before normalisation
precisely so the information normalisation destroys can be reinstated centrally.

But it is a deviation from each package's published behaviour, it is currently **not
declared** in `method_configs.json` alongside the cibersortx, dwls and quantiseq
deviations, and a reviewer who knows MuSiC would reasonably expect `music_prop` to be doing
its own cell-size estimation. It must be declared, with this measurement as the evidence.

That is the remaining work on D1, and it is documentation rather than a re-run.

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
rather than using `reference_frozen/cell_size_factors.csv`.

**STILL OPEN — a resolution was written and withdrawn the same day.** It said D1's
measurement settled this, because SCDC derives `ct.cell.size` from the cells it is given and
production normalises them, so what it derives would be flat. That reasoning holds only for
an export spanning the full gene space. Production exports a marker subset, over which
per-type library sizes are measured to span 2.17x on the real atlas, so what SCDC derives is
**not** flat and its conversion is not the identity.

What is measured: SCDC returns 0.7501 on a full-gene-space export and 0.5336 on a raw one,
so it does convert when its input gives it per-type variation to find. Production gives it
some. Whether passing `ct.cell.size` explicitly would be the second conversion or a
replacement for a partial one is exactly what D1's remaining experiment decides. Resolve D1
first — genuinely, this time.

**What remains is the declaration, not the code.** The substitution — our factors instead of
the package's own — is now listed in `DECLARED_DEVIATIONS` for `scdc` and `scdc_ensemble`
and appears in `method_configs.json`.

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


---

## D9 · The pre-registration ordering gate published a violation that never happened, and cannot be verified from a copy

**Status: the false verdict is FIXED. The gate's underlying weakness is DISCLOSED and
partly mitigated. Recorded 2026-09-12.**

### What was published, and why it was wrong

`RESULTS.md` §1a — the integrity section, above the leaderboard — carried this:

> REGISTERED, but 25 result file(s) are OLDER than the registration. The protocol's gate
> is that the registration precedes every result. Those files were produced before the
> constraints were registered and must be regenerated or labelled superseded.

**The true count was zero.** Measured: 45 result files, none with an mtime before
`2026-09-10T03:32:03Z`.

The cause was ordering, not arithmetic. `run_anatomic.run()` called
`registration.status()` partway through the run, and the gate works by comparing each
result file's mtime against the registration time. At that moment the run had not yet
overwritten most of its artefacts, so the mtimes compared were the **previous** run's. The
check therefore reported the previous run's files as violations on every run, and nothing
rewrote the verdict afterwards. Twenty-one result files ended up newer than the file
asserting they were too old.

A protocol violation asserted that did not happen is as damaging as one missed. A reader
cannot tell which kind of error they are looking at, and the natural remedy the text
itself suggests — "regenerate or label superseded" — would have thrown away sound results.

### Fixed

1. `registration_status.json` is now written **last**, after every other artefact, which
   is the only point at which every file carries the current run's mtime. The provisional
   copy written mid-run is documented as provisional at the call site.
2. `summarize_results.py` refuses to reprint an ordering verdict whose own file is older
   than the results it claims to have judged, and prints a STALE warning instead.
3. The verdict was recomputed: **REGISTERED, and every result file postdates the
   registration.** `release/` was updated to match.
4. `tests/test_registration_ordering_gate.py` pins all three, including a test that the
   repository's own tree is not carrying a stale verdict.

### The weakness that remains — and it matters for a reviewer

**mtimes do not survive a copy.** `git clone`, `rsync`, and unzipping an OSF download all
stamp every file with the moment of the copy. So for any tree a reviewer obtains, every
result file postdates the registration by construction, and the gate passes **vacuously**.
It cannot be the evidence that this project's results postdate its registration; it can
only catch the mistake locally, on the machine that produced them.

Partly mitigated: `anatomic_report.json` now records `run_started_utc` and
`run_finished_utc` inside the artefact, which travel with the file and can be checked
later. The gate does not yet read them.

**What actually establishes the ordering, and what the paper should cite:**

- the git history — the constraint file's commit and hash precede the result commits, and
  git timestamps are content-addressed rather than filesystem metadata;
- the OSF registration's own immutable timestamp, `2026-09-10T03:32:03Z`, which this
  project cannot alter;
- the constraint freeze hash `2d1fb47c…`, recorded in every artefact and matching the
  registered payload.

Those three are durable. The mtime gate is a local smoke test and should be described as
one, never as the provenance claim.

### The fix that would close it

Have the gate prefer `run_started_utc` from the run's own report when present, fall back
to mtimes with the vacuity stated in the verdict text, and refuse to report PASS from
mtimes alone on a tree whose files all share a copy timestamp.


---

## D10 · Both single-method re-measurements ran on a different gene space than the leaderboard they were compared against

**Status: the false comparisons are WITHDRAWN and the mechanism is FIXED. The
re-measurements themselves must be redone. Found 2026-09-12.**

### What happened

| | ACS | genes |
|---|---|---|
| every one of the 16 methods on the ACS leaderboard | — | **657** |
| genuine `R:DWLS` re-measurement | 0.7077 | **1,591** |
| genuine `R:BayesPrism` re-measurement | 0.8923 | **1,591** |

`scripts/remeasure_method.py` built its gene space with `frozen_gene_space()`. That is the
space `run_all.py` uses **only in its fallback branch**, when no cell-level reference is
available. The real path uses the benchmark's marker subset — `bench["genes"]` intersected
with the bulk — which is 657 genes. The two differ by a factor of 2.4.

Both re-measurements were then reported beside the leaderboard as though the only
difference were the implementation. The DWLS report printed "reimplementation in the
confirmatory run: ACS 0.7385" directly beneath its own 0.7077, and three documents drew
conclusions from the gap:

- `README.md`: the genuine package scores "*lower* than the reimplementation's 0.738, and
  still last", and "Running the real package shows the disagreement is genuine rather than
  an artefact of substituted software".
- `docs/RELATED_WORK.md`: "the disagreement with Avila Cobos survives the obvious
  explanation".
- `docs/ROAD_TO_PAPER.md` Tier 0.3, marked **DONE**.

**None of it follows.** All three are withdrawn in place.

### Why nothing caught it

Only the gene **count** was ever persisted — `n_signature_genes: 657` in
`method_selection_decision.json`, `n_genes: "657"` in the equal-footing certificate. The
gene *list* existed only in memory during a run. So a re-measurement had no artefact to
read, `frozen_gene_space()` looked like a reasonable way to rebuild "the" gene space, and
the mismatch was invisible to every check in the project. The equal-footing machinery
compares methods *within* a run and had no reason to look across runs.

This is the same shape as D9: a check that was structurally unable to see the thing it
was supposed to guarantee.

### Fixed

1. `run_all.py` writes `results/benchmark/signature_genes.json` — the full ordered gene
   list, its count, and an order-sensitive `sha256` (`config.sha256_strings`, added for
   this; order matters because two methods handed the same genes in a different column
   order have not received identical inputs).
2. `remeasure_method.py` reads that file and **aborts** if it is absent, if any recorded
   gene is missing from the bulk, or if the list does not reproduce the recorded hash. It
   never substitutes a gene space. The check runs before the multi-GB atlas load, so it
   fails in about two seconds rather than after an hour of CPU.
3. `results/dwls_remeasured.json` and `results/bayesprism_remeasured.json` carry a
   `COMPARABILITY: NOT COMPARABLE` header stating what each still establishes and what it
   does not.
4. The withdrawals are written into README, `RELATED_WORK.md` and Tier 0.3, which is
   reopened.

### What still stands from the re-measurements

Both published packages **run to completion** on this cohort, and both need roughly twice
the pipeline's 2,400 s budget — DWLS 3,873 s, BayesPrism 4,689 s. Those are real
measurements and they explain the fallbacks. Nothing about the *scores* is usable yet.

### What closes it

One full `run_all.py` (which writes the gene list), then re-run both re-measurements. The
two must not run concurrently on this machine: three BayesPrism Gibbs workers at ~900 MB
each already drive it into swap. Only then can the question the re-measurement was built
to answer — does genuine DWLS beat its reimplementation, and does the Avila Cobos
disagreement survive — be answered at all.


---

## D11 · EPIC returns a different estimand than every other method, and the cell-size conversion is applied to the wrong quantity

**Severity: high for EPIC's own numbers. EPIC is joint-2nd on the ACS leaderboard at 0.9846
and joint-first in one of the cell-size sensitivity orderings. Found 2026-09-12, measured.**

### What was measured

On the biased-subset probe (see D1 — per-type export library sizes spanning 1.92x, matching
the real atlas's 2.17x), the probe has two distinct correct answers available:

- true mRNA share of the big type over the **full** gene space: **0.7500**
- true mRNA share over the **exported subset**: **0.6094**

| package | measured | which one it returned |
|---|---|---|
| MuSiC | 0.7500 | full space |
| NNLS | 0.7506 | full space |
| SCDC | 0.7469 | full space |
| **EPIC** | **0.6104** | **the subset** |

Every other method returns the full-space mRNA share. EPIC returns the within-subset share.

### Why that is a defect rather than a preference

`to_cell_fractions` divides by cell-size factors derived from **full-transcriptome** library
sizes (`cell_totals`, captured in `build_from_h5ad` before normalisation). That conversion is
correct for a full-space mRNA share and wrong for a subset share, because the two differ by
each type's marker-space concentration.

Measured consequence on the probe, after the central conversion:

| package | after conversion | truth | error |
|---|---|---|---|
| MuSiC | 0.5000 | 0.5000 | +0.0000 |
| NNLS | 0.5008 | 0.5000 | +0.0008 |
| SCDC | 0.4959 | 0.5000 | −0.0041 |
| **EPIC** | **0.3431** | 0.5000 | **−0.1569** |

A 0.157 absolute error in the dominant cell type's fraction, where every other method is
within 0.005.

### The mechanism — ISOLATED AND CONFIRMED

`R/run_epic.R` calls `EPIC(..., scaleExprs = TRUE)`, EPIC's published default. `scaleExprs`
renormalises the reference profiles and the bulk onto their shared gene set, discarding the
column scaling that carries full-space information — the scaling `ref.subset_genes()`
preserves for every method that solves against the profile directly.

The distinguishing test was run, on the same biased-subset probe:

| `scaleExprs` | EPIC returns | which estimand |
|---|---|---|
| `TRUE` (published default, what production uses) | **0.6104** | within-subset mRNA share (true 0.6094) |
| `FALSE` (diagnostic only) | **0.7504** | full-space mRNA share (true 0.7500) |

**Confirmed.** `scaleExprs = TRUE` is the cause, and with it off EPIC agrees with MuSiC,
SCDC and NNLS to three decimals.

`R/run_epic.R` carries an `IVYGAP_EPIC_SCALE_EXPRS` environment hook solely so this test is
reproducible. **It defaults to TRUE and production behaviour is unchanged** — verified by
re-running with no variable set and with it set to `TRUE`, both giving 0.6104. The hook is
marked in the script as a diagnostic that must never be set for a real run.

### What must NOT be done

**Do not change `scaleExprs` to make EPIC's number look better.** It is EPIC's published
default, the protocol's fairness rule is defaults with no per-method tuning, and EPIC's ACS
is already known. Changing a parameter after seeing a score is the retune the protocol
forbids, and it would be indistinguishable from tuning EPIC up the leaderboard.

The legitimate routes are: declare the deviation in estimand and correct the conversion for
EPIC specifically (a documented, pre-stated rule about which estimand each method returns),
or report EPIC's row with the mis-scaling disclosed. Either requires knowing the mechanism
first.

### What this does and does not touch

**Does not touch** the other 15 methods. MuSiC, SCDC and the least-squares family are
measured correct to within 0.005 on this probe.

**Does touch** EPIC's accuracy numbers directly, and its ACS potentially — ACS is ordinal
within a cell type across structures, and a per-type scaling constant survives the division
but *not* the renormalisation that follows it, which is composition-dependent. That is the
same mechanism as D1's retracted-then-confirmed point about renormalisation, so an EPIC
mis-scaling can move ACS and cannot be assumed harmless.

### The size on real data is NOT established, and one line of evidence argues it is small

Two calculations disagree, and both belong here.

**An analytic propagation says the effect should be large.** EPIC's estimand is proportional
to `n_k · c_k · f_k`, where `f_k` is the share of type *k*'s transcriptome falling inside the
marker panel. The central conversion divides by `c_k`, leaving `n_k · f_k` instead of `n_k`,
so the residual distortion is exactly `f_k`. Measured on the real atlas export:

| cell type | `f_k` | relative to the mean |
|---|---|---|
| Tumor | 0.0300 | **0.594** |
| B_cell | 0.0409 | 0.810 |
| T_cell | 0.0437 | 0.866 |
| Endothelial | 0.0511 | 1.013 |
| NK_cell | 0.0516 | 1.023 |
| Macrophage_Microglia | 0.0603 | 1.193 |
| Astrocyte | 0.0613 | 1.213 |
| Oligodendrocyte | 0.0650 | 1.288 |

Propagating those through a plausible GBM composition gives Tumor 0.550 → 0.387
(**−0.163**) and Macrophage_Microglia 0.250 → 0.353 (**+0.103**), which matches the probe's
−0.157 closely.

**The pseudobulk benchmark argues against a distortion that large.** If EPIC's tumour
fraction were out by 0.16 on real mixtures, its error would be conspicuous. It is not:

| method | `mae_primary` |
|---|---|
| music | 0.0508 |
| elastic_net | 0.0551 |
| nnls | 0.0561 |
| bayesprism | 0.0585 |
| cibersortx | 0.0608 |
| svr | 0.0610 |
| **epic** | **0.0659** |
| cibersortx_smode | 0.0693 |

EPIC sits 7th of 15, modestly behind MuSiC, nowhere near where a 0.16 bias on the dominant
type would put it. **So the analytic propagation overstates the real-data effect**, and this
entry does not claim −0.163 for real data.

Candidate reasons the propagation is too crude, none tested: `withOtherCells = TRUE` gives
EPIC a sink column that can absorb part of the mis-scaling, and `write_proportions`
renormalises across the roster afterwards; the benchmark's own pseudobulk truth is built from
the same normalised cells, which may cancel part of the effect; the real 657-gene panel may
have a narrower spread than the 1,591-gene export these `f_k` come from; and EPIC's
constrained solve on a real, non-orthogonal signature does not reduce cleanly to "returns the
subset share" the way it does on the probe's block-diagonal one.

**What is established:** on a controlled probe EPIC returns the subset share where every
other method returns the full-space share, and `scaleExprs` is the cause.
**What is not:** how much that costs on Ivy GAP. The benchmark says less than the propagation
implies. Resolving it needs EPIC re-run on the real cohort with `scaleExprs = FALSE` and the
two compositions compared directly — which is a measurement, not an argument, and must be
recorded as a declared diagnostic rather than swapped into the leaderboard.


---

## D12 · The cell-size correction has never been applied. It is the identity on the reference every score was computed against.

**Severity: high, and it supersedes most of D1. Found 2026-09-14, measured.**

### What was measured

`results/bisque_remeasured.json` records the reference's `cell_size` vector for the first
time:

| cell type | cell_size |
|---|---|
| Tumor | 1,000,000.0000 |
| Macrophage_Microglia | 1,000,000.0000 |
| T_cell | 999,999.9375 |
| NK_cell | 1,000,000.0000 |
| B_cell | 1,000,000.0000 |
| Endothelial | 1,000,000.0000 |
| Oligodendrocyte | 1,000,000.0000 |
| Astrocyte | 1,000,000.0000 |

**Spread max/min = 1.0000.** `to_cell_fractions` divides by this vector and renormalises, so
it is **exactly the identity**. (The 999,999.9375 is float32 summation of values totalling
1e6 — it does not make the vector non-uniform in any meaningful sense.)

### The mechanism

`build_from_h5ad` does the right thing and says so in a comment at `reference.py:386`:

> Each cell's library size, captured BEFORE normalisation. Normalising first and measuring
> after would make every type's mean total 1e6 and turn the cell-size correction into a
> no-op.

It captures `raw_totals`, normalises each cell to 1e6, and passes `cell_totals=raw_totals`
into `build_reference`. That reference has real cell-size spread.

**`run_benchmark` then discards it.** At `run_benchmark.py:63-65` it rebuilds the reference
from the already-normalised matrix, training donors only, and does **not** pass
`cell_totals`:

```python
reference = build_reference(expression[train_cells], meta.loc[train_cells],
                            name=config.PRIMARY_REFERENCE)
```

`build_reference` falls back to `expression.sum(axis=0)` (`reference.py:152-153`), which is
1e6 for every cell because that matrix is the normalised one. So the reference handed to the
benchmark — and then, via `bench["reference"]`, to the anatomic leaderboard — carries a
uniform `cell_size`, and the conversion does nothing.

The comment warning against exactly this sits 230 lines away from the call that does it.

### Independent corroboration

Bisque was re-measured with the central conversion **skipped** (D1's fix). Its estimate table
came out **bit-identical** to the archived `ivygap_bisque.csv`, which was produced with the
conversion **applied** — max absolute difference 0.000000 across all 122 samples and 8 types,
and ACS 0.9231 with CI [0.8657, 0.9831] both ways.

Applying a transform and skipping it cannot give identical output unless the transform is the
identity. That is a second, independent proof, and it does not depend on reading any code.

### What this means

**It supersedes most of D1.** There is no double correction anywhere in the published
results, because there is no correction anywhere:

- Bisque's "confirmed double correction" is not double — it is zero. The fix committed for it
  is *correct in mechanism* and *inert in effect* on this reference, and the leaderboard row
  is unchanged.
- The "15 of 16 methods change ACS when the correction is applied a second time" experiment
  does **not** describe this pipeline. It must have used a cell-size vector with real spread
  (the frozen `cell_size_factors.csv` spans 737x), so it measured a counterfactual reference,
  not the one the leaderboard used. That experiment's numbers should not be cited as bearing
  on the published results.
- The `DECLARED_DEVIATIONS` entries added for MuSiC, SCDC, Bisque and BayesPrism say the
  harness "substitutes this project's cell-size factors" for the packages' own. **It
  substitutes nothing.** They must be corrected.

**It is a deviation from the registration**, which states the conversion is "applied once and
centrally". It is applied nowhere. That is a disclosure item, not a silent one.

**What it does NOT change.** Every leaderboard number, every benchmark number, ACS, rho, the
controls and the ISH validation all stand exactly as recorded — they were computed with an
identity transform in the chain, which is the same as without it. Nothing needs re-running to
*correct* anything. What needs deciding is whether the conversion *should* be applied, which
is a design question, not a bug fix.

### The decision this forces, and how to take it without breaking the protocol

Two defensible options, and the choice must be recorded before any number moves:

1. **Apply the conversion for real** — pass `cell_totals` through in `run_benchmark`, so the
   reference carries true per-type mRNA content. This changes every method's estimates and
   must be reported as a declared deviation with before/after per number. It also makes D1
   live again, and Bisque's fix load-bearing.
2. **Declare that this study reports mRNA proportions, not cell proportions** — which is what
   it has in fact been doing throughout — and remove the conversion rather than leave dead
   code implying otherwise. Cheaper, honest, and it changes no number; but it weakens the
   comparison with studies that report cell fractions.

**Option 1 must not be chosen because it might improve a score**, and option 2 must not be
chosen because it is less work. The measurement above is score-independent; the decision
should be too, and stated in the paper either way.


---

## D13 · EPIC runs with its defining feature disabled, and its output is mislabelled

**Severity: high. EPIC is joint-2nd on the ACS leaderboard at 0.9846. Found 2026-09-14 by
interrogating the installed package.**

Three separate breaks, all measured against `EPIC` as installed, all in `R/run_epic.R`.

### 1 · `refProfiles.var` is not supplied, so gene variance weighting is OFF

```
ref <- list(refProfiles = sig_mat, sigGenes = rownames(sig_mat))
```

EPIC emits a warning that this project has never recorded:

> `'refProfiles.var' not defined; using identical weights for all genes`

Weighting genes by their variability across reference samples is EPIC's published
contribution. Without it EPIC is a constrained least-squares solve with uniform weights.

**The variance matrix exists and is simply never passed.** `build_reference` computes it
(`reference.py:134-139`) and `ReferenceBundle` carries it as `sigma`. Nothing in `r_bridge` or
any script in `R/` mentions `sigma`, so it is never exported and EPIC never sees it.

This is the same shape as the project's own invariant about MuSiC — *"MuSiC without
cross-donor variance is NNLS"* — and it applies with equal force: **EPIC without
`refProfiles.var` is not EPIC.** By the invariant "never report a degenerate method under its
own name", the row should either supply the variance or be renamed and disclosed.

### 2 · `cellFractions` is taken as a cell fraction and is not one

`run_epic.R` reads `est$cellFractions`. EPIC converts mRNA share to cell share using
`mRNA_cell`, which this project leaves at `NULL`, so EPIC falls back to its shipped
`mRNA_cell_default`. That vector is keyed on EPIC's own roster:

```
Bcells 0.402  Macrophages 1.420  Monocytes 1.420  Neutrophils 0.130
NKcells 0.440  Tcells 0.395  ...  otherCells 0.400  default 0.400
```

**Not one of this project's eight roster names matches.** `Macrophage_Microglia` is not
`Macrophages`; `T_cell` is not `Tcells`; `B_cell` is not `Bcells`; `NK_cell` is not `NKcells`.
So EPIC applies `default = 0.400` to all eight, uniformly — and a uniform factor cancels in
the renormalisation.

Measured directly on the installed package with this project's roster names:

> `cellFractions` IDENTICAL to `mRNAProportions`, max difference **2.776e-17**.

So `cellFractions` here is an mRNA proportion wearing a cell-fraction label. Given **D12** —
the project's own conversion is also the identity — nothing downstream is inconsistent, but
the label is wrong and the paper must not call it a cell fraction.

This also explains D11 precisely: EPIC's divergence was never about `scaleExprs` alone.

### 3 · The `otherCells` compartment makes EPIC's estimand different, and its size is unrecorded

`withOtherCells = TRUE` gives EPIC a ninth column for signal the reference cannot explain,
which is defensible on this data — the roster drops 7.0% of the atlas. But the harness then
drops that column and renormalises across the eight roster types, so **EPIC's reported
fractions are conditional on the explained signal** while every other method's are not.

On a clean synthetic where the reference explains everything, `otherCells = 0.0000`, so the
effect is nil there. On the real cohort it is **not recorded at all**: `run_epic.R` writes the
sidecar next to its output inside a temporary directory that is deleted with the run, so no
artefact preserves it. The one number needed to judge how much this matters was thrown away.

### What must be done, and what must not

**Must be done, in this order:**

1. **Preserve the `otherCells` mass** as a first-class artefact. It is one line and it decides
   how serious break 3 is.
2. **Export `sigma` and pass it as `refProfiles.var`.** This restores EPIC to its published
   form. It will change EPIC's numbers, so it is a declared deviation with before/after
   published per number.
3. **Relabel.** EPIC's column is an mRNA proportion. Either pass a real `mRNA_cell` for this
   roster — which requires deciding what that vector should be, a design question — or state
   plainly that EPIC's row is an mRNA proportion.

**Must not be done:** none of these may be adopted because it improves EPIC's rank. Break 1
is a recording gap, break 2 is measured from a package warning, and break 3 is measured at
machine precision. All three are score-independent; the decisions must be too.

### MEASURED on the real cohort, 2026-09-14 — one break shrinks, one grows

Both diagnostics are now preserved. `run_epic.R` names them so `r_bridge` copies them out of
the temp directory into `results/diagnostics/`, and `config.DIAGNOSTICS_DIR` is re-pointed
under `results_synthetic/` for synthetic runs like every other output directory.

**Break 3 shrinks to almost nothing.** `otherCells` on the 122 real anatomic samples:

| mean | median | max | samples > 0.10 |
|---|---|---|---|
| 0.0013 | 0.0000 | 0.0427 | **0 of 122** |

EPIC's estimand is conditional on explained signal in principle; in practice the "other"
compartment absorbs essentially nothing, so dropping it and renormalising changes almost
nothing. Worth noting for its own sake: the roster discards 7.0% of the atlas and EPIC still
assigns ~0 to `otherCells` — it forces unexplained signal into the roster columns rather than
the compartment built for it.

**A fourth break appears, and it is the largest: EPIC did not converge on 52% of samples.**

| convergeCode | n | meaning |
|---|---|---|
| 0 | 59 | converged |
| 1 | 29 | iteration limit |
| 10 | 34 | degenerate simplex |
| | **63 of 122 did NOT converge** | |

This was emitted as an R warning on every run and never read, because nothing looked at
`fit.gof`. **EPIC is joint-2nd on the leaderboard at 0.9846, and more than half of its
estimates come from a failed optimisation.**

**Is it structure-dependent?** No: converged rates run 0.400 (CT) to 0.680 (MVP), chi-square
p = 0.275.

**Does it change the estimates?** Two different answers, and the difference matters.

*Endothelial: no — the marginal difference is a confound.* Converged samples show Endothelial
0.219 against 0.100 for non-converged, but MVP is 28.8% of the converged group and 12.7% of
the non-converged, and MVP carries Endothelial 0.724 against 0.003-0.024 elsewhere. Within
each structure the difference collapses to between -0.021 and +0.008. Simpson's paradox.

*Tumor: yes, and it survives the same test.* Within structure:

| structure | converged | not converged | delta | n (conv / not) |
|---|---|---|---|---|
| CT | 0.483 | 0.615 | **+0.132** | 12 / 18 |
| LE | 0.070 | 0.162 | **+0.092** | 8 / 11 |
| PAN | 0.796 | 0.638 | **−0.158** | 11 / 13 |
| IT | 0.340 | 0.337 | −0.003 | 11 / 13 |
| MVP | 0.151 | 0.150 | −0.002 | 17 / 8 |

Differences up to 0.16 in the dominant compartment, **and the sign flips between
structures**. Tumor carries **C1** (CT > LE) and **C7** (LE < IT < CT), so this lands on two
of the seven constraints — the same two the ISH check could not adjudicate.

**What is established and what is not.** Established: 52% non-convergence, not
structure-dependent, associated with tumour differences up to 0.16 that survive
stratification. **Not** established: that non-convergence *causes* them. Per-cell counts are
8-18, the comparison is observational, and a sample that is hard to fit may be genuinely
unusual rather than mis-fitted. The clean test is to re-run the non-converged samples with a
higher iteration limit and see whether their estimates move — which changes an EPIC parameter
and must be run as a declared diagnostic, reported beside the row and never substituted in.
