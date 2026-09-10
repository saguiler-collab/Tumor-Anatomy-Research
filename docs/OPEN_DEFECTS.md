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

### MuSiC, SCDC, Bisque, BayesPrism — SUSPECTED, not yet verified

Each returns something documented as *cell type proportions* and each models cell size
somewhere internally. None has been verified the way EPIC was, and none should be changed
until it has been.

| method | driver line | returns | why suspected |
|---|---|---|---|
| MuSiC | `R/run_music.R:61` | `est$Est.prop.weighted` | **Upgraded to near-confirmed 2026-09-10** — see below |
| SCDC | `R/run_scdc.R:46` | `est$prop.est.mvw` | `ct.cell.size` defaults to library size computed from the data |
| Bisque | `R/run_bisque.R:35` | `t(est$bulk.props)` | documented as cell proportions |
| BayesPrism | `R/run_bayesprism.R:72` | `get.fraction(..., "type")` | theta is a cell-type fraction |

#### MuSiC, and why it makes the fix harder than EPIC's

Source inspection of `MuSiC::music_prop` (no execution): it takes a `cell_size`
parameter, and when that is NULL it derives `M.S`, the mean cell size per type, from the
data. Its model is `Y_jg = sum_k S_k p_k theta_kg` — cell size is divided out internally,
so `Est.prop.weighted` is a **cell** proportion. Our correction is therefore a second one.

The difficulty: MuSiC returns only `Est.prop.weighted` and `Est.prop.allgene`
(`p.weight`, `p.nnls`). **It exposes no mRNA-proportion variant.** Unlike EPIC, you cannot
simply take the other output.

That splits the fix into two incompatible designs, and choosing between them is a real
decision rather than a patch:

**(A) Every method returns an mRNA share; the project converts once, centrally.**
This is what `METHODS.md` and the registration describe, and it keeps every method on
identical cell-size factors — which is the point of doing it centrally. But it is not
reachable for MuSiC without multiplying its output back by `M.S`, reconstructing a
quantity the package deliberately does not return.

**(B) Each method returns cell fractions using its own factors; the project applies
nothing to those methods.** Reachable for every package. But methods would then be using
*different* cell-size factors — EPIC's `mRNA_cell`, MuSiC's data-derived `M.S`, ours for
the least-squares family — which breaks the uniformity the central correction exists to
guarantee, and weakens the equal-footing claim.

There is a third option worth considering: pass the project's factors *into* each package
that accepts them (`MuSiC(cell_size=)`, `SCDC_ENSEMBLE(ct.cell.size=)`, `EPIC(mRNA_cell=)`)
and skip the central step for those methods. That preserves both uniformity and each
package's own machinery, at the cost of a per-method wiring change and a declared
deviation for each.

**No design should be adopted until every method is verified**, because (B) and the third
option require knowing exactly which methods already correct and which do not.

**How to verify each**, the same way EPIC was: run the package on its own bundled example,
obtain both an uncorrected and a corrected quantity if it exposes them, and compare. Where
a package exposes only one, read whether its vignette calls the output a transcript share
or a cell share. Do not guess — a wrong "fix" here silently biases every result.

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

The bias is **uniform within a method across all samples** — a per-type constant factor
applied twice. It therefore does not disturb any *within-method, between-structure*
ordering, which is what ACS scores. ACS is a rank statistic over structures inside a
tumour, and multiplying a cell type's column by a constant cannot reorder that column
across structures.

So: the ACS leaderboard is **not** invalidated. What moves is the **accuracy** side of the
agreement test, and therefore rho.

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
