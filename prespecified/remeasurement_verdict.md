# What the BayesPrism / DWLS re-measurement can and cannot conclude

**Written 2026-09-19, BEFORE the re-measurement finished.** The point of writing it first is that
the interpretation cannot then be chosen to suit the numbers.

## What is being measured

Two methods have been reported in every real run of this study as `python-reimplementation`, never
as the published R package:

- **`dwls`** — exceeded its 2,400 s budget. `scripts/remeasure_method.py` re-runs it with a
  14,400 s budget. Nothing about the method changes; only wall-clock is allowed to differ.
- **`bayesprism`** — failed with `Error in unserialize(node$con)` in 155.3 s, well inside budget.
  That is a `parallel` socket-cluster worker being killed for memory on an 8.6 GB machine
  (OPEN_DEFECTS D18). Re-run with `--cores 1`, which removes the cluster.

Both are scored on the **anatomic (Ivy GAP) cohort** with the leaderboard's own 657-gene space,
so the numbers are comparable to the existing ACS leaderboard.

## The four outcomes, and what each one means

| outcome | what it means | what goes in the paper |
|---|---|---|
| **Both run as the genuine R package** | Limitation 8 is substantially resolved. The panel becomes 13 of 15 genuine packages plus two declared reimplementations (`cibersortx`, `cibersortx_smode`, which implement a licence-gated published algorithm). | Report the genuine ACS for both **alongside** the original rows, not substituted for them. State that the reimplementation rows in the archived runs stand as recorded. |
| **BayesPrism runs, DWLS still times out** | The memory diagnosis is confirmed and the DWLS cost is genuinely beyond this hardware. | Report BayesPrism genuine; keep DWLS as a reimplementation with the measured runtime, and say that DWLS is the one method this study could not run as published. |
| **DWLS runs, BayesPrism still fails** | The `n.cores` fix is wrong or incomplete, and D18's root cause needs revisiting **again**. | Say so. Two wrong causes are already on the record; a third would be recorded the same way. |
| **Neither runs** | The hardware is the binding constraint for both, and that is a finding about reproducibility on commodity machines, not a gap in the study. | Report both as reimplementations with the measured failure modes, and state the RAM figure that caused it. |

## Committed in advance

1. **The result is reported whichever way it falls.** If genuine BayesPrism or genuine DWLS scores
   *worse* than its reimplementation, that is the finding. If it scores better, the original rows
   are still not deleted — they are labelled.
2. **No substitution into the archived runs.** `results_archive/` is hash-verified and stands as
   recorded. A re-measurement is an additional, separately labelled artefact.
3. **ACS must not select anything.** These two rows are being *measured*, not chosen between. The
   project's standing invariant applies: ACS is the object under test and may not be used to pick
   a method, a parameter or an implementation.
4. **A materially different result from the same model is itself reportable.** If single-threaded
   BayesPrism differs substantially from the three-worker run on the pseudobulk benchmark, that is
   a fact about the method's determinism and it gets stated, not smoothed.
5. **This does not change any headline.** The lymphoid failure, the model-fit bound, and
   anatomy-detects-but-cannot-rank rest on neither method. The most this can do is remove a
   limitation.

## What this re-measurement cannot do

It cannot make the panel fully comparable. `cibersortx` and `cibersortx_smode` remain this
project's implementation of a published algorithm whose hosted service is licence-gated, and that
is declared in `docs/METHODS.md`. It also cannot speak to the **purity** arm: `remeasure_method.py`
scores the anatomic cohort, so an improved implementation here does not update the DNA-purity
recovery figures, which continue to carry the reimplementation rows.
