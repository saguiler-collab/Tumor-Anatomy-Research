# Full-project sweep, 2026-09-19

Ran the project top to bottom looking for failures and defects, rather than re-running only what
had changed. **Three real defects found, all three fixed. One serious performance defect found and
still open.** Everything else verified clean.

The sweep harness and per-stage logs live in the session scratchpad; the checks are reproducible
from the commands recorded below.

---

## Defects found

### 1. FIXED — an unorderable sample was credited to `T_cell` by column order

A cross-artefact check flagged the GBM *T-ranked-first* figure as irreproducible from raw data:
the artefact said **0.6387**, recomputing gave **0.6323**.

`ordering_of` used `sort_values(ascending=False)`, and pandas sorts NaN **last**. A sample where
EpiDISH returns exactly 0 for T, B *and* NK renormalises to all-NaN, sorted to the frame's own
**column order**, and came back `"T>B>NK"` — counted as "T ranked first" purely because `T_cell`
is the first column. `TCGA-14-0862-01` is that sample.

One sample in 155, **and in the direction that flatters this project's own registered prediction**
(that T exceeds B), which is why it was fixed rather than noted. `ordering_of` now returns `None`
for any NaN row, callers exclude those rows and print how many, and the artefacts record
`n_orderable_samples` and `n_unorderable_zero_lymphoid` so the denominator is never implicit.

- **Corrected:** GBM T-ranked-first **0.6387 → 0.6364** (154 orderable of 155). LGG unaffected
  (zero unorderable, 0.7925 unchanged).
- **No claim changes.** T > B is a pairwise comparison on a different denominator: 94.8% / 93.2%,
  untouched.
- **The test gap is instructive.** `test_all_zero_lymphoid_row_becomes_nan_not_zero` verified that
  `renormalise` *produces* NaN, but nothing checked what `ordering_of` *did* with it — the suite
  covered half of a two-step path. Now covered for all-NaN and partial-NaN rows.

### 2. FIXED — two superseded scripts did not compile

`scripts/fetch_ivygap_ish.py` and `scripts/figure_data.py` both raised `SyntaxError`. Each is
correctly superseded, and the deprecation notice had been prepended as a **second docstring above
the original** — which makes the original a bare statement and the following
`from __future__ import annotations` illegal. Anyone running them got a confusing SyntaxError
instead of the "SUPERSEDED — do not use, and do not cite its output" message written for exactly
that moment.

**How it hid:** `ast.parse` accepts a misplaced `__future__` import; only `compile()` rejects it.
A parse-based audit reported all 51 scripts clean. The audit now compiles. Both docstrings are
merged, notice first, original kept verbatim.

### 3. FIXED — `RESULTS.md` went stale, and its own test caught it

Defect 1 changed the artefacts, so the shipped `RESULTS.md` no longer matched the generator:
`test_the_shipped_results_md_matches_the_generator` failed (1 failed, 282 passed). That test
exists because `RESULTS.md` claims in its own header that nothing in it is transcribed by hand,
and it is the only thing enforcing the claim. Regenerated.

### 4. OPEN, serious — `run_all.py --synthetic` does not complete

**This is the CLAUDE.md validation gate** (*"validate end-to-end with `python scripts/run_all.py
--synthetic` before touching the real-data path"*), and it is currently unusable.

Measured: **8 h 54 min elapsed, 11 m 56 s of CPU (2.2%), and `results_synthetic/` contained only
empty directories and `.run.lock`** — every one timestamped at the moment the run started. Nine
hours, zero artefacts. Killed.

Relaunched with `-u` to get streaming output, which located it: the run reaches **Stage 3, the
benchmark**, and stalls there. Per-method timings on 500 synthetic mixtures show the cost is real
and badly distributed — `nnls` **0.1 s**, `svr` **175.7 s**. With fifteen methods this stage alone
is tens of minutes before anything writes.

**The cause I first proposed was wrong, and is withdrawn.** I attributed it to OPEN_DEFECTS
**D18** — orphaned R cluster workers holding the inherited stdout pipe so `communicate()` never
returns. Tested: a child whose detached grandchild keeps the pipe open past the budget raises
`TimeoutExpired` at **3.0 s against a 3 s budget** on the *pre-fix* code. The pipe-EOF story does
not explain it. See D18's 2026-09-19 correction.

**A second candidate, now being tested rather than assumed: resource starvation of my own making.**
The 8 h 54 min run overlapped the rest of this sweep — a full `pytest` pass, two h5ad cohort legs,
several probe scripts — on a 4-core machine whose load average reached **50.8**. At the 48-minute
mark that process had 9 m 21 s of CPU; nearly six hours later it had 10 m 51 s. It was not spinning,
it was **starved**. A relaunch with streaming output and a quiet machine is progressing normally
through the method panel (`nnls` 0.1 s, `svr` 175.7 s, `cibersortx` 278.0 s, `cibersortx_smode`
239.1 s, `elastic_net` 37.4 s, `bayesian` 203.9 s), which is slow but finite.

**If that relaunch completes, the "hang" was my scheduling, not the pipeline**, and the honest
finding shrinks to: Stage 3 is expensive enough (minutes per method, ~15 methods) that the
validation gate takes hours and writes nothing until it finishes, which makes it look hung and
makes it easy to starve. That is still worth fixing — a gate that cannot be run casually is a
gate that is not run — but it is a performance and observability problem, not a deadlock.

**Why it matters beyond inconvenience:** a validation gate nobody can run is a gate that is not
being used, and every result in this project is supposed to sit behind it. It also means the
synthetic path's negative controls — the ones that matter more than the positive ones — have not
been exercised in this sweep.

---

## Verified clean

| check | result |
|---|---|
| `pytest tests/` | **283 passed, 1 skipped** after the fixes (1 failure before, defect 3) |
| every `results/*.json` parses | **58 of 58** |
| every script compiles | **51 of 51** (2 failures before, defect 2) |
| `results_archive/` hash integrity | **4 of 4 runs INTACT** — 96, 97, 102 files matching recorded hashes |
| **constraint freeze hash** | **identical across all four archived runs and live** — `2d1fb47c98832adf…`. The file has never changed. Only one commit has ever touched it. |
| pre-registration | `register.py`: registered at OSF `dm2t8`, 2026-09-10T03:32:03Z, **and every result file postdates it** |
| `RESULTS.md` vs generator | `CHECK OK` (565 lines) |
| ACS numbers in live docs vs artefacts | `check_doc_numbers.py`: **67 cells across 50 files, all match** |
| withdrawn wording in artefacts | none — no artefact still says "problem as posed", "identifiability failure" or "333 HM450" |
| cross-artefact recovery figures | `cross_cohort_validation` vs `equal_footing_ranking` agree within 0.03 for every method checked |
| the 18 analysis re-runs | **all exit 0** (see below) |

### The 18 analysis scripts, re-run from current artefacts

All completed; the two flagged lines were checked by hand and are legitimate — the comparability
module reporting `bayesian_hierarchical failed in this run` and `dwls failed in this run`, which is
the exclusion machinery working, not an error.

```
summarize 5s · doc_numbers 2s · manuscript 0s · identifiability 2s · model_fit 21s
lymph_gbm 2s · lymph_lgg 2s · lymph_gbm_h5ad 1s · lymph_lgg_h5ad 1s
ef_gbm 3s · ef_lgg 2s · agree_frozen 18s · agree_h5ad 16s
meth_gbm 9s · meth_lgg 12s · cross_cohort 1s · method_complete 3s · register 4s
```

---

## What this says about the project's readiness

The scientific artefacts are in good shape: every number in a live document traces to an artefact,
every archived run is hash-intact, the constraint file is provably unchanged since registration,
and the analysis scripts all reproduce their outputs in seconds.

The two classes of problem found are worth distinguishing:

- **Defect 1 is the one that mattered.** It was small, silent, and biased toward this project's
  own prediction — the exact shape of error that a sweep exists to catch and that no amount of
  re-reading would have found.
- **Defect 4 is the one still outstanding.** It does not affect any published number, because no
  result depends on the synthetic path. It affects the project's ability to *check itself* before
  the next change, which is worse in the long run and should be fixed before new datasets are
  added.

**Recommendation before adding more data**, revised now that the D18 cause is withdrawn:

1. **Let the quiet relaunch finish** and record whether the gate completes. That single fact
   decides whether this is a deadlock or a scheduling problem, and it is being measured rather
   than argued.
2. **Do not run the gate concurrently with anything else.** Load average 50.8 on four cores is
   what turned a slow stage into an apparent hang, and it was self-inflicted.
3. **Give Stage 3 progress output and a per-method budget**, so a slow method is visibly slow
   rather than indistinguishable from a hang. The run already prints per-method timings; the
   problem is that nothing reaches disk until the stage ends.
4. `_run_bounded` now captures to temporary files rather than pipes. That is defensive — no
   descendant can apply backpressure, and 2 MB of output is covered by test — but it is **not**
   presented as the cure, because the cause it was written for turned out not to be the cause.
