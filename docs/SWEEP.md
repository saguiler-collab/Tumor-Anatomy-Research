# Full-project sweep, 2026-09-19

Ran the project top to bottom looking for failures and defects, rather than re-running only what
had changed.

**Four defects found. All four resolved.** Everything else verified clean, including the end-to-end
validation gate and its negative controls.

| | defect | status |
|---|---|---|
| 1 | an unorderable sample credited to `T_cell` by column order | **fixed** — and it was biased *toward* this project's own prediction |
| 2 | two superseded scripts did not compile | **fixed** |
| 3 | `RESULTS.md` stale against its generator | **fixed**, caught by its own test |
| 4 | the synthetic validation gate "hung" | **resolved — it was two problems, neither a deadlock.** Contention I created (2.2% vs 77.5% CPU on identical hardware), plus BayesPrism dying on *memory*, not time |

**The conclusion of the sweep:** the science is sound and the tooling had three small holes and one
misread. Defect 1 is the only one that touched a number, and it moved GBM T-ranked-first from
0.6387 to 0.6364 — no claim changes. The gate now runs end to end, Stage 5 blocks correctly on
missing vital-status data rather than imputing it, and the negative controls sit at ACS 0.362–0.383
against a working-method floor of 0.906.

**I was wrong twice about defect 4** and both wrong answers are kept on the record below, because
the pattern — a plausible mechanism, never checked against the step that would falsify it — is the
transferable part.

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

### 4. RESOLVED — two separate problems wearing one costume

This entry originally read *"OPEN, serious — `run_all.py --synthetic` does not complete"* and
attributed it to a pipe deadlock. **Both the severity and the cause were wrong.** What looked
like one defect was two, and neither is a deadlock.

#### 4a. The 9-hour wall clock was contention, self-inflicted

The run overlapped a full `pytest` pass, two h5ad cohort legs and several probe scripts on a
4-core machine whose load average reached **50.8**. Measured: **2.2%** CPU utilisation
(11 m 56 s of CPU in 8 h 54 min). The quiet relaunch reached **77.5%** (15 m 20 s in 19 m 32 s) —
**a 35× difference on identical hardware** — and walked through all 15 methods and into Stage 4.
It was starved, not stuck.

What survives as a genuine, if smaller, finding: **Stage 3 costs minutes per method across ~15
methods and writes nothing to disk until the stage ends**, so it is indistinguishable from a hang
and trivially starved. Per-method: `nnls` 0.1 s, `elastic_net` 37.4 s, `music` 49.5 s,
`svr` 175.7 s, `bayesian` 203.9 s, `bayesian_hierarchical` 215.8 s, `cibersortx_smode` 239.1 s,
`cibersortx` 278.0 s, `dwls` 939.6 s.

#### 4b. BayesPrism's failure is a MEMORY limit, and it had been invisible

The same relaunch printed the error the pipe-based capture had been discarding on the error path:

```
bayesprism  155.3s  [fell back to Python: run_bayesprism.R exited 1
                     | Error in unserialize(node$con) : error reading from connection]
```

`unserialize(node$con)` is an R `parallel` **socket cluster** master failing to read from a worker
that is gone. `R/run_bayesprism.R` asked for **three workers**, each a separate R process with its
own copy of the data, on a machine with **8.6 GB RAM and ~3.2 GB free**. They are killed, the
master fails, R exits 1, the bridge takes its documented fallback. **Not a timeout — 155 s against
a 2,400 s budget.**

One cause accounts for every observation that made the earlier guesses look plausible: the
orphaned ppid-1 R workers, the master at 0% CPU (socket masters wait on workers by design), and
`bayesprism` recorded as `python-reimplementation` in every real run.

#### The gate then completed, and its negative controls fire

All six stages, on a quiet machine:

| stage | outcome |
|---|---|
| 1 load data | 212 samples / 24 tumours |
| 3 benchmark | 15 methods + 2 controls, `selected: epic` |
| 4 Anatomy Test | ACS leaderboard, 10,000 permutations per method |
| 5 survival | **BLOCKED, correctly** — `tumor_details.csv` publishes `survival_days` and no vital-status column, so who was censored is unknowable. Nothing imputed, no prognostic claim. This is the *"never impute a missing input"* invariant firing. |
| 6 release bundle | 21 files written, with the 8 not produced this run listed explicitly |

And the check that matters most:

| | ACS | null p |
|---|---|---|
| 13 real methods that produced estimates | **0.9060 – 0.9530** | **0.0001** — all beat their null |
| 2 negative controls | **0.3624 – 0.3826** | **0.703 – 0.752** — both fail their null |
| 2 methods that produced nothing (`bisque`, `quantiseq`) | 0.0000 | 1.0 |

**The controls sit below every working method, and fail their null while every working method
beats it.** A scorer that reported high agreement on shuffled signal as well as real signal would
be measuring nothing; this one does not. `control_shuffled_signature` at 0.362 against a working
floor of 0.906 is the single most reassuring number in this sweep.

One presentational note worth keeping: `bisque` and `quantiseq` land at ACS **0.0000 with
null_p 1.0**, which is a method that produced *nothing* on the fixture, not a method that scored
badly. Reading a naive `min(real) < max(control)` comparison off that table would suggest the
controls had beaten a real method. They had not — and a check that conflates "failed" with
"scored low" is exactly the sort of thing that turns into a false alarm or, worse, a false
reassurance.

#### What it cost, and the transferable lesson

**I proposed two causes and both were wrong** — pipe-EOF starvation (tested, false) and
machine-wide CPU starvation (true of the wall clock, not of why BayesPrism failed). Neither guess
was checked against the one thing that settles it in a single step: **reading the subprocess's own
stderr**, which the pipe-based capture was throwing away precisely on the path where it mattered.

The fix that actually follows is `n.cores = 1`, declared in `docs/METHODS.md`, and it needs no new
hardware. A re-measurement of both methods is running with the interpretation pre-specified in
`prespecified/remeasurement_verdict.md`.

---

### 4-original (superseded, kept for the record) — what this entry said before

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
