# Close-out report: the verification is complete

**2026-10-08.** Covers the work from the core close-out (`e908c55`, 2026-10-07 21:00) to `2c85320` (2026-10-08
08:28). Every number here is from `docs/VERIFICATION_RERUN.md`, `results/verification/state.json` or the commits
listed at the end.

## In one paragraph

Every step of the verification has now been re-run, 75 of 75, and none failed. Each output either reproduced the
registered result or differs for a cause that was found and checked difference by difference. The results freeze
held throughout: none of the 677 frozen files changed. The last runs turned up one finding that matters for the
write-up (D39, reference sensitivity) and one more method result that now reproduces (ReCIDE). They also filled a
known data gap in the extension panel. One decision is open, below.

## Status

| check | result |
|---|---|
| Verification steps re-run | **75 of 75** (43 fast, 21 medium, 11 heavy); none failed, none left |
| Core close-out (`scripts/check_core_closeout.py`) | **8 of 8 pass** |
| Results freeze (`--freeze`) | **0 of 677 files changed** |
| Full test suite | **483 passed, 1 skipped, 0 failed** |
| Claims ledger | **8 of 8 claims VERIFIED** (2 of them with explained differences) |
| Registration audit | 8 rules provably precede their data; 1 not provable here, with the reason recorded (GIMiCC) |

## Verification totals

| verdict | outputs |
|---|---|
| REPRODUCED (identical within 1e-9) | 74 |
| REPRODUCED + NEW ROWS (every registered row identical; rows added) | 4 |
| REPRODUCED + NEW FIELDS (every registered value identical; fields added) | 2 |
| NEW (created after the snapshot; each reproduced its own first run) | 5 |
| PASSED / AGREES (checks with a verdict, not an output) | 2 / 1 |
| DIFFERS, every difference explained and checked | 9 |

### The 9 outputs that differ, and why

Each cause is written into the harness, and `explanation_check` confirms every recorded difference matches it.

| step | outputs | differences | cause |
|---|---|---|---|
| `registered_pipeline` | 3 | 34 | D35: in the registered run, BayesPrism and DWLS overran their 2,400 s budget by seconds and fell back to this project's versions. Every other method reproduced exactly. |
| `tcga_gbm_frozen`, `tcga_lgg_frozen` | 2 | 13 + 13 | Additions only. `bayesian_hierarchical` now runs (D20) and two provenance fields are new. Every registered value is identical. |
| `anatomy_vs_biology` | 1 | 1 | A text correction: "94 minutes", where it said 95. |
| `reference_sensitivity_neftel`, `_darmanis` | 2 | 107 + 90 | **D39** (below): the registered results used a superseded marker set. |
| `recide_anatomic` | 1 | 2 | No result differs. The two fields are the estimates file's path (a re-run now keeps its estimates beside its own report) and the tail of ReCIDE's console output. |

## What the last runs showed

### The extension panel reproduces in all four arms

`extension_tcga` was re-fitted on the registered code path in 3 h 7 min. It covers GBM and LGG, each with the frozen
signature and the donor-level reference.
- All 4 summaries reproduced exactly.
- The two frozen-signature per-sample tables are the same rows in a different order. The merge appends re-fitted
  methods at the end, and the harness now matches rows by key.
- **It closes a gap the worklog had logged.** On the two donor-level arms the re-fit produced per-sample estimates
  for ARIC, FARDEEP, LinDeconSeq and RNA-Sieve: 616 GBM rows and 2,040 LGG rows. Until now only their summaries
  existed. Their summaries equal the registered ones.
- To respect the freeze, these rows are kept in `results/verification/rerun_outputs/extension_tcga/`, not merged
  into `results/`.

### The supplementary re-fits

These five had been skipped under "core steps only". All have now run.

| step | time | verdict |
|---|---|---|
| `cdseq_anatomic` (reference-free arm) | 8.7 min | all 6 outputs REPRODUCED |
| `unmix_s2_anatomy` (E2's secondary arm) | 5.1 min | REPRODUCED |
| `recide_anatomic` (ReCIDE on Ivy GAP) | 5.6 h | reproduces: ACS, CI, null p, pairs and all 1,098 estimate cells identical |
| `reference_sensitivity_neftel` | 4.9 min | DIFFERS: D39 |
| `reference_sensitivity_darmanis` | 5.9 min | DIFFERS: D39 |

ReCIDE took two attempts:
- **The first failed** at its last line after 93 minutes, from a path bug introduced the evening before in
  `run_recide.py`.
- **The fix** makes the path absolute (`ea07427`), and the test now exercises the call that failed.
- **The retry took 5.6 h** because the Mac slept for about an hour and Low Power Mode halved the CPU speed, until it
  was turned off.

## D39: the one finding that touches the write-up

**What happened.**
- The reference-sensitivity analysis asks how much the method ranking changes when the single-cell reference
  changes.
- Its registered results were computed on 2026-09-15, on the leaderboard's marker set of that day: 657 genes,
  selected on the log-layer build.
- On 2026-09-21 the registered leaderboard was rebuilt on raw counts, with 651 marker genes. This analysis was never
  recomputed afterwards.
- The re-run uses the registered 651-gene set.

| comparison | shared genes | ordering agreement, registered (old genes) | on the registered gene set |
|---|---|---|---|
| GBmap vs Neftel | 654 -> 648 | 0.5099 | **0.8867** |
| GBmap vs Darmanis | 651 -> 647 | 0.3655 | **0.0098** |

**What it means.**
- **The conclusion holds:** the ACS ranking of methods depends on which reference is used. On the current gene set
  it is nearly preserved with one alternative reference and destroyed with the other.
- **The magnitudes do not hold.** Changing about 6 of 650 genes moved them by 0.38 and 0.36, so neither number is a
  stable property. Describe reference sensitivity as large and unpredictable, not as a single correlation.
- Every difference falls in the method scores and the ordering statistics derived from them. Sample counts and
  implementations are unchanged. The full record is in `docs/OPEN_DEFECTS.md`, entry D39.

## Problems found and fixed since the close-out

Each is logged in `docs/OPEN_DEFECTS.md` or the worklog, with a test where code changed.

- **Harness safeguard after the freeze.** A re-run now copies each declared output before it runs and puts the
  frozen bytes back afterwards, keeping the re-run copy aside (`dc2365c`).
- **Two more writers that could overwrite registered files.** ReCIDE now keeps its estimates beside its own report,
  and CDSeq's files are declared outputs (`acd6b1d`).
- **Two false alarms in the comparison.** Identical tables in a different row order, and outputs written beside the
  registered file, are now compared correctly (`acd6b1d`, `2c85320`).
- **The ReCIDE path bug** described above (`ea07427`).

## Decision needed

**D39: which reference-sensitivity numbers should the write-up use?**

- **Keep the registered numbers** (0.5099 and 0.3655) and note D39. Nothing changes; the freeze stands as is.
- **Use the numbers on the registered gene set** (0.8867 and 0.0098), as a labelled update to the freeze (a
  "freeze v2" that records exactly which files changed and why).
  - This is the recommended option: the current numbers are computed on the same gene set as every other registered
    result.
- **Report both**, with D39 as the explanation.

### Optional, not needed for the write-up

- Merge the extension panel's new per-sample rows into `results/`. This is also a post-freeze change.
- SCDC ENSEMBLE with a second single-cell reference (memory-heavy; documented in `docs/METHOD_REPAIRS.md`).
- The CPTAC authors' curated single-nucleus labels (a data request; draft in `docs/IMPACT_KIT.md`).

## Commits since the close-out

| commit | time | what |
|---|---|---|
| `e908c55` | 10-07 21:00 | the core close-out: repairs by default, D36-D38 restored, results frozen |
| `dc2365c` | 10-07 21:15 | after the freeze, a re-run never changes a frozen file |
| `acd6b1d` | 10-07 22:23 | the remaining re-fits prepared (ReCIDE path rule, CDSeq outputs declared, row order) |
| `7cc0824` | 10-08 00:15 | the extension panel reproduces all four arms |
| `bdb215a` | 10-08 00:49 | D39 recorded and explained |
| `ea07427` | 10-08 03:35 | ReCIDE path fix; re-run |
| `2c85320` | 10-08 08:28 | verification complete: 75 of 75 steps, every difference explained |

All are on `anatomy-test-full-database`, and `main` has been fast-forwarded to each.

## Addendum: decisions taken (2026-10-08, results freeze v2)

| decision | what changed |
|---|---|
| **D39: "Freeze v2, current numbers"** | The two 651-gene reference-sensitivity results became the registered artefacts: ordering agreement **0.8867** (Neftel) and **0.0098** (Darmanis). Wording that followed the old numbers was corrected in the manuscript builder and the RESULTS.md renderer. Live documents updated; history documents marked superseded; the statistics checker pins both values. |
| **"Merge the new rows"** | The extension panel's per-sample estimates for ARIC, FARDEEP, LinDeconSeq and RNA-Sieve (616 GBM and 2,040 LGG rows) were appended to the donor-level tables. The registered rows are byte-identical. `lymphoid_pooled_tnk.json` was re-run on them: it gains the extension methods, each flagged, and every registered number is unchanged. |

**Freeze v2** (`docs/RESULTS_FREEZE.md`) records exactly 5 changed files against v1, each with its reason; v1 is kept.
The freeze check passes against v2.

