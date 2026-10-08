# Core close-out: what "done" means for this project

Defined 2026-10-07 at the user's request ("define core somewhere rather than existing in this plan"). The project
is closed when all eight items below pass. Each one names the evidence that shows it, and
`scripts/check_core_closeout.py` reads that evidence. It changes nothing.

    python3 scripts/check_core_closeout.py            # the eight items; exit 0 only if all pass
    python3 scripts/check_core_closeout.py --freeze   # re-hash every file the results freeze recorded

An exit code is not evidence on its own (CLAUDE.md), so the checker prints what it read for every item.

| # | core item | done when | evidence the checker reads |
|---|---|---|---|
| 1 | **Repaired DWLS / BayesPrism / quanTIseq** | The repairs are the default. The genuine DWLS, BayesPrism and quanTIseq packages run with no time budget. `IVYGAP_REPAIRED=0` restores the registered path and its budgets. The repaired results exist. | `config.REPAIRED_METHODS` and `r_bridge.timeout_for` in a fresh process, with the variable unset and with it set to 0; the default in `R/run_dwls.R`; `results/repaired/` |
| 2 | **Registered-vs-repaired verification** | Every core verification step has been re-run with the registered code path. Every output reproduced, or differs only for a cause stated in `EXPLAINED` and checked difference by difference. No step failed. Repaired results are reported beside the registered ones. | `results/verification/state.json`; `verify_rerun.explanation_check`; `docs/VERIFICATION_RERUN.md`; `docs/EVALUATION_MATRIX.md` |
| 3 | **D36 restoration** | The two shared extension per-sample files equal the 2026-10-06 snapshot byte for byte, and the guard against gene-subset writes is tested. | `cmp` against `results_snapshot_20261006/extension/`; the test in `tests/test_extension_panel.py` |
| 4 | **pytest** | The full suite passed with 0 failures and 0 errors, run after the last change to code or tests. | `results/verification/pytest_core.log`, its time against the newest file under `ivygap/`, `scripts/`, `R/`, `tests/` |
| 5 | **Claims ledger** | `docs/CLAIMS_LEDGER.md` was rebuilt after the verification report, and no claim is NOT VERIFIED or still pending a re-run. | the ledger's summary table; file times |
| 6 | **Registration audit** | `docs/REGISTRATION_AUDIT.md` was rebuilt after the verification report. Every rule PRECEDES the data it governs, except a rule whose note records why its time cannot be proven. | the audit table; file times |
| 7 | **Documentation audit** | Both document checkers are CLEAN, and every figure in `docs/figures/` has an entry in `docs/FIGURES.md`. | `check_doc_numbers.py --against-current`, `check_doc_statistics.py --strict`, figure file names against the index |
| 8 | **AI-use log** | `~/STS_AI_LOG/PROMPT_LOG.md` was exported after the latest prompt the student typed. | the newest typed prompt in the live Claude Code transcripts against the log's export time |

## Core verification steps

Every step in `scripts/verify_rerun.py` is core, except those in its `SUPPLEMENTARY` table:
- `extension_tcga`, `cdseq_anatomic`, `reference_sensitivity_neftel`, `reference_sensitivity_darmanis`,
  `unmix_s2_anatomy`, `recide_anatomic`;
- their downstream statistics reproduced from the saved fits in the fast and medium tiers;
- the user's decisions of 2026-10-07: "core steps only", then "you completing the project NOW";
- run any of them with `--only <step>`.
- **After the close-out they were run anyway** (user, 2026-10-07: "Still run what's left to do"), queued behind
  one another on the registered path. Their verdicts update `docs/VERIFICATION_RERUN.md`; no core item depends
  on them, and they cannot change a frozen file.

## The results freeze

- **When:** after the verification report is written, the results are frozen. Every figure, supplementary table
  and document is regenerated from that frozen run.
- **What is recorded:** `scripts/freeze_results.py` writes `docs/RESULTS_FREEZE.md` and
  `docs/results_freeze_manifest.tsv`, with a SHA-256 for every file under `results/`, `docs/figures/` and
  `docs/supplementary/`. Excluded: the verification folder, and R console logs (`*_r_stdout.log`), which every
  R call rewrites, tests included.
- **How it is checked:** `check_core_closeout.py --freeze` re-hashes them. Any later run that touches a frozen file
  shows up there.
- **Limitation:** `results/` is not in git, so the manifest is the only record of those files. On a fresh clone they
  cannot be re-checked.
