# Results freeze

*Written by `scripts/freeze_results.py` on 2026-10-07 20:39:10 EDT, from git `2415f6b` plus the uncommitted close-out changes, which are committed together with this manifest. Checked by `python3 scripts/check_core_closeout.py --freeze`.*

From this point every number, figure and table used for the write-up comes from this run. A file that changes after the freeze shows up in the check. The definition is in `docs/CORE_CLOSEOUT.md`.

| tree | files | size |
|---|---|---|
| `docs/figures/` | 32 | 7.9 MB |
| `docs/supplementary/` | 24 | 0.2 MB |
| `results/` | 621 | 356.1 MB |
| **total** | **677** | **364.2 MB** |

Every file's SHA-256 is in `docs/results_freeze_manifest.tsv`. Excluded: `results/verification/` (the re-run harness's own logs and state), R console logs (`*_r_stdout.log`, rewritten by every R call) and hidden files.

**Limit.** `results/` is not in git. The manifest records those files, but on a fresh clone they cannot be re-checked. The figures and tables in `docs/` are in git and can be.
