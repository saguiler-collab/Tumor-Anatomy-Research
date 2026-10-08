# Results freeze v2

*Written by `scripts/freeze_results.py` on 2026-10-08 08:43:48 EDT, from git `2a62ac5` plus the uncommitted close-out changes, which are committed together with this manifest. Checked by `python3 scripts/check_core_closeout.py --freeze`.*

From this point every number, figure and table used for the write-up comes from this run. A file that changes after the freeze shows up in the check. The definition is in `docs/CORE_CLOSEOUT.md`.

| tree | files | size |
|---|---|---|
| `docs/figures/` | 32 | 7.9 MB |
| `docs/supplementary/` | 24 | 0.2 MB |
| `results/` | 621 | 356.4 MB |
| **total** | **677** | **364.6 MB** |

Every file's SHA-256 is in `docs/results_freeze_manifest.tsv`. Excluded: `results/verification/` (the re-run harness's own logs and state), R console logs (`*_r_stdout.log`, rewritten by every R call) and hidden files.

**Limit.** `results/` is not in git. The manifest records those files, but on a fresh clone they cannot be re-checked. The figures and tables in `docs/` are in git and can be.

## Changes from freeze v1

Computed by comparing this manifest with `docs/results_freeze_manifest_v1.tsv`. Every change has a stated reason; the freeze refuses to write otherwise. The previous version is kept as `docs/results_freeze_manifest_v1.tsv` and `docs/RESULTS_FREEZE_v1.md`.

| file | change | reason |
|---|---|---|
| `results/extension/estimates_full_h5ad_extension.csv` | changed | merged: the extension re-fit's per-sample estimates for aric, fardeep, lindeconseq, rnasieve (616 rows) appended after the 462 registered rows, which are byte-identical |
| `results/extension/estimates_full_lgg_h5ad_extension.csv` | changed | merged: the extension re-fit's per-sample estimates for aric, fardeep, lindeconseq, rnasieve (2040 rows) appended after the 1530 registered rows, which are byte-identical |
| `results/lymphoid_pooled_tnk.json` | changed | re-run on the merged extension tables: adds the extension-panel methods (each flagged extension_panel = True; 80 added fields); every registered method's numbers are unchanged (B > T+NK in most samples: 6 of 8, 4 of 8, 9 of 12, 3 of 10) |
| `results/reference_sensitivity_gbmap_linear_vs_darmanis.json` | changed | D39: recomputed on the registered 651-gene marker space (the registered artefact used the superseded 657-gene log-build space) |
| `results/reference_sensitivity_gbmap_linear_vs_neftel.json` | changed | D39: recomputed on the registered 651-gene marker space (the registered artefact used the superseded 657-gene log-build space) |
