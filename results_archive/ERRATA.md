# Errata for the archived runs

Archives are immutable and verified by hash (`scripts/archive_run.py --verify <stamp>`).
Nothing inside an archive directory is ever edited or added — `verify` reports an added
file as ALTERED, and an archive that cannot be verified is worth less than one with a
known, documented fault. So corrections live here, outside them.

Citing an archive means citing what it actually says. Read this file alongside it.

---

## 2026-09-10T2039 — the archived `RESULTS.md` does not describe this run

**Found 2026-09-12.**

The archive is INTACT: all 102 files match their recorded hashes. The fault is not
corruption, it is that two of its files disagree with each other.

| | says |
|---|---|
| `RESULTS.md` (archived) | **UNREGISTERED** — "no public registration receipt exists … this project must not describe itself as pre-registered" |
| `results/anatomic/registration_status.json` (archived) | **REGISTERED**, at `2026-09-10T03:32:03Z`, "but 25 result file(s) are OLDER than the registration" |

Neither is the state of the project.

**Cause.** `scripts/run_all.py` archives the results tree but never regenerates
`RESULTS.md` first, and `summarize_results.py` only printed to stdout, so refreshing the
write-up depended on a human remembering to redirect it. The archived `RESULTS.md` has an
mtime of **2026-09-09 08:38**, more than a day before this run was archived
(2026-09-10 20:39). The archive therefore froze a write-up produced before the run it
accompanies — and before the registration was recorded, which is why it says UNREGISTERED.

**Also wrong, in the artefact rather than the write-up.** The "25 result file(s) are OLDER
than the registration" verdict is false; the true count is zero. `registration.status()`
was called partway through the run, before most artefacts had been rewritten, so the mtimes
it compared belonged to the *previous* run. See `docs/OPEN_DEFECTS.md` D9.

**What is sound in this archive.** Everything measured. The leaderboard, the per-constraint
and per-tumour tables, the agreement test, the equal-footing certificate, the
implementation report and the estimates are this run's and are unaffected — the fault is
confined to the pre-registration *narration*. The manifest's own summary fields are also
correct: `registration_state: REGISTERED`, `best_acs_method: music`, `best_acs: 1.0`,
`constraint_freeze_hash: 2d1fb47c…`, `git_commit: f22dfa7a`.

**The true registration state for this run:** REGISTERED with OSF at
`2026-09-10T03:32:03Z`, constraint file `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`,
and every result file postdates it.

**Fixed for future runs.** `summarize_results.py` gained `--write` and `--check`;
`registration_status.json` is now written last, after every artefact; the summariser
refuses to reprint a stale ordering verdict; and `run_all.py` regenerates the write-up
before archiving so an archive cannot again freeze a write-up older than its own results.

---

## 2026-09-05T2154, 2026-09-06T1103 — UNREGISTERED is correct for these

Both predate the OSF registration of 2026-09-10, and both say so. No erratum: they are
pilot runs, correctly labelled, and the leaderboards they carry are superseded rather than
wrong. Their DWLS rows are the Python reimplementation at ACS 0.877 — see
`docs/OPEN_DEFECTS.md` and `RESULTS.md` for why the genuine R package measures differently.

## `_safety_20260905T1647` — not citable

An unverified mid-run copy, kept for recovery only. Not committed, not a snapshot of any
completed run, and must not be cited.
