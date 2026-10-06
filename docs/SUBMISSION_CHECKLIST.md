# Submission checklist

> ## ⚠ Updated 2026-10-01 for Regeneron STS — read this before the rest
>
> The venue is now **Regeneron STS**, not CJSJ. Three entries below are superseded:
>
> | was | now | where |
> |---|---|---|
> | **headline: 0 of 12 (GBM) and 0 of 12 (LGG)** reproduce T > B | true on the **frozen** signature only. On the donor-level raw/X reference **3 of 14** and **2 of 13** do; Bisque's agreement is its reference's composition by construction. Report both references. | OPEN_DEFECTS **D23**; `docs/WHY_B_OVER_T.md` §7 |
> | DWLS "the package scores 0.7846" | that figure was on the **log** layer (D22). On raw/X the genuine package returns no estimate for 73 of 122 samples; ACS 0.7000 on 26 pairs is not comparable. | OPEN_DEFECTS **D22**; MANUSCRIPT §4.7 |
> | "ReCIDE absent" | ReCIDE is vendored (`vendor/ReCIDE`, commit 31bbd6b) and being gated; four more top-tier methods (FARDEEP, LinDeconSeq, ARIC, RNA-Sieve) form a labelled post-registration extension panel. | MANUSCRIPT §4.10; `ivygap/deconv/extension.py` |
>
> Also new and load-bearing: **methylation's per-sample lymphoid split failed a reference-free
> positive control**, so per-sample statements describe methylation and are not a validated
> truth (`results/lymphoid_tracking.json`). Cohort-level T > B stands.
>
> **New 2026-10-03, load-bearing: the truth itself was tested** (WHY_B_OVER_T §7n; MANUSCRIPT §4.4).
> - A glioma-specific methylation method (GIMiCC [52]) **agrees in GBM** (T above B, every control
>   passing) and **disagrees in LGG** (B above T). Its registered reading is INCONCLUSIVE.
> - Direct cell counts favour T above B: an independent atlas 18/18 patients (2 of them LGG); GBmap
>   96/98 donors (all GBM); flow cytometry [53].
> - ~~State the lymphoid headline for GBM; for LGG only with a caveat.~~ **Updated 2026-10-06:** flow
>   cytometry of 17 IDH-mutant gliomas (Klemm [53], Figure 1F, measured and validated) puts T 21-fold
>   above B. **State the headline in both cohorts.** Report that one of two methylation instruments
>   (GIMiCC) disagrees in LGG.
>
> The registered statistics — ACS leaderboard, control margin, ρ = +0.6372 and +0.0810 — are
> **unchanged**.


Everything the study produces is done and verified. What remains is writing and packaging,
and this is the list of it. Ticked items were checked by a command, not by recollection —
the command is named so you can re-run it.

---

## Verified and closed — do not redo

| | check | how |
|---|---|---|
| ✅ | Results frozen and hash-verified | `python3 scripts/archive_run.py --verify 2026-09-23T1221` → 242 files INTACT |
| ✅ | Pre-registration holds | OSF `dm2t8`, 2026-09-10; constraint hash `2d1fb47c…` matches the recorded hash |
| ✅ | Every ACS in every document matches the current run | `scripts/check_doc_numbers.py --against-current --strict` |
| ✅ | Every headline statistic matches its artefact field | `scripts/check_doc_statistics.py --strict` |
| ✅ | Every correlation, p-value and CI independently recomputed | `docs/STATISTICAL_VERIFICATION.md` |
| ✅ | p-values checked against exact/permutation alternatives | no difference exceeded 0.004 |
| ✅ | Data availability: every public dataset listed with accession, version and checksum | `docs/DATA_INVENTORY.md`, Table S7; 0 open problems on 2026-10-01 |
| ✅ | All 51 references verified | PDF on disk, or resolved against Crossref (35 on 2026-09-28; [36]–[44] and the nine former `cited` entries on 2026-10-01; [45]–[49] on 2026-10-02; [50]–[51] on 2026-10-03). [10] GBmap was wrong until 2026-10-01 — corrected, D25 |
| ✅ | No T/B label transposition behind the headline | 9 of 9 B markers in `B_cell`; 0 T markers there |
| ✅ | Figures render from artefacts, objective titles, published method names | `docs/figures/`, 11 figures in `docs/FIGURES.md` (Figure 11 added 2026-10-02 for Extension E2), each PDF + PNG |
| ✅ | Supplementary tables built | `docs/supplementary/`, 7 × CSV + HTML + `all_tables.tex` |

---

## The numbers, final

Quote these; they are the current run and nothing else is.

| | |
|---|---|
| cohort | 122 anatomic samples, 9 evaluable tumours (Ivy GAP) |
| panel | 15 estimators + 2 negative controls |
| best ACS | `music` 0.9692 |
| control separation | 0.7077 vs 0.4000 = **+0.3077** on fully-scored methods |
| registered outcome | ρ = **+0.6372**, p = 0.0143 — **meets the bar** |
| orthogonal truth | ρ = **+0.0810**, p = 0.80 — **INCONCLUSIVE**, not refuted |
| **headline** | **0 of 12 (GBM) and 0 of 12 (LGG)** reproduce T > B |
| truth strength | T > B in **94.8%** (GBM) / **93.2%** (LGG) of samples | *(EpiDISH's per-sample split; not a validated per-sample truth. A second methylation method agrees in GBM and puts B above T in LGG -- §4.4)*

---

## What is left, and it is all yours

### 1 · Write the prose

`docs/MANUSCRIPT.md` has every section, every number and every figure placement. The `> WRITE:`
blocks are beats and facts, not sentences to keep. Sections in order: Abstract (complete,
edit for voice only), I. Introduction, II. Methods, III. Statistical Analysis, IV. Results,
V. Discussion, Limitations, Acknowledgements, References.

### 2 · Three things to get right, because a reviewer will check them

- **Say INCONCLUSIVE and mean it.** The ACS-vs-accuracy result is underpowered at n = 12, not
  refuted. Claiming refutation is the easiest way to lose the paper.
- **quanTIseq could not be evaluated**, it did not perform badly — 15 of 57 constraint pairs,
  p = 0.31 against its own null. Saying otherwise understates your own control margin.
- **Do not claim a mechanism for the lymphoid failure.** Five were tested and rejected; two
  suspects are named in `docs/WHY_B_OVER_T.md` as properties of the reference, not causes.

### 3 · Cite a cost figure you can source

The Introduction's equity argument needs a per-sample single-cell sequencing cost. **This
project holds no costing artefact**, so that number must come from a source you cite. Do not
assert "expensive" without one.

### 4 · Acknowledgements

Supervisors; the Allen Institute for Brain Science (Ivy GAP); TCGA and its contributing
patients; the authors of every package in the panel — this is a test *of* their software and
they are owed the courtesy. **Disclose AI assistance in the form the venue requires.** Both
CJSJ and Regeneron STS ask, and it is far better volunteered than discovered.

### 5 · Decide the figure budget

Seven figures exist. If the venue caps it, Figures 4 and 5 (the purity scatters) are the
natural move to supplementary — they support Result 2, which is not the headline.

---

## Known limits, to state rather than be caught on

| limit | say it as |
|---|---|
| One tumour type | glioma only; GBM discovery, LGG registered replication |
| Twelve comparable methods | the power ceiling; no reanalysis widens it |
| ReCIDE absent | published during this study, not evaluated |
| No deep-learning methods | Scaden and DAISM-DNN train on simulated bulk, compounding a disclosed limitation |
| Two genuine R packages fell back | DWLS at 2,593 s vs a 2,400 s budget; the leaderboard carries the reimplementation's 0.7231 where the package scores 0.7846 |
| Full-database arm is on the older matrix | labelled in place; section 6 of `RESULTS.md` says so |

---

## If something needs re-running

```
python3 -m pytest tests/ -q                       # 334 tests, ~6 min on an idle machine
python3 scripts/check_doc_numbers.py --against-current --strict
python3 scripts/check_doc_statistics.py --strict
python3 scripts/summarize_results.py --check
python3 scripts/archive_run.py --verify 2026-09-23T1221
python3 scripts/independent_verification.py       # recomputes 4 headline quantities, no project code
```

Regenerate derived documents after any artefact change:

```
python3 scripts/summarize_results.py --write
python3 scripts/build_manuscript.py
python3 scripts/build_paper_figures.py && python3 scripts/build_stats_figures.py
python3 scripts/build_figure_index.py && python3 scripts/build_supplementary.py
```

**Do not hand-edit generated documents** — `MANUSCRIPT.md`, `RESULTS.md`, `FIGURES.md`,
`SUPPLEMENTARY.md` are all rebuilt from artefacts and your edits would be lost. Edit the
generator, or copy the text out into your manuscript file and edit it there.
