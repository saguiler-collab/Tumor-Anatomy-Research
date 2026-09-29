# Submission checklist

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
| ✅ | All 35 references verified | PDF on disk, or resolved against Crossref 2026-09-28 |
| ✅ | No T/B label transposition behind the headline | 9 of 9 B markers in `B_cell`; 0 T markers there |
| ✅ | Figures render from artefacts, objective titles, published method names | `docs/figures/`, 7 × PDF + PNG |
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
| truth strength | T > B in **94.8%** (GBM) / **93.2%** (LGG) of samples |

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
