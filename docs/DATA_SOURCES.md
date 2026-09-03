# Data sources

## Ivy GAP — the bulk cohort

**Ivy Glioblastoma Atlas Project**, Allen Institute for Brain Science.
<https://glioblastoma.alleninstitute.org/>

Openly downloadable. No application, no credential, no data-use agreement — a real
practical advantage over TCGA's controlled tiers, and one reason this cohort is a
better basis for a reproducible method comparison.

> Puchalski RB, Shah N, Miller J, et al. *An anatomic transcriptional atlas of human
> glioblastoma.* Science 360:660-663 (2018). doi:10.1126/science.aaf2666

### What gets downloaded

`python -m ivygap.data.download_ivygap`

| File | Contents | Size |
|---|---|---|
| `gene_expression_matrix_2014-11-25.zip` | `fpkm_table.csv`, `rows-genes.csv`, `columns-samples.csv` | ~50 MB |
| `tumor_details.csv` | per-tumour clinical table, carries `survival_days` | 2.9 KB |
| `rna_seq_samples_details.csv` | per-RNA-seq-sample metadata: study, structure, and both ids | 101 KB |

The downloader tries several URLs per file, because the Allen Institute has moved these
between well-known-file ids across releases. The expression archive is verified to
contain all three expected members before extraction.

**The working clinical URL, and how it was found (2026-09-03).** Every `tumor_details.csv`
URL previously recorded in this project is dead, and one fails in the worst possible way:

```
.../well_known_file_download/305873922   -> HTTP 404               (honest)
.../static/download/tumor_details.csv    -> HTTP 200, text/csv, HTML body
```

The second is what `_assert_looks_like_csv` exists to catch. The fix was to stop guessing
well-known-file ids and read the portal's own download page, which links an
`/api/v2/gbm/` namespace documented nowhere else:

```
https://glioblastoma.alleninstitute.org/api/v2/gbm/tumor_details.csv
https://glioblastoma.alleninstitute.org/api/v2/gbm/rna_seq_samples_details.csv
```

**`rna_seq_samples_details.csv` is not optional convenience.** It is the only published
file containing both `donor_id` and `tumor_id`:

```
tumor_details.csv    keyed on donor_id   (12111, 14734, ...)
columns-samples.csv  keyed on tumor_id   (711547, 705757, ...)
```

Neither contains the other's key. Without this third file the clinical table shares no
key with the expression matrix and cannot be joined to it at all — which is why survival
was previously BLOCKED in principle, not merely for want of a download.

### What makes this cohort different — and the trap inside it

**Verified against the 2014-11-25 release, not assumed.** The archive's 270 RNA-seq
samples are **two different studies sharing one file**, and telling them apart is the
single most important thing the loader does.

| Study | Samples | Tumors | How structure was assigned | Use here |
|---|---|---|---|---|
| **Anatomic structures** | **122** | **10** | H&E histology, by a neuropathologist, with no reference to RNA | **The instrument.** |
| Cancer-stem-cell clusters | 148 | 34 | ISH with an 18-probe reference set — i.e. *using expression* | Excluded. |

`structure_abbreviation` distinguishes them by suffix:

```
CT-reference-histology       30      <- anatomic study
CTmvp-reference-histology    25      <- anatomic study
CTpan-reference-histology    24      <- anatomic study
IT-reference-histology       24      <- anatomic study
LE-reference-histology       19      <- anatomic study      total 122
CT-CD44, CTpnz-PI3, CT-control-TGFBR2, CThbv-POSTN, ...     total 148
```

**The trap.** A prefix match on the acronym pulls all 148 cluster samples into the
anatomic cohort — `CT-CD44` starts with `CT` — roughly tripling the apparent sample size
with samples whose structure labels were assigned *from expression data*. That would
make the Anatomy Test circular in exactly the way it exists to avoid, while appearing to
strengthen it. `config.is_anatomic_sample` requires the `-reference-histology` suffix,
and `test_ish_cluster_samples_are_not_mistaken_for_anatomic_ones` guards it.

**Structure availability per tumor is uneven, and it bounds every constraint.** All 10
tumors contributed CT — which is why CT is the comparator — but:

| Structure | Tumors contributing |
|---|---|
| CT | 10 / 10 |
| MVP | 9 / 10 |
| PAN | 8 / 10 |
| IT | 8 / 10 |
| LE | 8 / 10 |

So C5 (PAN > LE, macrophage) is evaluable on only the tumors contributing *both*. The
per-constraint table reports `n_tumors_evaluable` for exactly this reason: a constraint
resting on six tumors is a weaker claim than one resting on ten, and the difference must
be visible rather than buried in a pooled score.

### This bounds ACS before any method runs

Maximum evaluable tumours per constraint, measured from the structure availability above:

| C1 | C2 | C3 | C4 | C5 | C6 | C7 |
|---|---|---|---|---|---|---|
| CT>LE | LE>CT | MVP>CT | MVP max | PAN>LE | MVP>CT | LE<IT<CT |
| 8 | 8 | 9 | 9 | **6** | 9 | 8 |

Maximum weighted denominator: **65.0**, across **9 distinct tumours** (tumour 710262
supplied only CT and PAN, so no constraint is evaluable for it at all).

No deconvolution method, hyperparameter, gene set or reimplementation raises any of these
numbers, and neither does deconvolving the other 148 samples — their structure labels are
not admissible evidence. Ivy GAP has ten H&E-annotated tumours and that is all it has.
Widening ACS's evidence base requires a second tissue, which is protocol step 7.

### What the other 148 samples are used for

They are deconvolved (`--full-database`) and never scored. Two things follow:

| | anatomic only | full database |
|---|---|---|
| samples deconvolved | 122 | **270** |
| tumours with composition | 10 | **37** |
| tumours with composition *and* a recorded survival time | 7 | **29** |
| samples contributing to ACS | 122 | 122 |
| tumours contributing to ACS | 10 | 10 |

The prognostic cohort more than quadruples; the ACS evidence base does not move by one
sample. That asymmetry is the whole reason the two sets are kept separate in code.

### A discrepancy worth flagging

The protocol document cites **41 tumors**. The RNA-seq release actually contains **37**
(10 anatomic + 34 cluster, with 7 tumors contributing to both). The 122 / 10 figure for
the anatomic study matches exactly. Verify the 41 against the current portal before
citing it — it may refer to the whole project including ISH-only tumors that contributed
no RNA-seq.

### Tolerant parsing, intolerant validation

Column names have shifted across releases (`structure_acronym` vs
`structure_abbreviation`, `tumor_id` vs `donor_id`, survival in days vs months). Each
field is resolved from a list of accepted names and fails with a message naming what was
actually present, rather than pinning one spelling or falling back to positional
indexing.

Validation is not tolerant:

- A `patient_id` column that is unique per sample is **rejected** — that means a
  sample-level identifier was resolved and the nested design would be silently lost.
- Duplicate sample or patient ids are rejected.
- Survival times ≤ 0 are dropped; an unrecognised vital-status value becomes NaN and the
  row is dropped, never defaulted to "censored" — defaulting would convert deaths into
  survivors and bias every estimate in one direction.
- A clinical table with zero events aborts.

Counts are **measured and recorded**, never asserted against hard-coded literals.

### Protocol step 3's gate is closed against the portal, not against ourselves

The loader derives its cohort description from `columns-samples.csv`, which ships inside
the archive it is describing. That is an internal count: a parsing error that
mis-assigned structures would produce a self-consistent and entirely wrong answer.

`ivygap/data/portal_metadata.py` checks it against `rna_seq_samples_details.csv`, which
the Allen Institute exports from their LIMS independently of this pipeline, carries the
study assignment as an explicit `study_name` column rather than as a suffix convention we
have to parse, and names every sample's structure. Measured, not asserted:

| | archive parse | portal table |
|---|---|---|
| anatomic samples | 122 | 122 |
| CT / IT / LE / MVP / PAN | 30 / 24 / 19 / 25 / 24 | 30 / 24 / 19 / 25 / 24 |
| study assignment | — | agrees on all 270 shared samples |
| samples in archive but not portal | **0** | — |

The portal currently lists **279** RNA-seq samples against the archive's 270. The nine
extras are all Cancer Stem Cells samples added after the 2014-11-25 freeze; that is
release skew and is reported as `portal_only`, not as a count mismatch. The reverse — a
sample in the archive the portal does not describe — is the case that can only be a
defect, and it hard-fails.

Two negative controls guard the reconciler itself
(`tests/test_full_database.py`): planting a corrupted structure label must produce
`MISMATCH`, and planting a corrupted study assignment must produce `DEFECT`. A
reconciler that always says RECONCILED is decoration.

### The clinical table is available. The event indicator is not.

`tumor_details.csv` downloads cleanly from the `/api/v2/gbm/` endpoint and joins to the
expression matrix through `rna_seq_samples_details.csv`. What it contains:

| | |
|---|---|
| tumours | 42 (37 with RNA-seq in the archive) |
| `survival_days` recorded | 32 of 42 · **29 of the 37 with RNA-seq** |
| vital-status / censoring column | **none — the Allen Institute publishes none** |
| other fields | `molecular_subtype`, `extent_of_resection`, `surgery`, `mgmt_methylation`, `egfr_amplification`, `initial_kps`, `age_in_years` |

**Survival is therefore still BLOCKED by default, for a different and more specific
reason than before.** A survival time without a censoring indicator cannot produce a
C-index: treating every recorded time as a death invents 100% mortality, and treating the
blanks as censored invents follow-up durations nobody published.
`clinical.load_clinical()` raises rather than choosing either.

**The blanks are not missing at random**, and this is measured rather than assumed:

| | recorded `survival_days` | blank |
|---|---|---|
| MGMT methylated | 9 / 31 | **8 / 9** |
| median age | 61 | 52 |

Fisher exact p = **0.0021** for MGMT; Mann-Whitney p = **0.015** for age. MGMT
methylation and younger age are the two strongest favourable prognostic factors in
glioblastoma, so the blanks are concentrated among the patients most likely to have been
alive at the data freeze. A complete-case analysis is biased toward short survival by
construction. `clinical.audit_missingness()` computes this and it is written to
`results/clinical_missingness_audit.json` on every run.

**The declared-assumption escape hatch.** `event_policy="observed-time-is-death"` reads a
recorded time as an observed death and drops the blanks — never censors them, never
imputes a time for them. It must be named by the caller
(`--assume-recorded-survival-is-death`), it is stamped into every artefact it touches,
and its output goes to `results/survival/declared_event_policy/` so it can never be
mistaken for the canonical result. It does not rescue the analysis: 29 events against the
80 the design needs leaves the comparison INCONCLUSIVE regardless.

**`age_in_years` ships as `"61 yrs"`.** A plain `to_numeric` yields all-NaN and silently
drops the only baseline covariate the survival model has. It is parsed with an explicit
digit extraction, and `test_age_survives_ivygap_string_formatting` guards it.

---

## Single-cell reference

Reference-based deconvolution needs cell-type expression profiles. This project builds
them from a GBM single-cell atlas rather than shipping a fixed signature matrix, because
MuSiC, Bisque and SCDC need statistics that a collapsed signature matrix has already
averaged away.

**Primary — GBmap.** Place at `data/reference/gbmap_core.h5ad`. Annotations are read
from `annotation_level_3` and collapsed onto the project roster via
`config.GBMAP_CELL_TYPE_MAP`, carried over verbatim from the TCGA project so the
reference is built identically across the two cohorts.

**Second reference — for SCDC ENSEMBLE.** ENSEMBLE requires more than one reference to
do anything; with one it degenerates to plain SCDC. A second independent GBM atlas
(e.g. Neftel et al. 2019) makes cross-reference weighting real and makes "robust to
reference choice" a *measurable* property rather than an assumption.

### What the reference bundle carries

Built once, in `data/reference.py`, so no two methods can receive subtly different
references:

| Field | Used by | Why it cannot be recomputed later |
|---|---|---|
| `profile` | all | the signature matrix S |
| `sigma` | MuSiC | **cross-donor** variance per gene per type — destroyed by averaging |
| `cell_size` | all | mean transcript content per type; the RNA→cell conversion |
| `donor_profiles` | Bisque, SCDC | per-donor profiles for the assay transform |

**Donor balance.** Cells are capped per donor per cell type before averaging, and the
signature is the mean *across donor means*, not across cells. Without the cap one deeply
sequenced donor supplies most cells of a common type and the "reference" becomes that
donor's profile — which then inflates the apparent cross-donor consistency MuSiC is
trying to measure, exactly backwards.

**Structural zeros are masked.** A donor contributing no cells of a type has a
structural zero there, not a measured zero. Including it would drag the mean toward zero
and inflate the variance.

### Cell-level export for the R packages

When real single-cell data is present, `export_cell_level` writes
`sc_counts_<name>.csv.gz` and `sc_meta_<name>.csv`. Without them `r_bridge` reports
itself unavailable and the pipeline falls back to the Python reimplementations,
recording the reason per tool.

---

## The synthetic fixture path

`build_synthetic` generates a reference with known generative structure: per-type gene
blocks, a shared housekeeping background, donor-level offsets (so cross-donor variance
is real and MuSiC's weighting has something genuine to find), and differing per-type
transcript content (so the cell-size correction is not a no-op).

`scripts/run_all.py --synthetic` runs the complete pipeline on it — including an
Ivy-GAP-*shaped* bulk cohort with uneven block counts and occasionally missing
structures, because a tidy balanced fixture would never exercise the patient-equal
aggregation the real data demands.

**Synthetic results are labelled as such everywhere they touch disk**, including
`release/manifest.json` (`run_type: SYNTHETIC FIXTURE`) and a warning banner in
`release/RELEASE_INDEX.md`. They validate that the pipeline runs and that its scoring
behaves correctly. They are not biological findings.

---

## What Ivy GAP does not have, and what replaces it

| TCGA had | Ivy GAP | Consequence |
|---|---|---|
| ABSOLUTE DNA-based tumour purity | — | orthogonal purity validation is not available |
| Methylation leukocyte fraction | — | orthogonal immune validation is not available |
| 115 patients / 93 events | ~40 tumours | survival is bounded and **not used to rank methods** |
| — | **anatomic structure labels** | histology becomes the ground truth on real tissue |

The trade is deliberate. TCGA's orthogonal assays validate *two* quantities (purity,
leukocyte fraction) on one sample per patient. Ivy GAP's histology validates the
*direction* of several cell types across five niches, paired within patient — and every
method is scored on the identical set of directional questions.
