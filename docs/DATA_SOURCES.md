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

---

## Datasets identified but not yet used

Researched 2026-09-10. None of these has been downloaded or read; each is recorded here
with what it would buy and what it would cost, so the decision is made deliberately
rather than by whichever file arrives first. Citations are in `sources/SOURCES.md`, all
DOIs Crossref-verified.

### Siletti et al. 2023 — Human Brain Cell Atlas v1.0

**Would supply:** the single-cell reference the second tissue needs.

The GBM roster has no neuron column and one cannot be added from GBmap, which contains
**22 neurons** (`results/anatomic/reference_coverage.json`). Normal cortex is mostly
neurons, so a normal-brain run is impossible without a different reference. This is that
reference.

| | |
|---|---|
| Source | CELLxGENE collection `283d65eb-dd53-496d-adb7-7570c7caa443` |
| Publication | Siletti K, et al. *Science* (2023) |
| All neurons | 2,480,956 cells · 32.9 GB h5ad |
| All non-neuronal | 888,263 cells · 4.7 GB |
| Per-supercluster | e.g. Oligodendrocyte 490,246 cells · 2.5 GB |
| Assay | **10x 3′ v3** |

The assay line matters more than the cell counts: it is the same platform family as
GBmap, so this project's CPM normalisation, its gene-space handling and its B-mode/S-mode
reasoning all carry over unchanged. The per-supercluster files also mean the 33 GB neuron
archive is not required — a donor-balanced reference can be built from the smaller ones.

### Allen Human Brain Atlas — bulk expression by region

**Would supply:** the second tissue's bulk arm. Fetcher written and tested:
`scripts/fetch_allen.py`. No API key required.

| product | donors | size | platform | weighted denominator at 7/donor |
|---|---|---|---|---|
| Microarray | **6** | 1.57 GB | microarray | **42** |
| RNA-seq | 2 | 44 MB | RNA-seq | 14 |

Ivy GAP's weighted denominator is 65, so the RNA-seq product would be *materially weaker
than the cohort that already exists* and would not answer the question a second tissue is
for. The microarray product is comparable in power but crosses a platform boundary. The
trade is a study-design decision and is recorded in `docs/SECOND_TISSUE.md` rather than
settled here.

**Structure ontology.** Region names come from `graph_id 10` ("Human Brain Atlas", 1,839
structures), whose depth-1 split is `Br / GM | WM | SS` — the grey/white distinction the
brain constraints turn on is the ontology's own, not one imposed on it. Note that
`ontology_id 1` is the **mouse** graph; querying it returns mouse cortical labels, which
is a mistake this project made once and fixed.

### Neftel et al. 2019 — GBM single-cell atlas

**Would supply:** a second GBM reference, and with it the single highest-value fix
available short of a second tissue.

`scdc_ensemble` currently runs **degenerate** — "one reference supplied ('gbmap'), so
there is nothing to weight across: SCDC ENSEMBLE reduces exactly to SCDC." A published
tool is running with its entire contribution switched off, disclosed but crippled.

A second reference would make SCDC ENSEMBLE a genuine method, and would additionally
support a robustness claim this project cannot currently make: **is the ACS ranking stable
across references, not merely across cohorts?** The cohort-sensitivity check already
exists (`acs_cohort_sensitivity.csv`); a reference-sensitivity check does not, for want of
a second reference.

doi:10.1016/j.cell.2019.06.024 — already cited in this file for cross-reference weighting.

### CELLxGENE

The platform both single-cell atlases are obtained through. Its `var` is indexed by
**Ensembl id** where this project uses HGNC symbols, which is why gene-symbol mapping is
explicit in `reference.py` rather than assumed.

Megill C, et al. *bioRxiv* (2021), doi:10.1101/2021.04.05.438318.

### Deliberately not pursued

**More least-squares deconvolution methods.** The 13 comparable methods already produce
only **8 distinct ACS values**, with four tied at 0.9846 (NNLS, SVR, Elastic Net, EPIC) on
a scale with 66 possible values. Another variant of the same family would land in that tie
group and reduce the rank correlation's resolution rather than improve it.

Methods that differ in *kind* would add information — a reference-free method such as
CDSeq, which uses no signature matrix at all, or a deep-learning method such as Scaden.
Those would spread across the accuracy range instead of clustering.

**The hosted CIBERSORTx service.** Web- and licence-gated. Both modes are implemented
here from the paper and its Supplementary Note 1, and no result in this project comes
from Stanford's servers.

**ABSOLUTE purity for TCGA-GBM.** Still wanted as an orthogonal yardstick that never
touches RNA (`docs/EXTERNAL_ACTIONS.md` item 3, Carter et al. 2012, doi:10.1038/nbt.2203).
Not obtained.

---

## Second and third GBM references — obtained 2026-09-10, not yet ingested

Both were downloaded to `pipeline packages / repos/SCDC/`. Inspected by metadata only; no
expression matrix has been read into a run.

### Neftel et al. 2019 (GSE131928) — the usable second reference

| file | size | content |
|---|---|---|
| `GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv` | 1.5 GB | **16,201 cells · 9 patients (21 samples) · 10x** |
| `GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv` | 944 MB | SMART-Seq2 version |
| `GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx` | 848 KB | sample metadata |

**Take the 10x file.** It is platform-matched to GBmap (10x), so pairing them in SCDC
ENSEMBLE tests *reference choice* rather than confounding it with *platform*. The
SMART-Seq2 file is the same study on a different chemistry and would introduce that
confound; it is worth keeping for a deliberate cross-platform experiment, not for D2.

**CORRECTED 2026-09-12 — it is 9 patients, not 21.** An earlier entry here said 21
donors, counted from the cell-ID prefix before the final underscore. That over-counts:
the authoritative `tumour name` column in the series xlsx gives **9 tumours**, and the
prefix splits repeat samples of one patient — `105A`, `105_B1`, `105_B2`, `105_C1`,
`105_C2`, `105_D1`, `105_D2` are all MGH105.

| tumour | cells | | tumour | cells |
|---|---|---|---|---|
| MGH105 | 5,513 | | MGH115 | 1,283 |
| MGH124 | 2,415 | | MGH118 | 539 |
| MGH143 | 2,314 | | MGH114 | 473 |
| MGH102 | 1,822 | | MGH126 | 229 |
| MGH125 | 1,613 | | | |

21 is the number of **samples**, which is what SCDC's `sample` argument wants; 9 is the
number of **patients**, which is what a "donors" figure means. Both must be reported with
the right label — cross-subject variance is estimated across patients, not across repeat
samples of the same patient, so calling this a 21-donor reference would overstate its
statistical footing by more than a factor of two.

All 16,201 10X cells are **adult**. The study also contains paediatric GBM, but none of it
is in the 10X subset, so mixing paediatric tissue into an adult reference is not a risk
here. It would be for the Smart-seq2 subset, which is not the one to use.

### Mossi Albiach et al. 2023 — NOT usable as a second reference

`d45b4ce6-9725-4d79-b97a-70a44158bdbf.h5ad`, 1.41 GB, from CELLxGENE collection
`113a558a-e96e-4643-81db-140e95c58578`.

135,482 cells and **one donor** (`SL040`, right temporal lobe, 10x 3' v3). Eight times
Neftel's cell count and none of its structure.

That rules it out for the purpose it was fetched for. SCDC and MuSiC exist to model
cross-subject variance; a single-donor reference cannot supply it, and SCDC_prop run on
this alone would be degenerate for the same reason `scdc_ensemble` currently is. Pairing
it with GBmap in an ENSEMBLE would weight across two references of which one has no
subject-level variance to weight with.

**What it is good for instead.** Its `obs` carries `Zone`, `Location` and `NeftelClass` —
the study sampled tumour core through to macroscopically normal cortex, so cells are
annotated by position along that gradient. That is a spatially-resolved single-cell view
of the same axis the Ivy GAP constraints describe, and it could support an independent
check on the constraints themselves rather than on the methods. Recorded here as an
opportunity, not a plan.

### The blocker on Neftel: GSE131928 ships no cell-type labels

Verified 2026-09-10. The GEO release contains expression matrices and a submission
metadata template. It contains **no cell-type annotation**. The bundled reanalysis repo
(`scRNA_GBM_Neftel2019-main`) confirms this by construction — its pipeline *derives* the
labels: QC and doublet filtering, Harmony batch correction, clustering, marker-gene
scoring, inferCNV, then malignant classification and Neftel-state scoring.

A deconvolution reference is a matrix of per-cell-type expression profiles. Without
labels there is nothing to profile, so the 1.5 GB download is not yet usable.

**Confirmed 2026-09-12 by checking all four GEO supplementary files.** The series ships
`GSE131928_RAW.tar`, the per-cell xlsx, and `filelist.txt`. The xlsx holds 24,131 per-cell
records with columns `Sample name, title, source name, organism, molecule, processed data
file, instrument model, tumour name, adult/paediatric` — **no cell-type column**. The
labels are produced by the authors' analysis, not distributed with the series.

**Two ways forward, and they are not equivalent.**

*Preferred — obtain the authors' annotations.* Neftel's annotated data is distributed
through the Broad Single Cell Portal rather than through GEO. Author-supplied labels keep
the reference independent of this project.

*Not preferred — annotate it here.* Running the reanalysis pipeline would work, but it
would put **this project's clustering and marker-threshold choices inside a reference
used to judge deconvolution methods**. That is not circular in the ACS sense — the
constraints are still independent — but it means the "second reference" would carry our
decisions rather than being an independent check on them, which is most of the value a
second reference was wanted for. If it is done, it must be declared prominently, and the
annotation choices recorded with the same discipline as the constraint file.

**Consequence for D2.** `scdc_ensemble` cannot be un-degenerated until annotated cell
types exist for a second multi-donor reference. Neither file currently on disk satisfies
that: Neftel has 9 patients (21 samples) and no labels; Albiach has labels and one
donor.

---

## Darmanis 2017 (GSE84465) — fetched 2026-09-14, two distinct uses

**Source.** GEO GSE84465, Darmanis et al. 2017, *Single-cell RNA-seq analysis of infiltrating
neoplastic cells at the migrating front of human glioblastoma*. `GSE84465_GBM_All_data.csv.gz`
(20.5 MB, 23,465 genes x 3,589 cells, raw counts, space-delimited despite the `.csv`) and
`GSE84465_series_matrix.txt.gz`. Both under `data/raw/darmanis_2017/`, with
`cell_metadata.csv` extracted from the series matrix: cell type, tissue, patient, plate, well
and **sorting gate**.

**What it contains.** 4 patients, tumour core and periphery, 7 author-annotated cell types.
Labels come straight from the GEO metadata (`characteristics_ch1.6`) — no derivation, no
clustering, no inferCNV.

### Use 1 — the only cross-patient constraint test in this project

`scripts/darmanis_constraint_check.py`, rendered in `RESULTS.md` §5c. **C1's direction only**,
tested within a fixed FACS gate so the selection bias is held constant. 4 of 4 scored gates
support it, two at p = 0.0001 with 3/3 patients each.

**It is NOT a composition test.** Only 665 of 3,589 cells are `Unpanned`, and the unpanned
periphery is **13 cells** — twelve from BT_S4, one from BT_S6, none from BT_S1 or BT_S2. A
composition over sorted cells measures the sort. This was my own initial recommendation and it
was wrong; see `docs/EXTERNAL_ACTIONS.md` item 15.

### Use 2 — a second deconvolution reference

`scripts/build_darmanis_reference.py` -> `data/reference/darmanis_2017/`
(profile.csv, sigma.csv, cell_size.csv, provenance.json). Panning ruins composition and is
**irrelevant to a reference**, which needs per-type profiles rather than per-region
abundances; sorting helps by enriching rare types.

| roster type | from Darmanis label | cells | donors |
|---|---|---|---|
| Tumor | Neoplastic | 1,091 | 4 |
| Macrophage_Microglia | Immune cell | 1,847 | 4 |
| Astrocyte | Astocyte *(GEO's spelling)* | 88 | 4 |
| Oligodendrocyte | Oligodendrocyte | 85 | 4 |
| Endothelial | Vascular | 51 | 4 |

22,799 genes after restricting to the genes the Ivy GAP bulk carries. `cell_size` spread
**1.63** — computed from RAW library sizes captured before normalisation, so unlike GBmap's
this reference can actually drive a cell-size conversion (see D12).

**Validated on marker placement:** PTPRC, CD68, MBP, PLP1, PECAM1, VWF and EGFR all peak in
the expected column. **GFAP peaks in Tumor, not Astrocyte** — which is not a defect but
independent corroboration of the invariant that Astrocyte is unidentifiable against Tumor in
glioblastoma, where the tumour cells are astrocytic in lineage.

**Limits, all structural.** `T_cell`, `NK_cell` and `B_cell` are **not recoverable** —
Darmanis pools every lymphoid and myeloid cell into one `Immune cell` label, so mapping it onto
Macrophage_Microglia is an approximation that includes T cells. `OPC` (406 cells) and `Neuron`
(21) are EXCLUDED rather than folded in; putting OPC into Oligodendrocyte would inflate exactly
the type C2 is about. Astrocyte, Oligodendrocyte and Vascular are thin (88, 85, 51) and thinner
per donor. It is Smart-seq2 where GBmap is 87% 10x — a platform confound, and also the regime
CIBERSORTx's S-mode exists for.

**So it supports a DECLARED 5-type sub-roster only** and cannot replace GBmap on the full
eight. What it makes possible: `SCDC ENSEMBLE` non-degenerate for the first time (D2), and a
reference-sensitivity check that separates "how these methods behave" from "how these methods
behave *on GBmap*".

---

## Independent confirmation of the bulk matrix — 2026-09-15

The Ivy GAP FPKM table is published twice, four years apart, by different routes. Both were
compared here rather than assumed equivalent:

| | lines | SHA-256 (first 16) |
|---|---|---|
| Allen Institute archive, `gene_expression_matrix_2014-11-25.zip` | 25,874 | `72a97ea792a81098` |
| **GEO GSE107559**, `GSE107559_ivygap_fpkm_table.csv.gz` | 25,874 | `72a97ea792a81098` |

**Bit-identical.** The input this study deconvolves is the same matrix the Ivy GAP authors
deposited alongside Puchalski et al., *Science* 2018;360(6389):660–663 (PMID 29748285), so the
provenance of the bulk does not rest on one download.

The same record also establishes what is **not** available: *"Raw data not provided for this
record"*, and *"The raw RNA-Seq and SNP array data will be submitted to dbGaP."* There are no
public read counts for Ivy GAP by either route. See `CORRECTIONS_REGISTRATION.md` **C9**.

---

## Read counts — obtained 2026-09-15, and they close the provenance chain

The published Ivy GAP matrix is normalised FPKM. The **counts underneath it are public**, on the
Allen Institute portal's download page rather than in GEO, as one RSEM `genes.results` file per
sample:

```
gene_id  transcript_id(s)  length  effective_length  expected_count  TPM  FPKM
1        NM_130786_3       1766.00  1604.54          45.56           3.41  2.28
```

25,873 genes, ~1.5 MB each, 270 samples. Fetched by `scripts/fetch_ivygap_counts.py` into
`data/raw/ivygap_counts/` with a per-file SHA-256 recorded. Assembled by
`scripts/build_ivygap_counts_matrix.py`.

The portal also publishes **anonymized BAMs**, unauthenticated, 414 MB each (~112 GB for 270),
manifest in `FPKM/bam_manifest.csv`. Not used: RSEM has already counted, and recounting would
substitute this pipeline's choices for the authors'. They remain available if a read-level
question ever needs them.

### What these files establish about the published matrix

Mapping RSEM's **Entrez** ids onto Ivy GAP's internal **`gene_id`** through `rows-genes.csv`:

| check | result |
|---|---|
| genes mapped / shared with the published matrix | 25,873 / 25,873 |
| `corr(published FPKM, per-sample raw FPKM)` | **1.00000000** |
| published ÷ raw, across genes within a sample | **one constant**, CV **1.4e-15** |

The published matrix is exactly **raw FPKM × one scale factor per sample** — the authors'
normalisation "based on genes not enriched in particular anatomic structures". So:

1. these per-sample files are provably the **source** of the matrix this study deconvolves;
2. that normalisation is a **no-op for this pipeline**, which rescales every bulk column to sum
   to 1e6 before deconvolving, and a per-sample constant is removed exactly by that step. This
   had been an untested assumption.

**The id spaces are both numeric, and a position-based join looks like it works.** It shares
7,741 of 25,873 ids by coincidence and gives `corr = 0.12` while producing a full-looking
matrix. `build_ivygap_counts_matrix.py` therefore gates on the two checks above for **every**
sample and writes nothing if any fails.

**`expected_count` is RSEM's posterior expectation and is fractional** (45.56, not 46). It is
not a raw integer count. Rounding, where a model needs integers, is declared at the point of use
rather than done silently in the loader.

---

# Vendored third-party resources

Everything under `pipeline_packages /` is someone else's work, kept locally so a run is
reproducible without a network. Nothing there is modified. This section records what each item
is, where it came from, what it is used for here, and — where it matters — what it does **not**
contain.

## GBMdeconvoluteR (Ajaib et al.) — marker sets, and an external ground truth

**Source.** `pipeline_packages / repos/GBMDeconvoluteR/GBMDeconvoluteR-main`, the authors'
Shiny application, GPL-3. Live instance: <https://gbmdeconvoluter.leeds.ac.uk>.

**Paper.** Ajaib S, Lodha D, Pollock S, Hemmings G, Finetti MA, Gusnanto A, Chakrabarty A,
Ismail A, Wilson E, Varn FS, Hunter B, Filby A, Brockman AA, McDonald D, Verhaak RGW, Ihrie RA,
Stead LF. *"GBMdeconvoluteR accurately infers proportions of neoplastic and immune cell
populations from bulk glioblastoma transcriptomic data."* **Neuro-Oncology**
2023;25(7):1236–1248. PDF vendored alongside the repo.

**What the tool is.** MCPcounter given GBM-tissue-specific marker genes. Not a new algorithm —
a new marker set for an existing one, which is precisely why it is useful here.

**What the `data/` directory contains** — marker gene sets, all `.rds`:

| file | contents |
|---|---|
| `Ajaib_et_al_2022_GBM_Immune_markers.rds` | 183 genes over 8 immune populations: TAM 39, Microglia 49, Mast 26, NK 26, T 16, Monocytes 14, B 8, DC 5 |
| `Moreno_et_al_2022_lvl3_immune_markers.rds` | 784 genes over 16 populations (50 each, Plasma B 34), including Astrocyte, Endothelial, Oligodendrocyte, OPC, Neuron, Mural, Radial glial |
| `Moreno_et_al_2022_lvl3_neoplastic_markers.rds` | AC, MES, NPC, OPC |
| `Neftel_et_al_2019_four_state_neoplastic_markers.rds` | AC, MES, NPC, OPC |
| `Neftel_et_al_2019_all_neoplastic_markers.rds` | AC, G1S, G2M, MES1, MES2, NPC1, NPC2, OPC |
| `Wang_et_al_2017_GBM_TI_markers.rds` | 11,529 tumour-intrinsic gene symbols, used by the tool to filter neoplastic markers |
| `TGCA_GBM_example.rds` | 19,938 genes x 18 TCGA samples, the app's worked example |
| `plot_colors.rds`, `plot_order.rds` | presentation only |

**Single-cell datasets the markers were derived from**, per the paper's Table 1: **GSE141383**
(~18k cells, 5 primary IDHwt GBM), **GSE163120** (~21k primary + ~43k recurrent),
**GSE135437**, **GSE138794**, and **GSE131928** (Neftel, ~23k cells, used for their CIBERSORTx
reference profile). None of these are downloaded here; only the derived marker lists are used.

**WHAT IT DOES NOT CONTAIN, and this matters.** The repo holds **no imaging mass cytometry data
and no matched bulk RNA-seq**. The IMC validation cohort — ten IDHwt GBM samples, five paired
primary/recurrent, a 33-antibody panel, with bulk RNA-seq on the same tissue — is the authors'
own (Leeds Institute of Medical Research). The vendored PDF carries no data-availability
statement, and the per-sample IMC counts are not public through this route. So the IMC data
itself is **not** available to re-analyse; what is available is the authors' **published
verdict**, which is what this project uses.

**How it is used here.** Two ways, both recorded rather than assumed:

1. **As marker sets** for the genuine MCPcounter, in `scripts/imc_anchored_test.py`. Populations
   are rewritten into this project's eight-type roster and any population without a roster
   equivalent — monocytes, DC, mast cells, neutrophils, fibroblasts, mural, radial glial — is
   **dropped, never folded into a neighbouring type**.
2. **As an external ground truth**, via the paper's reported IMC correlations. See
   `prespecified/imc_anchored_prediction.md` and `CORRECTIONS_REGISTRATION.md` C10.

## MCPcounter

**Source.** <https://github.com/ebecht/MCPcounter>, installed from source, version **1.2.0**.
Its default marker file (`Signatures/genes.txt`, 111 genes over 10 populations) is cached at
`data/raw/mcp_default_genes.txt` so a run does not depend on GitHub being reachable.

**Why it is in this project.** MCPcounter returns **abundance scores, not proportions.** They do
not sum to one and are not comparable across cell types — only within a cell type across
samples. That is a property of the method, and it has a consequence worth stating plainly:
**ACS can score a method of this class and the pseudobulk accuracy arm structurally cannot**, so
the registered primary outcome could never have included one.

## CDSeq

**Source.** `pipeline_packages / repos/CDSeq/CDSeq_R_Package-master`, version **1.0.9**
(2021-03-12), Kang et al. Installed from source.

**Why it is here.** Reference-free deconvolution: it estimates cell types de novo and uses a
reference only to *label* them. That makes it the one instrument able to test whether the ACS
ordering is a property of the reference atlas (C4) without using a reference at all. It was
blocked until Ivy GAP read counts were obtained — see the retraction in C9.

**Measured cost.** The Gibbs sampler walks every **read**, not every gene, so `dilution_factor`
is the lever: 300 genes x 12 samples at `dilution_factor=1` is 53M reads x 150 iterations = 8e9
draws and did not finish in four minutes, while the same run at `dilution_factor=100` took
**11.7 s**.


---

# Reference list

Compiled 2026-09-20. **Every DOI below was extracted from a document in this repository** — from
PDF metadata or from the first page of the paper itself — or from this project's own provenance
records. None is written from recollection.

Entries still marked **[VERIFY]** are sources genuinely used whose citation is *not* recoverable
from anything on disk. They are left incomplete on purpose: a plausible-looking DOI that resolves
to the wrong paper is worse than a visible gap.

---

## A · Data sources

### A1 · Primary cohort

**Ivy Glioblastoma Atlas Project (Ivy GAP)** — Allen Institute for Brain Science, 2014-11-25
release. Laser-capture-microdissected bulk RNA-seq with anatomic structure annotations (LE, IT,
CT, MVP, PAN), plus in-situ hybridisation used for a validation that involves no deconvolution.
Supplies every ACS measurement.

> Puchalski RB, Shah N, Miller J, *et al.* An anatomic transcriptional atlas of human
> glioblastoma. *Science* **360**:660–663 (2018). doi:10.1126/science.aaf2666

### A2 · Validation cohorts — TCGA

| dataset | provenance | supplies |
|---|---|---|
| TCGA-GBM bulk expression | vendored from the predecessor project; delivered as log2(x+1), linearised here | tumour-purity yardstick, 154 samples |
| TCGA-LGG `star_counts` | UCSC Xena / GDC | cross-cohort replication, 510 samples |
| TCGA-GBM HumanMethylation450 | UCSC Xena, **TCGA Hub** | per-cell-type lymphoid truth, 155 samples |
| TCGA-LGG HumanMethylation450 | UCSC Xena, **TCGA Hub** | per-cell-type lymphoid truth, 530 samples |
| MC3 gene-level mutation calls (LGG) | TCGA Unified Ensemble "MC3" | IDH1 status as a candidate failure factor |
| ABSOLUTE tumour purity | PanCanAtlas | the DNA ground truth for tumour content — never touches RNA, which is what makes it independent |
| Leukocyte fraction from methylation | `TCGA_all_leuk_estimate.masked.20170107.tsv` | aggregate immune truth |
| CIBERSORT relative fractions | `TCGA.Kallisto.fullIDs.cibersort.relative.tsv` | comparison against a published deconvolution of the same cohort |

> **Both methylation matrices come from the same Xena hub, deliberately.** The GDC Hub serves
> these cohorts under a different pipeline and genome build; mixing hubs across cohorts would
> make a GBM-vs-LGG difference indistinguishable from a pipeline difference.

- Carter SL, Cibulskis K, Helman E, McKenna A, *et al.* Absolute quantification of somatic DNA
  alterations in human cancer. *Nature Biotechnology* **30**:413–421 (2012). doi:10.1038/nbt.2203
- Thorsson V, *et al.* The Immune Landscape of Cancer. *Immunity* (2018).
  doi:10.1016/j.immuni.2018.03.023
- MC3 mutation calls — Ellrott K, *et al.* **[VERIFY]**
- UCSC Xena — Goldman M, *et al.* **[VERIFY]**

### A3 · Single-cell reference atlases

| atlas | identifier | role |
|---|---|---|
| **GBmap** (primary) | CELLxGENE collection `283d65eb-dd53-496d-adb7-7570c7caa443`, `annotation_level_3`, 100 donors | the signature matrix every method solves against, in both its vendored-frozen and `raw/X`-rebuilt forms |
| Neftel | GSE131928 | reference-sensitivity arm |
| Darmanis | GSE84465 | reference-sensitivity arm |
| Albiach | — | reference-sensitivity arm, mesenchymal phenotype |
| Further accessions referenced in this file | GSE107559, GSE135437, GSE138794, GSE141383, GSE163120 | candidate references and cross-checks |

- *Bidirectional tumor-host interdependence in glioblastoma* — the GBmap source publication.
  *Cancer Cell* **40**:639 (2022). doi:10.1016/j.ccell.2022.05.009
- Megill C, *et al.* CELLxGENE: a performant, scalable exploration platform for high dimensional
  sparse matrices. *bioRxiv* (2021). doi:10.1101/2021.04.05.438318
- Neftel C, *et al.* An Integrative Model of Cellular States, Plasticity, and Genetics for
  Glioblastoma. *Cell* (2019). doi:10.1016/j.cell.2019.06.024
- Mossi Albiach A, Janusauskas J, Kjaer J, *et al.* Futile wound healing drives mesenchymal-like
  cell phenotypes in human glioblastoma. *bioRxiv* (2023). doi:10.1101/2023.09.01.555882
- Siletti K, *et al.* Transcriptomic diversity of cell types across the adult human brain.
  *Science* (2023). **[VERIFY]**
- Darmanis S, *et al.* Single-cell RNA-seq analysis of infiltrating neoplastic cells at the
  migrating front of human glioblastoma (2017). **[VERIFY]**

### A4 · Methylation deconvolution reference

**EpiDISH**, `centDHSbloodDMC.m` — 333 HM450 CpGs across 7 blood cell types, RPC mode. Supplies
the per-cell-type lymphoid truth (T / NK / B) in both cohorts.

Two limits, both measured rather than assumed: it is a **blood** reference applied to **brain
tumour** tissue, so it is used only for lymphoid *sub-composition* and never absolute scale, and
`Macrophage_Microglia` is explicitly not measurable by this route. Only **255 of 333** CpGs
survive complete-case filtering in these matrices — 64 are all-NA.

- Teschendorff AE, Breeze CE, Zheng SC, Beck S. EpiDISH. **[VERIFY]**

---

## B · Literature

### B1 · Benchmarks and reviews this study is positioned against

- Avila Cobos F, Alquicira-Hernandez J, Powell JE, Mestdagh P, De Preter K. Benchmarking of cell
  type deconvolution pipelines for transcriptomics data. *Nature Communications* **11**:5650
  (2020). doi:10.1038/s41467-020-19015-1
- Sturm G, Finotello F, Petitprez F, *et al.* Comprehensive evaluation of transcriptome-based
  cell-type quantification methods for immuno-oncology. *Bioinformatics* **35**:436–445 (2019).
  doi:10.1093/bioinformatics/btz363
- Nguyen H, Nguyen H, Tran D, Draghici S, Nguyen T. Fourteen years of cellular deconvolution:
  methodology, applications, technical evaluation and outstanding challenges. *Nucleic Acids
  Research* **52**:4761 (2024). doi:10.1093/nar/gkae267
- Gaspard-Boulinc LC, *et al.* Cell-type deconvolution methods for spatial transcriptomics.
  *Nature Reviews Genetics* **26**:828 (2025). doi:10.1038/s41576-025-00845-y
- Liu F, *et al.* DNA Methylation-Based Cell Type Deconvolution Reveals the Distinct Cell
  Composition in Brain Tumor Microenvironment. *bioRxiv* (2025). doi:10.1101/2025.01.19.633794

### B2 · Methods in or adjacent to the 15-method panel

**Cited with a verified record:**

- **CIBERSORTx** — Newman AM, *et al.* Determining cell type abundance and expression from bulk
  tissues with digital cytometry. *Nature Biotechnology* **37**:773 (2019).
  doi:10.1038/s41587-019-0114-2 *(B-mode and S-mode; Supplementary Table 1d records the mode
  chosen per dataset)*
- **Bisque** — Jew B, *et al.* Accurate estimation of cell composition in bulk expression through
  robust integration of single-cell information. *Nature Communications* (2020).
  doi:10.1038/s41467-020-15816-6
- **DWLS** — Tsoucas D, Sistig A. DWLS: Gene Expression Deconvolution Using Dampened Weighted
  Least Squares. *(package documentation on disk; the primary paper's DOI is not recorded here)*
  **[VERIFY]**

**Assessed or considered, not in the final panel:**

- **CDSeq** — Kang K, Meng Q, Shats I, Umbach DM, Li M, Li Y, Li X, Li L. CDSeq: A novel complete
  deconvolution method for dissecting heterogeneous samples using gene expression data.
  *PLoS Computational Biology* (2019). doi:10.1371/journal.pcbi.1007510
  *(reference-free; blocked on a macOS Fortran toolchain, see `docs/ROAD_TO_PAPER.md`)*
- **Scaden** — Menden K, *et al.* Deep learning–based cell composition analysis from tissue
  expression profiles. *Science Advances* **6**:eaba2619 (2020). doi:10.1126/sciadv.aba2619
  *(assessed and dropped: a deep model needing training data this cohort cannot supply)*
- **GBMdeconvoluteR** — *Neuro-Oncology* **25**(7):1236–1248 (2023). doi:10.1093/neuonc/noad021
  *(GBM-specific deconvolution)*
- **EcoTyper** — Profiling Cellular Ecosystems at Single-Cell Resolution and at Scale with
  EcoTyper. *Methods in Molecular Biology* chapter (2023). **[VERIFY]**

**In the panel, with no citation recoverable from anything on disk. These must be added before
submission:**

- **MuSiC** — Wang X, *et al.* **[VERIFY]** *(package vendored; its DESCRIPTION carries no DOI)*
- **SCDC** — Dong M, *et al.* **[VERIFY]** *(package vendored; DESCRIPTION carries no DOI)*
- **EPIC** — Racle J, *et al.* **[VERIFY]** *(installed as an R package, not vendored)*
- **quanTIseq** — Finotello F, *et al.* **[VERIFY]** *(installed as an R package)*
- **BayesPrism** — Chu T, *et al.* **[VERIFY]** *(installed as an R package)*

> `docs/METHODS.md` documents what each of these does and every deviation from its published
> algorithm. What it does not carry is a bibliographic record, which is why five entries above
> are incomplete rather than filled in.

### B3 · Technical background

- Li X, Gibson G, Qiu P. Gene representation in scRNA-seq is correlated with common motifs at the
  3′ end of transcripts. *Frontiers in Bioinformatics* **3**:1120290 (2023).
  doi:10.3389/fbinf.2023.1120290

---

## C · How to read this list

**Section A is load-bearing.** Remove any entry and a result disappears.

**Section B1 is positioning** — benchmarks this study's findings agree or disagree with, set out
in `docs/RELATED_WORK.md`. Listing a benchmark is not a claim to have reproduced it.

**Section B2 is mixed**: three methods in the panel are cited properly, five are not, and four
further tools were assessed and excluded with the reason stated.

**Outstanding before submission**

1. Complete the **13 [VERIFY] entries** against the publisher's record.
2. In particular, the **five panel methods in B2** — MuSiC, SCDC, EPIC, quanTIseq, BayesPrism —
   are used in every result and currently uncited.
3. Choose a citation style. The exemplar manuscript this is modelled on uses IEEE numeric.
