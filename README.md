# The Anatomy Test

**Anatomic Concordance: Neuropathology as a Ground-Truth-Free Benchmark for Cell-Type
Deconvolution in Glioblastoma**

Can a tumor's own anatomy stand in for ground truth when choosing a cell-type
deconvolution method — and does a method that gets the anatomy right also get the
biology right?

Read **[Anatomy_Test.md](Anatomy_Test.md)** first. It is the protocol: the hypothesis,
the pre-registered constraint set, the negative controls, and what each possible outcome
would mean.

## The one idea

Deconvolution estimates cell fractions from bulk RNA. To trust an estimate you need
something to check it against, and for most tissues nothing exists. So the field falls
back on a consensus of methods, or on the informal check that the answer "matches known
biology."

That informal check has a structural problem: when the anatomic pattern *falls out of*
the deconvolution and is then reported as a finding, the method can never fail.

```
HOW IT IS USUALLY DONE
    method ──solves──> fractions per region ──read as──> anatomy      (a CONCLUSION)

THE ANATOMY TEST
    constraints written in advance ──> anatomy ──scores──> method     (a TEST)
```

Writing the ordering down first turns the same pattern into a prediction a method can
miss. That is the whole project.

## The two results

**1. The ACS leaderboard.** Twelve methods — ten real, two negative controls — scored on
how well they reproduce seven pre-registered ordinal constraints, paired within tumor,
with a bootstrap CI over tumors and a p-value against a within-tumor label-permutation
null.

**2. The experiment.** Does the ACS ranking agree with the ranking from *real* ground
truth? One number: Spearman ρ, judged against a bar fixed before it was computed
(ρ ≥ 0.60, bootstrap CI excluding zero).

Both directions are publishable, and the code treats them that way — `agreement.py`
returns `NULL RESULT: ... reproducing known biology is NOT evidence that a composition
estimate is numerically correct` as a first-class outcome, not an error.

### Where it currently stands

See [RESULTS.md](RESULTS.md). In short, from the 2026-09-03 full-database run:

- **The leaderboard exists and the controls discriminate**, but ACS's resolution on this
  cohort is low: the weighted denominator is 65, so ten real methods produce **five
  distinct scores** and five of them tie at exactly 0.800.
- **Two real methods are not distinguishable from a randomly permuted signature.**
  Against the shuffled-signature control run as a distribution rather than a single
  draw, `bayesian` and `bayesian_hierarchical` fall inside its 5–95% range.
- **The agreement test is still not computable**, and now for two independent reasons:
  no yardstick covers more than two methods, *and* the ACS side supplies only five
  distinct values. Fixing the first no longer suffices.
- **Deconvolving all 270 archive samples does not move the ranking** (Spearman 0.982).
  Only methods that borrow strength across samples move at all, and Bisque moves *down*.
- **Prognosis is BLOCKED** — Ivy GAP publishes no vital-status column — and would be
  INCONCLUSIVE even under the most generous reading of what it does publish.

## The data, verified

Ivy GAP's 270 RNA-seq samples are **two studies in one file**:

| Study | Samples | Tumors | Structure assigned by | ACS-scored | Deconvolved |
|---|---|---|---|---|---|
| Anatomic structures | **122** | **10** | H&E histology, no reference to RNA | **yes** | yes |
| Cancer-stem-cell clusters | 148 | 34 | ISH probes — i.e. *using expression* | **never** | yes, with `--full-database` |

Only the first is scorable: its labels are an independent physical fact about each
sample. The trap is that a prefix match on the acronym pulls all 148 cluster samples in —
`CT-CD44` starts with `CT` — tripling the apparent cohort with samples whose labels came
from expression data, making the test circular while appearing to strengthen it.
`config.is_anatomic_sample` requires the `-reference-histology` suffix, `run_anatomic`
raises if a cluster sample reaches the scoring set, and tests guard both.

Deconvolving the cluster samples is a different matter and is safe: composition is well
defined on all 270, and the wider set is what takes the prognostic cohort from 7 tumours
to 29. What may never widen is the *scoring*.

**Checked against the portal, not against ourselves.** The loader's counts come from
`columns-samples.csv`, which ships inside the archive it describes. They are reconciled
sample-by-sample against the Allen Institute's own live metadata export
(`rna_seq_samples_details.csv`), which is produced independently of this pipeline:
122 anatomic samples, CT/IT/LE/MVP/PAN = 30/24/19/25/24, study assignment agreeing on all
270. That closes protocol step 3's gate. See
[docs/DATA_SOURCES.md](docs/DATA_SOURCES.md).

## Run it

```bash
pip install -r requirements.txt

python -m ivygap.data.download_ivygap     # expression + clinical + portal metadata

python scripts/run_all.py                 # the pre-registered cohort: 122 anatomic samples
python scripts/run_all.py --full-database # also deconvolve all 270; ACS still scored on 122

python scripts/run_all.py --synthetic     # end-to-end on generated data, no downloads
pytest tests/ -q
```

`--full-database` widens the deconvolution, never the scoring — that separation is
enforced inside `run_anatomic.run`, not left to the caller. The two runs write to
separate directories and cannot overwrite each other.

Survival is BLOCKED by default: Ivy GAP publishes `survival_days` and no vital-status
column, so who was censored is unknowable. `--assume-recorded-survival-is-death` will
proceed under that named assumption, stamped into every artefact and written to a
separate `declared_event_policy/` directory. Read
`results/clinical_missingness_audit.json` first — the blanks are not missing at random.

`reference_frozen/` holds the signature matrix and cell-size factors, committed with
recorded SHA-256 hashes (protocol step 1), so no multi-gigabyte atlas download is needed
to reproduce a run. Its `PROVENANCE.json` also states what it *cannot* support — see
below.

Optional, for the genuine R packages rather than this project's reimplementations:
`Rscript R/install_deps.R`.

## What this repo is careful about

**Selection never touches an outcome, and ACS never selects.** The frozen method is
chosen on donor-held-out synthetic mixtures where composition is known. ACS is the
*object under test* — using it to select before knowing whether it is safe to select on
would assume the conclusion. Survival is computed last and
`assert_selection_frozen` aborts if the recorded selection criterion so much as mentions
an outcome.

**Equal footing is verified, not asserted.** Ivy GAP's samples are nested in tumors with
uneven block counts (6 to 15). Everything collapses to one value per (tumor, structure)
before scoring, so a heavily dissected tumor cannot outvote a sparsely dissected one —
`test_uneven_block_counts_do_not_move_the_score` asserts the score is *identical*.
`bench/equal_footing.py` hashes every shared input and **raises** if two methods' inputs
disagree.

**Degenerate methods are labelled degenerate.** The vendored signature is collapsed, so
it carries no cross-donor variance. MuSiC's entire contribution is weighting genes by
cross-donor consistency, so on this reference it is arithmetically NNLS. It reports that
rather than appearing on a leaderboard under a name it hasn't earned. Same for SCDC
ENSEMBLE with one reference, and Bisque in no-overlap mode.

**Missing inputs are reported, never imputed.** Every documented `tumor_details.csv` URL
now returns an HTML page with a 200 status. The downloader rejects it; the survival stage
reports BLOCKED. No prognostic claim is made, and TCGA outcomes are never carried across.

## Layout

```
Anatomy_Test.md              the protocol — read first
reference_frozen/            vendored signature + cell sizes, hashed (step 1)
ivygap/
  anatomic/
    constraints.py           C1-C7, frozen and hashed for pre-registration
    acs.py                   the score, permutation null, bootstrap CI
    agreement.py             the primary result
    run_anatomic.py          leaderboard + experiment
  deconv/                    10 methods + 2 negative controls, one contract
  bench/                     pseudobulk ground truth, metrics, equal-footing guard
  survival/                  out-of-fold prognosis, power stated before fitting
docs/METHODS.md              every method, and every deviation from its paper
tests/                       101 tests, negative controls included
```

## Data

Ivy GAP is openly downloadable — no application, no data-use agreement.
Puchalski RB et al., *An anatomic transcriptional atlas of human glioblastoma*,
Science 360:660-663 (2018). <https://glioblastoma.alleninstitute.org/>
