# Road to the paper

Everything that must happen before writing, in order, with dependencies. Status 2026-09-11, after the confirmatory run (`results_archive/2026-09-10T2039`).

`docs/EXTERNAL_CHECKLIST.md` ranks items by value. **This file sequences them**, and
separates what *blocks* a credible paper from what *strengthens* one.

---

## The honest headline first

**You can write a defensible paper today**, on one tissue, if the scope claim is
"a case study demonstrating the method on glioblastoma" rather than "a general method".
Everything in Tier 0 below is about not publishing something false; nothing in it is
about getting a better number.

What you cannot yet claim is that anatomic concordance *is a method* for choosing
deconvolution tools. `Anatomy_Test.md` says why: one tissue is a case study, two is a
method. That is Tier 2.

---

## Tier 0 · Blocking — do not write the paper without these

These are not improvements. Each one prevents a statement in the paper from being wrong.

### 0.1 The cell-size correction · **RESOLVED 2026-09-12 for MuSiC, SCDC and least squares**

Four revisions in one day, all recorded in `docs/OPEN_DEFECTS.md` D1. The conclusion rests on
a probe that mirrors production's actual condition — a marker-subset export whose per-type
library sizes span 1.92x with the tumour type lowest, matching the 2.17x measured on the real
atlas (Tumor 29,999, Oligodendrocyte 65,003).

The probe has two correct answers available, which is what makes it decisive: the true mRNA
share is **0.7500** over the full gene space and **0.6094** over the exported subset, against
a true cell fraction of **0.5000**.

| package | measured | estimand | after the central conversion | error |
|---|---|---|---|---|
| MuSiC | 0.7500 | full space | **0.5000** | +0.0000 |
| NNLS (control) | 0.7506 | full space | 0.5008 | +0.0008 |
| SCDC | 0.7469 | full space | 0.4959 | −0.0041 |
| **EPIC** | **0.6104** | **subset** | **0.3431** | **−0.1569** |

**Not a double correction.** MuSiC's own `M.S` is non-flat and does apply a conversion, but
it corrects the *subset bias*, not cell size — it cannot correct cell size, because
normalising each cell to 1e6 destroyed that information before MuSiC saw the cells. What it
leaves is exactly the full-transcriptome mRNA share, which is what this project's conversion
expects. The two compose correctly.

**So the earlier "every leaderboard number is provisional" is withdrawn, and so is the
opposite over-claim made and retracted the same day.** For MuSiC, SCDC and the least-squares
family the conversion is applied once and lands within 0.005 of truth.

**Still open, and moved to its own entry:**

1. **EPIC — now Tier 0.4, and a new defect (D11).** `scaleExprs = TRUE`, its published
   default, renormalises onto the shared gene set and makes EPIC's estimand the within-subset
   mRNA share. The central conversion then divides a subset share by full-transcriptome
   factors. Mechanism confirmed: with `scaleExprs = FALSE` EPIC returns 0.7504. EPIC is
   joint-2nd on the leaderboard at 0.9846, so this is not a footnote.
2. **Bisque and BayesPrism remain unmeasured.** The probe cannot read either — ≥31% of their
   mass lands on roster types absent from the mixture. A probe with no absent types would
   settle it. Their conventions must not be asserted meanwhile, and the old claim that Bisque
   was double-corrected is withdrawn as unmeasured.
3. **Declare the substitution**, which is done: `DECLARED_DEVIATIONS` now carries it for
   music, scdc, scdc_ensemble, bisque and bayesprism, each stating what is measured and what
   is not.

### 0.2 Publish the registration corrections · **EXTERNAL, 20 minutes**

The registration is permanent and correct not to edit. Two statements in it are wrong:

- *"Two published packages exceeded their wall-clock budget in one stage"* — true of
  `2026-09-05T2154`, understates `2026-09-06T1103` (two stages) and the confirmatory run.
- *"Fifteen deconvolution methods"* — the plan; both pilot runs it cites ran fourteen.

**Steps.** Write the correction into the paper's methods section and into `RESULTS.md`,
phrased as a correction *beside* the registration. Do not amend the OSF record. A
registration with a published correction is stronger than one nobody checked.

### 0.3 Settle how DWLS and BayesPrism are reported · **REOPENED 2026-09-12 — was marked DONE in error**

Both genuine R packages were measured, and both completed:

| method | genuine R package | elapsed | pipeline budget | genes |
|---|---|---|---|---|
| `R:DWLS` | ACS 0.7077 [0.576, 0.864] | 3,873 s | 2,400 s | 1,591 |
| `R:BayesPrism` | ACS 0.8923 [0.800, 0.972] | 4,689 s | 2,400 s | 1,591 |

**Why this is reopened.** Every method on the leaderboard ran on the benchmark's **657**
gene marker subset. Both re-measurements ran on **1,591** —
`scripts/remeasure_method.py` was reading the gene space the pipeline uses only as a
*fallback*, and only the gene COUNT had ever been persisted, so nothing caught it. The
entry previously read "It scores *lower* than the reimplementation (0.7385) and stays last.
The disagreement with Avila Cobos is real, not an artefact of substituted software."
**Neither sentence follows from measurements on gene spaces differing 2.4-fold, and both
are withdrawn.** See `docs/OPEN_DEFECTS.md` D10.

**Already done, and it is what makes the rest cheap:**
- `run_all.py` persists the exact gene list and its hash to
  `results/benchmark/signature_genes.json`.
- `remeasure_method.py` aborts unless that file exists and its hash reproduces, instead of
  silently substituting a gene space. It fails in ~2 s rather than after an hour of CPU.
- Both report files carry a `COMPARABILITY: NOT COMPARABLE` header.
- README, `RELATED_WORK.md` and this entry carry the withdrawal.

**Steps to close it.**
1. One full `python scripts/run_all.py` — writes `signature_genes.json`.
2. `python scripts/remeasure_method.py --method dwls --budget 14400`
3. `python scripts/remeasure_method.py --method bayesprism --budget 14400`
   (~1.1 h and ~1.3 h respectively, measured, and they must not run concurrently on this
   machine — three BayesPrism Gibbs workers at ~900 MB each already put it into swap.)
4. Then, and only then, state whether genuine DWLS scores above or below its
   reimplementation, and whether the Avila Cobos disagreement survives.

**Until then:** every DWLS and BayesPrism statement in the paper must name the software
**and** the gene space, and must not compare across the two. The archived leaderboard rows
are the reimplementations and stay as recorded.

---

### 0.4 EPIC's estimand · **BLOCKING, new 2026-09-12**

`docs/OPEN_DEFECTS.md` D11. EPIC returns the mRNA share **within the exported gene subset**;
every other method returns the full-space share. The cell-size conversion divides by
full-transcriptome factors, so for EPIC it is applied to the wrong quantity — measured error
**−0.157** in the dominant type on the probe, where every other method is within 0.005.

Mechanism confirmed: `scaleExprs = TRUE` (EPIC's published default). With it off EPIC returns
0.7504 and agrees with the others.

**Do not fix this by changing `scaleExprs`.** It is the published default, the fairness rule
is defaults with no per-method tuning, and EPIC's ACS is already known — changing a parameter
after seeing a score is the retune the protocol forbids and would be indistinguishable from
tuning EPIC up the leaderboard. EPIC is joint-2nd at 0.9846.

**Legitimate routes, and both need a decision recorded before any number moves:**

1. **State per method which estimand it returns**, and convert accordingly — a documented
   rule, pre-stated, applied to every method rather than to EPIC alone. This is the honest
   fix and it requires re-running the affected rows.
2. **Report EPIC's row with the mis-scaling disclosed** and make no correction, treating it
   as a property of running EPIC on a marker subset.

**First, measure the effect on the real cohort — and note the evidence is mixed.** An
analytic propagation of the per-type marker-space concentrations measured on the real atlas
predicts a large distortion (Tumor 0.550 → 0.387). The pseudobulk benchmark contradicts it:
EPIC's `mae_primary` is **0.0659**, 7th of 15 and only modestly behind MuSiC's 0.0508, which
is not where a 0.16 bias on the dominant type would put it. So the probe result is solid and
the real-data magnitude is genuinely unresolved, probably smaller than the propagation
implies — `withOtherCells = TRUE` gives EPIC a sink column, the output is renormalised across
the roster afterwards, and the benchmark's own truth is built from the same normalised cells.

The measurement that settles it: EPIC re-run on the real cohort with `scaleExprs = FALSE` as
a **declared diagnostic**, compared composition-by-composition against the production row.
Recorded alongside, never swapped into the leaderboard.


## Tier 1 · Strongly recommended — the paper is materially weaker without these

### 1.1 Neftel cell-type annotations · **EXTERNAL**

`OPEN_DEFECTS.md` D2. GEO ships expression only; the 1.5 GB file cannot become a reference
without labels, so `scdc_ensemble` runs degenerate — a published tool with its entire
contribution switched off.

**Steps.**
1. Obtain the **authors'** annotations for GSE131928 (Broad Single Cell Portal hosts the
   annotated version). Do not annotate it here: that puts this project's clustering and
   marker-threshold choices inside a reference used to judge methods.
2. Map its cell types onto the eight-type roster, as `GBMAP_CELL_TYPE_MAP` does.
3. Register it as a second cell source and pass both to `SCDC_ENSEMBLE` via
   `sc.eset.list`. The R driver already branches on `length(esets) > 1`.
4. Record its SHA-256 in a provenance JSON as GBmap has.

**Buys:** a real ENSEMBLE row, the `scdc = scdc_ensemble` tie broken, and a
**reference-stability check** — is the ACS ranking stable across references, not just
across cohorts? That check does not exist and cannot be built without this.

### 1.2 H&E figure panels · **EXTERNAL, one afternoon**

Two or three panels showing LE / CT / MVP / PAN.

**Why it matters more than it sounds.** The paper's central premise is that certain
orderings are *near-definitional* — that MVP **is** proliferating endothelium. A reader
who has never seen a GBM slide must take that on faith. One figure converts the premise
from an assertion into something checkable. Allen images require attribution.

### 1.3 Decide what the ISH result supports · **INTERNAL**

Done and committed: ESM1 supports C3 (p = 0.005) and C4 (p = 0.021), CD163 supports C6
(p = 0.0085). Three limits must travel with those numbers wherever they appear — C1/C7
are marker-dependent, C5 has one evaluable block, C2 is unmeasurable with that panel, and
the endothelial/myeloid results rest on 6-8 sub-blocks.

**Step:** write the limits into the figure captions, not only the supplement.

---

## Tier 2 · Changes what the paper can claim

### 2.1 A second tissue

The only thing that lifts the ceiling. MuSiC leads by **one constraint–tumour pair out of
57**, and six methods' intervals contain its score.

Two routes, and the decision is yours:

| | Allen Human Brain Atlas | Visium GBM (Ravi 2022) |
|---|---|---|
| Constraint file | **new one needed** — `constraints_brain.py`, drafted, not frozen | **existing one applies unchanged** |
| Registration | a second, separate one | none needed |
| Reference | Siletti 2023 (normal brain) | existing GBmap |
| Unit | microdissected region | 55 µm spot — a handful of cells |
| Power | 6 donors, denominator 42 (microarray) | unknown until examined |

**Steps for the Allen route**, in strict order — the ordering is the scientific claim and
cannot be repaired afterwards:
1. Settle the claims in `constraints_brain.py`, set `FROZEN = True`.
2. `freeze_hash()` and register publicly — a **second, separate** registration.
3. Only then download Allen bulk and build the Siletti reference.
4. Only then deconvolve and score.

`freeze_hash()` raises while `FROZEN = False`, so step 4 cannot silently happen first.

### 2.2 Methods that differ in kind — CDSeq and Scaden

**Status 2026-09-11: sources obtained, neither installed, and each is blocked by
something different.**

| | source | blocker | verdict |
|---|---|---|---|
| **CDSeq** | R package v1.0.9, in `pipeline packages / repos/CDSeq/` | **macOS Fortran toolchain missing.** `ld: library 'emutls_w' not found`; R is configured to use `/opt/gfortran/bin/gfortran`, which does not exist. Its R dependencies (`dirmult`, `RcppThread`, `ggpubr`, `harmony`) installed fine. | **Worth unblocking** |
| **Scaden** | Python source, in `~/Downloads/scaden-master` | **`tensorflow>=2.0` has no wheel for Python 3.14.5.** `pip` reports "No matching distribution found". | **Recommend dropping** |

**To unblock CDSeq — EXTERNAL, one installer.** Install the official R macOS toolchain
from <https://mac.r-project.org/tools/> (the gfortran package matching R 4.6, x86_64).
Then `install.packages("<path to CDSeq_R_Package-master>", repos=NULL, type="source")`.
Nothing else is missing; the failure is purely a linker looking for Fortran runtime
libraries that are not on this machine.

**Why drop Scaden rather than build a second Python environment.** Two reasons, and the
second is the real one:

1. TensorFlow would need a separate Python 3.11 or 3.12 environment maintained alongside
   the project's 3.14 — a standing maintenance cost for one method.
2. **Scaden trains on simulated bulk**, which is the construction Nguyen et al. (2024)
   single out: *"any algorithm would be the best, when applied to data that was simulated
   based on the same set of assumptions."* Adding a method trained on simulated mixtures
   to a benchmark already criticised for leaning on simulated mixtures buys less than it
   costs. CDSeq, being reference-free, buys something no current entry provides.

**Scaden is dropped, and the paper will say so.** Suggested wording:

> Scaden (Menden et al. 2020) was obtained and deliberately not evaluated. Two reasons.
> Practically, it requires TensorFlow, for which no wheel exists on this project's Python
> version, so including it would mean maintaining a parallel interpreter for one method.
> Substantively, Scaden is trained on simulated bulk mixtures — the construction Nguyen
> et al. (2024) identify as biasing benchmarks toward methods that share the simulator's
> assumptions. Since this study's accuracy arm is itself built on simulated mixtures,
> adding a method trained on them would compound a limitation the study already discloses
> rather than test against it. CDSeq was pursued instead, being reference-free and
> therefore independent of both the simulator and the signature matrix.

A deliberate exclusion with a stated rationale reads better than silence, and it is the
same discipline the constraint file already applies to T cells.

The 14 comparable methods produce only **8 distinct ACS values**, four tied at 0.9846.
Another least-squares variant joins that tie and *lowers* the rank correlation's
resolution. A method built on different assumptions spreads the range instead.

| | what it is | why it adds something |
|---|---|---|
| **CDSeq** | **Reference-free.** Uses no signature matrix at all; infers cell-type profiles and proportions jointly from the bulk. | Tests something no current entry does: can anatomy rank a method that never saw your reference? Every one of the 15 shares the same GBmap signature, so they share its biases. CDSeq does not. |
| **Scaden** | Deep learning, trained on simulated bulk. | **DROPPED — decision taken 2026-09-12.** Reported in the paper as a deliberate exclusion with its reason, not omitted silently. |

**CDSeq is the more valuable of the two**, for the reason in its row: it is the only
candidate that breaks the shared-reference dependency.

**The registration consequence.** The registration says *"Fifteen deconvolution methods
and two negative controls."* Adding either makes it sixteen or seventeen — a change to
the registered analysis plan, decided after the leaderboard was seen. Three honest options, in the **Registration** section below. Do not edit the existing registration.

### 2.3 ABSOLUTE purity — **data obtained 2026-09-11, not yet computed**

This is the protocol's **yardstick 2**, and it is closer than it looked.

**What already exists.** `agreement.py` evaluates four yardsticks, not one — the
confirmatory run reports `absolute_purity` as UNAVAILABLE rather than absent. So the
machinery is wired; it needs data, not code.

**What was obtained today.** `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` from the
PanCanAtlas, GDC file UUID `4f277128-f793-4354-a13d-30cc7fe9f6b5`, 901,812 bytes,
sha256 `f430a975433d82e0…`. 10,786 rows with both `purity` and `Cancer DNA fraction`.
Provenance in `reference_frozen/tcga_benchmark/ABSOLUTE_PROVENANCE.json`; the file itself
is gitignored, as data is.

**Overlap, measured:** the predecessor project's `TCGA-GBM.star_counts.tsv` has 175 bulk
samples; **156 of them have an ABSOLUTE purity call** — more than the 120 the predecessor
reported.

**Why it is still not computable.** The existing `orthogonal_validation.csv` records the
purity correlation for **one** method (`v2_STAR_SVR_clean`, rho = 0.761, n = 120). One
point cannot rank anything. To use this as a yardstick, all 15 methods must be run on the
TCGA-GBM cohort and each method's Tumor fraction correlated against purity.

**Steps.**
1. Load `TCGA-GBM.star_counts.tsv` (predecessor project, read-only) and normalise to CPM.
2. Deconvolve all 15 methods against the frozen signature — the same one, for equal footing.
3. Per method, Spearman-correlate the Tumor column against `Cancer DNA fraction` across
   the 156 shared samples.
4. Feed those correlations to `load_absolute_purity()`; the agreement test picks them up
   with no further change.

**Why this is worth more than a sixteenth method.** It is a *second yardstick*, measured
from **DNA**, sharing no failure mode with any RNA-based method under test. Nguyen et al.
(2024) criticise simulation-based benchmarking for favouring methods that share the
simulator's assumptions; a DNA-derived yardstick is immune to that criticism in a way no
RNA method can be. It also gives the agreement test a second, independent rho — turning a
single correlation into a replication.

**Also worth considering:** Liu, Qian & Ma 2025 — DNA-methylation deconvolution of the
brain-tumour microenvironment. Same orthogonality argument, and specifically about brain
tumours.

---

## Tier 2.5 · Registration: what you may and may not do

**You cannot edit or update a published OSF registration. That is the point of one.**
A record revisable after seeing results proves nothing, which is exactly why yours is
worth something.

Withdrawing is worse than doing nothing: it leaves a public tombstone naming you and
saying the registration was withdrawn, and a reader will assume something serious was
wrong rather than "the stage count was off by one".

So for anything that changes the registered plan — adding CDSeq or Scaden, adding the
ABSOLUTE yardstick, changing the cell-size mechanism — there are three honest routes:

| route | when it fits | cost |
|---|---|---|
| **Declared deviation** | Small additions reported separately and kept out of the primary correlation. | A paragraph in the methods. Nothing filed. |
| **Exploratory** | Your registration already reserves this: *"Everything not named above is exploratory and labelled as such."* | Nothing. Already covered. |
| **A second registration** | A materially different analysis — a second tissue, or a substantially expanded method set. | A new OSF record citing the first. The original stays. |

**What I would do for each pending item:**

- **CDSeq / Scaden** — declared deviation, reported outside the primary correlation, which
  stays on the fifteen registered methods. Or fold into a second registration if the
  second tissue happens anyway, since that needs one regardless.
- **ABSOLUTE purity** — **no registration change needed.** It is already a registered
  yardstick; `agreement.py` evaluates it and currently reports UNAVAILABLE. Supplying the
  data executes the registered plan rather than departing from it.
- **The cell-size fix** — declared deviation. The registered *intent* ("applied once") is
  preserved; only the mechanism changes.
- **The two 0.2 corrections** — published beside the record. Not a deviation at all, just
  errata.

**The order that matters:** do not file a second registration until the plan is settled.
A second registration filed now, then amended again in a week, is worth less than one
filed once when the method set and yardsticks are final.

---

## Tier 3 · Publication hygiene · **EXTERNAL**

- [ ] **Link ORCID** to the OSF registration so the record outlives an institutional email.
- [ ] **Confirm `main` is the GitHub default branch** (Settings → Branches).
- [ ] **Zenodo deposit** for the 3.8 GB vendored package tree, if it is to be preserved
      with a DOI. Git cannot hold it and should not.
- [ ] **Decide on a preprint.** bioRxiv is how methods work reaches the field; it also
      timestamps priority, which matters for a student-led project.
- [ ] **Ivy GAP vital status** — almost certainly unobtainable. Worth one message to
      <https://community.brain-map.org/>, but plan on survival staying BLOCKED. That is a
      result, not a gap.

---

## What "done" looks like

**Minimum defensible paper:** Tier 0 complete. Scope claimed as a single-tissue case
study. Reports rho = 0.7501, the controls behaving, the ISH support for C3/C4/C6, the
null concordance result, and every limitation in `OPEN_DEFECTS.md`.

**Strong paper:** Tier 0 + Tier 1. Adds a working SCDC ENSEMBLE, reference stability, and
figures a reader can check the premise against.

**The paper the protocol set out to write:** Tier 0 + 1 + 2.1. Two tissues. Only then is
the claim "anatomic concordance is a method for choosing a deconvolution tool" supported
rather than suggested.
