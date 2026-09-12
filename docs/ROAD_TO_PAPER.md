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

### 0.1 Verify the cell-size double-correction · **INTERNAL + package docs**

`OPEN_DEFECTS.md` D1. **Verification is COMPLETE as of 2026-09-10** — all five packages
were tested on a reference where two cell types differ 3x in mRNA content, mixed 50/50 by
cell count, so 0.50 means the package already divided out cell size and 0.75 means it did
not:

| method | returned | reports | our correction is |
|---|---|---|---|
| **MuSiC** | 0.500 | cell share | **a second correction — defect** |
| **Bisque** | 0.500 | cell share | **defect** |
| **EPIC** | offers both | cell share, as taken | **defect** |
| BayesPrism | 0.750 | mRNA share | correct, the first |
| SCDC | 0.653 | between the two | **partial** |

**Three of the top methods are double-corrected, including MuSiC**, which leads ACS at
1.000 and the accuracy benchmark at MAE 0.0508.

**Why blocking — and this is worse than stated here until 2026-09-11.** An earlier
version of this file said the defect "does not change ACS". That was asserted from the
shape of the formula, never tested, and is **wrong**. `to_cell_fractions` renormalises
each sample after dividing by cell size, and the renormalisation divisor depends on that
sample's composition, so the effective scaling varies sample to sample and *can* reorder a
cell type across structures.

Measured by applying the correction a second time to every archived estimate table:
**15 of 16 methods change ACS**, the largest change is 0.0462, and **7 of 14 ranks move** —
DWLS rises from 14th to 12th, EPIC and SVR rise to joint first.

So it moves **ACS, the leaderboard ordering, the accuracy arm, and rho.** The design is
untouched — controls, permutation null, constraint file and ISH validation do not depend
on the cell-size conversion — but every leaderboard number is provisional until this is
fixed.

**What remains is a decision, not an investigation.** SCDC's 0.653 settles it: it is
neither answer, so neither "apply ours" nor "skip ours" is right for it. Only the third
design in D1 handles all five — pass this project's cell-size factors *into* each package
that accepts one (`MuSiC(cell_size=)`, `SCDC_ENSEMBLE(ct.cell.size=)`, `EPIC(mRNA_cell=)`)
and skip the central step for those.

**Steps.**
1. Adopt the third design.
2. Apply it to all five at once. **Do not** fix them one at a time — a half-corrected
   leaderboard is worse than a uniformly wrong one.
3. Re-run, and publish before/after numbers per method so the effect is visible.
4. Record it as a declared deviation: the registration says the conversion is "applied
   once and centrally", and passing factors into each package is a different mechanism for
   the same intent.

### 0.2 Publish the registration corrections · **EXTERNAL, 20 minutes**

The registration is permanent and correct not to edit. Two statements in it are wrong:

- *"Two published packages exceeded their wall-clock budget in one stage"* — true of
  `2026-09-05T2154`, understates `2026-09-06T1103` (two stages) and the confirmatory run.
- *"Fifteen deconvolution methods"* — the plan; both pilot runs it cites ran fourteen.

**Steps.** Write the correction into the paper's methods section and into `RESULTS.md`,
phrased as a correction *beside* the registration. Do not amend the OSF record. A
registration with a published correction is stronger than one nobody checked.

### 0.3 Settle how DWLS is reported · **DONE 2026-09-10**

Resolved. The genuine R package was measured at **ACS 0.7077** [0.576, 0.864], 57 pairs,
null_p 0.0001, in 3,873 s. It scores *lower* than the reimplementation (0.7385) and stays
last. The disagreement with Avila Cobos is real, not an artefact of substituted software.

**Remaining step:** every DWLS statement in the paper must quote 0.7077 and name the
software. The archived leaderboard row is the reimplementation and stays as recorded.

---

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
