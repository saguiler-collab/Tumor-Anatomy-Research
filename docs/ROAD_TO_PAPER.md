# Road to the paper

Everything that must happen before writing, in order, with dependencies. Status
2026-09-10, after the confirmatory run (`results_archive/2026-09-10T2039`).

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

`OPEN_DEFECTS.md` D1. EPIC is **confirmed** corrected twice — the package returns
`mRNAProportions` and `cellFractions` separately and the driver takes the corrected one,
after which `fit_predict` corrects again. MuSiC is near-confirmed. SCDC, Bisque and
BayesPrism are suspected and unverified.

**Why blocking:** it changes every accuracy number, therefore the accuracy ranking,
therefore rho — the paper's headline. It does **not** change ACS, because a constant
per-type factor cannot reorder a cell type across structures.

**Steps.** For each of MuSiC, SCDC, Bisque, BayesPrism: run the package on its own bundled
example; establish from its vignette or return structure whether the output is a
transcript share or a cell share; record the evidence. Then choose one of the three
designs in D1 and apply it uniformly. Publish before/after numbers for every affected
method so the effect on the leaderboard is visible, not absorbed.

**Do not** fix it method-by-method as each is verified — a leaderboard half-corrected is
worse than one uniformly wrong.

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

### 2.2 A method that differs in kind

CDSeq (reference-free, no signature matrix) or Scaden (deep learning). The 14 comparable
methods produce only **8 distinct ACS values**, four tied at 0.9846. Another least-squares
variant joins the tie and *lowers* the rank correlation's resolution; a method with
different assumptions spreads the range.

### 2.3 An orthogonal yardstick

- **ABSOLUTE purity** for TCGA-GBM (Carter 2012) — DNA-based, never touches RNA.
- **Liu, Qian & Ma 2025** — DNA-methylation deconvolution of the brain-tumour
  microenvironment. Orthogonal modality, and specifically about brain tumours.

Either shares no failure mode with RNA-based deconvolution, which is what makes it worth
more than a fourth RNA method.

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
