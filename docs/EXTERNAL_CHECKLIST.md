# External checklist

Everything that has to happen outside this repository for the deconvolution comparison to
be as accurate and as complete as it can be. Ordered within each section by what it buys,
not by effort.

`docs/EXTERNAL_ACTIONS.md` holds step-by-step instructions for the older items.
`docs/OPEN_DEFECTS.md` holds what is known to be wrong internally.

Status as of 2026-09-10.

---

## A · Data — essentially done

- [x] **Ivy GAP bulk RNA-seq** — `gene_expression_matrix_2014-11-25.zip`, 270 samples / 37 tumours
- [x] **Ivy GAP clinical + sample metadata** — `tumor_details.csv`, `rna_seq_samples_details.csv`
- [x] **Ivy GAP ISH quantification** — `gene_expression_details.csv`, 18,778 rows of
      expression energy per (gene × sub-block × structure)
- [x] **Ivy GAP structure areas** — `sub_block_details.csv`, normalised area per structure
- [x] **GBmap Core** single-cell atlas — 338,564 cells
- [x] **Neftel GSE131928** — 10x (16,201 cells / 21 donors) and SMART-Seq2
- [x] **Mossi Albiach** — 135,482 cells, 1 donor
- [ ] **Nothing further from Ivy GAP.** BAM files and the reference genome are not needed —
      nothing in this pipeline consumes reads or re-aligns.

---

## B · Reference quality — the largest accuracy gap

- [ ] **Cell-type annotations for Neftel GSE131928.** ← *highest value item on this page*

  GEO ships expression only. Without labels the 1.5 GB download cannot become a
  reference, so `scdc_ensemble` stays degenerate — a published tool running with its
  entire contribution switched off.

  Get the **authors'** annotations (Broad Single Cell Portal hosts the annotated version).
  Do **not** settle for annotating it here: that would put this project's clustering and
  marker-threshold choices inside a reference used to judge deconvolution methods, which
  forfeits most of what an independent second reference is for.

  *Buys:* a real SCDC ENSEMBLE row, and a **reference-stability check** — is the ACS
  ranking stable across references, not just across cohorts? That check does not exist
  today and cannot be built without this.

- [ ] **Decide what to do with Mossi Albiach.** One donor, so it cannot serve as a second
      reference. Its `Zone` / `Location` annotations along a tumour-core-to-cortex gradient
      make it interesting for checking the constraints instead. Decide: use for that, or
      set aside.

---

## C · Method fidelity — affects every accuracy number

These are internal defects, listed because resolving them may need external material.

- [ ] **Verify the cell-size double-correction per method** (`OPEN_DEFECTS.md` D1).
      EPIC confirmed corrected twice; MuSiC near-confirmed; SCDC, Bisque, BayesPrism
      suspected. Needs each package's documentation or vignette to establish whether its
      output is a transcript share or a cell share.
      *Buys:* correctness of the accuracy arm, hence of rho. Does **not** affect ACS.

- [ ] **Decide the correction design** once verification is done — convert centrally
      (registered, unreachable for MuSiC), let each package use its own factors (breaks
      uniformity), or pass this project's factors into each package that accepts one.
      A study-design decision, not a patch.

---

## D · Independent validation — makes the yardstick itself checkable

- [ ] **Rebuild the ISH constraint check** against `gene_expression_details.csv`, paired
      within `sub_block_id`, all rows, per marker, with a within-unit permutation null.

  The first attempt used 11% of the data, pooled across specimens, and pooled markers that
  disagree with each other. Its output is in `results/ish_constraint_check.json` and
  **must not be cited**. The correct input was on disk the whole time.

  *Buys:* constraints resting on measurement rather than expert expectation. Note the
  panel is a cancer-biology screen, not a cell-type panel: C3/C4 (endothelial) are well
  covered, C5/C6 (myeloid) weakly, C1/C7 (tumour) ambiguously, and **C2 not at all** —
  no myelin marker exists in it, and OLIG2 must not be substituted.

- [ ] **ABSOLUTE purity for TCGA-GBM** (`EXTERNAL_ACTIONS.md` item 3).
      An orthogonal yardstick derived from DNA that never touches RNA, so it cannot share
      a failure mode with any method under test.
      *Buys:* a second accuracy axis independent of the pseudobulk construction.

- [ ] **Two or three H&E panels** showing LE / CT / MVP / PAN, for the paper.
      *Buys:* a reader can verify that "MVP is defined by proliferating endothelium" is
      near-definitional instead of taking it on faith. Allen images require attribution.

---

## E · Statistical power — the ceiling nothing else lifts

- [ ] **Freeze `constraints_brain.py`** (`FROZEN = True`) once the claims are settled.
- [ ] **Register the brain constraint hash** as a *second, separate* registration.
- [ ] **Choose the Allen platform** — microarray (6 donors, weighted denominator 42) or
      RNA-seq (2 donors, denominator 14). The RNA-seq product is weaker than the cohort
      that already exists. Recorded in `docs/SECOND_TISSUE.md`.
- [ ] **Siletti 2023 Human Brain Cell Atlas** — the normal-brain reference. GBmap has 22
      neurons and normal cortex is mostly neurons.
- [ ] **Run the second tissue.**

  *Buys:* the difference between a case study and a method. MuSiC currently leads by
  **one constraint–tumour pair out of 57**, with six methods' intervals containing its
  score. Nine tumours cannot resolve the top of the leaderboard and no further analysis
  inside Ivy GAP will.

---

## F · Methods worth adding, and one that is not

- [ ] **CDSeq** — reference-free deconvolution, no signature matrix at all.
- [ ] **Scaden** — deep learning on simulated bulk.

  *Buys:* the 14 comparable methods currently produce only **8 distinct ACS values**, four
  tied at 0.9846. A method that differs in *kind* spreads the range; another least-squares
  variant joins the tie and lowers the rank correlation's resolution.

- [ ] ~~More least-squares or Bayesian variants~~ — actively unhelpful, for the reason above.
- [ ] ~~Hosted CIBERSORTx~~ — licence-gated, and both modes are implemented here from the
      paper and its Supplementary Note 1.

---

## G · Publication hygiene

- [x] **OSF registration live and public** — <https://osf.io/dm2t8>, hash `2d1fb47c…`
- [ ] **Publish the corrections** the registration needs: it says two packages fell back in
      *one* stage, which understates `2026-09-06T1103` (two stages) and the confirmatory
      run (three). It says fifteen methods, which is the plan; the two pilot runs it names
      ran fourteen. Corrections go **beside** the registration, never into it.
- [ ] **Link ORCID** to the registration so the record outlives an institutional email.
- [ ] **Confirm `main` is the GitHub default branch** (Settings → Branches).
- [ ] **Zenodo deposit** for the 3.8 GB vendored package tree, if it is to be preserved
      with a DOI. Git cannot hold it and should not.
- [ ] **Ivy GAP vital status** (`EXTERNAL_ACTIONS.md` item 4). Almost certainly
      unobtainable — the release does not publish it. Worth one message to the Allen
      community forum, but plan on survival staying BLOCKED. That is a result, not a gap.

---

## If only three things get done

1. **Neftel annotations** (B) — un-degenerates a published method and unlocks a robustness
   check that does not currently exist.
2. **The ISH rebuild** (D) — validates the yardstick itself, and the data is already local.
3. **The second tissue** (E) — the only item that lifts the statistical ceiling.

Everything else is improvement. These three change what the study can claim.
