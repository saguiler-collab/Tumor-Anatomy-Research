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
- [x] **Neftel GSE131928** — 10x (16,201 cells / **9 patients**, 21 samples) and
      SMART-Seq2. **No cell-type labels** — confirmed against all four GEO
      supplementary files, so it cannot yet serve as a reference.
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

- [x] **ISH constraint check — DONE 2026-09-10.** `scripts/ish_constraint_check.py`,
      paired within `sub_block_id`, all 18,778 rows, per marker, 2,000-draw within-block
      permutation null. **8 of 17 marker-constraint tests reach p < 0.05.**

  ESM1 supports C3 (MVP > CT) at 1.000 over 8 blocks / 6 donors, p = 0.005, and C4
  (MVP is the maximum) at 0.750, p = 0.021. CD163 supports C6 (MVP > CT) at 1.000 over
  7 blocks, p = 0.0085. So the endothelial and perivascular-myeloid claims — the
  strongest in the set on neuropathology — now have independent measured support from
  data that never touched a deconvolution method.

  Three limits, all reported rather than smoothed: C1 and C7 are **marker-dependent**
  (CD44 agrees at 0.900 over 137 blocks, SOX2 and PTPRZ1 disagree at 0.000 and 0.056 —
  all three are real tumour markers and none measures tumour cell *density*); C5 has one
  evaluable block; C2 is unmeasurable, no myelin marker exists in the panel and OLIG2 was
  refused because in glioma it marks tumour cells and OPCs. Sample sizes for the
  endothelial and myeloid results are 6-8 blocks, which is thin and must be said.

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

## F2 · Two validations that need no new tissue — added 2026-09-10

- [ ] **Concordance with published independent benchmarks.** ← *cheapest high-impact item on this page*

  Two peer-reviewed benchmarks already ranked many of these same methods against real
  ground truth, on different tissues, by different authors:

  > Avila Cobos F, Alquicira-Hernandez J, Powell JE, Mestdagh P, De Preter K.
  > Benchmarking of cell type deconvolution pipelines for transcriptomics data.
  > *Nat Commun* 11:5650 (2020). doi:[10.1038/s41467-020-19015-1](https://doi.org/10.1038/s41467-020-19015-1)

  > Sturm G, Finotello F, Petitprez F, et al. Comprehensive evaluation of
  > transcriptome-based cell-type quantification methods for immuno-oncology.
  > *Bioinformatics* 35:i436 (2019). doi:[10.1093/bioinformatics/btz363](https://doi.org/10.1093/bioinformatics/btz363)

  The question this answers is the study's own question, asked once more with someone
  else's data: **does the ACS ranking agree with a ranking produced independently, by
  other people, on other tissue, against real ground truth?** Agreement is external
  validation of ACS that costs one afternoon and no new data. Disagreement is equally
  informative and equally publishable.

  Method overlap is partial — the published benchmarks cover CIBERSORT, EPIC, quanTIseq,
  MuSiC, DWLS, Bisque and SCDC to varying degrees — so this is a rank correlation over
  whatever overlaps, with the overlap stated. It must be pre-specified which methods and
  which published metric are used, before looking, or it becomes cherry-picking.

- [ ] **Spatial transcriptomics as a third tissue.** Visium spots are small enough that a
      pathologist annotates regions on the same section, so anatomic labels and expression
      come from the same physical tissue rather than from adjacent blocks.

  > Ravi VM, Will P, Kueckelhaus J, et al. Spatially resolved multi-omics deciphers
  > bidirectional tumor-host interdependence in glioblastoma. *Cancer Cell* 40:639 (2022).
  > doi:[10.1016/j.ccell.2022.05.009](https://doi.org/10.1016/j.ccell.2022.05.009)

  *Buys:* a second GBM cohort without leaving glioblastoma, so the GBM constraint file
  applies unchanged — no new constraint set, no new registration, no roster change. That
  makes it cheaper than the Allen brain route, though it is a different modality and the
  deconvolution unit becomes a spot rather than a microdissected block.

  *Costs:* spot-level deconvolution has its own literature and its own failure modes, and
  a Visium spot is 55 microns — a handful of cells, not a tissue block. Whether ACS is
  even well defined at that scale needs thinking through before any data is downloaded.

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
