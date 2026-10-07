# Per-sample truth from the same tissue: CPTAC glioblastoma, bulk RNA against single-nucleus RNA

**Written 2026-10-06 20:45:29 EDT** (the file's birth time). Exploratory, post-registration.

**At writing, no CPTAC data value had been seen.**
- The files were still downloading (`scripts/fetch_cptac_gbm.py`).
- Only GDC metadata was known: 18 cases; 17 with snRNA-seq filtered count matrices; 18 bulk STAR
  count files; 18 methylation beta files; GDC Seurat cluster tables for 17 cases.
- The paper's methods text had been read: "We performed snRNA-seq on 18 GBM samples using the same
  cryopulverized material from previous analyses."

## Why

Every truth in the study so far is either cohort-level (methylation, flow cytometry, atlases) or one
quantity per sample (ABSOLUTE purity). The user's question is whether computational methods give
*actual* results where direct measurement is unaffordable.

The decisive test is per-sample: bulk RNA and a direct count of cells, from the same piece of tissue.
CPTAC provides that for 17 glioblastomas. It also gives a second, independent test of the paper's
central question with a per-sample truth: does the method that gets the anatomy right (ACS, from Ivy
GAP) also get the composition right?

## The truth: single-nucleus composition per sample

- **Nuclei:** the deposited CellRanger QC-filtered barcodes, as GDC serves them.
- **Clusters:** GDC's own Seurat clustering (`seurat.analysis.tsv`), per sample, as deposited.
  - If that table has no cluster assignment for the barcodes, this arm is BLOCKED.
  - The data are never re-clustered here; a clustering chosen after seeing the data would be a free
    parameter.
- **Naming:** the rule declared on 2026-10-01 for the Abdelfattah atlas
  (`prespecified/abdelfattah_cluster_mapping.md`; `PANELS` in `scripts/build_abdelfattah_reference.py`),
  unchanged.
  - Per cluster, the mean log1p-CPM of each marker panel.
  - A cluster takes the roster type whose panel scores highest, only if it beats the runner-up by
    >= 0.5. Otherwise, or if the best panel is the excluded one (pericyte, plasma, mast, DC), the
    cluster is unmapped.
  - A T/NK cluster within 0.5 is split per nucleus by the sign of (T panel - NK panel).
  - Counts per nucleus come from the filtered matrix.
- **Per-sample truth:** each roster type's share of the **mapped** nuclei, over the 8 roster types.
  Unmapped nuclei are counted and reported, never assigned.
- **A sample is excluded (counted) if fewer than 70% of its nuclei map.** This threshold is fixed
  now, not tuned.

## The estimates

- **Bulk:** each case's GDC STAR counts, column `unstranded`.
  - The four `N_` summary rows are dropped; gene symbols are GDC's `gene_name` (GENCODE v36, the
    annotation of the TCGA bulk); duplicates are collapsed by mean; then CPM.
  - This is the TCGA bulk's construction (`build_tcga_bulk.py`), minus its log2 inversion, which these
    files do not need.
- **Methods:** the registered panel (`build_methods(prefer_r=True)`), under both reference builds,
  exactly as the TCGA arms run them (`absolute_purity_yardstick.py`):
  - the frozen signature, and the donor-level reference rebuilt from raw/X counts;
  - marker genes chosen by `select_signature_genes` from the reference alone.
- **Small-cohort change:** the TCGA arms skip a method with fewer than 50 finite estimates, and the
  cohort here is 17. The threshold becomes "every sample finite". Nothing else changes.

## Controls (read before any accuracy result)

1. **The truth must agree with DNA.**
   - Per sample, the snRNA tumour share against an independent DNA-based purity: GIMiCC's Layer-0
     purity (InfiniumPurify) from the CPTAC methylation of the same case.
   - **Pass if Spearman >= 0.40.**
   - If it fails, nuclei capture does not reflect tumour content here, and every accuracy reading
     below is INCONCLUSIVE.
2. **Mapping coverage:** the median share of nuclei mapped must be >= 70%.
3. **Negative control:** every per-method Spearman is given a permutation p (10,000 shuffles of
   which bulk sample is paired with which snRNA sample, seed `config.RANDOM_SEED`).

## Questions and readings, fixed now

- **C1, which compartments are recovered, per sample.**
  - For each method and roster type, and for leukocytes (Mac + T + NK + B) and lymphoid (T + NK + B):
    Spearman across samples between the estimate and the snRNA share.
  - Per compartment, the median across methods.
  - **Tracks** if the median is >= 0.40 (E2's bar); otherwise **not recovered**.
  - With n = 17 a single rho has a wide interval, so every value is reported with its permutation p.
- **C2, is the anatomy-chosen method the accurate one?**
  - Across the registered methods, Spearman between the anatomic ACS (the registered leaderboard) and
    each method's tumour-share accuracy (C1's rho for Tumor).
  - **Reading by the registered bar:** rho >= 0.60 means anatomy ranks methods by accuracy here;
    otherwise it does not.
  - This is a second test of the registered question, with a per-sample truth.
- **C3, does the reference build change per-sample accuracy?** Per method, C1's tumour rho under the
  frozen signature against the donor-level reference. Descriptive.
- **C4, lymphoid composition.**
  - Pooled and per sample: snRNA T against B against NK nuclei. The methods' cohort-mean
    within-lymphoid shares are set beside them.
  - **Reading:** the methods reproduce the snRNA direction if their T exceeds B where snRNA's does.
  - If fewer than 5 samples have any B nuclei, the per-sample B comparison is UNRESOLVED and only the
    pooled direction is read.
- **C5 (secondary), the per-sample methylation truth.** EpiDISH and GIMiCC on the same cases'
  methylation, per-sample T/(T+B) against snRNA's, where both are defined. This is the per-sample
  question OPEN_DEFECTS D23 could not answer. Descriptive.

## Limits, stated now

- These are nuclei, not cells: nuclear capture differs by cell type.
- 17 tumours, all glioblastoma.
- Bulk and nuclei come from the same cryopulverised material, but from different aliquots of it.
- Clusters are GDC's automated ones; the naming rule is this project's.

## Addendum 1 -- the declared T/NK split, implemented as declared (written 2026-10-06 20:49:57 EDT, the file's modification time after appending, copied from stat; the CPTAC download was still running and no CPTAC file had been opened)

A planted-data test (`tests/test_cptac_per_sample.py`) showed the T/NK split could never trigger.
- Where it came from: the code copied from `build_abdelfattah_reference.py` applies the 0.5 margin
  between the best panel and the runner-up **first**.
- When T and NK are within 0.5 of each other, the runner-up *is* the other lymphoid panel, so the
  cluster is discarded as ambiguous before the split is reached.
- The rule's text (2026-10-01) says otherwise: "If T and NK fall in one cluster ... cells of that
  cluster are split per cell by the sign of (mean T panel - mean NK panel)".

**Implemented as declared:**
- **Clusters whose best panel is T or NK:** the margin is between max(T, NK) and the best panel
  outside {T, NK}.
  - Below 0.5: unmapped.
  - Otherwise, if |T - NK| >= 0.5, the larger of the two.
  - Otherwise, split per nucleus by the sign rule (>= 0 is T).
- **Every other cluster:** unchanged.

**The Abdelfattah atlas is unaffected,** as recomputed from the panel scores recorded in its provenance.
- Its only T/NK cluster, C3, has T 4.88 against NK 1.24, so |T - NK| = 3.64.
- No other cluster's best panel is T or NK.
- The builder is corrected the same way and logged as OPEN_DEFECTS D28 (latent, no effect).

## Addendum 2 -- two flaws in the truth, found by inspecting it before any accuracy was computed (written 2026-10-06 20:58:55 EDT, the file's modification time after appending, copied from stat; the methods stage had not finished, and no method's estimate had been compared with the truth)

Inspection of the registered truth (`results/cptac/snrna_truth_detail.json`) showed two problems.

**1. The declared NK panel labels non-immune nuclei NK in brain tissue.**
- In C3N-03186 (clusters 6, 7, 9; 794 nuclei) and C3N-02783 (cluster 13; 104 nuclei), the "NK" call
  rests entirely on NCAM1 (mean log1p-CPM 5.5-6.1).
  - NKG7, GNLY and KLRF1 are about 0.
  - PTPRC (CD45) is about 0.2: these nuclei are **not leukocytes**.
- Genuine NK clusters elsewhere (C3N-01798 cluster 15, C3N-03188 cluster 15) show NKG7 3-7, GNLY 4-6
  and PTPRC 7.9.
- NCAM1 (CD56) is expressed by neurons and glioma cells. The panel was declared for an
  immune-enriched atlas (Abdelfattah), where this never arose: no cluster there was NK-labelled.
- Logged as OPEN_DEFECTS D29.

**2. Myeloid nuclei are probably absorbed into tumour-labelled clusters.** Several samples have no
myeloid-labelled nuclei at all, which is implausible for glioblastoma. GDC's automated clusters
appear to merge them.

**Decisions:**
- **The registered truth stays primary.** Every reading is made on it first.
- **Sensitivity truth S1** (declared now): the NK panel without NCAM1 (NKG7, GNLY, KLRD1, KLRF1); every
  other part of the rule unchanged. C1 and C4 are reported under both.
  - Where S1 and the registered truth disagree on NK, leukocytes or lymphoid, the S1 value is the
    biologically valid one, and the text says so.
- **Diagnostic M1** (declared now; it does not change the truth): per sample, the share of nuclei in
  tumour-labelled clusters whose own myeloid panel exceeds their tumour panel by >= 0.5.
  - **If M1 exceeds 10% in any scored sample, the Macrophage_Microglia and leukocyte truths are
    flagged as underestimated,** and their C1 readings are reported as unreliable.
