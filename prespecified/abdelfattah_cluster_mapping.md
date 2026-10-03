# PRE-SPECIFIED: naming the Abdelfattah 2022 (GSE182109) clusters, and the test they serve

**Written 2026-10-01 at ~10:30 EDT, before any GSE182109 expression value was on this machine.**
At the time of writing only `GEO_GSE182109/Cluster_GBM.txt` (cluster ids C1-C12 per cell, no
cell-type names) and `genes.tsv` (gene names, no values) were present. Nothing below was chosen
after seeing expression.

## Precedence

1. **The authors' own cell-type labels**, if the study's metadata file is obtained (Single Cell
   Portal metadata or the paper's Supplementary Data). Used as published; their names are mapped
   onto this project's roster exactly as GBmap's are (`ivygap.config.GBMAP_CELL_TYPE_MAP` logic).
2. **Only if (1) is unavailable**, the marker rule below.

## The marker rule (used only without the authors' labels)

For each cluster, the mean log1p-CPM of each panel over all its cells; the cluster takes the roster
type whose panel scores highest, **only if** that score exceeds the runner-up by at least 0.5
(log scale). A cluster failing the margin, or whose best panel is "excluded", is left unmapped and
dropped -- never forced into a roster column.

| roster type | panel (fixed) |
|---|---|
| Tumor | SOX2, PTPRZ1, EGFR, NES, PDGFRA |
| Macrophage_Microglia | CD14, CD68, AIF1, C1QA, CSF1R |
| T_cell | CD3D, CD3E, CD3G, TRAC, CD2 |
| NK_cell | NKG7, GNLY, KLRD1, KLRF1, NCAM1 |
| B_cell | MS4A1, CD79A, CD79B, CD19, BANK1 |
| Endothelial | PECAM1, VWF, CLDN5, FLT1 |
| Oligodendrocyte | MBP, PLP1, MOG, MOBP, MAG |
| Astrocyte | AQP4, GJA1, ALDH1L1 |
| *excluded (no roster home)* | pericyte RGS5, PDGFRB; plasma JCHAIN, MZB1; mast TPSAB1, CPA3; DC CLEC9A, FCER1A |

**If T and NK fall in one cluster** (common at this resolution), cells of that cluster are split
per cell by the sign of (mean T panel - mean NK panel), declared here, not chosen after.

## The question

Does the lymphoid inversion persist when the reference's lymphoid profiles come from an atlas that
is **not** one of GBmap's 16 source studies? The reference is built with this project's own
builder rules (CPM per cell, capped per donor and type, seeded), patients as donors. TCGA GBM and
LGG are refit with the registered NNLS and SVR on the frozen-arm gene-selection rule, and the
lymphoid ordering is compared with DNA methylation exactly as in `scripts/lymphoid_ordering.py`.

**Reading, fixed now:** if T > B appears under the independent atlas where GBmap gives B > T, the
inversion is a property of GBmap's B (or T) profile. If B > T persists, it is not specific to
GBmap's annotation or donors -- it lies in how bulk tumour signal meets any droplet-derived
lymphoid profile.

## Addendum, written 2026-10-01 11:53 EDT -- before the independent-atlas test's result was seen

The built reference's B column (cluster C11, B cells together with plasma cells, mapped by the rule
above) has **62 immunoglobulin genes admitted to its marker space, carrying 86.7% of the column's
mass**, and no ribosomal genes. Plasma-cell antibody transcripts are abundant in bulk tissue, so
under this atlas a B estimate partly measures immunoglobulin output.

Declared now, before the main result is visible: a **sensitivity refit with immunoglobulin genes
removed** (the rule ^IG[HKL][VDJC] | ^IGH[GAMDE] | ^IGKC | ^IGLC | ^IGLL | ^JCHAIN, the same rule used
for GBmap's signature) on the same samples, methods and comparison. Reading: if the main result and the
Ig-removed result agree, immunoglobulin content is not driving it; if they differ, the Ig-removed
result is the one that speaks to B cells, and the difference is reported as an immunoglobulin effect.

## Addendum 2, written 2026-10-01 15:58 EDT -- AFTER the main and immunoglobulin-removed results were seen

**Seed sensitivity of the independent reference.** The build samples at most 50 cells per
(patient, type) with one seed. The main result is MIXED. Its most consequential cell -- LGG
NNLS ordering T above B where GBmap gives B above T -- could be a property of that one draw.
Rebuild with cell-sampling seeds 1, 2 and 3 (everything else unchanged; one shared streaming pass:
`build_abdelfattah_reference.py --seeds 1 2 3`). Then refit the same samples and methods
(`independent_atlas_test.py --seed k`).

**Reading, fixed before any seed is run.** A cohort-method cell is **seed-stable** if its
cohort-level ordering matches the registered build's under all three seeds, and **seed-dependent**
otherwise. Seed-dependent cells are reported as such. This was declared after the main result was
seen, so it can only qualify the main reading, never upgrade it: a seed-stable cell keeps its
reading, and a seed-dependent one loses it.
