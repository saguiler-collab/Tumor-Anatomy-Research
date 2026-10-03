# T against B in an independent single-cell atlas (Abdelfattah 2022, GSE182109) -- a non-methylation check of the lymphoid truth

**Written 2026-10-03 13:14:27 EDT** (the file's birth time; `date` in the same command printed
13:14:27). Exploratory, post-registration.

**What had been seen when this was written:** only the per-patient cell counts of two clusters, from
`Cluster_GBM.txt` and the existing cluster mapping (`data/reference/abdelfattah_2022/provenance.json`):
- C3 (T and NK together) against C11 (B): C3 is larger in 18 of 18 patients; pooled 19,570 against
  1,205 cells.
- No per-cell split of C3 into T and NK had been computed.

## Why

The headline's truth is methylation's T > B, measured by EpiDISH with a blood reference. GIMiCC
(`gimicc_truth_confirmation.md`) is a second methylation method. A single-cell atlas counts cells
directly: no methylation and no deconvolution. Its patients are a different cohort, and dissociation
can bias which cells survive. But T and B cells are both lymphocytes of similar size and fragility,
and both are CD45+, so CD45 enrichment cannot create a T : B ratio.

## The split, already declared

`prespecified/abdelfattah_cluster_mapping.md` (2026-10-01) declared the per-cell split of a mixed
T/NK cluster: the sign of (mean T panel - mean NK panel), in log1p-CPM per cell.
- T panel: CD3D, CD3E, CD3G, TRAC, CD2.
- NK panel: NKG7, GNLY, KLRD1, KLRF1, NCAM1.
- Counts are recovered from the processed matrix exactly as `scripts/build_abdelfattah_reference.py`
  recovers them (that identity check must pass again).

**Ties.** A cell with no detected marker of either panel has a difference of exactly 0. The declared
rule (`>= 0` -> T) would credit such cells to T.
- Here they are counted separately, as unresolved.
- The comparison uses strict T (difference > 0), which can only work against T > B.
- The declared-rule version is reported beside it.

B = every C11 cell. C12 (unmapped) is reported, never assigned.

## Readout and reading

Per patient and pooled: T (strict), NK, unresolved, B, and within-lymphoid shares (T, NK, B
renormalised), beside EpiDISH's TCGA within-lymphoid means.
- **Supports T > B** if pooled strict T > pooled B **and** strict T > B in a majority of the 18
  patients.
- **Contradicts** if pooled B >= T, or B >= T in a majority of patients.
- **Mixed** otherwise.
