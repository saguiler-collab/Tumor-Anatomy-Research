# PRE-SPECIFIED: does ACS reproduce an orthogonal protein ground truth?

**Written and committed 2026-09-15, BEFORE the experiment was run.** Nothing below was
adjusted after seeing a result. The commit date is the evidence.

## Why this is the most direct test the project has of its own question

The registered question is whether a tumour's own anatomy can substitute for ground truth when
choosing a deconvolution method. Until now the project has had **no external ground truth at
all** — the accuracy arm is synthetic pseudobulk built from the same atlas the ACS arm
deconvolves against, which is the shared-dependence problem recorded as C4.

Ajaib et al., *"GBMdeconvoluteR accurately infers proportions of neoplastic and immune cell
populations from bulk glioblastoma transcriptomic data"*, **Neuro-Oncology** 2023;25(7):1236–1248,
supplies one. They applied **single-cell-resolution imaging mass cytometry** — protein, 33
antibodies — to **ten IDHwt GBM samples (five paired primary/recurrent)** with **matched bulk
RNA-seq on the same tissue**, then scored deconvolution approaches against it by Pearson
correlation.

Their reported result, averaged across all samples (their Supplementary Tables S8/S9):

| approach | r vs IMC, immune | r vs IMC, neoplastic |
|---|---|---|
| MCPcounter + **GBM-specific** markers (MCP_GBM) | **0.37** | **0.43** |
| MCPcounter + **default** markers (MCP_default) | 0.27 | — |
| MCPcounter + **GBmap-derived** markers (MCP_GBmap) | **0.06** | 0.22 |
| CIBERSORTx | 0.05 | 0.02 |

**This is a near-controlled contrast of exactly the kind this project uses.** The top three rows
are **one algorithm** — MCPcounter — differing only in **which marker set** it is given. The
algorithm, the samples, the bulk data and the ground truth are all held fixed. What varies is
the source of the biological knowledge.

## The prediction

I will run **MCPcounter 1.2.0** (the genuine package) on the Ivy GAP anatomic cohort with the
same three marker sources, and score each with ACS.

> **PREDICTION: ACS will order them MCP_GBM > MCP_default > MCP_GBmap.**

**If ACS reproduces that ordering**, then anatomic concordance on one cohort has recovered a
ranking established by protein ground truth on a different cohort, by a different group, in a
different modality. That is the strongest evidence the anatomy test could receive, and it is
evidence no amount of internal benchmarking can supply.

**If ACS does not**, that is evidence against the registered hypothesis, and it will be reported
as such. Specifically:

- **MCP_GBmap scoring at or near the top would be the most damaging outcome.** It is the marker
  set IMC ranks near chance (0.06), and it is derived from the atlas this whole project uses. ACS
  preferring it would show ACS rewarding agreement with its own reference rather than with
  tissue.
- Any ordering that inverts MCP_GBM and MCP_GBmap falsifies the prediction outright.

## Exactly what will be run, fixed now

| | |
|---|---|
| algorithm | `MCPcounter::MCPcounter.estimate` v1.2.0, `featuresType = "HUGO_symbols"` |
| bulk | Ivy GAP anatomic samples only (`-reference-histology` suffix), the same 122 every other ACS number in this project uses |
| preprocessing | log2(x + 1) on the FPKM matrix, matching GBMdeconvoluteR's own log-scale input; declared, not tuned |
| MCP_GBM markers | Ajaib et al. 2022 GBM immune markers (183 genes) + Neftel four-state neoplastic markers, both as shipped in the GBMdeconvoluteR repo |
| MCP_default markers | MCPcounter's own published `genes.txt` |
| MCP_GBmap markers | this project's `select_signature_genes` applied to its own GBmap reference |
| scoring | `ivygap.anatomic.acs.score`, 10,000 permutations, 2,000 bootstraps, the frozen constraint file, unchanged |
| types compared | see below |

**No parameter will be tuned after seeing an ACS.** If a run fails, the failure is reported.

## The honest limits, stated in advance

1. **Different cohort and different question.** IMC r is measured across samples within a cell
   type on Leeds tissue; ACS is measured across anatomic structures within a cell type on Ivy
   GAP. Both are *within-type* comparisons, which is why they are comparable at all — but they
   are not the same measurement.
2. **The types the two can both see are Tumor and Macrophage_Microglia.** ACS constrains Tumor,
   Macrophage_Microglia, Endothelial and Oligodendrocyte. The IMC panel measured Tumor
   (AC/MES/NPC/OPC), Macrophages, Microglia, Monocytes, NK and T cells. **B cells were not in
   the IMC panel**, and Endothelial and Oligodendrocyte were not either. Per the instruction
   that prompted this work, the extra IMC populations (monocyte, DC, mast) are **not** forced
   into the eight-type roster.
3. **My MCP_GBmap is not exactly theirs.** Ajaib used marker genes from the GBmap preprint's
   supplementary data; I derive mine from this project's own GBmap reference with
   `select_signature_genes`. Same atlas, different derivation. The correspondence is real but
   looser than the other two rows, and no conclusion will rest on MCP_GBmap alone.
4. **MCPcounter returns abundance scores, not proportions.** They do not sum to one and are not
   comparable across cell types — only within a cell type across samples. This is a property of
   the method, not a defect, and it has a consequence worth stating: **ACS can score this method
   and the accuracy arm structurally cannot.** So this comparison cannot be cross-checked
   against the pseudobulk benchmark, and the registered primary outcome could never have
   included a method of this class.
5. **One ordering of three is weak evidence.** It is a single ordinal prediction with three
   entries. It is worth running because it is external, pre-specified and falsifiable — not
   because it is decisive.

## What this is not

It is **not** a use of ACS to select anything, and **not** a tuning of ACS. The constraint file
is frozen and hashed; no constraint, weight or sample rule is touched. ACS is the object being
tested here, exactly as the protocol intends.

---

# ADDENDUM, written before the experiment was run

While preparing the run I established a structural fact about the marker sets that changes how
the comparison must be scored. It is recorded here, before any ACS was computed, because
discovering it afterwards and quietly adjusting would destroy the point of pre-specifying.

## The fact

The seven frozen constraints need exactly four cell types: **Tumor** (C1, C7),
**Oligodendrocyte** (C2), **Endothelial** (C3, C4), **Macrophage_Microglia** (C5, C6). The four
marker sets do not cover the same four:

| marker set | Tumor | Oligodendrocyte | Endothelial | Macrophage_Microglia | constraints expressible |
|---|---|---|---|---|---|
| **MCP_GBM** (Ajaib immune + Neftel neoplastic) | yes | **no** | **no** | yes (TAM + Microglia) | C1, C5, C6, C7 |
| **MCP_GBM_moreno** (Moreno lvl3 + Neftel) | yes | yes | yes | yes | all 7 |
| **MCP_default** (MCPcounter's own 111 genes) | **no** | **no** | yes | approx. (Monocytic lineage) | C3, C4, C5, C6 |
| **MCP_GBmap** (this project's GBmap markers) | yes | yes | yes | yes | all 7 |

Scoring each arm on whatever it happens to cover would compare ACS values computed over
**different constraint subsets**, and D8 already established that C6 alone carries most of the
ranking. That comparison would be meaningless, and it would flatter whichever arm happened to
draw the easier constraints.

## How it is handled, fixed now

**Every pairwise comparison is scored on the constraints BOTH arms can express**, by setting the
non-shared cell types to NaN so the scorer drops exactly those pairs from the weighted
denominator. That uses the existing, tested partial-coverage path — the same one quanTIseq goes
through — and touches neither `constraints.py` nor its hash.

**The primary test becomes MCP_GBM vs MCP_GBmap on C1, C5, C6, C7** — Tumor and
Macrophage_Microglia, which are exactly the two cell types the IMC panel also measured, and the
pair the paper's contrast maps onto most directly (immune 0.37 vs 0.06; neoplastic 0.43 vs 0.22).

> **PRIMARY PREDICTION, restated and unchanged in direction: ACS(MCP_GBM) > ACS(MCP_GBmap) on
> the shared four constraints.**

**MCP_default is demoted to a secondary arm.** It cannot express Tumor at all, so it cannot join
the primary comparison; it is scored against the others on C5 and C6 only, and the three-way
ordering in the original prediction is therefore reported as two pairwise tests rather than one
ranking. This is a loss of power that I am declaring rather than papering over.

**MCP_GBM_moreno is added as a secondary arm** because it is the only GBM-specific set with full
constraint coverage, so it is the one that can be compared with MCP_GBmap on all seven. Its
result is secondary because the paper's 0.37 refers to the Ajaib markers, not Moreno's.

## What would still falsify the prediction

Unchanged: MCP_GBmap scoring at or above MCP_GBM on the shared constraints. Naming the adverse
outcome in advance is the whole point, and the coverage fix does not soften it.
