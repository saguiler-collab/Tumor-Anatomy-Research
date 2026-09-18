# PRE-SPECIFIED: the immune arm — a second, independent ground truth

**Written and committed 2026-09-17, before any immune error was computed.** What had been
examined when this was written: the leukocyte-fraction file's coverage and value distributions,
its join to the expression matrices, and the paper's methods section defining the measurement.
None of that involves an error. The commit date is the evidence.

## Why this arm exists

The tumour arm compares an RNA-derived tumour proportion against **DNA copy number** (ABSOLUTE
purity). It answers *when does tumour deconvolution fail*. It says **nothing** about immune
accuracy — and immune infiltration is what deconvolution is most often deployed to estimate.

This project's pseudobulk arm reports that every method over-calls T cells by **4.6–7.9×** in
high-purity tumour. That is currently a **synthetic** observation: the mixtures were built from
GBmap's own cells, so it could be an artefact of the simulator. It has never been tested on
tissue, because there was no immune ground truth.

**There is now.** Thorsson et al. 2018, *Immunity* 48(4), estimated **leukocyte fraction** in
10,817 TCGA aliquots from **DNA methylation** — a third molecule, independent of both the RNA the
methods consume and the copy number the tumour arm uses.

> Their method, quoted: *"Overall leukocyte content … was assessed by identifying DNA methylation
> probes with the greatest differences between pure leukocyte cells and normal tissue, then
> estimating leukocyte content using a mixture model"* — 2,000 loci on HumanMethylation450,
> 1,000 in each direction.

## Coverage, measured

| cohort | expression samples | with a leukocyte fraction | LF median | LF range |
|---|---|---|---|---|
| TCGA-GBM | 174 | **141** | 0.1398 | 0.008–0.868 |
| TCGA-LGG | 530 | **530** | 0.0827 | 0.005–0.576 |

GBM is more immune-infiltrated than LGG, which is correct and is a free sanity check on a file
nothing has been tuned against.

## The outcome

    immune error(i, m) = [ Macrophage_Microglia + T_cell + NK_cell + B_cell ](i, m)
                         - leukocyte_fraction(i)

Roster immune types summed because leukocyte fraction is a **single aggregate** — it does not
resolve cell types, so neither may the comparison.

## THE DENOMINATOR PROBLEM, declared before any result

This is the part that could invalidate the whole arm if left implicit, so it is stated first.

**The two quantities do not have identical denominators.**

- **Leukocyte fraction** = leukocyte share of *all* nucleated cells contributing DNA.
- **Our estimate** = immune share of the *eight roster types*, and by **D12** it is an **mRNA**
  proportion rather than a cell proportion, because the conversion is the identity here.

So this compares **an RNA-derived proportion of a 8-type model** with **a DNA-derived cell
fraction of the whole sample**. Three specific gaps:

1. **The roster is incomplete.** It drops monocytes, DC, mast, neutrophils, neurons and others —
   **7.0%** of the GBmap atlas (`reference_coverage.json`). Some of those are leukocytes, so our
   immune sum is missing populations the methylation estimate counts. **This biases our estimate
   DOWNWARD relative to LF**, which is the opposite direction to the predicted over-call, making
   the prediction **conservative**.
2. **RNA versus cells.** Immune cells carry less mRNA than tumour cells (measured here: 10x
   Tumor/T_cell ≈ **2.245×**, `prespecified/mrna_content_gbmap_10x.csv`). An mRNA proportion
   therefore **understates** immune cell fraction — again conservative for the prediction.
3. **Methylation is itself an estimate** from a mixture model, with its own error, and it is not
   a gold standard.

**The tumour arm has the identical RNA-versus-DNA mismatch**, which is what makes the two arms
comparable to each other: any systematic RNA-proportion-versus-cell-fraction offset applies to
both. That is the basis on which a *difference in sign* between the arms can be interpreted at
all, and it is the reason this arm is worth running.

## The predictions, fixed now

**P1 — direction.** Methods **over-call** immune content on tissue: median error **> 0**.
Grounded in the pseudobulk 4.6–7.9× T-cell over-call. **Both denominator gaps above push the
other way, so a positive error would have to overcome them.**

**P2 — the asymmetry, and this is the one that matters.** The **signs are opposite**: tumour
under-called (established, all 14 methods), immune over-called.

> If P2 holds, deconvolution does not merely add noise — it **systematically redistributes signal
> from the malignant compartment into the immune compartment**, measured against two different
> molecules in the same samples.

**P3 — shared factors.** The biological factors predicting tumour error also predict immune
error, tested with the **same** five-factor model and the same Holm correction.

**P4 — recovery.** Immune recovery fraction, defined exactly as for tumour (`1 + slope`), is
**lower** than tumour recovery, because immune populations are a smaller share of signal and the
pseudobulk arm shows them at the detection floor.

## What falsifies what

- **P1 fails** if median error ≤ 0 → the pseudobulk over-call is a simulator artefact. That is a
  major correction to this project's own most-quoted clinical claim, and it would be published.
- **P2 fails** if both arms have the same sign → errors are a common scaling problem, not a
  redistribution. Simpler, and equally reportable.
- **P3 fails** if no factor is shared → failure mechanisms are compartment-specific.

## What this arm cannot do

Leukocyte fraction is **one aggregate number**. It cannot test whether *T cells specifically* are
over-called — only the immune compartment as a whole. The 4.6–7.9× figure is T-cell-specific and
this arm can only corroborate its direction at the aggregate level, never its magnitude per type.

## The CIBERSORT file is NOT ground truth, and will not be used as one

`TCGA.Kallisto.fullIDs.cibersort.relative.tsv` is CIBERSORT LM22 run on Kallisto-quantified TCGA
RNA. It is **another RNA deconvolution**, sharing every failure mode of the methods under test.
Using it as truth would be circular.

Its one legitimate use, declared here: a **comparison arm** — how does this project's panel
relate to the canonical published pan-cancer run on the same samples? Reported as
method-versus-method agreement, never as accuracy.

---

# ADDENDUM — per-cell-type immune truth from DNA methylation, prediction fixed before computing

Leukocyte fraction is one aggregate number, so the immune arm above could corroborate the
*direction* of an over-call and never a per-type magnitude. **EpiDISH** (already installed) closes
that with a DNA-methylation reference resolving **B, NK, CD4T, CD8T, Mono, Neutro, Eosino** from
333 HM450 CpGs. All 333 are present in the TCGA-LGG methylation matrix with no missing values.

Mapping to the roster: **CD4T + CD8T → T_cell**, **B → B_cell**, **NK → NK_cell**.

**`Macrophage_Microglia` is NOT mapped and will not be compared.** Microglia are brain-resident
and appear in no blood reference; `Mono` is a different population and substituting it would
invent a correspondence the reference does not have. That type is reported **not measurable by
this route**, which is a limitation of the method and not a result about it.

## The prediction

The immune arm found an unpredicted anomaly: **methods place more B cells than T cells in GBM**
(CIBERSORTx 0.046 vs 0.003) while GBmap itself holds **5× more T than B**. That observation had no
way to be checked, because no per-type truth existed. It does now.

> **PREDICTION: DNA methylation will show T cells > B cells.**

**Grounded, not guessed.** T cells are the dominant lymphocyte in glioma infiltrate; B cells are
rare. The reference atlas agrees (T 0.165, B 0.032). If methylation also says T > B while the
methods say B > T, the methods are placing lymphocyte signal in the wrong column — an
**identifiability failure**, not a proportion error.

**What would falsify it:** methylation showing **B ≥ T**. That would mean the methods may be
right and my "anomaly" was an artefact of comparing estimates against a reference composition
rather than against tissue. **That outcome would be reported and the anomaly withdrawn.**

## Stated limits

- EpiDISH's reference is a **blood** reference applied to **brain tumour** tissue. It estimates
  the leukocyte sub-composition; it does not know about microglia, and its absolute scale is a
  fraction of leukocytes rather than of all cells.
- The comparison is therefore **within the immune compartment** — relative composition — and not
  an absolute per-type fraction of the sample.
- This is LGG, not GBM. The anomaly was measured in GBM. **Whether it exists in LGG at all is
  checked first**; if it does not, this tests nothing about the GBM observation and will say so.
