# The two problems: RNA contribution, then cellular abundance

**Decision recorded 2026-09-15.** This is the answer to the question `OPEN_DEFECTS.md` D12
forced — *does this study report mRNA proportions or cell proportions?* — and it is a better
answer than either option D12 offered.

> **Methodological deconvolution is reported in mRNA proportions. Clinical translation happens
> in cell proportions.** These are two nested problems, and the study is clearer when they are
> separated than when one is silently substituted for the other.
>
> **Problem 1.** Can we accurately infer the RNA contributions?
> ↓
> **Problem 2.** Can we convert those contributions into actual cellular abundance?
>
> **The research question:** *which factors govern the robustness of bulk RNA deconvolution, and
> under what conditions can RNA-level estimates be translated into cellular abundance?*

---

## Why this is the right decision and not merely a convenient one

D12 found the mRNA-to-cell conversion had never been applied — `cell_size` is 1e6 for every type
on the reference every number was computed against, so `to_cell_fractions` is the identity. The
two options I put forward were to start applying it, or to declare that the study reports mRNA
proportions.

Both were wrong in the same way: they treated one quantity as the answer. The methods in this
panel **do not agree about which quantity they report** — measured, not assumed:

| package | what it returns | measured on the probe |
|---|---|---|
| MuSiC, SCDC, EPIC, BayesPrism | mRNA share | 0.7502–0.7503 |
| **Bisque** | **cell share** | **0.5000** |
| NNLS, SVR, CIBERSORTx, Elastic Net | mRNA share | 0.7500–0.7743 |

*(true mRNA share 0.7500, true cell fraction 0.5000; `results/cell_size_semantics_*.json`)*

So "the estimate" is not one thing, and a benchmark that scores all of them against one truth is
measuring units as much as accuracy. **D15** is that error in its concrete form: mRNA-share
estimates scored against cell-fraction truth, which is exactly what `pseudobulk.py`'s own
docstring warns against.

Separating the problems makes the disagreement visible instead of absorbing it.

---

## What each problem is, operationally

### Problem 1 — can the RNA contributions be inferred?

**Estimand.** Each cell type's share of the bulk mRNA pool.

**Truth.** `PseudobulkSet.truth_mrna`, added 2026-09-15: the same mixtures' composition expressed
as each type's share of the mRNA, computed from the cells actually drawn.

**This is what the study has in fact been measuring all along**, because the conversion was the
identity. Reporting it as Problem 1 is therefore not a change of result; it is a correct
description of the existing one.

**What the project already knows about Problem 1**, and it is the substantive part:

- **The ordering is a property of the reference atlas, not of the methods.** Holding the atlas and
  changing the sequencing platform preserves it (rho 0.817); holding the platform and changing the
  atlas destroys it (0.221). `RESULTS.md` §5d, D14.
- **Separation from the negative controls is robust** under every reference tested. No control
  approaches a real method in any arm.
- **One design property measurably helps: weighting genes by cross-donor consistency.**
  Regularisation, batch correction and a donor-level hierarchy buy nothing on this cohort.
  `docs/CLINICAL_READINESS.md` §2.
- **Accuracy collapses where tissue actually sits.** Error grows monotonically with tumour
  content, and at high purity every method under-calls tumour by 0.33–0.51 — with the missing
  mass appearing as **T cells**, over-called 4.6–7.9×.

### Problem 2 — can those contributions be converted into cellular abundance?

**Estimand.** Each cell type's share of the cells.

**Truth.** `PseudobulkSet.truth`, the cell-fraction truth that has been there all along.

**The conversion.** Divide each type's mRNA share by its mRNA content per cell, then renormalise.
That needs a per-type mRNA-content vector, and **this is where the project's real difficulty
lives**:

- **The conversion has never been applied** (D12), because `run_benchmark` rebuilds the reference
  from an already-normalised matrix without passing `cell_totals`.
- **When it is available, it is not small.** The frozen factors span **737×** (B cell 3,624 to
  Astrocyte 197,731). The Darmanis reference gives a 1.63× spread; GBmap's assay subsets 1.58–1.60×.
  Three sources, three answers.
- **It is unrecoverable from some data entirely.** Neftel's published matrix is TPM, already
  normalised per cell, so per-cell mRNA content is gone at source. That reference cannot support
  Problem 2 at all, and says so in its provenance.
- **The renormalisation is composition-dependent**, so the conversion is not a per-type constant
  in its effect — it can reorder a cell type across samples. This was retracted once and then
  confirmed; see D1.

**So Problem 2 is not a post-processing step. It is a measurement problem with its own error
budget, its own missing data, and no agreed source of truth.** That is the finding, and it is why
"under what conditions" belongs in the research question rather than being assumed away.

---

## What this framing changes, concretely

| | before | after |
|---|---|---|
| what the leaderboard reports | "cell fractions" | **mRNA proportions** — Problem 1 |
| what the accuracy arm scores against | cell-fraction truth | **`truth_mrna`** for Problem 1; `truth` for Problem 2 |
| the cell-size conversion | a silent, inert step | **the object of study in Problem 2** |
| D12 | a defect to fix | a **finding**: the conversion was never applied, and that is informative |
| Bisque's anomaly | unexplained | the one method already in cell-fraction units (D15) |

**It also reframes D14 without softening it.** "The ordering is a property of the atlas" is an
answer to *which factors govern robustness* — the first half of the research question. It is a
result, not only a limitation.

---

## What must happen, and the rules that apply

1. **Re-run the benchmark scoring both arms** — Problem 1 against `truth_mrna`, Problem 2 against
   `truth` after a declared conversion. Publish before-and-after per method; the accuracy ranking
   may move and that must be visible.
2. **Handle Bisque explicitly.** It returns cell fractions, so under Problem 1 it is the one
   method needing conversion *backwards* to be comparable, or it must be reported separately with
   the reason. Silence advantages it under Problem 2 and penalises it under Problem 1, both times
   for a reason unrelated to its accuracy.
3. **Pick the mRNA-content vector before seeing what it does to the ranking**, and record the
   choice and its provenance. Three sources give three spreads; choosing among them by outcome
   would be the selection this project forbids.
4. **Do not retrofit the registration.** It says the conversion is "applied once and centrally".
   It was not. That correction is `docs/CORRECTIONS_REGISTRATION.md` C3, and this framing is an
   addition beside it, not an amendment to it.

---

## Why this is a stronger paper than the one originally planned

The registered study asked whether anatomy can choose a deconvolution method. On its own terms the
answer came out qualified: anatomy separates good methods from bad robustly, and orders them only
within one atlas.

The two-problem framing asks something that the evidence actually supports and that matters more:
**bulk deconvolution has two distinct failure modes, and the field routinely reports the second as
though the first were solved.** This project can show, with measurements rather than argument,
that:

- the RNA-level ordering depends on the reference atlas (D14);
- the RNA-to-cell conversion is frequently absent, inconsistent across sources by up to 737×, and
  sometimes unrecoverable (D12, D15);
- published packages silently disagree about which quantity they return (D1);
- and the anatomic constraints that validate the RNA level hold on measured composition in a
  second and third specimen (§5b, §5c), while being blind to the largest error mode (T cells).

None of those depends on the leaderboard's ordering being stable, which is the part that did not
survive.
