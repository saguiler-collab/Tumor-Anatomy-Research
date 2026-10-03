# Does the scale of the fitting loss decide the lymphoid ordering? -- a post hoc mechanism test

**Written 2026-10-02 07:05, after DESeq2 `unmix` returned T > B on TCGA-GBM with the frozen
signature** (`results/extension/tcga_gbm_frozen.json`: purity rho 0.780, B > T in 23.8%). Every
registered method returns B > T there. The hypothesis below was formed after seeing that, so this
test is **post hoc**. Its reading rule is fixed here, before any ablation run.

## Hypothesis

Least squares in linear space (NNLS, SVR, CIBERSORT-type methods) is dominated by the most highly
expressed genes -- in a tumour, the tumour's own genes -- so a poorly fitted tumour signal spills into
a minor column (B under SVR, §7h of WHY_B_OVER_T). `unmix` computes its loss on log(x + shift), which
evens out the genes' weights and lets low-abundance lymphoid markers count.

## Design

The genuine package throughout. Only `shift` and `power` change; same samples, genes and reference
as the registered TCGA-GBM frozen arm.

- **shift in {1 (the pre-declared choice), 10, 100, 1000, 10000}.** As shift grows, log(x + shift)
  becomes affine in x, so the loss approaches linear-space L1.
- **power in {1, 2} at shift 1.** Is it the scale of the loss, or L1 versus L2?
- GBM first; LGG if GBM reads either way.

## Reading rule (fixed now)

- **Loss scale implicated** if T > B holds at shift 1 and the cohort-level ordering turns to B > T at
  shift >= 1000 -- the near-linear end -- in GBM. The share of samples with B > T should rise
  monotonically with shift.
- **Not implicated** if T > B persists at shift 10000.
- **L1 versus L2:** if power 2 at shift 1 also gives T > B, the loss scale, not robustness to
  outliers, carries it.
- **Reported alongside:** purity rho at each setting, since a loss that gets the lymphoid ordering
  right must not be bought with worse tumour-content tracking.
- Exploratory; selects nothing; no setting replaces the pre-declared shift in any reported `unmix` row.
