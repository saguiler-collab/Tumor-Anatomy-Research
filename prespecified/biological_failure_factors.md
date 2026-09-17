# PRE-SPECIFIED: which biological properties of a glioblastoma predict deconvolution error

**Written and committed 2026-09-17, BEFORE any relationship between error and any factor below
was examined.** What had been looked at when this was written, and nothing else: the marginal
distributions of the candidate factors (ranges, medians) and the join coverage between tables.
Neither involves the outcome. The commit date is the evidence.

## The hypothesis

> **Biological properties of GBM that alter the relationship between malignant-cell abundance
> and transcriptomic composition will predict systematic errors in bulk RNA deconvolution.**

This names a mechanism rather than a list. Any individual factor below may die without taking
the hypothesis with it, and the hypothesis may die while individual factors survive. That is the
point of stating it this way.

## Why this question and not the registered one

The registered question — can anatomy substitute for ground truth when choosing a method — has
been answered: it detects but does not rank (`docs/ENDPOINT.md`). That answer is a **null**, and
it motivates this one. If the evaluation cannot tell us *which* method to trust, the useful
question becomes *when* any of them can be trusted. That is answerable here because ABSOLUTE
purity gives a per-sample, DNA-measured ground truth for the one column that matters clinically.

## Design

**Discovery cohort:** TCGA-GBM, **154 samples** with an ABSOLUTE DNA purity call, deconvolved by
**12 comparable methods** against the frozen signature (`results/absolute_purity_yardstick.json`).

**Outcome, primary:** the signed per-sample tumour-fraction error

    e(i, m) = estimated Tumor fraction (sample i, method m) - ABSOLUTE purity (sample i)

**Outcome, secondary:** |e|, the error magnitude, for factors expected to degrade accuracy
without a directional bias.

## Candidate factors, with a direction fixed in advance

Direction is stated so the test can fail. A factor whose effect is significant in the OPPOSITE
direction counts as a **failed** prediction, not a discovery.

| # | factor | source | pre-specified direction | mechanism |
|---|---|---|---|---|
| F1 | **tumour purity** | ABSOLUTE | error becomes **more negative** as purity rises | the pseudobulk arm already shows a 0.33–0.51 tumour under-call at high purity; if real, it reproduces against DNA |
| F2 | **ploidy** | ABSOLUTE | **\|e\| increases** | aneuploidy scales transcript output away from the diploid reference profile |
| F3 | **genome doublings** | ABSOLUTE | **\|e\| increases** | whole-genome doubling is the discrete version of F2 |
| F4 | **subclonal genome fraction** | ABSOLUTE | **\|e\| increases** | one profile per cell type is the model's core assumption; subclonal heterogeneity violates it |
| F5 | **Verhaak expression subtype = Mesenchymal** | cBioPortal `gbm_tcga_pub2013` | error becomes **more negative** | MES tumours are transcriptionally myeloid-like, so malignant signal should leak into the macrophage column |
| F6 | **IDH1 mutation** | cBioPortal | **\|e\| increases** | IDH-mutant GBM is a biologically distinct entity under-represented in the reference atlas |
| F7 | **G-CIMP methylation** | cBioPortal | **\|e\| increases** | global methylation reprogramming shifts expression away from the reference |

**MGMT status is deliberately EXCLUDED** from the primary set: it is a treatment-response marker
with no mechanistic path to a composition-estimation error, and it has the weakest coverage
(115 of 154). Testing it would add a multiple-comparison penalty for a factor we have no reason
to expect. It is reported descriptively only.

## The statistical model, fixed now

Per method, an ordinary least squares fit of the outcome on **all seven factors jointly**,
because they are correlated with each other and a marginal test would attribute shared variance
arbitrarily:

    e ~ F1 + F2 + F3 + F4 + F5 + F6 + F7

**A factor counts as supported only if its coefficient is consistently signed in the
pre-specified direction across independent methods.** Twelve methods are twelve semi-independent
replications; a real biological effect should not depend on which solver was used.

**Primary statistic per factor:** a two-sided sign test on the 12 method-level coefficients
against the null that positive and negative signs are equally likely, plus the median
coefficient and its bootstrap interval.

**Multiple testing:** Holm correction across the **seven** primary tests. alpha = 0.05.

**Degenerate methods** (MuSiC, SCDC ENSEMBLE, Bisque on this reference) are retained in the
primary analysis and the whole analysis is repeated without them, because a method that has
silently reduced to another contributes a duplicate replication.

## What would falsify the hypothesis

**No factor reaching Holm-corrected significance with a consistent sign across methods.** That
outcome is publishable and will be published: it would mean deconvolution error in GBM is not
predictable from the tumour's measured biology, which is itself worth knowing and would close
the question rather than leave it open.

**Partial support is the most likely outcome** and is not a licence to re-frame. If F1 is strong
and F2–F7 are null, the finding is "error tracks purity and nothing else we measured" — stated
plainly, with the nulls reported at the same prominence as the hit.

## External validation, specified before discovery so it cannot be tuned

If and only if at least one factor survives Holm correction with a consistent sign, the
**identical model** is applied to **TCGA-LGG**. The factors, directions, model and statistic
are those fixed above. LGG is not used to discover anything and no factor may be added,
dropped or re-signed after seeing GBM.

Reproduction in LGG would show the mechanism is a property of glioma deconvolution. Failure to
reproduce would localise it to GBM. **Both are reportable and neither is a reason to change
this document.**

## What this analysis cannot do, stated in advance

- **Only the Tumor column has DNA ground truth.** Full-composition error remains measurable only
  on synthetic pseudobulk. Nothing here speaks to immune-cell accuracy.
- **Purity confounds with almost everything** — with subtype, with ploidy, with sample handling.
  The joint model addresses this only as far as linear adjustment can.
- **ABSOLUTE purity is itself an estimate**, from copy number, with its own error. It is an
  orthogonal measurement, not a gold standard.
- **The reference is GBmap-derived for every method**, so this measures methods-on-a-GBmap-
  reference, and cannot separate "the biology breaks deconvolution" from "the biology is
  under-represented in GBmap".
