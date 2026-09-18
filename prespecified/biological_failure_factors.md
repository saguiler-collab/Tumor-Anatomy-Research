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

---

# ADDENDUM, written before the model was fitted

Building the covariate table surfaced three structural facts. None involves the outcome — they
are properties of the predictors alone — and all three are recorded here, before any fit, because
finding them afterwards and adjusting would destroy the point of pre-specifying.

## 1. F6 and F7 are severely underpowered, and nearly the same variable

| factor | non-null | positives |
|---|---|---|
| is_mesenchymal (F5) | 150 | **52** (34.7%) |
| **idh1_mutant (F6)** | 148 | **8** (5.4%) |
| **gcimp (F7)** | 151 | **8** (5.3%) |

And they are not independent: **Spearman +0.868**, with 7 of 8 samples carrying both.

That is biologically correct — IDH1 mutation *causes* the G-CIMP phenotype — but it means F6 and
F7 are close to a single variable measured twice, and eight positives cannot support separating
them. **Neither the model nor the factor list is changed**; what is fixed here, in advance, is
the rule for reading them:

> **A null result for F6 or F7 is uninformative** and will be reported as *underpowered*, not as
> evidence of absence. A positive result for either will be reported as a single
> IDH-mutant/G-CIMP signal rather than as two independent findings, because these data cannot
> distinguish them.

The honest version of this test needs an IDH-mutant-enriched cohort. **That is exactly what
TCGA-LGG is**, where IDH mutation is the majority class rather than 5% — so F6/F7 move from
"cannot be tested here" to "the natural question for the validation cohort", which strengthens
the case for LGG rather than weakening it.

## 2. Ploidy and genome doublings are correlated (+0.619), as expected

F3 is close to the discrete version of F2. The joint model handles this — they share variance and
the coefficients split it — so **neither will be interpreted alone**; the pair is read together.

## 3. The design is nonetheless well conditioned

Condition number of the standardised design matrix with all seven factors: **6.6**. Collinearity
is present and is not severe enough to destabilise the fit; the concern above is power, not
numerical conditioning.

## What does not change

The seven factors, their directions, the joint OLS, the sign test across twelve methods, Holm
across seven, alpha = 0.05, and the falsification criterion. All as written above.

---

# ADDENDUM 2 — how F5 is carried into LGG, declared before it is related to any error

F5 was pre-specified as the **Verhaak categorical call** `EXPRESSION_SUBTYPE == Mesenchymal`.
That classification **does not exist for TCGA-LGG** — lower-grade glioma is classified by IDH
status and 1p/19q codeletion — so F5 cannot be validated in LGG with the same variable.

Substituting a different variable and presenting it as the same pre-specified test would be the
quiet swap this project refuses. So:

**F5's validation in LGG is a LABELLED SECONDARY ANALYSIS**, using a continuous mesenchymal
program score computed identically in both cohorts from the Neftel four-state neoplastic markers
already vendored in the GBMdeconvoluteR drop (`scripts/mes_score.py`).

**Score definition, fixed here:** per sample, the mean rank-percentile of a state's marker genes
within that sample's own expression profile, then `MES − mean(AC, NPC, OPC)`. Within-sample
percentiles make it independent of library size and comparable across cohorts with no
cross-cohort normalisation.

**The hypothesis and the direction are unchanged.** Mesenchymal character should depress the
tumour estimate. Only the measurement changes, and the change is declared rather than absorbed.

## The score was validated against the Verhaak call, NOT against error

On the 150 GBM samples carrying a Verhaak class:

| variable | MES-classified median | other subtypes | one-sided p | AUC |
|---|---|---|---|---|
| `MES` raw | 0.9671 | 0.9531 | 6.5e-08 | 0.762 |
| **`MES − mean(AC,NPC,OPC)`** | **0.0316** | **0.0040** | **1.3e-11** | **0.831** |

**The contrast form is chosen, and the basis for choosing it is stated so it cannot be
mistaken for outcome-selection.** It was selected because (a) it agrees better with the
Verhaak call that F5 was originally defined on, and (b) subtracting the other three states
controls for MES markers simply being highly expressed genes, which is a reason that exists
independently of any result. **It was not compared against deconvolution error before being
chosen.** The outcome of this study is error; agreement with a subtype label is variable
definition, not outcome selection.

AUC 0.831 is good agreement and is not identity, which is expected — a continuous program score
and a categorical classifier are different constructs. That imperfection is why this is
secondary.

---

# ADDENDUM 3 — the spatial layer, with its direction fixed before computing

If mesenchymal character predicts deconvolution error (GBM discovery), the natural next question
is **where in a tumour that character lives**. Ivy GAP is the only cohort here with H&E-guided
anatomic labels, so it can answer it — and it is used for *interpretation*, not for statistical
validation, because nine tumours cannot carry that weight.

**Prediction, fixed before the score is related to structure:**

> **The mesenchymal program score is HIGHEST in pseudopalisading cells around necrosis (PAN).**

**Mechanism, stated so the prediction can fail.** The MES-like malignant state is established as
hypoxia-associated (Neftel et al. 2019). PAN is, by definition, the rim of viable tumour cells
palisading around a necrotic core — the most hypoxic niche the Ivy GAP roster contains. If MES
character is hypoxia-driven it should peak there. Microvascular proliferation (MVP) is the
plausible runner-up, being a hypoxia-driven angiogenic response.

**A failure of this prediction is reportable and does not touch the GBM finding**, which stands
on its own cohort with its own ground truth.

## Design, and it obeys the same invariants as everything else here

- The score is the **identical** `MES − mean(AC, NPC, OPC)` contrast used in GBM and LGG
  (`scripts/mes_score.py`), computed on the Ivy GAP anatomic bulk.
- **Samples are nested in tumours**, so estimates are collapsed to one value per
  (tumour, structure) before anything is aggregated — the same rule ACS uses.
- The null is **within-tumour permutation of structure labels**, 10,000 draws. Not a t-test
  across samples, which would treat multiple blocks from one tumour as independent.
- Only the five pre-registered primary structures are scored, and only H&E-labelled samples;
  the 148 expression-labelled ISH samples are excluded as circular, exactly as in ACS.

## What it can and cannot show

**Can:** whether the biological property associated with error is spatially concentrated, which
would give the failure mechanism an anatomical location and a reason a clinician would care.

**Cannot:** establish that deconvolution error itself is higher in that niche. Ivy GAP has no
per-sample ground truth, so error cannot be measured there at all. The link between niche and
error is **inferential, through the score**, and will be stated that way rather than implied.
