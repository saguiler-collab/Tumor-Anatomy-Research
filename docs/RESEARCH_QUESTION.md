# The research question, and whether this study answers it

Verified 2026-09-24, before the manuscript was written. The purpose is to check that the
question asked, the hypothesis registered, the test run and the conclusion reported are the
same object — because a mismatch between a registered question and a reported answer is the
first thing a reviewer or a judge looks for, and the easiest thing to fix before writing and
the hardest after.

---

## What was registered, verbatim

From `Anatomy_Test.md`, the protocol:

> **The question.** Can a tumor's own anatomy stand in for ground truth when choosing a
> cell-type deconvolution method — and does a method that gets the anatomy right also get the
> biology right?

> **Hypothesis.** Deconvolution methods that better reproduce the cell-composition gradients a
> pathologist would predict across a tumor's anatomic structures are also more accurate
> against independent ground truth — making anatomic concordance a usable substitute for
> ground truth in tissues that have none.

> **Null it can prove.** That "the answer matches known biology" is not evidence that the
> answer is right. Which is a finding, because that check is exactly what most papers
> currently rely on.

The question has **two clauses**, and they are not the same question. Clause 1 is about
*selection*: can this be used to pick a method? Clause 2 is about *implication*: does getting
the anatomy right entail getting the biology right? A study can answer one and not the other.

---

## Clause 1 — can anatomy stand in for ground truth when choosing a method?

**Answered: partly, and the useful half is the negative half.**

| what was asked | evidence | verdict |
|---|---|---|
| Does anatomic concordance distinguish real methods from noise? | Worst fully-scored method 0.708 vs best control 0.400; both controls fail their own permutation nulls; every real method clears its null at p < 1e-4 | **Yes — it detects** |
| Does the ordering it produces predict accuracy against independent truth? | Spearman +0.081, n = 12, p = 0.80, CI [−0.629, +0.687]; on matched reference builds +0.314, p = 0.27 | **No — and underpowered** |

The ranking claim is reported **INCONCLUSIVE**, not refuted. Twelve comparable methods cannot
separate "no relationship" from "a relationship this test is too small to see", and the
confidence interval is wide enough to contain the other arm's value. Saying so is the
protocol's own instruction: an underpowered result is an answer.

**The honest one-line answer to clause 1: anatomy detects, it does not rank.**

---

## Clause 2 — does a method that gets the anatomy right also get the biology right?

**Answered: no, decisively, and this is the clause the study answers best.**

This clause was previously addressed only through the clause-1 correlation, which is the weak
route. It is answered directly by the study's own **pre-registered biological criterion**:

> *"PREDICTION: DNA methylation will show T cells > B cells."*
> — `prespecified/immune_failure_factors.md`, committed **2026-09-17 23:37:27**

The LGG methylation measurement was produced at **2026-09-18 01:12:23** — ninety-five minutes
later. The prediction genuinely precedes the data, and the document names its own falsifier:
*"methylation showing B ≥ T … the anomaly withdrawn."*

Measured, on that criterion (`results/anatomy_vs_biology.json`):

| | GBM | LGG |
|---|---|---|
| methylation places T above B in | 94.8% of 155 samples | 93.2% of 530 samples |
| methods reproducing T > B | **0 of 12** | **0 of 12** |
| method ranked **first** by anatomic concordance | `music`, ACS 0.9692 | `music`, ACS 0.9692 |
| …does it get T > B right? | **no** | **no** |
| …fraction of samples it places B above T | **100%** | **98.5%** |

**The method that best reproduces tumour anatomy places B cells above T cells in every single
glioblastoma sample it scores.** That is clause 2 answered categorically: reproducing known
anatomy does not entail getting the immune biology right, and the study has a case where it
entails the opposite.

The categorical claim needs no correlation and no power argument. A weaker supporting
observation — that the rank correlation between ACS and lymphoid failure runs in the
*opposite* direction to the hypothesis (+0.45 GBM, +0.48 LGG) — is **not significant at n = 12**
(p = 0.14, p = 0.12) and must be reported as directional only, if at all.

---

## Does the study prove its registered null?

> *"That 'the answer matches known biology' is not evidence that the answer is right."*

**Yes — by clause 2, not by clause 1.** This distinction matters for how the paper is written.

- Clause 1's correlation is underpowered, so it cannot carry the null on its own.
- Clause 2 carries it outright: methods that reproduce anatomy to 0.97 concordance place the
  lymphoid compartment backwards in ~100% of samples, against a pre-registered criterion on an
  orthogonal instrument.

A paper that rests the null on the correlation is resting it on the weakest evidence it has.
A paper that rests it on the lymphoid result rests it on the strongest.

---

## Two things the manuscript currently gets wrong

### 1. The title answers a different question from the one registered

Current draft title: *"Bulk RNA deconvolution cannot report the lymphoid compartment of a
glioma."* Registered title: *"Anatomic Concordance: Neuropathology as a Ground-Truth-Free
Benchmark for Cell-Type Deconvolution in Glioblastoma."*

The draft title is a claim about deconvolution's capability. The registered question is about
whether the Anatomy Test works. Both are supported by the data, but a reader comparing the
registration against the report sees a study that asked one thing and reported another — which
invites the reading that the registered question gave an underpowered answer and a stronger
secondary finding was promoted to the headline.

**It was not promoted post hoc** — the T > B criterion is itself pre-registered, with a
timestamp. But the title should *show* that, by naming both halves. A title that connects them
is both more honest and stronger, because it states the null the protocol set out to prove:

> **Anatomic concordance detects broken deconvolution but cannot choose a working one: the
> best-scoring method inverts the lymphoid compartment in both glioma cohorts**

### 2. GBM and LGG are presented as equivalent cohorts; they are not

The manuscript's Result 3 table sets GBM beside LGG as two confirmations. The
pre-specification is explicit that they play different roles:

> *"This is LGG, not GBM. The anomaly was measured in GBM. Whether it exists in LGG at all is
> checked first."*

So **GBM is discovery and LGG is the pre-registered replication.** Presenting them as two
equal cohorts overstates GBM and, more importantly, *throws away the strongest structural
feature of the design* — an anomaly observed, a prediction registered with a named falsifier,
and an independent confirmation in a different tumour type at 3.4× the sample size. That is a
better story than "we saw it twice", and it is the one that actually happened.

---

## What the paper may and may not claim

| claim | standing |
|---|---|
| Anatomic concordance separates real methods from negative controls | **May state plainly.** Margin +0.308 over the best control on fully-scored methods. |
| No method reproduces the pre-registered T > B ordering | **May state plainly.** 0 of 12 in both cohorts; truth at p = 1.1e-24 and 1.1e-76. |
| The anatomy-best method inverts the lymphoid compartment | **May state plainly.** `music`, ACS 0.9692, places B over T in 100% (GBM) and 98.5% (LGG). |
| Reproducing known biology is not evidence of correctness | **May state plainly**, sourced to the clause-2 result, not the correlation. |
| Anatomic concordance does not predict numerical accuracy | **Must be qualified as INCONCLUSIVE.** ρ = +0.081, p = 0.80, n = 12. Underpowered. |
| Higher anatomy scores go with worse lymphoid accuracy | **Must not be claimed.** Direction is opposite to the hypothesis but p = 0.14 / 0.12 at n = 12. Report as directional only. |
| Deconvolution cannot report the lymphoid compartment *in general* | **Must be bounded** to glioma, to this reference, and to the methods tested. |
