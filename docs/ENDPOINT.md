# The endpoint: what this project found, and what the paper says

> **Every correlation, p-value and confidence interval in this study was independently recomputed and, where it rests on an approximation, checked against an exact or permutation alternative** — see [STATISTICAL_VERIFICATION.md](STATISTICAL_VERIFICATION.md).


> **Limits of this study's measurements** — quanTIseq is INCONCLUSIVE rather than bad (scored on 15 of 57 constraint pairs, p = 0.31 against its own null), two genuine R packages fell back to reimplementations on this hardware, and two more could not be built at all. All of it is in [METHOD_LIMITATIONS.md](METHOD_LIMITATIONS.md).


> # ⚠ NUMBERS UPDATED 2026-09-22 — the reference was rebuilt from genuine counts
>
> Every figure in this document was originally computed against a reference built from GBmap's
> **`X`**, which is `log1p(counts × size factor)` and which OPEN_DEFECTS **D16** identifies as
> indefensible when solved against linear bulk. `run_all.py --matrix raw/X` has now rebuilt the
> ACS leaderboard, the benchmark and the agreement test from the **genuine counts**.
>
> | | log `X` (superseded) | **`raw/X` (current)** |
> |---|---|---|
> | primary outcome, pseudobulk yardstick | +0.7501 | **+0.6372** — still clears the 0.60 bar |
> | ACS vs DNA purity, frozen signature | +0.0698 | **+0.0810** (n=12) |
> | ACS vs DNA purity, excl. degenerate | −0.1044 | **-0.1255** (n=9) |
> | **ACS vs purity, TRUE parity (both arms `raw/X`)** | −0.0156 | **+0.3142** (n=14), p = 0.274 |
> | best ACS (`music`) | 1.0000 | **0.9692** |
> | `control_shuffled_signature` | 0.1385 | **0.3692** |
>
> **The conclusion is unchanged and the numbers are now defensible.** Every external yardstick
> still fails the pre-registered 0.60 bar. The ordering is largely preserved (Kendall tau
> **+0.6719** between the two leaderboards, Spearman **+0.8203**, median |Δ| **0.0308**), no
> method changed implementation between the runs, and the control separation survives — worst
> real method **0.600** against best control **0.400**. That margin is the weakest form of the
> claim, because the method at the bottom of the range is `quantiseq`, which is scored on only
> **15 of 57** constraint-tumour pairs and does not clear its own permutation null (p = 0.31) —
> *not evaluable* rather than *bad*. Excluding it, the worst **fully-scored** method is
> `bisque` at **0.708**, so the true margin over the best control is **+0.308**, not +0.200.
>
> **One statement genuinely weakens, and it is the parity row.** Under true parity the
> correlation is **+0.3142**, not ≈ 0. It still fails the bar and is
> not significant (p = 0.274, CI [-0.294,
> +0.951] spans zero), so "anatomy does not rank" holds — but the
> honest phrasing is now *"a weak positive that does not reach significance and does not clear
> the registered bar"*, not *"indistinguishable from zero"*. Passages below that still say
> ≈ 0 are superseded by this table.
>
> Superseded numbers are left in place below rather than overwritten, so the change is visible.
> `results_archive/` holds both runs; `results/matrix_arm_comparison.json` holds the per-method
> comparison.


**Status 2026-09-15. Everything below is measured and has an artefact behind it.** This is the
synthesis document — the claim, the evidence for it, the evidence against the original
hypothesis, and what is deliberately left undone.

---

## The closing measurement: the registered outcome, computed against BOTH yardsticks

The registration asks whether the ACS ranking agrees with a ranking from real ground truth. On
the `raw/X` build it reports **rho = 0.6372**, above the pre-registered bar of 0.60
(0.7501 in the superseded log-`X` arm). Correction **C4** observed that both arms of that
comparison use GBmap. That was an argument. This is the measurement:

| yardstick | shares with the ACS arm | rho | 95% CI | meets the registered bar? |
|---|---|---|---|---|
| synthetic pseudobulk | **GBmap** — mixtures built from its cells | **+0.6372** (p = 0.0143) | [+0.1040, +0.9430] | **yes** |
| **ABSOLUTE purity, from DNA** | **nothing** | **+0.0810** (p = 0.80) | **[-0.6286, +0.6866]** | **no** |
| ... excluding 3 degenerate methods | | **-0.1255** (p = 0.75) | | no |

> **These are the `raw/X` values**, read from `results/anatomic/agreement_report.json`. The
> log-`X` arm reported +0.7501 [0.328, 0.957], +0.0698 [−0.510, +0.769] and −0.1044
> [−0.753, +0.743]. Every verdict in the right-hand column is the same under both builds — the
> pseudobulk arm clears the bar, the independent arm does not — which is why the rebuild
> changed no conclusion. `scripts/check_doc_statistics.py` pins each of these to its artefact
> field so a future rebuild cannot leave this table behind again.

**The registered criterion is definitively not met against an independent yardstick** — it
requires rho >= 0.60 *and* a CI excluding zero, and this fails both.

**But it is a failure to reproduce, not a refutation, and the difference matters.** With twelve
methods the interval [-0.6286, +0.6866] is wide enough to contain the
pseudobulk arm's 0.6372 itself. (On the superseded log-`X` arm the same sentence named
0.7501, which the current interval would *not* reach — the argument survives the rebuild but its
numbers do not carry over.) What twelve points cannot do is
distinguish "no relationship" from "a relationship this test is too small to see". Per this
project's own rule, an underpowered result is reported as **INCONCLUSIVE** — which is an answer.

What is *not* inconclusive: the pre-registered bar was cleared against a yardstick sharing the
reference and missed against one that does not.

**Descriptively, the ordering is close to inverted.** The two methods ACS ranks *last* have the
*highest* correlation with DNA purity:

| method | ACS | rho vs DNA purity |
|---|---|---|
| MuSiC | **0.9692** | 0.4271 |
| SVR | 0.9538 | 0.7184 |
| EPIC | 0.7692 | 0.3732 |
| Bisque | 0.7077 | **0.0284** |
| Bayesian | **0.7538** | **0.7421** |
| DWLS | **0.7385** *(last)* | **0.6865** |

**What differs between the arms, because this is not a like-for-like re-run.** ACS is scored on
Ivy GAP against the h5ad-built GBmap reference; the purity correlations are scored on TCGA-GBM
(154 samples with a DNA call) against the vendored frozen signature. Cohort *and* reference
build differ. What is constant is the set of named methods being ranked — the same form the
registered outcome takes, since its own two arms are also different cohorts.

Artefacts: `results/absolute_purity_yardstick.json`, `results/yardstick_agreement.json`.

---

## The one-sentence finding

> **Anatomic concordance is a valid detector of whether a deconvolution method is working, and
> an invalid ranker of which method is working best.**

The tumour's anatomy does carry real, independently verifiable information about cell-type
composition, and methods that recover it are doing something genuine. But the *ordering* it
produces among working methods is not stable enough, and not externally valid enough, to choose
a method with.

**That is a negative result on the registered hypothesis, and it is the more useful finding.**
"Use biological plausibility instead of ground truth" is an increasingly popular idea across
deconvolution, spatial transcriptomics and single-cell annotation. This project tested it
properly — pre-registered, with negative controls, against three independent external yardsticks
— and found where it breaks. A field-wide caution supported by measurement is worth more than
one more benchmark that agrees with the others.

---

## Part 1 — What works. The detector claim, and it is strong.

### The constraints are real biology, verified in tissue that this project did not deconvolve

The seven pre-registered anatomic constraints were checked directly against **measured**
single-cell composition, with no deconvolution method involved:

| source | what it is | result |
|---|---|---|
| **Albiach et al. 2023** | 135,482 annotated cells with anatomic zones | **3 of 4 testable constraints satisfied.** C1 (Tumor CT > LE) **13.34-fold**, p < 5e-5. C2 (Oligodendrocyte LE > CT) **3.53-fold**, p = 0.00075. C7 by ordering. C5 violated and marked WEAK. |
| **Darmanis et al. 2017** | 4 patients, FACS-gated neoplastic fraction, core vs periphery | **4 of 4 gates support C1**, two at p < 5e-5 with 3 of 3 patients. The only cross-patient constraint test in the project. |
| **Ivy GAP ISH** | 18,778 rows, 480 genes, 899 sub-blocks, 42 donors | Markers declared before values were read; C2 explicitly **not checkable** because the panel has no oligodendrocyte-lineage marker, and OLIG2 was refused on the grounds that in glioma it measures tumour. |

**The constraint file was frozen and hashed before any of this** (`2d1fb47c…5807a`) and was not
edited when a constraint failed. C5's violation is reported, not repaired.

### ACS separates real methods from controls decisively

| | ACS | beats permutation null? |
|---|---|---|
| 14 real methods | 0.7385 – 1.0000 | **all yes, p < 1e-4** |
| `control_random` | **0.3385** | no, p = 0.78 |
| `control_shuffled_signature` | **0.1385** | no, p = 0.999 |

This separation **holds under every reference atlas tested** — GBmap, GBmap-linear, Neftel,
Darmanis, and both GBmap assay subsets. It is the most robust result in the project.

### The registered primary outcome passed its pre-registered bar

**Spearman rho = 0.6372** between the ACS ranking and the accuracy ranking, CI
[+0.1040, +0.9430], p = 0.0143, on
14 methods. The bar, fixed in advance, was rho ≥ 0.60 with a CI excluding zero.
(The log-`X` arm reported 0.7501, CI [0.328, 0.957], p = 0.0020. Both clear the bar.)

**It passed. And Part 2 is about what it actually measures.**

---

## Part 2 — What fails. The ranker claim, refuted four independent ways.

### 2a. It failed against every external ground truth tested

**Three independent anchors, none of which this project controls:**

| anchor | what it is | result |
|---|---|---|
| **Avila Cobos 2020** (*Nat Commun* 11:5650) | benchmark against known composition | **Mann-Whitney p = 0.50.** No evidence their top-tier methods score higher on ACS. **DWLS is in their top tier and ranks 14th of 14 here.** |
| **Ajaib 2023** (*Neuro-Oncology* 25:1236) | **imaging mass cytometry**, 33 antibodies, protein, 10 GBM with matched bulk | **Pre-specified prediction FAILED**, p = 0.0090 in the wrong direction. ACS ranked GBmap-derived markers (0.6667) above the GBM-specific markers IMC ranks first (0.5128). |
| **Sturm 2019** (*Bioinformatics* 35:i436) | flow cytometry + immune mixtures | Descriptive only — 3 overlapping methods. EPIC ranks 2nd here, consistent with their recommendation. |

The IMC test is the important one because it was **pre-registered**: the prediction, the arms,
the preprocessing and the constraint subsets were committed before the run, and the outcome that
occurred was named in advance as the most damaging one.

### 2b. The ordering is substantially a property of the reference atlas, not the methods

Rebuilding the same atlas from its own counts to break a confound (D16):

| comparison | atlas | expression space | Spearman |
|---|---|---|---|
| GBmap_10x vs GBmap_SS2 | same | same | +0.8170 |
| GBmap_log vs GBmap_linear | same | different | **+0.9161** (shared genes) / **+0.6852** (per-arm) |
| GBmap_linear vs Neftel | different | same | **+0.5099** |
| GBmap_linear vs Darmanis | different | same | **+0.3655** |

Changing platform costs little. Changing expression space costs little **when both arms share a
gene space** and appreciably more when each reference nominates its own (0.9161 → 0.6852), which
is itself a symptom of §2d. **Changing the atlas leaves about half the ordering** — and the
atlas rows have not been rerun per-arm, so they are compared like-for-like only on the shared
space. And the registered primary outcome compares two rankings that
**both** use GBmap — the ACS arm deconvolves against it, the accuracy arm builds its mixtures
from its cells (C4). Part of rho = 0.6372 is two GBmap-based rankings agreeing about
GBmap.

**And it is not only the atlas — it is what the reference CARRIES.** Measured 2026-09-18
against DNA purity on TCGA-GBM, holding the atlas fixed and changing only whether the reference
supplies cross-donor variance and donor profiles:

> **Kendall tau between the frozen-signature ranking and the sigma-carrying ranking = +0.214**
> on the 8 methods rankable under both, scored on the same 154 samples.

Near-random. **MuSiC — whose entire published contribution is variance weighting — moves from
excluded-as-arithmetically-NNLS (17.5% recovery) to FIRST at 60.2%**, and DWLS falls from rank 3
to rank 12. The obvious confound was tested: only 1 of the 8 changed implementation
(reimplementation → R package), and excluding it tau is **+0.143**, lower still.
`docs/EQUAL_FOOTING.md`.

This is a stronger statement than §2b's. A ranking that reshuffles when the *atlas* changes might
be blamed on biology. A ranking that reshuffles when the **same atlas is merely supplied in the
form the methods were designed to consume** is measuring the harness, not the methods.

### 2f. The last objection to the ranker refutation is now closed

Every result above compared an ACS ranking scored on Ivy GAP **against the h5ad-built GBmap
reference** with purity correlations scored on TCGA **against the vendored frozen signature**.
That left one honest objection available: the two arms disagree because `music` names **R:MuSiC**
in the ACS arm and, with an all-zero sigma, **arithmetically NNLS** in the purity arm. The label
was constant; the algorithm was not. A rank correlation between two rankings whose rows mean
different things is not evidence about anything.

With the h5ad purity run in hand, the *algorithm-identity* half of that objection can be removed
rather than argued about — both arms score a reference that carries sigma, so `music` is R:MuSiC
on both sides:

> **CORRECTED 2026-09-20.** This passage originally claimed the two arms were scored on "the same
> reference build". They are not, and the mismatch was introduced by this project. `run_all.py`
> builds the ACS arm's reference with `--matrix` defaulting to **`X`** — GBmap's
> `log1p(counts × size factor)`, the matrix OPEN_DEFECTS **D16** identifies as indefensible
> against linear bulk. The purity h5ad arm was later patched to `raw/X`. So before that patch
> both arms shared the *same* (wrong) matrix; afterwards they differ, with the ACS side left on
> the log matrix.
>
> **What is genuinely closed:** the objection that `music` named a different algorithm in each
> arm. Both arms now carry cross-donor variance, and the correlation is ≈ 0 either way.
>
> **What is not:** full parity. Closing it requires re-running `run_all.py --matrix raw/X`, which
> would rebuild the ACS leaderboard and the primary outcome on the defensible matrix. That has
> not been done, and until it is, the tau and rho figures below compare a log-space ACS ranking
> against a counts-space purity ranking.
>
> **Why this does not rescue the ranker claim:** the correlation is ≈ 0 under *every* combination
> tested — both arms on the frozen signature (+0.070), both excluding degenerate methods
> (−0.104), and the mismatched pair below (−0.016, +0.007). A mismatch that flipped a null result
> to a positive one would have to be extraordinarily well-aimed. It is recorded as a real
> imperfection in the comparison, not as a reason to doubt the direction.



**SUPERSEDED TABLE — values below are from the log-`X` arm.** Kept as the record of what was published; the `raw/X` values and the corrected conclusion follow immediately after it.

| purity arm scored on | n methods | Spearman(ACS, purity recovery) | 95% CI | meets bar (≥0.60) |
|---|---|---|---|---|
| frozen signature (as published) | 12 | **+0.0698** | [-0.510, +0.769] | no |
| frozen, excluding degenerate | 9 | **-0.1044** | [-0.753, +0.743] | no |
| **h5ad — same build as the ACS arm** | **14** | **-0.0156** | [-0.723, +0.660] | **no** |
| **h5ad, excluding degenerate** | **12** | **+0.0071** | [-0.760, +0.708] | **no** |

> ### ⚠ SUPERSEDED — the four values in the table above are the log-`X` arm's
>
> They are kept because they are the record of what was published, not because they are
> current. On `raw/X` the same quantities are **+0.0810**
> (n = 12),
> **-0.1255** (n = 9),
> **+0.3142** (n = 14,
> p = 0.274) and
> **+0.2817** (n = 12,
> p = 0.3751). The paragraph below has been rewritten to the
> `raw/X` numbers. The pre-`raw/X` wording is preserved in git history — the commit
> *"Close the CJSJ version on raw/X; fix a sign error in the orthogonal yardstick"* is the
> last one to carry it. (It is **not** in `results_superseded/`; that directory holds
> superseded result artefacts, not superseded prose, and an earlier draft of this note said
> otherwise.)

**Reference parity does not rescue it — but it does not leave the correlation at zero either,
and that distinction is the honest one.** On matched references the correlation is
**+0.3142** across 14 methods
and **+0.2817** across the
12 that are comparable. Neither is significant
(p = 0.274 and p = 0.3751), neither comes near
the pre-registered 0.60, and both confidence intervals span zero — so the claim *anatomy does
not rank methods* stands. What does **not** stand is the older phrasing *"indistinguishable from
zero"*: a weak positive that fails to reach significance is a different statement, and the
earlier version of this paragraph overstated the result in the study's own favour.

So the refutation no longer rests on a comparison anyone can object to on those grounds. **ACS
ranks methods; DNA purity ranks methods; on the same cohort-independent footing the two rankings
show no relationship that survives its own confidence interval.**

*(The bar was met only against the GBmap-derived pseudobulk yardstick, rho = **0.6372** on
`raw/X` — 0.7501 in the superseded log-`X` arm — the arm that shares its atlas with ACS, which is
correction C4's point. Every yardstick that shares nothing with ACS fails, and that is the
pattern, not one bad run.)*

### 2g. The most promising replacement for ACS was tested, and it also fails

The residual measurement above invites an obvious hope. Per-sample model fit needs **no ground
truth** — which is the whole appeal of ACS — and a sample the model cannot fit is a sample whose
estimate should be distrusted. If per-sample R² predicted per-sample error, this project would
have found the ground-truth-free trust signal it set out to look for, just not the one it
registered.

**It does not.** Spearman(R², |estimate − purity|) per method; the usable direction is *negative*
(better fit, smaller error):

| | GBM | LGG |
|---|---|---|
| samples | 154 | 510 |
| control: Spearman(R², purity) | -0.1276 (p = 0.115) | -0.0182 (p = 0.681) |
| methods where better fit → smaller error (p<0.05) | **1 of 12** | **6 of 12** |
| methods significant in the **wrong** direction | 0 | **2** |
| verdict | **NOT a usable trust signal** | **NOT a usable trust signal** |

The control matters and it passes: R² is **not** a proxy for purity (p = 0.115 and 0.681), so the
test is not circular — a positive result would have been real.

But the result is not positive. It is **1 of 12 in GBM and 6 of 12 in LGG** — not a replication, a
coin flip that landed differently — and in LGG **2 methods run significantly backwards**, where a
*better* model fit predicts a *larger* error. Even where the correlation is significant its
magnitude is |rho| ≈ 0.1–0.2, which accounts for 1–4% of the variance in error. Nothing in that
table would justify telling a clinician which of their samples to distrust.

**So the model-fit residual bounds the problem without diagnosing any individual case.** It says
the mixing model is violated for most of the signal; it does not say which samples are the bad
ones. That distinction is worth stating precisely, because the first claim is the useful finding
and the second is the one a reader will assume follows from it.

**This is the fourth candidate explanation or diagnostic rejected by measurement**, after
macrophage spillover, high purity, and signature non-separability. The project's central negative
claim therefore stands in its strongest form: *on this panel, on two cohorts, against three
independent molecules, no ground-truth-free quantity tested here tells you when a deconvolution
estimate can be believed.*

### 2c. It cannot resolve the top of its own leaderboard, and it is blind to magnitude

| rank | method | ACS |
|---|---|---|
| 1 | MuSiC | 1.0000 |
| 2= | NNLS, SVR, elastic net, EPIC | **0.9846 — a four-way tie** |
| 6= | CIBERSORTx, CIBERSORTx S-mode, Bisque | **0.9692 — a three-way tie** |

**MuSiC's win is one satisfied pair out of 57.** The question "why did MuSiC win?" is
unanswerable at this resolution, and saying so is more honest than constructing a story around a
single pair.

### 2d. It is movable by choices the method never made

- **One constraint carries the ranking.** Drop C6 and the ordering changes: rho 0.678, and 13 of
  14 methods move (D8).
- **Marker-set size moves ACS by 0.31.** Same atlas, same algorithm, same bulk, same constraints,
  same tumours — only markers-per-type changed, and ACS went **0.3846 → 0.6923**, non-overlapping
  CIs, from indistinguishable-from-null to significant, **non-monotonically**. The whole
  leaderboard spans ~0.60 to 1.00. A nuisance parameter moves ACS across three-quarters of the
  range it is being used to rank fifteen methods across.
- **ACS is invisible to a change that moves tumour content by 8 points.** EPIC's published
  contribution is variance weighting; the harness never supplied it, so what ran was a
  constrained least squares with uniform weights. Restoring it changes the estimates
  substantially — mean absolute difference **0.0829** on Tumor, 0.0529 on Astrocyte, and
  **0.3539** at maximum on a single value — and changes EPIC's ACS by **exactly 0.0000**
  (0.9846 both ways, identical CIs). The mechanism is plain: **ACS is a rank statistic**, so any
  change preserving the ordering across structures is invisible to it by construction. That is
  fine for a detector and disqualifying for certifying accuracy. `ROAD_TO_PAPER.md` 0.4.
- **ACS prefers the reference that violates the mixing model.** The registration states the model
  is *"additive on the linear scale"* and that data is *"never log-transformed for
  deconvolution"*. The reference was in fact built from log values (**C7**). Rebuilding it from
  the atlas's own counts — correcting the defect — makes **every one of the 14 methods score
  lower**, by 0.02 to 0.26, mean about 0.12. That is not a gene-space artefact: it holds when
  each reference nominates its own markers at matched density, and quanTIseq, which ignores the
  supplied reference, is unchanged to exactly 0.0000. **Fixing a model violation costs anatomic
  concordance.** A selector that rewards a preprocessing error is not measuring what it is being
  asked to measure.

### 2e. The reference-free test — and it cuts BOTH ways, which is the useful part

Run 2026-09-15. Every other method in the panel is handed GBmap, so none of them can separate
"the methods recover the anatomy" from "GBmap recovers the anatomy". **CDSeq estimates cell types
de novo from the bulk alone**, which makes it the only instrument that can. It was blocked until
the Ivy GAP read counts were obtained (**C9**, a retraction).

Genes were chosen by **variance in the bulk**, not from the leaderboard's 657 markers — those
were selected on GBmap, and using them would reintroduce the dependence this arm exists to
remove. Every method is scored on the **same 41 constraint–tumour pairs**.

| | ACS | p |
|---|---|---|
| MuSiC, CIBERSORTx S-mode | 1.0000 | < 1e-4 |
| elastic net, EPIC, NNLS, SVR | 0.9756 | < 1e-4 |
| BayesPrism, SCDC | 0.9268 | < 1e-4 |
| Bisque | 0.7077 | < 1e-4 |
| **DWLS** | **0.7231** | 0.0002 |
| **CDSeq — atlas for labelling only** | **0.7561** | 0.0002 |
| **CDSeq — no atlas at all** | **0.6829** | 0.0009 |
| `control_random` | 0.5122 | **0.2233 — fails its null** |
| `control_shuffled_signature` | 0.1951 | 0.9977 |

**Reading one: the constraints are recoverable without any atlas.** A method that never saw
GBmap — not for estimation, not even for naming its own output — reaches **0.6829 at
p = 0.0009**, decisively above a random control that does not beat its own null. The anatomic
signal is in the bulk, not manufactured by the reference. This is the strongest support the
**detector** claim has.

**Reading two, and it is a genuine counterweight to this project's own thesis: the reference
earns its place.** CDSeq lands at the *bottom* of the real methods, tying DWLS exactly. Thirteen
of fourteen reference-using methods beat it. So the reference is not merely a source of
arbitrary ordering — it contributes real signal that bulk-only estimation does not recover.

**That matters for honesty about C4.** "The ordering is substantially a property of GBmap" is
established. "GBmap contributes nothing but its own bias" is **not**, and this measurement is
what rules it out. Both sentences belong in the paper.

**The caveat is severe and is the reason the two CDSeq rows differ.** CDSeq returns anonymous
cell types, and mapping them onto the roster changed the answer: correlating with the GBmap
profile called four of eight estimated types **Tumor**; marker enrichment called **none** of them
Tumor and three of them Astrocyte. The de novo types are not cleanly identifiable against this
eight-type roster, and a reference-free method that cannot name its own output is of limited
clinical use whatever it scores.

Artefacts: `results/cdseq_anatomic.json`, `results/cdseq_matched_comparison.json`.

---

## Part 3 — The clinical findings, which stand on their own

These do not depend on the ranking question and are directly reportable.

- **Every method under-calls tumour by 0.33–0.51 at high purity** — exactly where real tissue
  sits. Error grows monotonically with tumour content. **This is not an artefact of denying the
  methods their inputs.** Re-measured 2026-09-18 on a reference carrying cross-donor variance and
  donor profiles, **12 of 12 methods under-call, by a median of 0.47 against a mean true purity of
  0.75** — worse than the 7 of 12 on the frozen signature. Equal footing changes *which method
  wins*; it does not fix *the level*. `docs/EQUAL_FOOTING.md`.
- **The missing mass becomes T cells**, over-called **4.6–7.9×** — *in the pseudobulk
  benchmark only.* **WITHDRAWN 2026-09-18 as a claim about tissue: it is a simulator artefact.**
  On TCGA-GBM tissue, measured against methylation-derived leukocyte fraction, the same methods
  place T cells at **0.0002–0.0085**, two orders of magnitude below their pseudobulk estimates,
  and the immune compartment as a whole is **UNDER**-called (median 0.661× after the
  pre-specified mRNA correction). `docs/IMMUNE_ARM.md`.
- **Methods cannot report the lymphoid compartment of a glioma at all.** Pre-registered with
  an explicit falsifier, then tested against DNA methylation in **both** cohorts (EpiDISH RPC;
  255 of 333 reference CpGs in LGG, 258 in GBM — the rest are all-NA). Methylation puts T above B
  in **94.8%** of GBM and **93.2%** of LGG samples; the true ratio is **4.76:1** in GBM and
  **2.34:1** in LGG. Two distinct failures follow, and conflating them overstates how tidy the
  result is: **four of twelve methods return *exactly zero* T, B and NK** in most samples
  (`music`/`nnls` in 86.9% of LGG and 98.2% of GBM samples), and of the eight that do return
  lymphocytes, **8 of 8 in LGG and 6 of 8 in GBM put B above T** per-sample. **Zero of twelve
  reproduce the true ordering in either cohort.** Its cause is not identified — the natural
  explanation, that the signature cannot separate T from B, was tested and refuted (NNLS recovers
  the planted ratio exactly on mixtures built from the reference, survives 100% noise and the
  deletion of a whole cell type; condition number 5.3). `docs/IMMUNE_ARM.md` §6.
- **The mixing model itself explains a minority of the signal, and that is measurable without
  ground truth.** Fitting the best non-negative combination of reference profiles to each sample
  leaves a median **63.8%** (GBM) and **76.4%** (LGG) of marker-space variance
  unexplained, with **83%** and **89%** of samples below R² = 0.5. Bracketed by controls:
  shuffling which gene belongs to which cell type drops R² to **-0.0077**, so the reference
  does carry real structure — but eight dimensions *could* explain **98%**, and the reference
  reaches only **39%** (GBM) and **25%** (LGG) of that. It does not identify which
  assumption fails, but it bounds every explanation that assumes the model holds.
  `results/model_fit_residual.json`.
- **ACS is structurally blind to it.** `T_cell` carries no constraint, by pre-registered design,
  because no anatomic fact about T-cell distribution was defensible enough to register.
- **1 of 15 methods reports uncertainty at all**, and its intervals cover **9.25%** against a
  nominal 95%.
- **Split-conformal intervals on tumour content span ±0.49 (Bisque) to ±0.67** — the latter is
  the width of the whole simplex, i.e. no information.
- **Only cross-donor variance weighting buys accuracy** (5 of 6 outcomes). Regularisation, batch
  correction and donor hierarchy buy nothing on this cohort.

**Clinical bottom line:** no method in this panel is fit to report a patient's immune content,
and the evaluation framework most likely to be used to choose between them cannot see the error.

---

## Part 4 — The methodological contribution

Seventeen defects were found, measured and recorded, several of which invalidate work that
looked finished. The ones a reader should see:

| | what it was | why it matters beyond this project |
|---|---|---|
| **D12** | the central mRNA→cell conversion was **the identity**, applied nowhere | proven twice, the second time without reading any code: applying the transform and skipping it gave bit-identical output |
| **D16** | the reference was built from **log-transformed** values treated as linear | the accuracy arm is internally consistent and survives; the **anatomic arm is not**. One defect explains why the benchmark looked healthy while tissue accuracy was poor |
| **D14** | the ordering does not survive a change of atlas | the first attribution ("it is the atlas") was confounded with expression space; corrected after measurement, magnitude withdrawn |
| **D17** | variant builds silently overwrote the primary reference's provenance record | a figure script read that exact path |
| **C9** | claimed Ivy GAP counts were not public | **retracted the same day**: GEO says "raw data not provided", the Allen portal serves per-sample RSEM counts. Wrong because one authoritative source was mistaken for an exhaustive search |

**Three claims in this project were published and then withdrawn after measurement** (D1 three
times, D14's magnitude, C9). Each retraction is kept in place rather than edited away. That
record is part of the result.

---

## What the paper claims, precisely

1. Anatomic constraints in glioblastoma are **real and independently verifiable** — confirmed in
   two external single-cell datasets that were never deconvolved.
2. Deconvolution methods **satisfy them far above chance**, and negative controls do not. ACS is
   a working **detector**. The constraints are recoverable **without any atlas at all**: CDSeq,
   estimating cell types de novo from the bulk and never touching GBmap even to name its own
   output, reaches **0.6829 at p = 0.0009** against a random control that fails its own null.
3. **ACS is not a valid method selector.** It fails against published ground-truth benchmarks
   (p = 0.50), fails a pre-registered test against protein ground truth (p = 0.0090, wrong
   direction), cannot resolve the top of its own leaderboard (a four-way tie), rests largely on
   one constraint, and moves by 0.31 on marker-set size alone.
4. The apparent success of the registered outcome (rho = 0.6372) is **substantially a shared
   dependence on one reference atlas**, and every analysis that breaks that dependence weakens
   it. **But the atlas is not merely a source of bias.** Reference-free estimation lands at the
   bottom of the real methods, tying DWLS, with 13 of 14 reference-using methods above it — so
   the reference contributes real signal that bulk-only estimation does not recover. "The
   ordering is substantially a property of GBmap" is supported; "GBmap contributes nothing but
   its own bias" is refuted.
5. Separately and robustly: **current deconvolution methods are not clinically usable for immune
   content in GBM**, and anatomy-based evaluation is structurally blind to that failure.

**The honest headline is not "anatomy works" and not "anatomy fails". It is: anatomy tells you
whether a method is working, not which one is best — and the difference between those two
questions is where the field has been assuming rather than measuring.**

**Why the negative result is credible, and this belongs in the paper rather than being left for
a reader to notice.** The hypothesis was registered with its constraint file frozen and hashed
before any deconvolution output was seen; the decision rule was fixed in code; the adverse
outcome of the external test was named in advance as adverse; and when the test failed, the
constraint file was not edited — nor was it edited when one of its own constraints failed an
external check. Three claims made during this work were published and then **withdrawn after
measurement**, and the retractions are kept in place. A result produced under those conditions
is worth more than a positive one produced without them.

---

## What is deliberately NOT done, and why that is fine

- **The full linear re-run (B1).** D16 means the anatomic arm should ideally be recomputed on a
  linear reference. It is declared as a limitation with its direction measured (+0.9161 — the
  expression space barely moves the ordering), so the paper does not depend on it.
- **CDSeq on the anatomic arm.** Unblocked and feasibility-measured (~30 min), not run. It would
  strengthen Part 2b; it is not required for any claim made.
- **A second tissue.** Nine evaluable tumours is the statistical ceiling and is stated as one.
- **Retuning ACS so it performs better.** Forbidden by the project's own invariants, and the
  reason the negative result is credible.
