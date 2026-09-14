# What makes a deconvolution method clinically usable

This project's purpose is not to crown a winner. It is to find out whether a tumour's own
anatomy can pick a method, and — the question behind that — whether any of these methods can
carry a clinical decision. This file answers the second question from the measurements
already on disk.

Everything here is rendered from artefacts: `results/design_factors.json`,
`results/uncertainty_conformal.json`, `results/method_completeness.json`,
`results/benchmark/benchmark_per_type.csv`, `benchmark_per_niche.csv`.

---

## 1. "Why did MuSiC win?" is the wrong question, and it cannot be answered

MuSiC leads ACS at 1.000 against 0.9846 for four other methods. That margin is **one
constraint–tumour pair out of 57**. Six methods have confidence intervals containing MuSiC's
score. ACS takes only 66 distinct values on a weighted denominator of 65, and
`constraint_sensitivity.json` shows that dropping one constraint (C6) moves 13 of the 14
ranked methods.

Attributing a one-pair margin to a design feature would be storytelling. The protocol exists
to stop exactly that, and the honest answer is that **the top of this leaderboard is not
resolved.**

What *is* answerable: **which design property, isolated, changes accuracy?**

---

## 2. Near-controlled contrasts — the only causal-ish evidence available

Several pairs in this panel differ in exactly one property while sharing the reference, the
gene space, the harness and the cell-size handling.

| property isolated | ACS | MAE | tumour MAE | tumour \|bias\| | tumour r | worst-niche MAE |
|---|---|---|---|---|---|---|
| **cross-donor variance weighting** (MuSiC − NNLS) | **+0.0154** | **−0.0053** | **−0.0226** | +0.0347 | **+0.0843** | **−0.0043** |
| L1/L2 regularisation (ElasticNet − NNLS) | 0.0000 | −0.0010 | −0.0054 | +0.0190 | +0.0066 | +0.0009 |
| batch correction (CIBERSORTx B − SVR) | −0.0154 | −0.0001 | +0.0018 | −0.0132 | −0.0244 | +0.0011 |
| S-mode vs B-mode | 0.0000 | +0.0084 | +0.0029 | +0.0427 | −0.0295 | +0.0088 |
| donor-level hierarchy (Hier − Bayes) | 0.0000 | +0.0009 | +0.0001 | +0.0003 | −0.0067 | +0.0003 |
| **ensembling across references** (SCDC ENS − SCDC) | **0.0000** | **0.0000** | **0.0000** | **0.0000** | **0.0000** | **0.0000** |

The last row is a **built-in null control**: with one reference SCDC ENSEMBLE reduces exactly
to SCDC, so this contrast must be zero on every outcome. It is 0.0000 on all six. Any
non-zero value there would have meant the machinery was reading noise.

**One property buys anything: weighting genes by their consistency across donors.** It
improves five of six outcomes, most strongly the tumour correlation (+0.084) and tumour MAE
(−0.023). It is also the only property with a mechanism that fits this cohort: the reference
has **88 training donors**, so cross-donor variance is estimable. In a single-donor reference
MuSiC degenerates to NNLS, which this project already treats as an invariant.

**Three properties buy nothing here.** Regularisation is within noise. Batch correction makes
ACS *worse* and tumour correlation worse — CIBERSORTx's base solver is this project's SVR, so
that contrast is the correction layer alone. A donor-level hierarchy changes nothing at four
decimal places.

That is worth stating plainly because it runs against expectation: on a 10x-derived reference
against 2014 laser-capture bulk — the exact cross-platform gap batch correction exists for —
batch correction did not help.

---

## 3. The finding that matters clinically: every method fails where it counts

Mean absolute error hides this completely. Broken out by the *true* tumour content of the
mixture:

**mean |error| on tumour, by true tumour fraction**

| method | 0.0–0.1 | 0.1–0.3 | 0.3–0.5 | 0.5–0.7 | **0.7–1.0** |
|---|---|---|---|---|---|
| music | 0.032 | 0.089 | 0.109 | 0.187 | **0.351** |
| nnls | 0.049 | 0.120 | 0.138 | 0.198 | **0.383** |
| bayesprism | 0.051 | 0.123 | 0.135 | 0.203 | **0.362** |
| cibersortx | 0.057 | 0.113 | 0.117 | 0.184 | **0.465** |
| epic | 0.030 | 0.110 | 0.173 | 0.242 | **0.497** |
| bisque | 0.157 | 0.081 | 0.087 | 0.202 | **0.391** |
| scdc | 0.018 | 0.134 | 0.233 | 0.271 | **0.515** |
| dwls | 0.045 | 0.095 | 0.143 | 0.207 | **0.511** |

**signed bias at high purity (0.7–1.0): every method under-calls tumour**

| music | bayesprism | nnls | elastic_net | bisque | cibersortx | dwls | epic | scdc |
|---|---|---|---|---|---|---|---|---|
| −0.339 | −0.325 | −0.354 | −0.376 | −0.391 | −0.457 | −0.492 | −0.489 | −0.511 |

Error rises monotonically with tumour content for every method without exception, and at high
purity every method under-estimates tumour by **0.33 to 0.51**.

**This is where real tissue lives.** MuSiC's tumour estimates on the 122 real Ivy GAP samples
have median **0.68** and 75th percentile **0.86**. The benchmark's mixtures average 0.33, with
142 of 500 below 0.10. So the regime the benchmark is dominated by is the regime real tissue
is *not* in, and the aggregate MAE is therefore flattering.

### Where the missing tumour mass goes: it becomes T cells

Mean (predicted − true) per type, high-purity mixtures only:

| method | Tumor | Macrophage | **T_cell** | NK | Endothelial | Oligo | Astrocyte |
|---|---|---|---|---|---|---|---|
| music | −0.339 | +0.059 | **+0.157** | +0.010 | +0.025 | +0.032 | +0.039 |
| nnls | −0.354 | +0.047 | **+0.219** | +0.034 | +0.006 | +0.007 | +0.038 |
| bayesprism | −0.325 | +0.030 | **+0.190** | +0.037 | +0.017 | +0.007 | +0.033 |
| cibersortx | −0.457 | +0.062 | **+0.295** | +0.037 | +0.009 | +0.016 | +0.032 |
| scdc | −0.511 | +0.114 | **+0.306** | +0.025 | +0.008 | +0.031 | +0.027 |

True mean T-cell content in those mixtures is **0.044**. Predicted:

| music | bayesprism | nnls | cibersortx | scdc |
|---|---|---|---|---|
| 0.202 (**4.6×**) | 0.234 (5.3×) | 0.263 (6.0×) | 0.339 (7.7×) | 0.350 (**7.9×**) |

**Every method invents T cells in high-purity tumour, by a factor of 4.6 to 7.9.** The
Astrocyte sidecar — the population this project already flags as unidentifiable against Tumor
— absorbs only a third to a fifth as much.

For a clinical purpose this is the worst available failure mode. Immune infiltration is what
deconvolution of a tumour is most often used to estimate, and an immune-cold high-purity
glioblastoma would be reported as T-cell infiltrated by every method in this panel.

### And ACS is structurally blind to it

The frozen constraint file carries constraints on **Tumor, Macrophage_Microglia, Endothelial
and Oligodendrocyte**. `T_cell` carries none — it is *explicitly excluded*, for a defensible
reason recorded at registration (no anatomic fact about T-cell distribution across these five
structures is near-definitional enough to pre-register).

The consequence is not defensible to leave unsaid: **the anatomy test cannot see the single
largest error mode in the panel.** A method could invent T cells without limit and score 1.000
on ACS. That is a limit on what the method this study proposes can certify, and it belongs in
the paper's limitations rather than its discussion.

---

## 4. No method has usable uncertainty, and that alone disqualifies all of them

Of 15 methods, **one** reports per-sample uncertainty (`bayesian`), and its 95% intervals
achieve **9.25% coverage**. An interval a clinician would read as near-certain contains the
truth about one time in eleven. That is worse than no interval.

`scripts/uncertainty_conformal.py` fixes the availability problem without re-running anything:
the donor-held-out benchmark already has 500 mixtures with known truth, and split-conformal
prediction turns residuals into intervals with a coverage guarantee. All 17 methods reach
~95% measured coverage on held-out data.

The widths are the finding.

**95% interval half-width on the tumour compartment**

| method | ± | vs `control_random` |
|---|---|---|
| **bisque** | **0.489** | 0.67 |
| music | 0.557 | |
| bayesprism | 0.583 | |
| nnls | 0.630 | |
| elastic_net | 0.634 | |
| cibersortx | 0.634 | |
| epic | 0.642 | |
| scdc | 0.658 | |
| dwls | 0.675 | |

**Most methods cannot bound tumour purity to better than ±0.63 at 95% confidence, against
±0.67 for a Dirichlet random draw.** Restricting to the clinically plausible band
(0.30–0.80 true tumour) does not narrow them — it widens them slightly, because it removes the
easy near-zero mixtures.

Only **Bisque (±0.49)** and **MuSiC (±0.56)** are meaningfully better than random on tumour
purity, and neither is close to clinically actionable. A decision turning on 10 percentage
points of purity is not supportable by any method measured here.

---

## 5. The Bisque paradox: ACS rank and clinical reliability are different axes

Bisque ranks **10th of 14** on ACS (0.9231) and 9th on MAE (0.0741). On the two properties a
clinical use actually needs it is **first**:

| | Bisque | MuSiC | best? |
|---|---|---|---|
| tumour bias | **−0.0256** | −0.0757 | Bisque |
| tumour 95% half-width | **0.489** | 0.557 | Bisque |
| ACS | 0.9231 | 1.000 | MuSiC |
| MAE (primary types) | 0.0741 | 0.0508 | MuSiC |
| tumour MAE | 0.1550 | 0.1200 | MuSiC |

Bisque has *higher mean error* but *less bias and a lighter tail*. For a clinical assay the
tail is what matters: mean error describes the typical case, and a clinical failure is by
definition not the typical case.

Two caveats keep this from being a recommendation. Bisque runs in its documented **no-overlap
mode** here, because Ivy GAP and GBmap share no subjects — so this is not Bisque at full
strength. And Bisque is the panel's **only confirmed double cell-size correction** (D1), so
these numbers will move when that is fixed. Which direction is unknown, and guessing it would
be the same error this project has already made twice.

---

## 6. What a clinically usable method would have to demonstrate

Derived from the failures above, not from a wish list. None of the six is satisfied by any
method in this panel.

1. **Calibrated per-sample uncertainty.** Not a global error estimate — an interval for *this*
   patient's sample, with measured coverage. Conformal calibration makes this achievable for
   any method today; what is missing is that nobody ships it.
2. **Accuracy that does not collapse at high purity.** Every method here degrades
   monotonically with tumour content and under-calls by a third or more where real tissue
   sits. A clinical method must be validated *in the regime of use*, not on mixtures spanning
   the whole simplex.
3. **Bias small enough to state.** A systematic −0.34 on the dominant compartment is
   correctable if it is characterised; it is disqualifying if it is not.
4. **No invented populations.** A 4.6–7.9× over-call of T cells in immune-cold tissue would
   mislead exactly the decision such an assay is used for.
5. **Worst-case, not average, robustness.** Report the worst anatomic niche. MuSiC's spread
   across niches is 0.0106 and the Bayesian models' is 0.044 — a four-fold difference invisible
   in any mean.
6. **The genuine published software, on a recorded input.** Two of the panel's rows are
   reimplementations because the published packages exceeded their budgets, and their
   re-measurements used a gene space nobody had recorded (D10). A clinical claim cannot rest
   on software that was not the software.

---

## 7. What this says about the Anatomy Test itself

The honest summary is two-sided, and both sides belong in the paper.

**In favour.** ACS ranks methods approximately as ground truth does (rho = 0.75, CI
[0.33, 0.96]) using no ground truth and no outcomes. The one design property that measurably
improves accuracy — cross-donor variance weighting — is also the one that improves ACS. That
is the agreement the study was built to test, and it holds.

**Against.** ACS cannot resolve the top of its own leaderboard, its ordering below the top
rests largely on one constraint, and it is **structurally blind to the largest error mode in
the panel** because T_cell carries no constraint. A method selected by anatomy could be one
that invents T cells.

So anatomy is usable for *separating good methods from bad*, which is what the leaderboard's
controls establish. It is not usable for *certifying a method for clinical use*, and this
study should not be read as doing so. Whether a constraint set could be built that does cover
the failure modes is a real question for a second tissue — and it would need constraints on
the immune compartment, which is precisely where near-definitional anatomic facts are hardest
to pre-register.
