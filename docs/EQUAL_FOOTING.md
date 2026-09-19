# Equal footing: which methods were handicapped, and what it costs the ranking

**Written 2026-09-18 because the question was asked directly and the answer is uncomfortable.**

A leaderboard is only meaningful if every method got the same problem. On the two TCGA arms it did
not. **Five of fourteen methods ran with a defining feature disabled or failed outright**, and the
handicap is **not random with respect to the result** — it fell on exactly the methods whose
contribution is variance weighting.

## What the reference gave every method

The TCGA arms solve against the vendored frozen signature. Measured:

| property | frozen signature *(TCGA arms)* | h5ad-built reference *(Ivy GAP arm)* |
|---|---|---|
| profile | 17,850 × 8 | 16,758 × 8 |
| **sigma (cross-donor variance)** | **all zero** | **non-zero** |
| donor profiles | **None** | present |
| `has_cross_donor_variance` | **False** | True |

## The five compromised methods, and why

| method | state on TCGA | cause | fixable without new data? |
|---|---|---|---|
| **MuSiC** | **degenerate → NNLS** | sigma is all-zero, so its gene weighting is constant | **yes** — the h5ad reference has sigma |
| **EPIC** | uniform gene weights | `refProfiles.var` comes from sigma, which is zero (D13) | **yes** — same |
| **SCDC ENSEMBLE** | degenerate → SCDC | one reference supplied | **yes** — Neftel, Darmanis and gbmap_linear are all built |
| **CIBERSORTx S-mode** | **FAILED, absent** | needs the single cells the signature came from; none registered | **yes** — `sc_counts_gbmap.csv.gz` exists |
| **hierarchical Bayesian** | **FAILED, absent** | `manifest["patient_id"]` missing | **yes — FIXED**, TCGA patient = first three barcode fields |
| **Bisque** | degenerate, no-overlap mode | needs subjects assayed as **both** bulk and single cells | **NO — this needs external data** |

## What this does to the cross-validated ranking, stated plainly

The immune arm reported CIBERSORTx (1.03×) and SVR (0.87×) as best calibrated against both
molecules, and that finding is **confounded**:

> **MuSiC and EPIC — the two methods whose published contribution is variance weighting — were
> competing with that feature switched off. CIBERSORTx and SVR have no such feature to lose.**

So part of "CIBERSORTx and SVR win" is "CIBERSORTx and SVR were not handicapped." The ranking is
**not** safe to report as a method comparison until the reference carries sigma. The *measurements*
— 21.5% / 17.2% tumour recovery, the immune under-call, the withdrawn T-cell over-call — do **not**
depend on the ranking and stand as they are.

MuSiC returning byte-identical numbers to NNLS in both cohorts is the proof rather than an
inference: two different published algorithms cannot agree exactly unless one has collapsed into
the other.

## The fix, and its cost

Rebuild the TCGA reference from `gbmap_core.h5ad` instead of the frozen signature. That single
change restores sigma and donor profiles, which un-degenerates MuSiC, restores EPIC's variance
weighting, and lets S-mode run. Roughly **10 minutes to build plus ~3 hours to re-run both
cohorts**.

It is worth doing before the method ranking is published. It is **not** required for anything else.

---

## MEASURED 2026-09-18 — equal footing changes the ranking almost completely

The rebuild was run (`--reference h5ad`, `matrix="raw/X"`). It delivered what this document
predicted: **sigma all-zero: False; cross-donor variance: True**, 15,311 cells across 110 donors.
MuSiC, EPIC and CIBERSORTx S-mode became runnable in their intended modes, and
`bayesian_hierarchical` — which had been failing outright on a missing `patient_id` — joined the
panel. The GBM run scored **14 real methods against 12** on the frozen signature.

### Recovery of true tumour-content variation, GBM, on the 154 samples scored by BOTH runs

| method | frozen | h5ad (sigma) | comparable frozen | comparable h5ad |
|---|---|---|---|---|
| music | 17.5% | 60.2% | excluded | ✓ |
| cibersortx | 63.7% | 58.4% | ✓ | ✓ |
| svr | 59.1% | 55.0% | ✓ | ✓ |
| scdc | 20.5% | 53.3% | ✓ | ✓ |
| scdc_ensemble | 20.5% | 53.3% | excluded | excluded |
| bisque | -1.5% | 35.9% | excluded | excluded |
| bayesprism | 9.8% | 30.0% | ✓ | ✓ |
| nnls | 17.5% | 23.7% | ✓ | ✓ |
| bayesian_hierarchical | — | 22.7% | excluded | ✓ |
| bayesian | 43.1% | 21.3% | ✓ | ✓ |
| cibersortx_smode | — | 18.5% | excluded | ✓ |
| epic | 29.1% | 16.5% | excluded | ✓ |
| elastic_net | 20.1% | 15.2% | ✓ | ✓ |
| dwls | 43.5% | 14.7% | ✓ | ✓ |

### The ranking does not survive

> **Kendall tau between the two rankings = +0.214** on the 8 methods rankable under both.

That is near-random agreement. The frozen-reference ranking was **substantially an artefact of the
inputs those methods were denied**, not a measurement of the methods.

**The single largest movement is the one this document predicted.** MuSiC — whose entire published
contribution is cross-donor variance weighting — was excluded on the frozen signature as
arithmetically NNLS (17.5%). Given the variance it is designed to use, it recovers **60.2% and
ranks first.** DWLS falls furthest in the other direction, **43.5% → 14.7%, rank 3 → 12.**

**What survives:** CIBERSORTx and SVR stay near the top (63.7% → 58.4%, 59.1% → 55.0%, ranks 1→2
and 2→3). The immune arm's "CIBERSORTx and SVR are best calibrated" claim is **not overturned** —
but it is no longer the whole story, because the method that now beats them could not previously
be measured at all.

### The confound was checked, not assumed away

Registering the cells does two things at once: it supplies sigma, **and** it makes several R
packages runnable where the frozen path fell back to a Python reimplementation. A rank change
caused by swapping a reimplementation for a package would say nothing about equal footing.

Of the 8 methods ranked under both, **1 changed implementation** (`scdc`,
python-reimplementation → R:SCDC). Excluding it, **tau = +0.143** — lower than the
+0.214 overall. The rank change is therefore **not** an artefact of the implementation
swap; every large mover (dwls −9, elastic_net −5, bayesian −4, bayesprism +3) ran the *same*
implementation in both runs.

### Equal footing improves the RANKING and worsens the CALIBRATION — these are different questions

Measured on the same 154 GBM samples, mean true purity 0.752:

| | frozen signature | h5ad (sigma + raw counts) |
|---|---|---|
| median bias (estimate − purity) | **−0.0277** | **−0.4720** |
| methods under-calling tumour | 7 of 12 | **12 of 12** |
| Kendall tau agreement with the other ranking | — | **+0.214** |

**The rebuilt reference tracks variation better and sits further from the truth.** That is not a
contradiction, and it is not a reason to prefer one reference as simply "the right one":

- **Recovery** is `1 + slope` of (estimate − truth) regressed on truth. It measures whether a
  method *tracks* tumour content. It is **invariant to a constant offset** — verified directly:
  shifting every MuSiC estimate by −0.30 or +0.50 leaves recovery at 60.2% to one decimal, and
  `tests/test_recovery_is_bias_invariant` locks that in.
- **Bias** measures whether a method lands on the right *level*, and is exactly what an offset
  moves.

So the tau = +0.214 reshuffle is **not** the calibration shift in disguise. The two findings are
independent, and both belong in the paper:

> **Giving every method its intended inputs changes which method wins, and does not fix the
> systematic under-call of tumour content. The under-call gets worse.**

**Why the calibration moves is NOT established.** The plausible cause is that the frozen
signature is built from log-transformed values (D16) while the rebuilt one uses `raw/X` counts, so
the Tumor profile's magnitude relative to the other seven types differs, shifting every absolute
fraction. That is a hypothesis consistent with the numbers, not a measurement, and it is recorded
as such — the frozen signature's near-zero median bias may equally be a coincidence of the log
transform rather than a virtue. **Neither reference is endorsed as correctly calibrated here.**

The project's standing clinical finding — *no method in this panel is fit to report a patient's
tumour or immune content* — is **strengthened**, not weakened, by this. Under the reference that
gives every method what it was designed to consume, **all twelve under-call tumour, by a median
of 0.47 against a mean true purity of 0.75.**

### REPLICATION IN LGG — and it does NOT clearly hold

The GBM reshuffle above was reported before the LGG leg finished. It has now finished, and it
qualifies the claim rather than confirming it.

| | GBM | LGG |
|---|---|---|
| methods rankable under both references | 8 | 7 |
| Kendall tau, all rankable methods | **+0.214** | **+0.333** |
| methods whose *implementation* also changed | 1 (`scdc`) | 1 (`scdc`) |
| **tau excluding the switched method** | **+0.143** | **+0.733** |

**The confound check points opposite ways in the two cohorts.** In GBM, removing the one method
whose implementation changed *lowers* tau (+0.214 → +0.143) — the reshuffle is real and if anything
understated. In LGG, removing it *raises* tau sharply (+0.333 → **+0.733**) — meaning on the six
methods whose implementation was constant, **the two rankings largely agree.**

> **So LGG does not support the strong claim.** "Equal footing reshuffles the ranking" is
> demonstrated in GBM and **not replicated** in LGG once the implementation confound is removed.

**A defect in this project's own tooling is what surfaced it.** `equal_footing_ranking.py`
originally printed *"so the rank change is NOT an artefact of swapping a reimplementation for a
package"* **unconditionally** — a conclusion that survived either outcome, which is not a check at
all. It would have stamped that sentence onto the LGG result, which contradicts it. The script now
reads the direction and prints a caution when the switched method carried the effect.

**Both taus rest on six or seven methods**, where Kendall tau moves in steps of ~0.1–0.13 per
discordant pair. The honest reading is that the cohorts **disagree about whether the ranking is
stable**, not that one of them has measured a precise effect size. The script prints this warning
whenever fewer than eight methods are rankable.

**What survives the qualification, in both cohorts:**

- **MuSiC moves from excluded-as-degenerate to rank 1** given the variance it is designed to use
  (GBM 17.5% → 60.2%; LGG 20.0% → 52.6%). That is not a tau artefact — it is a method that could
  not previously be measured at all, and it is the single clearest consequence of equal footing.
- **The absolute under-call gets worse, not better** — GBM median bias −0.028 → −0.472; 7 of 12
  methods under-calling → 12 of 12.
- **`dwls` fell out of the LGG h5ad panel entirely** (1 finite estimate) and `epic` collapsed from
  36.5% to 3.4%, so the LGG comparison rests on fewer methods than GBM's and is the weaker arm on
  that ground too.

### Three honest caveats

1. **`bisque` and `scdc_ensemble` remain excluded under both references.** Their missing inputs
   are properties of the cohort and the pipeline, not of the reference build: TCGA has no
   subjects assayed as both bulk and single cells, and `DeconvolutionInput` is constructed with
   exactly one reference. ENSEMBLE returned numbers identical to SCDC to four decimals on the
   h5ad run, which is what a correctly-reported degeneracy looks like.
2. **`bayesprism` is the Python reimplementation in BOTH runs**, labelled as such and never
   reported as the published package. In the h5ad run its R attempt was terminated by me after I
   misread an R `parallel` socket-cluster master — which sits at 0% CPU by design while its
   workers compute — as a stalled process. The frozen run had already fallen back to the same
   reimplementation, so the two sides are consistent; but the h5ad R attempt was cut short rather
   than allowed to fail on its own, and that is recorded rather than smoothed over.
3. **This is GBM.** The LGG leg is the replication, and a ranking that reshuffles under a
   reference change is exactly the kind of result that needs one.

---

## The one genuine external need

> **Paired bulk RNA-seq and single-cell RNA-seq from the same subjects.**

Bisque's method *is* the assay transform learned from subjects measured both ways. Without
overlapping subjects it estimates that transform from marginal distributions instead — which is its
documented degraded mode, and it is why Bisque carries **no tumour information at all**
(r = −0.006) while over-calling immune content 2.27×.

**TCGA has no single-cell data, so this cannot be fixed with anything on disk.** It has been open
since the project started (`EXTERNAL_ACTIONS` item 11) and it is the only equal-footing gap that
external data would close.

**Until then Bisque must be reported as running in its degraded mode, never as "Bisque performed
poorly."** That distinction is this project's own invariant — *never report a degenerate method
under its own name* — and it applies to Bisque, MuSiC and SCDC ENSEMBLE alike.

## What is NOT an equal-footing problem

- **quanTIseq** brings its own TIL10 signature and receives the full gene space. That is declared
  (`OWN_SIGNATURE_METHODS`) and correct: ranking genes against a reference it does not use would be
  the wrong genes for it. It returns 0 finite estimates on both TCGA cohorts and is excluded from
  every statistic.
- **CDSeq** is reference-free by construction, so it cannot be handicapped by a reference.
