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
