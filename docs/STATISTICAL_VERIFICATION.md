# Statistical verification: every correlation, p-value and confidence interval, rechecked

Performed 2026-09-23, before any manuscript text was written. This is a different exercise
from `scripts/check_doc_numbers.py` and `scripts/check_doc_statistics.py`, which ask whether
a **document quotes an artefact correctly**. This asks whether the **artefact's number is
itself correct** — recomputed from the underlying data, and where the reported statistic
rests on an approximation, checked against an exact or simulation-based alternative.

Everything below was recomputed with `scipy` and `numpy` directly, from the source artefacts,
without calling the project's analysis code.

---

## 1. The rank correlations reproduce exactly

Each was rebuilt from the two artefacts that feed it — per-method ACS from
`acs_leaderboard.csv`, per-method accuracy from the yardstick JSON — restricted to comparable
methods with controls excluded, exactly as `agreement.test_agreement` does.

| statistic | n | recomputed | artefact | agrees |
|---|---|---|---|---|
| ACS vs synthetic-mixture accuracy | 14 | +0.637174 | +0.637174 | to 1e-9 |
| ACS vs DNA tumour purity | 12 | +0.080986 | +0.080986 | to 1e-9 |
| ACS vs DNA purity, excluding degenerate | 9 | −0.125524 | −0.125524 | to 1e-9 |
| ACS vs purity, true parity (both arms `raw/X`) | 14 | +0.3142 | +0.3142 | to 1e-4 |
| ACS vs purity, true parity, excluding degenerate | 12 | +0.2817 | +0.2817 | to 1e-4 |

**The inputs were verified too**, not just the correlation of them. Each method's
`spearman_vs_purity` was recomputed from the per-sample estimate table against the DNA purity
column: largest disagreement **0.000048** across 12 methods in GBM (154 samples) and
**0.000049** in LGG (510 samples). The per-method sample counts match the artefact in every
case.

**Direction was checked, not assumed.** `synthetic_mixtures` carries mean absolute error, so
the truth vector is negated; `absolute_purity` carries a correlation, so it is not. Feeding
the wrong flag inverts the sign — that is OPEN_DEFECTS **D21**. Recomputing with the flag
flipped reproduces the retracted value, which confirms the current one is the right way up.

---

## 2. The p-values are approximations, and the approximation was tested

`scipy.stats.spearmanr` returns a p-value from an asymptotic *t* approximation,
t = ρ·√((n−2)/(1−ρ²)) on n−2 degrees of freedom. At n between 9 and 14, with ties present,
that approximation is worth checking rather than trusting. Each was therefore compared with a
permutation p-value — exact by full enumeration where the factorial permits, Monte-Carlo
otherwise. Spearman was computed as Pearson on average ranks, so ties are handled correctly
and the permutation is exact in the rank domain.

| statistic | n | asymptotic p | permutation p | method | difference |
|---|---|---|---|---|---|
| ACS vs synthetic mixtures | 14 | 0.014251 | 0.016286 | 2,000,000 draws | 0.002035 |
| ACS vs DNA purity | 12 | 0.802431 | 0.801291 | 300,000 draws | 0.001141 |
| ACS vs purity, excl. degenerate | 9 | 0.747618 | **0.751201** | **exact, 362,880 permutations** | 0.003583 |

No difference exceeds 0.004. **The registered outcome survives either choice**: ρ = 0.6372
clears the 0.60 bar, and p is below 0.05 under both the asymptotic (0.0143) and the
permutation (0.0163) calculation. The two null results stay null under both.

The true-parity arm's p-values (0.274 at n = 14, 0.3751 at n = 12) were reproduced to four
decimals from ρ and n through the same *t* relation, and independently from that artefact's
own stored per-method points.

---

## 3. The confidence intervals reproduce bit-for-bit

`_bootstrap_rho` is a percentile bootstrap resampling **methods**, which is the right unit:
the sampling variability that matters is which methods happened to be in the comparison. It
was re-run outside the project code with the recorded seed.

| interval | recomputed | artefact | resamples skipped |
|---|---|---|---|
| synthetic mixtures | [+0.104038, +0.943029] | [+0.104038, +0.943029] | 0 of 5,000 |
| DNA purity | [−0.628586, +0.686574] | [−0.628586, +0.686574] | 0 of 5,000 |

Identical to every digit at `RANDOM_SEED = 0` with 5,000 draws. **Method order matters** and
was matched: the bootstrap indexes into `sorted(set(acs) & set(truth))`, and resampling a
differently-ordered array with the same seed gives a different interval — a first attempt
that used dictionary order produced [−0.631, +0.713] and did not match, which is how the
ordering dependence was confirmed rather than assumed.

Two properties of the procedure are stated rather than left to be discovered:

- **`ci_excludes_zero` tests only the lower bound** (`lo > 0`). It is one-sided by design,
  because the pre-registered criterion is directional — ρ ≥ 0.60 — and a negative correlation
  fails the threshold regardless. It is not a two-sided significance test and is not reported
  as one.
- **Resamples with no spread are skipped.** Zero were skipped in either arm here, so the
  intervals above condition on nothing. (A comment in `_bootstrap_rho` claimed the skipped
  count "is not hidden"; the function does not in fact return it. The count is reported here
  instead.)

---

## 4. A reporting error found and corrected: permutation p-values at their floor

`null_p = (hits + 1) / (draws + 1)` over 10,000 draws. With zero hits it equals
**1/10001 = 9.999e-05** and cannot go lower. **14 of 17 methods sit at exactly that value.**

Several documents printed it as `p = 0.0001`, which reads as a quantity measured to four
decimal places. It is not: it means *no permutation out of 10,000 reached the observed ACS*,
and the resolvable statement is `p < 1e-4`. The 20,000-draw constraint checks
(Albiach, Darmanis) have the same shape at 1/20000, printed as `0.0001` after rounding.

Corrected everywhere to `< 1e-4` and `< 5e-5`. `scripts/summarize_results.py:fmt_perm_p`
now formats these, handling both floor conventions in use, and
`tests/test_permutation_p_floor.py` fails if any document reintroduces a point value —
verified by planting one.

**No conclusion changes.** Every affected method was already far past significance; what
changes is a claim of precision the test cannot support.

---

## 5. The headline result's independent truth, recomputed

The claim that no method reproduces the methylation-measured lymphoid ordering rests on the
methylation side being right. It was recomputed from `methylation_celltypes{,_lgg}.csv`,
building T as CD4T + CD8T and renormalising within {T, NK, B}:

| cohort | n | mean T | mean B | samples with T > B | Wilcoxon T vs B |
|---|---|---|---|---|---|
| GBM | 154 | 0.5224 | 0.1098 | **95.5%** | p = 1.1 × 10⁻²⁴ |
| LGG | 530 | 0.4758 | 0.2034 | **93.2%** | p = 1.1 × 10⁻⁷⁶ |

The full T > NK > B ordering holds on the means in both cohorts. The result is unaffected by
the tie-handling convention: `wilcox`, `pratt` and `zsplit` all give p < 10⁻²⁴ in GBM and
identical p in LGG. The sample counts match `n_orderable_samples` in both artefacts.

So the truth side of the headline claim is not marginal — it is overwhelming, in both
cohorts, on an instrument that shares nothing with the deconvolution arm.

---

## What this exercise does not establish

It verifies that each reported number is the correct output of its stated procedure, and that
the procedures' approximations hold at these sample sizes. It cannot establish that the
procedure is the right one to have chosen, or that the input data mean what the source
describes them as meaning. Those are judgements, and they are argued in `docs/METHODS.md` and
bounded in `docs/METHOD_LIMITATIONS.md` rather than settled here.

The study's power ceiling also remains what it was: twelve comparable methods. No amount of
recomputation widens it, which is why the ACS-vs-accuracy result is reported as
**inconclusive** rather than as a refutation.
