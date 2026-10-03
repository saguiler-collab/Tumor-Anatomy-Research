# DESeq2 `unmix` -- configuration fixed before any output

**Written 2026-10-02, before R/run_deseq2_unmix.R existed and before any estimate.** Post-
registration extension panel; never enters a registered statistic.

**Why this method.** DESeq2's `unmix` is one of the seven methods Nguyen et al. 2024 rank "in top
10 in all three scenarios", and one of the two of those seven this study had not run.

**The package and the call.** DESeq2 1.52.0 (Bioconductor, installed), `unmix(x, pure, alpha,
shift, power = 1, format = "matrix", quiet = TRUE)`. Love MI, Huber W, Anders S, *Genome Biol*
15:550 (2014).

**Inputs.** The same as every signature-based method in this study (FARDEEP, LinDeconSeq, EPIC):
- `pure` = the reference profile matrix on the run's gene space;
- `x` = the bulk matrix on that gene space.

Both are on a linear CPM scale. That is TPM-like, so `unmix`'s documented TPM argument, `shift`,
applies, not `alpha`, which is for negative-binomial counts.

**The one choice the documentation leaves to the analyst: `shift`.** The help page says "the shift
which approximately stabilizes the variance of log shifted TPMs. Can be assessed with
vsn::meanSdPlot". That is a visual judgement. Its automatable form, fixed here:
- compute log(x + s) for s in {0.1, 0.5, 1, 2, 5, 10, 20, 50};
- for each s, take each gene's mean and SD across bulk samples;
- choose the s whose |Spearman correlation of SD with mean| is smallest (the flattest mean-SD
  trend -- what the meanSdPlot inspection looks for);
- break ties toward the smaller s.

It uses the bulk matrix alone -- no reference, no estimate, no outcome -- and the chosen s is
recorded.

**Everything else is default.** `power = 1` (L1 loss, the default), `format = "matrix"`.

**Output and scale.** `unmix` returns mixing proportions of the pure profiles, each row summing to
1. Like FARDEEP's relative beta, that is an mRNA proportion on the signature's scale, so the
central cell-size conversion applies (not in `R_RETURNS_CELL_FRACTIONS`).

**Gate before any real data** (`tests/test_deseq2_unmix_gate.py`). Mixtures built from the
reference profiles with planted proportions must be recovered with per-type Pearson >= 0.90 and
max |error| <= 0.10. The same mixtures with the reference's gene labels shuffled must fail it.
