#!/usr/bin/env Rscript
# run_fardeep.R — the genuine FARDEEP package (CRAN 1.0.1).
#
# FARDEEP (Hao, Yan, Heath, Lei & Xie, PLoS Comput Biol 2019;15:e1006976; doi:10.1371/journal.pcbi.1006976) is least-trimmed-squares regression
# with an adaptive outlier threshold: genes whose residuals are outliers are dropped from the
# fit, per sample, before the coefficients are estimated. It is in the top tier of reference-
# based methods in Nguyen et al. 2024, the 53-method review this project cites.
#
# POST-REGISTRATION EXTENSION PANEL (added 2026-09-30). FARDEEP is not one of the methods
# named in Anatomy_Test.md and was not part of the registered panel. It is run and reported
# under its own name in a separate, labelled extension panel and never enters the registered
# leaderboard or any registered statistic.
#
# It consumes a signature matrix, not cells, and is given this project's reference profile on
# the same gene space as every other method -- the EPIC arrangement.
#
# Every tuning argument is the package default (alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1,
# nn = TRUE, intercept = TRUE, lognorm = TRUE, QN = FALSE). Nothing is tuned against ACS.
# The ONE departure is permn = 0: the permutation loop computes only the goodness-of-fit
# p-values, which nothing here reads, and the estimates do not depend on it (FARDEEP source,
# lines 54-75; verified: max |diff| 0 between permn = 0 and permn = 25).
#
# FARDEEP returns `relative.beta`, abs.beta normalised to sum to one per sample -- an mRNA
# proportion on the signature's scale, so the central cell-size conversion applies to it as
# it does to every other signature-based method (it is NOT in R_RETURNS_CELL_FRACTIONS).

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(FARDEEP))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sig_mat  <- load_signature(args$signature)

shared <- intersect(rownames(sig_mat), rownames(bulk_mat))
if (length(shared) < 20) {
  stop(sprintf("only %d genes shared between the FARDEEP reference and the bulk", length(shared)))
}
sig_mat  <- sig_mat[shared, , drop = FALSE]
bulk_mat <- bulk_mat[shared, , drop = FALSE]

# PARALLEL OVER BLOCKS OF SAMPLES. FARDEEP fits each sample independently (its own tuningBIC
# and trimming; QN = FALSE, so nothing is shared across samples), so splitting the samples
# cannot change an estimate -- verified 2026-09-30: max |joint - blocked parallel| = 0. Blocks,
# not single samples, because fardeep() drops a dimension on a one-column Y and errors
# ("Y[YinX, ]: incorrect number of dimensions"). Workers from IVYGAP_FARDEEP_CORES, default 1.
.cores <- max(1L, as.integer(Sys.getenv("IVYGAP_FARDEEP_CORES", "1")))
.n <- ncol(bulk_mat)
.nb <- max(1L, min(.cores, .n %/% 2L))
# One block is the joint fit. cut() rejects a single interval ("invalid number of intervals"), which
# made the default (1 worker) fail outright until 2026-10-02; every earlier run had set 2 workers.
.blocks <- if (.nb == 1L) list(seq_len(.n)) else split(seq_len(.n), cut(seq_len(.n), .nb, labels = FALSE))
.fits <- parallel::mclapply(.blocks, function(j)
  fardeep(X = sig_mat, Y = bulk_mat[, j, drop = FALSE], permn = 0), mc.cores = .nb)
.err <- vapply(.fits, function(f) inherits(f, "try-error") || is.null(f$relative.beta), logical(1))
if (any(.err)) stop(sprintf("FARDEEP failed in %d of %d sample blocks", sum(.err), length(.err)))
props <- do.call(rbind, lapply(.fits, `[[`, "relative.beta"))[colnames(bulk_mat), , drop = FALSE]
est <- list(k.value = unlist(lapply(.fits, `[[`, "k.value")))
cat(sprintf("FARDEEP: %d samples in %d block(s) on %d worker(s)\n", .n, .nb, .nb))

# A sample whose every coefficient is zero divides 0/0 into NaN. Report it rather than let a
# downstream simplex projection turn it into a plausible-looking uniform composition.
bad <- rowSums(!is.finite(props)) > 0
if (any(bad)) {
  cat(sprintf("FARDEEP: %d of %d samples returned no finite composition (all-zero fit)\n",
              sum(bad), nrow(props)))
}
cat(sprintf("FARDEEP: outlier-trimming k per sample, median %.1f\n", stats::median(est$k.value)))

write_proportions(props, args$cell_types, args$out)
cat("FARDEEP complete\n")
