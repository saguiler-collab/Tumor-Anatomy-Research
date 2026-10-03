#!/usr/bin/env Rscript
# run_linseed.R -- linseed 0.99.3 (genuine, vendored from github.com/ctlab/LinSeed @ 2435a8ce,
# installed from vendor/linseed; see vendor/linseed/PROVENANCE.md).
#
# Zaitsev K, Bambouskova M, Swain A, Artyomov MN. Nat Commun 2019;10:2209;
# doi:10.1038/s41467-019-09990-5. Reference-free complete deconvolution: genes whose expression
# is mutually linear across samples sit near the corners of a simplex whose vertices are the cell
# types; corners give signatures, and proportions follow.
#
# POST-REGISTRATION, EXPLORATORY: the reference-free arm's second instrument, beside CDSeq.
# Configuration fixed before any output: prespecified/linseed_config.md -- the package README's
# tutorial step for step, with k fixed a priori (the tutorial's visual SVD choice cannot be
# pre-declared).
#
# args.json: bulk (genes x samples, linear), k, top_genes, sig_iters, pval, seed, out_prop, out_sig.
# Writes proportions (samples x k, as returned) and signatures (genes x k) on the gene-expression
# scale: Linseed's row-normalised W times each gene's total in its input (the addendum to the
# configuration file).

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(linseed))

args <- read_args()
set.seed(args$seed)
bulk <- as.matrix(load_bulk(args$bulk))

lo <- LinseedObject$new(bulk, topGenes = args$top_genes)
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(args$sig_iters)
lo$filterDatasetByPval(args$pval)
lo$setCellTypeNumber(args$k)
lo$project("filtered")
lo$smartSearchCorners(dataset = "filtered", error = "norm")
lo$deconvolveByEndpoints()

H <- lo$proportions                      # k x samples
W <- lo$signatures                       # filtered genes x k, row-normalised scale
if (!identical(colnames(H), colnames(bulk))) {
  stop("Linseed's sample names/order differ from the input's; refusing to map them by position")
}
totals <- rowSums(lo$exp$full$raw[rownames(W), , drop = FALSE])
sig <- W * totals                        # back to the gene-expression scale

write.csv(t(H), args$out_prop)
write.csv(sig, args$out_sig)
cat(sprintf("linseed complete: %d of %d genes passed the significance filter; k = %d; %d samples\n",
            nrow(W), nrow(lo$exp$full$raw), args$k, ncol(H)))
