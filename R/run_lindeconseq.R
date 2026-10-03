#!/usr/bin/env Rscript
# run_lindeconseq.R -- the genuine LinDeconSeq package (GitHub lihuamei/LinDeconSeq @ 20f1aec75012).
#
# LinDeconSeq (Li H, Sharma A, Ming W, et al. BMC Genomics 2020;21:652;
# doi:10.1186/s12864-020-06888-1, verified against Crossref 2026-10-01) solves a weighted
# robust regression of each bulk sample on a signature matrix. Top tier of reference-based
# methods in Nguyen et al. 2024.
#
# POST-REGISTRATION EXTENSION PANEL (added 2026-10-01). Not in the registered panel; reported
# under its own name in its own panel, never in a registered statistic.
#
# The package has two stages: findMarkers() derives a signature from PURE bulk samples, and
# deconSeq() solves bulk against any signature ("signature marker genes must be known in
# advance", README). There are no pure bulk samples of glioma cell types, so stage one does not
# apply; deconSeq() is given this project's reference profile on the shared gene space, the same
# arrangement as EPIC and FARDEEP. All deconSeq() arguments are package defaults (weight = TRUE,
# intercept = TRUE, scale = FALSE, QN = FALSE). Planted-truth gate 2026-10-01: max |est - truth|
# 0.023, per-type r >= 0.998.
#
# Returns proportions on the signature's scale (mRNA shares); the central cell-size conversion
# applies, as for every signature-based method.

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(LinDeconSeq))

args <- read_args()
set.seed(args$seed)
bulk_mat <- load_bulk(args$bulk)
sig_mat  <- load_signature(args$signature)
shared <- intersect(rownames(sig_mat), rownames(bulk_mat))
if (length(shared) < 20) stop(sprintf("only %d genes shared between signature and bulk", length(shared)))

props <- as.matrix(deconSeq(bulk_mat[shared, , drop = FALSE], sig_mat[shared, , drop = FALSE],
                            verbose = FALSE))
bad <- rowSums(!is.finite(props)) > 0
if (any(bad)) cat(sprintf("LinDeconSeq: %d of %d samples returned no finite composition\n",
                          sum(bad), nrow(props)))
# SAMPLE NAMES. deconSeq returns rows in the bulk's column order, but its data.frame step mangles
# names that look numeric (Ivy GAP ids such as "300173630" come back as "X300173630"), so indexing
# by name failed on real data ("subscript out of bounds", 2026-10-01). Names are restored by
# POSITION only after proving the order: exact match, or exact match once R's "X" prefix is
# stripped. Anything else stops -- a reordered result must never be relabelled silently.
rn <- rownames(props); cn <- colnames(bulk_mat)
if (nrow(props) != length(cn)) stop(sprintf("deconSeq returned %d rows for %d samples", nrow(props), length(cn)))
if (!identical(rn, cn) && !identical(sub("^X", "", rn), cn) && !identical(make.names(cn), rn)) {
  stop("deconSeq row names do not correspond to the bulk's samples in order; refusing to relabel")
}
rownames(props) <- cn
write_proportions(props, args$cell_types, args$out)
cat("LinDeconSeq complete\n")
