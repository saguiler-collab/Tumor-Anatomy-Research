#!/usr/bin/env Rscript
# run_dwls.R — the genuine DWLS package.
#
# DWLS builds its own signature by differential expression (buildSignatureMatrixMAST),
# then solves with dampened weighted least squares one sample at a time. The
# dampening constant is re-selected per sample by findDampeningConstant, which is
# what lets a rare cell type register at all.

# Locate common.R relative to this script. `sys.frame()$ofile` is only set under
# source(); under Rscript the path arrives in commandArgs() as --file=.
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(DWLS))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sc_eset  <- load_sc_eset(args$sc_counts, args$sc_meta)

sc_counts <- exprs(sc_eset)
labels    <- as.character(pData(sc_eset)$cell_type)

# CACHE THE SIGNATURE ACROSS CALLS.
#
# buildSignatureMatrixMAST is the dominant cost in this driver — minutes, and it does
# not depend on the bulk at all, only on the cells and their labels. The pipeline calls
# DWLS three times per run (the benchmark, then each ACS cohort) with the SAME
# reference, so a per-call tempdir recomputes an identical signature three times.
#
# Keyed on the export path and the cell count, beside the export itself, so a different
# reference or a different donor split gets its own signature rather than silently
# reusing one built from other cells.
sig_key <- substr(tools::md5sum(args$sc_counts)[[1]], 1, 16)
sig_dir <- file.path(dirname(args$sc_counts), paste0("dwls_signature_", sig_key))
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)
if (length(list.files(sig_dir))) {
  cat(sprintf("DWLS: reusing cached signature build in %s\n", basename(sig_dir)))
}

# DOUBLE FILTERING, AND WHY THE THRESHOLDS ARE RELAXED HERE
# ---------------------------------------------------------
# The published DWLS workflow selects its signature genes from the whole
# transcriptome. This pipeline hands every method the SAME pre-filtered gene space
# (config.USE_SIGNATURE_GENE_SUBSET — top-N per cell type plus markers), because equal
# footing requires it. DWLS's MAST step then selects again, from an input that has
# already been reduced to the informative genes.
#
# At the published defaults that second pass left 31 genes for 8 cell types — about four
# each — and this driver stopped rather than deconvolve on them. The failure looked like
# DWLS being unable to run, when it is an artefact of asking it to re-do a selection
# that was already made.
#
# So the cutoffs are relaxed progressively until enough genes survive, and WHICH
# thresholds were used is printed into the run log. This is a documented deviation from
# the published defaults, made because of the pre-filtering this pipeline imposes, not
# because DWLS scored badly — its score is not known at this point.
# ONE BUILD, PERMISSIVE CUTOFFS, AND WHY THAT IS THE RIGHT SHAPE HERE
# -------------------------------------------------------------------
# The published DWLS workflow selects signature genes from the whole transcriptome. This
# pipeline hands every method the SAME pre-filtered gene space, because equal footing
# requires it — so by the time DWLS sees the data, the "find the informative genes" step
# has already been done, by `select_signature_genes`, identically for every method.
#
# Running MAST's differential-expression filter again on top of that is not a second
# opinion, it is a second filter. At the published cutoffs it left 31 genes for 8 cell
# types and this driver refused to deconvolve on them.
#
# An earlier attempt relaxed the cutoffs progressively, rebuilding the signature at each
# setting. buildSignatureMatrixMAST fits a hurdle model per gene across every cell, so
# each rebuild costs minutes and the loop ran for over an hour without finishing.
#
# So: ONE build, with cutoffs permissive enough to keep the gene space it was given.
# What remains of DWLS is its dampened weighted least squares solver — which is the part
# the method is named for, and the part this benchmark is comparing. Recorded as a
# deviation in docs/METHODS.md.
MIN_SIGNATURE_GENES <- 50
MAX_SIGNATURE_CELLS <- 3000

if (ncol(sc_counts) > MAX_SIGNATURE_CELLS) {
  set.seed(args$seed)
  keep_frac <- MAX_SIGNATURE_CELLS / ncol(sc_counts)
  idx <- unlist(lapply(split(seq_along(labels), labels), function(ix) {
    k <- max(1L, min(length(ix), as.integer(round(length(ix) * keep_frac))))
    if (k < length(ix)) sample(ix, k) else ix
  }), use.names = FALSE)
  idx <- sort(idx)
  cat(sprintf("DWLS: subsampling %d of %d cells for the signature build (within cell type)\n",
              length(idx), ncol(sc_counts)))
  sc_counts <- sc_counts[, idx, drop = FALSE]
  labels    <- labels[idx]
}

signature <- buildSignatureMatrixMAST(
  scdata = sc_counts, id = labels, path = sig_dir,
  diff.cutoff = 0, pval.cutoff = 1
)
shared <- intersect(rownames(signature), rownames(bulk_mat))
cat(sprintf("DWLS signature: %d genes, %d shared with the bulk (of %d supplied)\n",
            nrow(signature), length(shared), nrow(bulk_mat)))

if (length(shared) < MIN_SIGNATURE_GENES) {
  stop(sprintf(paste0("DWLS signature shares only %d genes with the bulk even with its ",
                      "differential-expression filter fully open. The %d-gene input is ",
                      "too small for DWLS to select from."),
               length(shared), nrow(bulk_mat)))
}
signature <- signature[shared, , drop = FALSE]

props <- t(sapply(colnames(bulk_mat), function(s) {
  trimmed <- trimData(signature, bulk_mat[shared, s])
  solveDampenedWLS(trimmed$sig, trimmed$bulk)
}))
rownames(props) <- colnames(bulk_mat)

write_proportions(props, args$cell_types, args$out)
cat("DWLS complete\n")
