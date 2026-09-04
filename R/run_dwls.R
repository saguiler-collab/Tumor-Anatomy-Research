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

# CELLS PER TYPE, NOT A SHARE OF THE WHOLE.
#
# buildSignatureMatrixMAST fits a hurdle model per gene across every cell, and its cost
# tracks the cell count: 758 s on 11,739 cells, and ~30 minutes again whenever the gene
# set changes and the cache misses.
#
# An earlier version took a proportional subsample of the whole reference. That is the
# wrong shape twice over. A signature is a per-cell-type summary, so the type
# PROPORTIONS are irrelevant to it — what matters is how many cells estimate each type.
# Sampling proportionally keeps the abundant types abundant and leaves the rare ones
# thin, which is exactly backwards: Astrocyte, at ~1% of cells, would have been
# estimated from a handful.
#
# Capping PER TYPE spends the budget where it buys accuracy, keeps every cell of a type
# that has fewer than the cap, and costs less.
MAX_CELLS_PER_TYPE <- 250

if (ncol(sc_counts) > MAX_CELLS_PER_TYPE * length(unique(labels))) {
  set.seed(args$seed)
  idx <- unlist(lapply(split(seq_along(labels), labels), function(ix) {
    if (length(ix) > MAX_CELLS_PER_TYPE) sample(ix, MAX_CELLS_PER_TYPE) else ix
  }), use.names = FALSE)
  idx <- sort(idx)
  kept <- table(labels[idx])
  cat(sprintf("DWLS: %d of %d cells for the signature build, <= %d per type (%s)\n",
              length(idx), ncol(sc_counts), MAX_CELLS_PER_TYPE,
              paste(sprintf("%s=%d", names(kept), as.integer(kept)), collapse = ", ")))
  sc_counts <- sc_counts[, idx, drop = FALSE]
  labels    <- labels[idx]
}

# WHICH OF DWLS'S TWO SIGNATURE BUILDERS
# --------------------------------------
# The package ships buildSignatureMatrixMAST and buildSignatureMatrixUsingSeurat. They
# are alternatives offered by the authors, not a published method and a workaround.
#
# MAST fits a hurdle model per gene per cell type and is the expensive one: 758 s on
# 11,739 cells, and still ~30 minutes after capping to 2,000 — the cost is dominated by
# the model fitting, not the cell count, so subsampling does not rescue it. The pipeline
# needs a signature per gene set and there are two gene sets per run, so MAST costs
# roughly an hour of every run for one method's preprocessing.
#
# The Seurat builder was tried as a cheaper alternative and is NOT faster: the cost is
# dominated by the condition-number search over gene counts, which both builders share,
# not by the differential-expression step. MAST is kept because it is the one verified
# end to end here (exit 0, valid proportions, on the real 11,739-cell export).
#
# The cutoffs are left fully open for the reason above: this pipeline has already
# selected the informative genes, identically for every method.
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

# PER-SAMPLE ISOLATION.
#
# solveDampenedWLS succeeds on most mixtures and fails on some — the signature and the
# solver are fine, but a particular sample can trim down to a degenerate system. Letting
# one such sample abort the whole method sent DWLS to the Python fallback for the entire
# cohort, and reported it as though the package could not run at all.
#
# A failed sample becomes NA, which is what the rest of this pipeline already means by a
# failed sample: project_to_simplex turns it into NaN and the benchmark counts it in
# n_failed_samples rather than scoring it. The count and the first error are printed, so
# "DWLS failed on 3 of 200 samples" is visible rather than inferred.
n_types <- ncol(signature)
first_error <- NULL
n_failed <- 0L

props <- t(sapply(colnames(bulk_mat), function(s) {
  out <- tryCatch({
    trimmed <- trimData(signature, bulk_mat[shared, s])
    solveDampenedWLS(trimmed$sig, trimmed$bulk)
  }, error = function(e) {
    if (is.null(first_error)) first_error <<- conditionMessage(e)
    n_failed <<- n_failed + 1L
    setNames(rep(NA_real_, n_types), colnames(signature))
  })
  # a solver that returns the wrong shape is also a failure, not something to reshape
  if (length(out) != n_types) {
    n_failed <<- n_failed + 1L
    out <- setNames(rep(NA_real_, n_types), colnames(signature))
  }
  out
}))
rownames(props) <- colnames(bulk_mat)

if (n_failed > 0) {
  cat(sprintf("DWLS: %d of %d samples failed to solve and are NA. First error: %s\n",
              n_failed, ncol(bulk_mat), first_error))
  if (n_failed == ncol(bulk_mat)) {
    stop(sprintf("DWLS failed on every sample. First error: %s", first_error))
  }
}

write_proportions(props, args$cell_types, args$out)
cat("DWLS complete\n")
