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

sig_dir <- file.path(tempdir(), "dwls_signature")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

signature <- buildSignatureMatrixMAST(
  scdata = sc_counts, id = labels, path = sig_dir, diff.cutoff = 0.5, pval.cutoff = 0.01
)

shared <- intersect(rownames(signature), rownames(bulk_mat))
if (length(shared) < 50) {
  stop(sprintf("only %d genes shared between the DWLS signature and the bulk matrix", length(shared)))
}
signature <- signature[shared, , drop = FALSE]

props <- t(sapply(colnames(bulk_mat), function(s) {
  trimmed <- trimData(signature, bulk_mat[shared, s])
  solveDampenedWLS(trimmed$sig, trimmed$bulk)
}))
rownames(props) <- colnames(bulk_mat)

write_proportions(props, args$cell_types, args$out)
cat("DWLS complete\n")
