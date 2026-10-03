#!/usr/bin/env Rscript
# run_mixture.R -- MIXTURE 0.0.1 (genuine, vendored from github.com/elmerfer/MIXTURE @ 6332e160,
# installed from vendor/MIXTURE; see vendor/MIXTURE/PROVENANCE.md).
#
# Fernandez EA, Mahmoud YD, Veigas F, Rocha D, et al. Brief Bioinform 2021;22(4):bbaa317;
# doi:10.1093/bib/bbaa317. nu-SVR with a linear kernel (nu tuned over 0.25/0.5/0.75), negative
# coefficients clipped, cell types under a 0.1% noise bound removed recursively; the signature is
# z-scored globally as in CIBERSORT. Nguyen et al. 2024 rank it in the top 10 of reference-based
# methods in all three of their scenarios.
#
# POST-REGISTRATION EXTENSION PANEL (added 2026-10-02). Not part of the registered panel; run and
# reported under its own name, never in a registered statistic.
#
# Configuration fixed before any output: prespecified/mixture_config.md. Package defaults
# throughout (functionMixture = nu.svm.robust.RFE, nullDist = "none"); one core. Inputs are the
# arm's linear profiles (`signature`) and bulk on one gene space.
#
# Returns proportions of the signature profiles (rows sum to 1): an mRNA proportion on the
# signature's scale, so the central cell-size conversion applies (NOT in R_RETURNS_CELL_FRACTIONS).

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(MIXTURE))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sig_mat  <- load_signature(args$signature)

shared <- intersect(rownames(sig_mat), rownames(bulk_mat))
if (length(shared) < 20) {
  stop(sprintf("only %d genes shared between the MIXTURE signature and the bulk", length(shared)))
}
sig_mat  <- sig_mat[shared, , drop = FALSE]
bulk_mat <- bulk_mat[shared, , drop = FALSE]

res <- MIXTURE::MIXTURE(expressionMatrix = as.matrix(bulk_mat), signatureMatrix = as.matrix(sig_mat),
                        useCores = 1L, verbose = FALSE, nullDist = "none")
props <- as.matrix(MIXTURE::GetMixture(res, type = "proportion"))

# MIXTURE names rows by the expression matrix's columns and columns by the signature's. Check
# rather than assume: a renamed or reordered sample must stop the run, never be mapped silently.
if (!identical(rownames(props), colnames(bulk_mat))) {
  stop("MIXTURE's sample names/order differ from the bulk's; refusing to map them by position")
}
if (!identical(colnames(props), colnames(sig_mat))) {
  stop("MIXTURE's cell-type names/order differ from the signature's")
}

bad <- rowSums(!is.finite(props)) > 0
if (any(bad)) {
  cat(sprintf("MIXTURE: %d of %d samples returned no finite composition\n", sum(bad), nrow(props)))
}
write_proportions(props, args$cell_types, args$out)
cat(sprintf("MIXTURE complete: %d samples x %d types on %d genes\n",
            nrow(props), ncol(props), length(shared)))
