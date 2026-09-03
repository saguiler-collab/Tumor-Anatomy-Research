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
MIN_SIGNATURE_GENES <- 50
grid <- list(
  list(diff = 0.5, pval = 0.01),   # published defaults, tried first
  list(diff = 0.3, pval = 0.05),
  list(diff = 0.1, pval = 0.10),
  list(diff = 0.0, pval = 0.50)
)

signature <- NULL
shared <- character(0)
for (g in grid) {
  sig <- tryCatch(
    buildSignatureMatrixMAST(scdata = sc_counts, id = labels, path = sig_dir,
                             diff.cutoff = g$diff, pval.cutoff = g$pval),
    error = function(e) NULL
  )
  if (is.null(sig)) next
  sh <- intersect(rownames(sig), rownames(bulk_mat))
  cat(sprintf("DWLS signature at diff.cutoff=%.2f pval.cutoff=%.2f: %d genes, %d shared with bulk\n",
              g$diff, g$pval, nrow(sig), length(sh)))
  if (length(sh) > length(shared)) { signature <- sig; shared <- sh }
  if (length(shared) >= MIN_SIGNATURE_GENES) break
  # buildSignatureMatrixMAST caches per-cutoff files; clear so the next pass recomputes.
  unlink(list.files(sig_dir, full.names = TRUE), recursive = TRUE)
}

if (is.null(signature) || length(shared) < MIN_SIGNATURE_GENES) {
  stop(sprintf(paste0("DWLS signature reached only %d genes shared with the bulk after ",
                      "relaxing to diff.cutoff=0, pval.cutoff=0.5. The input gene space ",
                      "(%d genes) is too small for DWLS to select from on top of it."),
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
