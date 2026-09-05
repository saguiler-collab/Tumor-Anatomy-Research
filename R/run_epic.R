#!/usr/bin/env Rscript
# run_epic.R — the genuine EPIC package.
#
# EPIC solves a constrained least squares against a reference profile, weighting genes
# by their variability, and estimates an "otherCells" fraction for populations the
# reference does not contain. That last part matters here: the roster drops ~7% of the
# atlas's cells (Mono, DC, RG and others), so a method that can say "some of this is
# something I don't model" is measuring a real property of this data rather than being
# forced to attribute it to the nearest column.
#
# EPIC takes a signature matrix, not cells, so unlike MuSiC/Bisque/SCDC it can run from
# the collapsed reference. It is given the same gene space as every other method.

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(EPIC))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sig_mat  <- load_signature(args$signature)

shared <- intersect(rownames(sig_mat), rownames(bulk_mat))
if (length(shared) < 20) {
  stop(sprintf("only %d genes shared between the EPIC reference and the bulk", length(shared)))
}
sig_mat  <- sig_mat[shared, , drop = FALSE]
bulk_mat <- bulk_mat[shared, , drop = FALSE]

# sigGenes: EPIC's variable-gene set. The pipeline has already chosen an informative
# gene space for every method, so all of it is offered rather than filtering twice —
# the same decision, and for the same reason, as DWLS's opened cutoffs.
ref <- list(refProfiles = sig_mat, sigGenes = rownames(sig_mat))

est <- EPIC(bulk = bulk_mat, reference = ref, withOtherCells = TRUE,
            constrainedSum = TRUE, scaleExprs = TRUE)

props <- est$cellFractions
# EPIC appends an "otherCells" column for material the reference cannot explain. It is
# reported rather than dropped: write_proportions renormalises across the roster, and
# the otherCells share is written beside the output so it is not silently absorbed.
if ("otherCells" %in% colnames(props)) {
  other <- props[, "otherCells", drop = FALSE]
  utils::write.csv(other, sub("\\.csv$", "_othercells.csv", args$out))
  cat(sprintf("EPIC otherCells fraction: mean %.4f, max %.4f\n",
              mean(other[, 1]), max(other[, 1])))
  props <- props[, setdiff(colnames(props), "otherCells"), drop = FALSE]
}

write_proportions(props, args$cell_types, args$out)
cat("EPIC complete\n")
