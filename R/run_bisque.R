#!/usr/bin/env Rscript
# run_bisque.R — the genuine BisqueRNA package.
#
# use.overlap = FALSE because Ivy GAP shares no subjects with any public GBM
# single-cell atlas. Bisque then estimates its gene-wise assay transform from the
# marginal distributions rather than from paired subjects. That is Bisque's
# documented degraded mode and the result must be reported as such.

# Locate common.R relative to this script. `sys.frame()$ofile` is only set under
# source(); under Rscript the path arrives in commandArgs() as --file=.
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(BisqueRNA))

args <- read_args()
set.seed(args$seed)

bulk_eset <- ExpressionSet(assayData = load_bulk(args$bulk))
sc_eset   <- load_sc_eset(args$sc_counts, args$sc_meta)

# Bisque expects its subject/cell-type columns under fixed names.
pData(sc_eset)$SubjectName <- pData(sc_eset)$donor
pData(sc_eset)$cellType    <- pData(sc_eset)$cell_type

est <- ReferenceBasedDecomposition(
  bulk.eset = bulk_eset, sc.eset = sc_eset,
  markers = NULL, use.overlap = FALSE
)

write_proportions(t(est$bulk.props), args$cell_types, args$out)
cat("Bisque complete (no-overlap mode)\n")
