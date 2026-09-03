#!/usr/bin/env Rscript
# run_scdc.R — the genuine SCDC package, single-reference or ENSEMBLE.
#
# ENSEMBLE mode expects several reference ExpressionSets and derives its combining
# weights from them. Given one reference it degenerates to SCDC_prop, and this script
# says so rather than reporting a one-member "ensemble".

# Locate common.R relative to this script. `sys.frame()$ofile` is only set under
# source(); under Rscript the path arrives in commandArgs() as --file=.
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(SCDC))

args <- read_args()
set.seed(args$seed)

bulk_eset <- ExpressionSet(assayData = load_bulk(args$bulk))

# The Python side may pass several comma-separated reference exports for ENSEMBLE.
counts_paths <- strsplit(args$sc_counts, ",")[[1]]
meta_paths   <- strsplit(args$sc_meta, ",")[[1]]
esets <- Map(load_sc_eset, counts_paths, meta_paths)

if (isTRUE(args$ensemble) && length(esets) > 1) {
  est <- SCDC_ENSEMBLE(
    bulk.eset = bulk_eset, sc.eset.list = esets,
    ct.varname = "cell_type", sample = "donor", ct.sub = args$cell_types,
    search.length = 0.01, grid.search = TRUE
  )
  props <- wt_prop(est$w_table[1, seq_along(esets)], est$prop.only)
  cat("SCDC ENSEMBLE weights: ", paste(round(est$w_table[1, seq_along(esets)], 3),
                                       collapse = ", "), "\n")
} else {
  if (isTRUE(args$ensemble)) {
    cat("ENSEMBLE requested but only one reference supplied - running SCDC_prop\n")
  }
  est <- SCDC_prop(
    bulk.eset = bulk_eset, sc.eset = esets[[1]],
    ct.varname = "cell_type", sample = "donor", ct.sub = args$cell_types
  )
  props <- est$prop.est.mvw
}

write_proportions(props, args$cell_types, args$out)
cat("SCDC complete\n")
