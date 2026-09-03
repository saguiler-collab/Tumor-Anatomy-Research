#!/usr/bin/env Rscript
# run_music.R — the genuine MuSiC package.
#
# MuSiC consumes single cells, not a signature matrix: music_prop() builds its own
# design matrix and, crucially, its own cross-subject variance Sigma from the donor
# labels. That Sigma is the entire method, so passing pre-averaged profiles would
# run the function while removing the thing that makes it MuSiC.

# Locate common.R relative to this script. `sys.frame()$ofile` is only set under
# source(); under Rscript the path arrives in commandArgs() as --file=.
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(MuSiC))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
bulk_eset <- ExpressionSet(assayData = bulk_mat)
sc_eset  <- load_sc_eset(args$sc_counts, args$sc_meta)

est <- music_prop(
  bulk.eset    = bulk_eset,
  sc.eset      = sc_eset,
  clusters     = "cell_type",
  samples      = "donor",
  select.ct    = args$cell_types,
  verbose      = TRUE
)

# Est.prop.weighted is the weighted (true MuSiC) estimate; Est.prop.allgene is the
# unweighted NNLS baseline MuSiC also returns. Taking the weighted one is the point.
write_proportions(est$Est.prop.weighted, args$cell_types, args$out)
cat("MuSiC complete\n")
