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
suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
})

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sc_eset  <- load_sc_eset(args$sc_counts, args$sc_meta)

# MuSiC 1.x changed its interface: music_prop() takes a bulk MATRIX and a
# SingleCellExperiment, where 0.x took two ExpressionSets. Calling it the old way fails
# with `argument "bulk.mtx" is missing, with no default`, the driver falls back to this
# project's Python reimplementation, and the run reports "music" as though MuSiC had
# produced it. Detect which interface is present rather than pinning a version.
sc_sce <- SingleCellExperiment(
  assays  = list(counts = Biobase::exprs(sc_eset)),
  colData = Biobase::pData(sc_eset)
)

if ("bulk.mtx" %in% names(formals(music_prop))) {
  est <- music_prop(
    bulk.mtx     = bulk_mat,
    sc.sce       = sc_sce,
    clusters     = "cell_type",
    samples      = "donor",
    select.ct    = args$cell_types,
    verbose      = TRUE
  )
} else {
  est <- music_prop(
    bulk.eset    = ExpressionSet(assayData = bulk_mat),
    sc.eset      = sc_eset,
    clusters     = "cell_type",
    samples      = "donor",
    select.ct    = args$cell_types,
    verbose      = TRUE
  )
}

# Est.prop.weighted is the weighted (true MuSiC) estimate; Est.prop.allgene is the
# unweighted NNLS baseline MuSiC also returns. Taking the weighted one is the point.
write_proportions(est$Est.prop.weighted, args$cell_types, args$out)
cat("MuSiC complete\n")
