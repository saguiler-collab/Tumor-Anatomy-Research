#!/usr/bin/env Rscript
# run_bayesprism.R — the genuine BayesPrism package.
#
# BayesPrism is the third named tool in Anatomy_Test.md's method list. It differs from
# every other entry here in kind, not degree: instead of solving a linear system it
# fits a Bayesian model in which the bulk is generated from cell-type-specific
# expression, and it updates those reference profiles toward the bulk rather than
# treating them as fixed. That is the property worth having in this benchmark — the
# reference is a modern droplet atlas and the bulk is 2014 laser-capture microdissection,
# so reference mismatch is the dominant error, and BayesPrism is the one method that
# models it instead of absorbing it.
#
# It consumes CELLS, like MuSiC/Bisque/SCDC, so it needs the cell-level export.

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(BayesPrism))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)              # genes x samples
sc_eset  <- load_sc_eset(args$sc_counts, args$sc_meta)
sc_counts <- exprs(sc_eset)                   # genes x cells
labels    <- as.character(pData(sc_eset)$cell_type)

shared <- intersect(rownames(sc_counts), rownames(bulk_mat))
if (length(shared) < 50) {
  stop(sprintf("only %d genes shared between the reference and the bulk", length(shared)))
}

# BayesPrism wants samples x genes and cells x genes.
bulk_t <- t(bulk_mat[shared, , drop = FALSE])
sc_t   <- t(sc_counts[shared, , drop = FALSE])

# Cost control. BayesPrism's Gibbs sampler scales with cells x genes x types, and this
# pipeline calls it three times per run. Cells are capped PER TYPE — a reference is a
# per-type summary, so the type proportions are irrelevant to it and a proportional
# subsample would starve the rare types. Same reasoning, and same seed discipline, as
# the DWLS driver. Recorded in docs/METHODS.md.
MAX_CELLS_PER_TYPE <- 300
if (nrow(sc_t) > MAX_CELLS_PER_TYPE * length(unique(labels))) {
  idx <- unlist(lapply(split(seq_along(labels), labels), function(ix)
    if (length(ix) > MAX_CELLS_PER_TYPE) sample(ix, MAX_CELLS_PER_TYPE) else ix),
    use.names = FALSE)
  idx <- sort(idx)
  cat(sprintf("BayesPrism: %d of %d cells, <= %d per type\n",
              length(idx), nrow(sc_t), MAX_CELLS_PER_TYPE))
  sc_t <- sc_t[idx, , drop = FALSE]
  labels <- labels[idx]
}

prism <- new.prism(
  reference   = sc_t,
  mixture     = bulk_t,
  input.type  = "count.matrix",
  cell.type.labels = labels,
  cell.state.labels = labels,     # no finer state layer in this roster
  key         = NULL,             # no single malignant key: Tumor is one roster column
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

n_cores <- max(1, min(4, parallel::detectCores() - 1))
res <- run.prism(prism = prism, n.cores = n_cores)

theta <- get.fraction(bp = res, which.theta = "final", state.or.type = "type")
write_proportions(theta, args$cell_types, args$out)
cat("BayesPrism complete\n")
