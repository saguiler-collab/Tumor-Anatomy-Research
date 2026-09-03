# common.R — shared I/O for the R driver scripts.
#
# Each driver receives one JSON file naming its inputs and its output path. Keeping
# the contract in one place means the Python side has exactly one thing to satisfy.

suppressPackageStartupMessages({
  library(jsonlite)
  library(Biobase)
})

read_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) < 1) stop("usage: Rscript <driver>.R <args.json>")
  fromJSON(argv[[1]])
}

# Bulk: genes x samples, as written by pandas .to_csv() with the gene symbol index.
load_bulk <- function(path) {
  m <- read.csv(path, row.names = 1, check.names = FALSE)
  as.matrix(m)
}

# Single-cell reference: genes x cells, plus per-cell donor and cell-type labels.
# Cells present in one file and not the other are dropped, with a message — a silent
# mismatch here would misattribute cells to the wrong donor.
load_sc_eset <- function(counts_path, meta_path) {
  counts <- as.matrix(read.csv(gzfile(counts_path), row.names = 1, check.names = FALSE))
  meta   <- read.csv(meta_path, row.names = 1, check.names = FALSE)

  shared <- intersect(colnames(counts), rownames(meta))
  if (length(shared) == 0) stop("no shared cell ids between counts and metadata")
  if (length(shared) < ncol(counts)) {
    message(sprintf("dropping %d cells absent from metadata", ncol(counts) - length(shared)))
  }
  counts <- counts[, shared, drop = FALSE]
  meta   <- meta[shared, , drop = FALSE]

  for (col in c("donor", "cell_type")) {
    if (!col %in% colnames(meta)) stop(sprintf("metadata is missing required column '%s'", col))
  }
  ExpressionSet(assayData = counts, phenoData = AnnotatedDataFrame(meta))
}

# Proportions out: samples x cell types, columns forced into the project's roster
# order. A type the method never estimated is written as 0 and named in a message,
# so a structurally absent type is visible rather than inferred from a short table.
write_proportions <- function(prop, cell_types, out_path) {
  prop <- as.data.frame(prop)
  missing <- setdiff(cell_types, colnames(prop))
  if (length(missing) > 0) {
    message("method returned no column for: ", paste(missing, collapse = ", "))
    for (m in missing) prop[[m]] <- 0
  }
  prop <- prop[, cell_types, drop = FALSE]
  write.csv(prop, out_path)
  invisible(prop)
}
