#!/usr/bin/env Rscript
# run_recide.R -- ReCIDE, vendored unmodified at vendor/ReCIDE (TianLab-Bioinfo/ReCIDE @ 31bbd6b, MIT).
#
# Li M, Su Y, Gao Y, Tian W. ReCIDE: robust estimation of cell type proportions by integrating
# single-reference-based deconvolutions. Brief Bioinform 2024;25:bbae422. doi:10.1093/bib/bbae422
#
# POST-REGISTRATION EXTENSION PANEL. Called by scripts/run_recide.py, never by the registered
# pipeline. TWO DISCLOSED COMPATIBILITY MEASURES, neither of which edits ReCIDE's code:
#   1. `summarize_each` -- defunct in dplyr >= 1.2 -- is defined before ReCIDE is sourced, with
#      its old semantics (apply .funs to every non-grouping column). ReCIDE pins dplyr 1.1.3; the
#      1.1.x series no longer compiles on R 4.6 (it uses non-API entry points R removed).
#   2. The reference is passed as a SPARSE matrix with Seurat.object.assay.version = "v3", so
#      Seurat 5 builds the v3 assays whose @counts/@data slots ReCIDE reads (Seurat 4 API).
# Kernel: DWLS, ReCIDE's default. All other arguments are ReCIDE's defaults.
options(Seurat.object.assay.version = "v3")
summarize_each <- function(.tbl, .funs, ...) dplyr::summarise(.tbl, dplyr::across(dplyr::everything(), .funs, ...))
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
suppressPackageStartupMessages({
  library(jsonlite); library(Matrix)
  source(file.path(.script_dir(), "..", "vendor", "ReCIDE", "ReCIDE_all_function.R"))
})
args <- fromJSON(commandArgs(trailingOnly = TRUE)[[1]])
set.seed(args$seed)
stamp <- function(x) { cat(sprintf("STAGE %s %s\n", format(Sys.time(), "%H:%M:%S"), x)); flush.console() }
X <- readMM(args$counts); X <- as(X, "CsparseMatrix")
rownames(X) <- readLines(args$genes); meta <- read.csv(args$meta, stringsAsFactors = FALSE)
colnames(X) <- meta$cell
bulk <- read.csv(args$bulk, row.names = 1, check.names = FALSE)
stamp(sprintf("reference %d genes x %d cells, %d donors; bulk %d samples", nrow(X), ncol(X),
              length(unique(meta$donor)), ncol(bulk)))
mgm <- MGM_build(SC_ref = X, EXP_df = bulk, label_celltype = meta$cell_type,
                 label_subject = meta$donor, n_cores = args$n_cores)
stamp(sprintf("MGM_build done: %d per-donor signatures", length(mgm)))
res <- ReCIDE_deconvolution(Sig_list = mgm, EXP_df = bulk, Method = "DWLS",
                            n_cores = args$n_cores, n_celltype = 0)
stamp("ReCIDE_deconvolution done")
E <- t(as.matrix(res[["results_final_df"]]))
write.csv(E, args$out)
cat(sprintf("ReCIDE types returned: %s\n", paste(colnames(E), collapse = ", ")))
cat("ReCIDE complete\n")
