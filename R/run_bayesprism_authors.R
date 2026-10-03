#!/usr/bin/env Rscript
# run_bayesprism_authors.R -- BayesPrism configured as its AUTHORS prescribe.
#
# Source of every setting: the package's own tutorial ("Tutorial: bulk RNA-seq deconvolution
# using BayesPrism", Chu & Danko; pipeline_packages / repos/BayesPrism-main/
# tutorial_deconvolution.pdf), whose worked example is TCGA-GBM:
#   * cleanup.genes(..., gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),
#     exp.cells = 5) -- "Gene expressed at high magnitude, such as ribosomal protein genes and
#     mitochondrial genes, may dominate the distribution and bias the inference ... We recommend
#     the removal of these genes." chrX/chrY because reference and mixture sexes differ.
#   * select.gene.type(..., gene.type = "protein_coding")
#   * new.prism(..., input.type = "count.matrix", key = <the malignant type>, outlier.cut = 0.01,
#     outlier.fraction = 0.1) -- key makes BayesPrism learn each mixture's OWN malignant expression
#   * run.prism with default gibbs.control / opt.control ("We recommend the use of default").
#
# DEVIATIONS, declared: (1) genes are restricted to this project's marker space before cleanup --
# the tutorial's optional select.marker step on ~16k protein-coding genes is ~10x the compute and
# is not feasible on a 2-core, 8 GB machine; (2) cell.state.labels = cell.type.labels -- the
# tutorial sub-clusters malignant cells per patient, which multiplies the Gibbs state space.
# The sample-specific malignant profile that `key` provides is kept.
#
# NO TIME BUDGET (user directive 2026-10-01). Called by scripts/bayesprism_tcga.py.
suppressPackageStartupMessages({ library(BayesPrism); library(jsonlite); library(Matrix) })
args <- fromJSON(commandArgs(trailingOnly = TRUE)[[1]])
set.seed(args$seed)
stamp <- function(x) { cat(sprintf("STAGE %s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), x)); flush.console() }

sc <- as.matrix(t(as(readMM(args$counts), "CsparseMatrix")))      # cells x genes, dense (tutorial)
colnames(sc) <- readLines(args$genes); meta <- read.csv(args$meta, stringsAsFactors = FALSE)
rownames(sc) <- meta$cell
bk <- as.matrix(read.csv(args$bulk, row.names = 1, check.names = FALSE))   # samples x genes
ct <- meta$cell_type
stamp(sprintf("reference %d cells x %d genes; mixture %d samples", nrow(sc), ncol(sc), nrow(bk)))

if (isTRUE(args$cleanup)) {
  sc <- cleanup.genes(input = sc, input.type = "count.matrix", species = "hs",
                      gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
                      exp.cells = 5)
  sc <- select.gene.type(sc, gene.type = "protein_coding")
  stamp(sprintf("after cleanup.genes + protein_coding: %d genes", ncol(sc)))
}
shared <- intersect(colnames(sc), colnames(bk))
writeLines(shared, sub("\\.csv$", "_genes_used.txt", args$out))
key <- if (is.null(args$key) || identical(args$key, "")) NULL else args$key
prism <- new.prism(reference = sc[, shared], mixture = bk[, shared, drop = FALSE],
                   input.type = "count.matrix", cell.type.labels = ct, cell.state.labels = ct,
                   key = key, outlier.cut = 0.01, outlier.fraction = 0.1)
stamp(sprintf("prism built on %d genes; key = %s; running Gibbs (defaults) on %d core(s)",
              length(shared), if (is.null(key)) "NULL" else key, args$n_cores))
res <- run.prism(prism = prism, n.cores = args$n_cores)
stamp("run.prism complete")
write.csv(get.fraction(bp = res, which.theta = "final", state.or.type = "type"), args$out)
write.csv(get.fraction(bp = res, which.theta = "first", state.or.type = "type"),
          sub("\\.csv$", "_first.csv", args$out))
cat("BayesPrism (authors' configuration) complete\n")
