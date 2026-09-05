#!/usr/bin/env Rscript
# run_quantiseq.R — the genuine quanTIseq package (quantiseqr).
#
# quanTIseq is built for immune deconvolution: it ships its own TIL10 signature and
# estimates ten immune populations plus an "Other" compartment. That makes it a
# different kind of entry in this benchmark, and the difference is the point.
#
# WHAT IT CAN AND CANNOT SPEAK TO HERE
# ------------------------------------
# The roster this project scores has eight columns, and quanTIseq's signature covers
# only part of it: T cells, NK, B, and the monocyte/macrophage compartment map onto the
# roster; Tumor, Endothelial, Oligodendrocyte and Astrocyte have no quanTIseq
# equivalent and land in its "Other" compartment.
#
# So quanTIseq cannot be scored on C1, C2, C3, C4 or C7 — five of the seven constraints,
# including both Tumor claims and both Endothelial ones. It CAN speak to C5 and C6,
# which are the macrophage constraints. Its row is therefore reported with the
# constraints it could evaluate, not padded to look comparable, and the ACS scorer
# already excludes unevaluable pairs from numerator and denominator alike.
#
# Reporting it any other way — zero-filling the types it does not model — would invent
# a composition it never estimated and make it fail constraints it was never shown.

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(quantiseqr))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)

# quanTIseq expects a gene-symbol x sample matrix on a linear scale, which is what this
# pipeline carries throughout (config.NORMALIZATION = "cpm", never log).
est <- run_quantiseq(bulk_mat, signame = "TIL10", is_arraydata = FALSE,
                     is_tumordata = TRUE, scale_mrna = TRUE)

# quantiseqr returns a data.frame with a Sample column plus one column per population.
rownames(est) <- est$Sample
est$Sample <- NULL

# Map TIL10 populations onto the roster. Anything without a roster equivalent stays in
# the unmapped remainder, which is written out rather than discarded.
map <- list(
  T_cell               = c("T.cells.CD4", "T.cells.CD8", "Tregs"),
  NK_cell              = c("NK.cells"),
  B_cell               = c("B.cells"),
  Macrophage_Microglia = c("Macrophages.M1", "Macrophages.M2", "Monocytes",
                           "Dendritic.cells", "Neutrophils")
)

out <- matrix(NA_real_, nrow = nrow(est), ncol = length(args$cell_types),
              dimnames = list(rownames(est), args$cell_types))
for (target in names(map)) {
  cols <- intersect(map[[target]], colnames(est))
  if (length(cols)) out[, target] <- rowSums(est[, cols, drop = FALSE])
}

covered <- names(map)
cat(sprintf("quanTIseq covers %d of %d roster types: %s\n",
            length(covered), length(args$cell_types), paste(covered, collapse = ", ")))
cat("Types it does not model are left NA, not zero-filled.\n")

utils::write.csv(as.data.frame(out), args$out)
cat("quanTIseq complete\n")
