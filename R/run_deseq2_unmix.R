#!/usr/bin/env Rscript
# run_deseq2_unmix.R -- DESeq2's `unmix` (Bioconductor DESeq2 1.52.0).
#
# Love MI, Huber W, Anders S. Genome Biol 2014;15:550; doi:10.1186/s13059-014-0550-8. `unmix`
# fits non-negative mixing proportions of "pure" profiles with the loss computed in a variance-
# stabilised space. Nguyen et al. 2024 rank it in the top 10 of reference-based methods in all
# three of their scenarios.
#
# POST-REGISTRATION EXTENSION PANEL (added 2026-10-02). Not part of the registered panel; run and
# reported under its own name, never in a registered statistic.
#
# Configuration fixed before any output: prespecified/deseq2_unmix_config.md. Inputs are the
# reference profile matrix (`pure`) and the bulk (`x`) on one gene space, both linear CPM
# (TPM-like), so the documented TPM argument `shift` applies. The documentation leaves `shift`
# to a visual meanSdPlot check; its automatable form is used: the s in a fixed grid minimising
# |Spearman(SD, mean)| of log(x + s) across genes, from the bulk alone. power = 1 (default).
#
# Returns proportions of the pure profiles (rows sum to 1): an mRNA proportion on the signature's
# scale, so the central cell-size conversion applies (NOT in R_RETURNS_CELL_FRACTIONS).

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 0) return(".")
  dirname(normalizePath(f[[1]]))
}
source(file.path(.script_dir(), "common.R"))
suppressPackageStartupMessages(library(DESeq2))

args <- read_args()
set.seed(args$seed)

bulk_mat <- load_bulk(args$bulk)
sig_mat  <- load_signature(args$signature)

shared <- intersect(rownames(sig_mat), rownames(bulk_mat))
if (length(shared) < 20) {
  stop(sprintf("only %d genes shared between the unmix reference and the bulk", length(shared)))
}
sig_mat  <- sig_mat[shared, , drop = FALSE]
bulk_mat <- bulk_mat[shared, , drop = FALSE]

# shift: the flattest mean-SD trend of log(x + s) over the bulk (prespecified grid and rule)
grid <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50)
flatness <- vapply(grid, function(s) {
  l <- log(bulk_mat + s)
  abs(suppressWarnings(stats::cor(rowMeans(l), apply(l, 1, stats::sd), method = "spearman")))
}, numeric(1))
shift <- grid[which.min(flatness)]           # which.min breaks ties toward the smaller s
cat(sprintf("unmix: shift = %g (|rho(SD, mean)| %s)\n", shift,
            paste(sprintf("%g:%.3f", grid, flatness), collapse = " ")))

# ABLATION ONLY (prespecified/unmix_loss_scale_ablation.md): when these environment variables are
# set, they override shift and power. Unset -- every reported unmix row -- nothing changes.
power <- 1
.env_shift <- Sys.getenv("IVYGAP_UNMIX_SHIFT", "")
.env_power <- Sys.getenv("IVYGAP_UNMIX_POWER", "")
if (nzchar(.env_shift)) { shift <- as.numeric(.env_shift); cat(sprintf("unmix: ABLATION shift override = %g\n", shift)) }
if (nzchar(.env_power)) { power <- as.numeric(.env_power); cat(sprintf("unmix: ABLATION power override = %g\n", power)) }

props <- DESeq2::unmix(as.matrix(bulk_mat), as.matrix(sig_mat), shift = shift, power = power,
                       format = "matrix", quiet = TRUE)
props <- as.matrix(props)
rownames(props) <- colnames(bulk_mat)
colnames(props) <- colnames(sig_mat)

bad <- rowSums(!is.finite(props)) > 0
if (any(bad)) {
  cat(sprintf("unmix: %d of %d samples returned no finite composition\n", sum(bad), nrow(props)))
}
write_proportions(props, args$cell_types, args$out)
cat(sprintf("DESeq2 unmix complete: %d samples x %d types on %d genes\n",
            nrow(props), ncol(props), length(shared)))
