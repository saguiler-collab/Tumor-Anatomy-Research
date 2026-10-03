#!/usr/bin/env Rscript
# run_gimicc.R -- GIMiCC 0.99.1 (genuine; the user's copy in `pipeline_packages / repos/GIMiCC-main`,
# byte-identical to github.com/SalasLab/GIMiCC @ 26cb8a15; library ExperimentHub EH9483).
#
# Pike SC, Wiencke JK, Zhang Z, Molinaro AM, Hansen HM, Koestler DC, Christensen BC, Kelsey KT,
# Salas LA. Acta Neuropathol Commun 2024;12. Glioma-specific hierarchical methylation deconvolution.
#
# Rules fixed before any TCGA run: prespecified/gimicc_truth_confirmation.md.
#   args: in (CpG x samples TSV, Xena layout), tumor_type, h, out (CSV, samples x types, percent),
#         shuffle_seed (optional: permute CpG labels -- the negative control)
# Missing values are NEVER imputed: a CpG missing in any sample is dropped (complete cases), and the
# count is printed alongside GIMiCC's own per-layer missing-probe messages.

# GIMiCC calls query() without importing it (it lives in AnnotationHub, attached with ExperimentHub);
# the authors' vignette attaches ExperimentHub first, so the documented usage is followed here.
suppressPackageStartupMessages({library(jsonlite); library(ExperimentHub); library(GIMiCC)})
args <- fromJSON(commandArgs(trailingOnly = TRUE)[[1]])

b <- as.matrix(read.delim(args$`in`, row.names = 1, check.names = FALSE, na.strings = c("", "NA", "NaN")))
n0 <- nrow(b)
b <- b[complete.cases(b), , drop = FALSE]
cat(sprintf("gimicc input: %d CpGs, %d dropped as incomplete (never imputed), %d kept; %d samples\n",
            n0, n0 - nrow(b), nrow(b), ncol(b)))
if (!is.null(args$shuffle_seed)) {
  set.seed(args$shuffle_seed)
  rownames(b) <- sample(rownames(b))
  cat(sprintf("NEGATIVE CONTROL: CpG labels permuted (seed %d)\n", args$shuffle_seed))
}

res <- GIMiCC_Deconvo(b, h = args$h, tumor.type = args$tumor_type)
if (!identical(rownames(res), colnames(b))) {
  stop("GIMiCC's sample names/order differ from the input's; refusing to map them by position")
}
write.csv(res, args$out)
cat(sprintf("gimicc complete: %d samples x %d types (tumor.type %s, h %s)\n",
            nrow(res), ncol(res), args$tumor_type, args$h))
