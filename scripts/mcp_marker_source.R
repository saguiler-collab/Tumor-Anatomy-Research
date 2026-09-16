#!/usr/bin/env Rscript
# mcp_marker_source.R — run the GENUINE MCPcounter with one marker set at a time.
#
# GBMdeconvoluteR (Ajaib et al., Neuro-Oncology 2023;25(7):1236-1248) is MCPcounter given
# GBM-tissue-specific markers. Its IMC validation is a near-controlled contrast: one
# algorithm, several marker sources, the same ten tumours and the same protein ground truth.
# This reproduces that contrast inside this project's pipeline. See
# prespecified/imc_anchored_prediction.md, fixed before any of this was run.
#
#   Rscript scripts/mcp_marker_source.R <markers.tsv> <bulk.csv> <out.csv>
#
# markers.tsv: two columns, "HUGO symbols" and "Cell population" — MCPcounter's own format.
# bulk.csv:    genes x samples, LINEAR expression; logged here, once, explicitly.
suppressMessages(library(MCPcounter))
a <- commandArgs(TRUE)
stopifnot(length(a) == 3)

markers <- read.delim(a[1], sep = "\t", check.names = FALSE, colClasses = "character")
stopifnot(all(c("HUGO symbols", "Cell population") %in% colnames(markers)))

bulk <- read.csv(a[2], row.names = 1, check.names = FALSE)
m <- as.matrix(bulk); storage.mode(m) <- "double"

# LOG SCALE, applied here and nowhere else. GBMdeconvoluteR logs its input
# (tinyscalop::exprs_levels(log_scale = TRUE)) because MCPcounter's estimate is a mean over
# marker genes and is otherwise dominated by a handful of high expressers. Declared, fixed
# before the run, and identical across every arm so it cannot favour one marker set.
m <- log2(m + 1)
m <- m[rowSums(m) > 0, , drop = FALSE]

cat(sprintf("markers: %d genes over %d populations; %d present in the bulk\n",
            nrow(markers), length(unique(markers[["Cell population"]])),
            sum(markers[["HUGO symbols"]] %in% rownames(m))))

s <- MCPcounter.estimate(m, featuresType = "HUGO_symbols", genes = markers)
cat(sprintf("scores: %d populations x %d samples\n", nrow(s), ncol(s)))
write.csv(as.data.frame(s), a[3])
