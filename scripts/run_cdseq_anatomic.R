#!/usr/bin/env Rscript
# run_cdseq_anatomic.R — reference-FREE deconvolution of the Ivy GAP anatomic cohort.
#
# Why this is the one method that matters for C4
# ----------------------------------------------
# The project's central open question is whether the ACS ordering is a property of the
# METHODS or of the reference atlas. Every other method in the panel is handed GBmap, so
# none of them can answer it. CDSeq estimates cell types de novo from the bulk alone.
#
# `reference_gep` is deliberately NOT passed. With CDSeq's default beta it affects only the
# LABELLING of estimated cell types, never the estimation — but leaving it out entirely
# removes the question, and labelling is done downstream in Python where it can be done two
# independent ways and compared.
#
#   Rscript scripts/run_cdseq_anatomic.R <counts.csv> <T> <iters> <dilution> <outprefix>
suppressMessages(library(CDSeq))
a <- commandArgs(TRUE); stopifnot(length(a) == 5)
Tn <- as.integer(a[2]); iters <- as.integer(a[3]); dil <- as.numeric(a[4])

m <- read.csv(a[1], row.names = 1, check.names = FALSE)
m <- as.matrix(m); storage.mode(m) <- "double"

# INTEGER COUNTS, declared here and nowhere else. RSEM's expected_count is a posterior
# expectation and is fractional (45.56, not 46); CDSeq is a multinomial read-count model and
# needs integers. Rounding is applied once, visibly, rather than left to the package.
m <- round(m)
m <- m[rowSums(m) > 0, , drop = FALSE]
cat(sprintf("bulk: %d genes x %d samples, %.0f total reads (dilution_factor=%g)\n",
            nrow(m), ncol(m), sum(m), dil))

t0 <- Sys.time()
r <- CDSeq(bulk_data = m, cell_type_number = Tn, mcmc_iterations = iters,
           dilution_factor = dil, block_number = 1, cpu_number = 1)
cat(sprintf("elapsed %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

prop <- r$estProp                       # cell types x samples
gep  <- r$estGEP                        # genes x cell types
if (is.null(rownames(gep))) rownames(gep) <- rownames(m)
cat(sprintf("estProp %s, estGEP %s\n", paste(dim(prop), collapse = "x"),
            paste(dim(gep), collapse = "x")))
write.csv(as.data.frame(prop), paste0(a[5], "_prop.csv"))
write.csv(as.data.frame(gep),  paste0(a[5], "_gep.csv"))
