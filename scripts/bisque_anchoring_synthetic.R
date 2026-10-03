#!/usr/bin/env Rscript
# Planted-truth test of Bisque no-overlap anchoring (OPEN_DEFECTS D23).
#
# Claim under test: with use.overlap = FALSE, BisqueRNA returns, as the COHORT MEAN, the
# single-cell reference's donor-mean composition -- whatever the bulk holds.
#
# Design. Synthetic cells for 4 types and 8 donors whose donor-mean composition (the PRIOR) is
# A .55 B .25 C .15 D .05. Bulk cohorts of 60 samples are mixed from the SAME type profiles,
# with composition drawn around a planted mean:
#   arm "far"     planted mean A .10 B .15 C .25 D .50   -- far from the prior  (the test)
#   arm "matched" planted mean = the prior                -- the control: Bisque should be right
# If the claim holds, "far" returns the prior, not the plant; "matched" is accurate. Per-sample
# Spearman with the plant is reported too: z-scoring keeps deviations, so tracking may survive
# even where the level is replaced.
suppressPackageStartupMessages({ library(BisqueRNA); library(Biobase); library(jsonlite) })
set.seed(0)
G <- 400; types <- c("A", "B", "C", "D"); nd <- 8
prof <- matrix(rgamma(G * 4, shape = 0.5, scale = 40), G, 4, dimnames = list(paste0("g", 1:G), types))
for (k in 1:4) prof[sample(G, 40), k] <- prof[sample(G, 40), k] * 8          # type markers
prior <- c(A = .55, B = .25, C = .15, D = .05)
cells <- list(); meta <- list()
for (d in 1:nd) {
  p <- as.numeric(rgamma(4, prior * 60)); p <- p / sum(p)
  n <- round(400 * p); n[n < 3] <- 3
  for (k in 1:4) {
    lam <- prof[, k] * 0.05
    m <- matrix(rpois(G * n[k], rep(lam, n[k])), G, n[k])
    cells[[length(cells) + 1]] <- m
    meta[[length(meta) + 1]] <- data.frame(cellType = types[k], SubjectName = paste0("d", d), stringsAsFactors = FALSE)[rep(1, n[k]), ]
  }
}
X <- do.call(cbind, cells); md <- do.call(rbind, meta)
colnames(X) <- rownames(md) <- paste0("c", seq_len(ncol(X))); rownames(X) <- rownames(prof)
sc <- ExpressionSet(X, phenoData = AnnotatedDataFrame(md))
# The prior Bisque sees: per-donor cell proportions, averaged over donors.
tab <- prop.table(table(md$SubjectName, md$cellType), 1); prior_obs <- colMeans(tab)[types]

arm <- function(mu) {
  P <- t(sapply(1:60, function(i) { p <- as.numeric(rgamma(4, mu * 30)); p / sum(p) }))
  colnames(P) <- types; rownames(P) <- paste0("s", 1:60)
  B <- prof %*% t(P); B <- B * matrix(rlnorm(length(B), 0, 0.15), nrow(B))
  colnames(B) <- rownames(P)
  be <- ExpressionSet(B)
  r <- ReferenceBasedDecomposition(be, sc, use.overlap = FALSE, verbose = FALSE)
  E <- t(r$bulk.props)[rownames(P), types]
  list(planted_mean = as.list(round(colMeans(P), 4)), estimated_mean = as.list(round(colMeans(E), 4)),
       L1_estimate_to_plant = round(sum(abs(colMeans(E) - colMeans(P))), 4),
       L1_estimate_to_prior = round(sum(abs(colMeans(E) - prior_obs)), 4),
       per_sample_spearman = as.list(round(sapply(types, function(k) cor(E[, k], P[, k], method = "spearman")), 4)))
}
res <- list(
  what_this_is = "Planted-truth test of BisqueRNA no-overlap anchoring (OPEN_DEFECTS D23). Genuine BisqueRNA, synthetic data.",
  bisque_version = as.character(packageVersion("BisqueRNA")),
  reference_prior_donor_mean = as.list(round(prior_obs, 4)),
  far = arm(c(A = .10, B = .15, C = .25, D = .50)),
  matched = arm(prior))
writeLines(toJSON(res, auto_unbox = TRUE, pretty = TRUE), "results/bisque_anchoring_synthetic.json")
for (a in c("far", "matched")) {
  x <- res[[a]]
  cat(sprintf("\n[%s] planted %s\n        estimated %s\n        L1 to plant %.3f | L1 to prior %.3f | per-sample rho %s\n",
      a, paste(unlist(x$planted_mean), collapse = "/"), paste(unlist(x$estimated_mean), collapse = "/"),
      x$L1_estimate_to_plant, x$L1_estimate_to_prior, paste(unlist(x$per_sample_spearman), collapse = "/")))
}
cat(sprintf("\nreference prior (donor mean): %s\n", paste(round(prior_obs, 3), collapse = "/")))
