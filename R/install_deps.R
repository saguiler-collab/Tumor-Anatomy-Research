#!/usr/bin/env Rscript
# install_deps.R — install the four published deconvolution packages.
#
# None of them are on CRAN, so each needs its own source. Run once:
#     Rscript R/install_deps.R
#
# If a package fails to install the pipeline still runs — ivygap.deconv.r_bridge
# detects the absence and the registry falls back to the Python reimplementation,
# recording that it did so. A partial install is a supported state, not a broken one.

options(repos = c(CRAN = "https://cloud.r-project.org"))

need <- function(pkg) !requireNamespace(pkg, quietly = TRUE)

if (need("devtools"))  install.packages("devtools")
if (need("BiocManager")) install.packages("BiocManager")

# Biobase / ExpressionSet underpins MuSiC, Bisque and SCDC.
if (need("Biobase")) BiocManager::install("Biobase", ask = FALSE, update = FALSE)
if (need("xbioc"))   devtools::install_github("renozao/xbioc")

if (need("MuSiC"))     devtools::install_github("xuranw/MuSiC")
if (need("BisqueRNA")) install.packages("BisqueRNA")          # this one IS on CRAN
if (need("DWLS"))      devtools::install_github("dtsoucas/DWLS")
if (need("SCDC"))      devtools::install_github("meichendong/SCDC")

cat("\n--- installation status ---\n")
for (p in c("MuSiC", "DWLS", "BisqueRNA", "SCDC")) {
  cat(sprintf("%-10s %s\n", p,
      if (requireNamespace(p, quietly = TRUE)) "OK" else "NOT INSTALLED"))
}
