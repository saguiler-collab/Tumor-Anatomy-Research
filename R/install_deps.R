#!/usr/bin/env Rscript
# install_deps.R — install the four published deconvolution packages.
#
# Run once:
#     Rscript R/install_deps.R
#
# If a package fails the pipeline still runs — ivygap.deconv.r_bridge detects the
# absence and the registry falls back to the Python reimplementation, recording that
# it did so. A partial install is a supported state, not a broken one.
#
# WHAT CHANGED, AND WHY IT MATTERS  (verified 2026-09-03)
# ------------------------------------------------------
# 1. **BisqueRNA is no longer on CRAN.** This script used to say "this one IS on CRAN"
#    and call install.packages("BisqueRNA"), which now fails silently into the Python
#    fallback. It is installed from the authors' repository instead.
# 2. **DWLS now IS on CRAN** (0.1.0), so it no longer needs a GitHub build — but it
#    depends on MAST, which is a Bioconductor package. install.packages() cannot see
#    Bioconductor, so DWLS fails with "dependency 'MAST' is not available", which reads
#    like a network fault. MAST is now installed first.
# 2b. **pkgmaker has been ARCHIVED from CRAN.** xbioc depends on it and SCDC depends on
#    xbioc, so the whole SCDC chain fails at the first link. Both now come from the
#    maintainer's GitHub.
# 3. **Every step retries.** Roughly half of a first attempt here failed with
#    "cannot download any files" — transient network errors, not missing packages.
#    A single attempt reports a package as unavailable when it is merely unlucky, and
#    the pipeline then quietly runs a reimplementation under the published tool's name.
# 4. **Local sources are preferred when present.** MuSiC and SCDC are GitHub-only; if
#    you have already downloaded them, point LOCAL_SOURCES at the unpacked directories
#    and they install from disk without needing GitHub at all.
#
# INSTALLING THESE IS NECESSARY BUT NOT SUFFICIENT
# ------------------------------------------------
# All four drivers in R/ consume CELL-LEVEL single-cell data, not a signature matrix:
# MuSiC builds its cross-subject variance from donor labels, Bisque estimates its assay
# transform across donors, SCDC weights across subjects. With only the collapsed
# `reference_frozen/` signature they cannot run even when installed, and the bridge
# will say so. Supply a single-cell atlas at data/reference/gbmap_core.h5ad to make
# these packages actually usable.

options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4, timeout = 600)

#: Unpacked package sources to prefer over a network fetch, in priority order.
#:
#: MuSiC and SCDC are GitHub-only, and GitHub was the least reliable host during this
#: install. If the sources are already on disk, use them: it is faster, reproducible,
#: and it works with no network at all.
LOCAL_SOURCE_ROOTS <- c(
  "pipeline packages / repos/%s",          # vendored alongside the repo
  "../pipeline packages / repos/%s",
  "~/Downloads/%s-master",
  "~/Downloads/%s"
)

#: package -> the directory name it unpacks to, when it differs from the package name.
LOCAL_DIR_HINTS <- list(
  MuSiC = c("MuSiC/MuSiC-master", "MuSiC-master"),
  SCDC  = c("SCDC/SCDC-master",  "SCDC-master")
)

find_local_source <- function(pkg) {
  cands <- character(0)
  for (hint in c(LOCAL_DIR_HINTS[[pkg]], pkg)) {
    cands <- c(cands,
               path.expand(sprintf("pipeline packages / repos/%s", hint)),
               path.expand(sprintf("../pipeline packages / repos/%s", hint)),
               path.expand(sprintf("~/Downloads/%s", hint)))
  }
  cands <- c(cands, path.expand(sprintf("~/Downloads/%s-master", pkg)))
  for (d in cands) if (dir.exists(d) && file.exists(file.path(d, "DESCRIPTION"))) return(d)
  NULL
}

need <- function(p) !requireNamespace(p, quietly = TRUE)
say  <- function(...) cat(sprintf(...), "\n")

#: Retry, because a first failure here is usually the network rather than the package.
attempt <- function(pkg, expr, tries = 4) {
  for (i in seq_len(tries)) {
    if (!need(pkg)) return(invisible(TRUE))
    say(">> %s (attempt %d/%d)", pkg, i, tries)
    try(eval(expr), silent = TRUE)
    if (!need(pkg)) return(invisible(TRUE))
    Sys.sleep(5)
  }
  invisible(!need(pkg))
}

attempt("remotes",     quote(install.packages("remotes")))
attempt("BiocManager", quote(install.packages("BiocManager")))

# Shared dependencies. xbioc (SCDC) needs pkgmaker, which needs these.
for (p in c("nnls", "ggplot2", "Matrix", "MCMCpack", "Rcpp", "L1pack", "reshape",
            "cowplot", "pheatmap", "assertthat", "registry", "checkmate", "digest",
            "stringr", "pkgmaker")) {
  attempt(p, bquote(install.packages(.(p))))
}

# Biobase / SingleCellExperiment underpin MuSiC, Bisque and SCDC. TOAST is a MuSiC
# dependency; MAST is a DWLS dependency. All four live on Bioconductor, which
# install.packages() cannot see — DWLS otherwise fails with
# "dependency 'MAST' is not available", which reads like a network fault and is not.
for (p in c("Biobase", "SingleCellExperiment", "TOAST", "MAST")) {
  attempt(p, bquote(BiocManager::install(.(p), ask = FALSE, update = FALSE)))
}

# xbioc needs pkgmaker, which has been ARCHIVED from CRAN. The maintainer's GitHub is
# the live source for both, and without them SCDC cannot install at all.
attempt("pkgmaker", quote(remotes::install_github("renozao/pkgmaker", upgrade = "never")))
attempt("xbioc",    quote(remotes::install_github("renozao/xbioc", upgrade = "never")))

# DWLS: on CRAN since 0.1.0, but only after MAST is present.
attempt("DWLS", quote(install.packages("DWLS")))

# BisqueRNA: archived from CRAN, so fetch from the authors.
attempt("BisqueRNA", quote(remotes::install_github("cozygene/bisque", upgrade = "never")))

# MuSiC and SCDC: local source if available, GitHub otherwise.
install_pkg <- function(pkg, gh) {
  local <- find_local_source(pkg)
  if (!is.null(local)) {
    say("   using local source: %s", local)
    attempt(pkg, bquote(remotes::install_local(.(local), upgrade = "never",
                                               force = TRUE)), tries = 2)
  }
  attempt(pkg, bquote(remotes::install_github(.(gh), upgrade = "never")))
}
install_pkg("MuSiC", "xuranw/MuSiC")
install_pkg("SCDC",  "meichendong/SCDC")

cat("\n--- installation status ---\n")
missing <- character(0)
for (p in c("MuSiC", "DWLS", "BisqueRNA", "SCDC")) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (!ok) missing <- c(missing, p)
  cat(sprintf("%-10s %s%s\n", p, if (ok) "OK" else "NOT INSTALLED",
      if (ok) paste0("  ", as.character(packageVersion(p))) else ""))
}
if (length(missing)) {
  cat(sprintf("\n%d package(s) missing: %s\n", length(missing),
              paste(missing, collapse = ", ")))
  cat("The pipeline still runs; those methods fall back to this project's Python\n",
      "reimplementations and every artefact records that they did.\n", sep = "")
} else {
  cat("\nAll four published packages are installed. They still require a cell-level\n",
      "single-cell reference to run — see the header of this file.\n", sep = "")
}
