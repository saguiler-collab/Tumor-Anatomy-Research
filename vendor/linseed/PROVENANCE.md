# linseed, vendored

- **Source:** https://github.com/ctlab/LinSeed, cloned 2026-10-03. The commit is in `COMMIT`
  (2435a8ce..., 2025-08-22). Package version 0.99.3.
- **Reference:** Zaitsev K, Bambouskova M, Swain A, Artyomov MN. Nat Commun 2019;10:2209;
  doi:10.1038/s41467-019-09990-5 (Crossref-verified).
- **Licence:** MIT + LICENSE file (included).
- **Contents:** DESCRIPTION, NAMESPACE, LICENSE, Readme.md, R/, src/, man/, data/, **unmodified**.
  `SHA256SUMS` lists them. Excluded: vignettes/ and Readme_files/ (documentation images).
- **Installed** with `R CMD INSTALL` from this directory, which compiles src/ (Rcpp, RcppArmadillo,
  BH) with the Xcode command-line clang++ and /opt/gfortran.
- **Missing imports installed first** (BiocManager, `update = FALSE`, binaries): GEOquery 2.80.0 and
  combinat 0.0.9. The install also added, as dependencies: clipr 0.8.1, vroom 1.7.1, tzdb 0.5.0,
  selectr 0.8-0, rlang 1.3.0, readr 2.2.0, rentrez 1.2.4, rvest 1.0.5, httr2 1.3.0.
  - **rlang may have replaced an older installed version.** It was listed as a dependency to install
    although ggplot2 and dplyr already import it, and the prior version was not recorded.
  - Every method package used in this study was checked to load afterwards.
