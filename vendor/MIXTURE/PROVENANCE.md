# MIXTURE, vendored

- **Source:** https://github.com/elmerfer/MIXTURE, cloned 2026-10-02. The commit is in `COMMIT`
  (6332e160..., 2024-08-19, "Bug fix"). Package version 0.0.1 (DESCRIPTION).
- **Reference:** Fernández EA et al. (MIXTURE, a nu-SVR recursive feature extraction model with
  noise constraints). Ranked in the top 10 of reference-based methods in all three scenarios by
  Nguyen et al. 2024 [3].
- **Licence:** DESCRIPTION states "MIT licence". The repository carries no LICENSE file at this
  commit.
- **Contents:** DESCRIPTION, NAMESPACE, README.md, R/, man/, data/ -- the files the package needs to
  install -- **unmodified**. `SHA256SUMS` lists them.
- **Excluded:** wiki/, vignettes/, Paper.R, BRCA_TCGA_MIXTURE_paper.R, the top-level MIXTURE.R and
  LM22.RData. They are documentation and paper scripts; LM22.RData duplicates data/.
- **Installed from this directory** with `R CMD INSTALL`. Its declared plotting and IO dependencies
  (ComplexHeatmap 2.28.0, ggExtra 0.11.0, openxlsx 4.2.9) were installed first; the deconvolution
  path does not use them.
