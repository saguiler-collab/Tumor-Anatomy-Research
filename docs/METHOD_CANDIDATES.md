# Deconvolution methods not yet applied -- researched candidates for a wider comparison

*Written 2026-10-02 from the reference papers in `celldecov_reference_papers/` (text searched, not
recalled). Availability checked on this machine the same day; repository addresses checked for
existence (HTTP 200).*

## What the sources say

- **Nguyen et al. 2024** (*Nucleic Acids Res* 52:4761) benchmarked reference-based bulk methods,
  verbatim: "Seven methods (DWLS, DESeq2's unmix, MuSiC, CIBERSORT, ARIC, LinDeconSeq and MIXTURE)
  are ranked in top 10 in all three scenarios, while 13 methods (DWLS, DESeq2's unmix, MuSiC,
  CIBERSORT, MethylResolver, ARIC, LinDeconSeq, FARDEEP, MIXTURE, MySort, AdRoit, ImmuCellAI and
  Scaden) are consistently among the top 15." And: "Linseed and Deblender consistently have the
  highest accuracy scores in their respective category (reference-free and semi-reference-free)."
  And: "In general, reference-based methods outperform reference-free, and semi-reference-free
  methods." They "recommend users analyze their data using multiple methods."
- **Gaspard-Boulinc et al. 2025** (*Nat Rev Genet* 26:828) review deconvolution for **spatial**
  transcriptomics, not bulk: "Cell2location and RCTD are top performing deconvolution tools."
  Ivy GAP samples are laser-microdissected regions -- in effect large spots -- so these could be
  run on them, off-label.
- Avila Cobos 2020 and Sturm 2019 benchmark the classic bulk set (DSA, dtangle, DeconRNASeq,
  EPIC, quanTIseq, xCell, MCP-counter, TIMER). Li 2026 (ReCIDE's own lab) benchmarks on real bulk.

## Already applied in this study

- **Registered panel (15):** NNLS, SVR (CIBERSORT, this project's implementation), CIBERSORTx B-
  and S-mode (implementations), elastic net, Bayesian, hierarchical Bayesian, MuSiC, DWLS, Bisque,
  SCDC, SCDC ENSEMBLE, EPIC, quanTIseq, BayesPrism.
- **Post-registration extension, genuine packages:** FARDEEP, LinDeconSeq, ARIC, RNA-Sieve,
  ReCIDE (Ivy GAP, TCGA-GBM and TCGA-LGG done), DESeq2 `unmix`, MIXTURE (2026-10-02), BayesPrism
  in its authors' configuration.
- **Other arms:** CDSeq (reference-free); Linseed (reference-free, all three arms done 2026-10-03:
  anatomic ACS 0.677 / 0.659, TCGA tumour rho 0.37 / 0.21 profile-named, no lymphoid ordering);
  MCP-counter with three marker sets (the IMC test).
- **Second methylation truth:** GIMiCC [52] (glioma-specific hierarchical methylation deconvolution;
  `prespecified/gimicc_truth_confirmation.md`), run 2026-10-03 to check the lymphoid truth itself.
- **All seven of Nguyen's "top 10 in all scenarios" methods are now covered** (DESeq2's unmix and
  MIXTURE were added on 2026-10-02).

## Candidates

| method | why (source) | needs | gives | here | plan |
|---|---|---|---|---|---|
| **DESeq2 `unmix`** | top 10 in all three scenarios (Nguyen) | signature profiles | proportions | **done** (DESeq2 1.52.0) | ran 2026-10-02 after its gate; MANUSCRIPT §4.10, E2 |
| **MIXTURE** | top 10 in all three scenarios (Nguyen) | signature matrix | proportions | **done**: vendored @6332e160 [50] | gate passed; Ivy GAP ACS 0.954; TCGA frozen + raw/X (MANUSCRIPT §4.10) |
| **Linseed** | best reference-free (Nguyen); requested by the user earlier | bulk only | proportions + profiles; components need naming by markers, as for CDSeq | **installed**: vendored @2435a8ce [51] | gate passed (r 1.00 / control 0.14); **done 2026-10-03**: anatomic ACS 0.677 (profile) / 0.659 (marker), beating the null; TCGA tumour vs ABSOLUTE 0.365 / 0.089 (GBM), 0.210 / -0.226 (LGG); no labelling claims T, NK and B, so no lymphoid ordering (`prespecified/linseed_config.md`) |
| **AdRoit** | top 15 (Nguyen) | single-cell reference | proportions | not installed; github.com/TaoYang-dev/AdRoit | ask |
| **InstaPrism** | fast re-implementation of BayesPrism's model (Hu & Chikina 2024) | single-cell reference | proportions | not installed; github.com/humengying0907/InstaPrism | ask: would make LGG feasible in hours where genuine BayesPrism needs days. Reported as InstaPrism, never as BayesPrism |
| **MethylResolver** | top 15 (Nguyen); **methylation** | HM450 betas | leukocyte fractions | not installed; github.com/darneson/MethylResolver | ask: a second, independent methylation estimate to cross-check EpiDISH's T > B truth |
| **RCTD** (spacexr) | spatial top performer (Gaspard-Boulinc) | single-cell reference | proportions | not installed; github.com/dmcable/spacexr | optional: off-label on microdissected regions |
| dtangle | Avila Cobos benchmark | markers | proportions | not installed; CRAN | low priority |
| DeconRNASeq | Avila Cobos benchmark | signature | proportions | not installed; Bioconductor | low priority |
| Scaden | top 15 (Nguyen) | trains on simulated pseudobulk from the single-cell reference | proportions | vendored | re-assess: the old note ("needs training data this cohort cannot supply") is doubtful, since Scaden simulates its own training data. The real cost is CPU training on this machine |
| Deblender | best semi-reference-free (Nguyen) | MATLAB | -- | -- | not feasible (no MATLAB) |
| MySort, ImmuCellAI | top 15 (Nguyen) | immune-only signatures | immune fractions | -- | not applicable: no tumour or glial columns, so the roster cannot be matched |

**A property worth stating in the paper:** score-based tools (MCP-counter, xCell, ESTIMATE,
GBMDeconvoluteR) can be scored by ACS exactly like proportion methods. ACS only uses each cell
type's order across regions, so it cannot tell a score from a proportion (`docs/WHAT_ACS_IS.md`).
