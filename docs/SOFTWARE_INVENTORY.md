# Software inventory

*Generated 2026-10-07 20:16 EDT by `scripts/build_software_inventory.py`, on macOS-14.8.9-x86_64-i386-64bit-Mach-O with R 4.6.0 and Python 3.14.5. Every version is read from the installed package. Files marked SUPERSEDED in their docstring are not counted.*

**What this is for.** It records which software produced the results, at which version. **It is not a reference list.** The authors of each package say how it should be cited: in R, run `citation("<package>")`; for a Python package, see its documentation. Read the source and build the reference list from it (STS 2027 rules, Appendix 2: "Entrants should generate their reference lists without the use of AI").

## R packages

| package | version | installed from | loaded by | runs |
|---|---|---|---|---|
| `BayesPrism` | 2.2.3 | GitHub Danko-Lab/BayesPrism@19052052a6f3 | `R/run_bayesprism.R`, `R/run_bayesprism_authors.R` | bayesprism, bayesprism_authors |
| `Biobase` | 2.72.0 | https://bioc-release.r-universe.dev | `R/common.R`, `R/run_music.R` | music |
| `BiocManager` | 1.30.27 | CRAN | `R/install_deps.R` | support |
| `BisqueRNA` | 1.0.5 | GitHub cozygene/bisque@8640250dff0b | `R/run_bisque.R` | bisque |
| `DESeq2` | 1.52.0 | https://bioc-release.r-universe.dev | `R/run_deseq2_unmix.R` | deseq2_unmix |
| `dplyr` | 1.2.1 | CRAN | `R/run_recide.R` | recide |
| `DWLS` | 0.1.0 | CRAN | `R/run_dwls.R` | dwls |
| `EPIC` | 1.1.7 | GitHub GfellerLab/EPIC@50a4f404f96c | `R/run_epic.R` | epic |
| `EpiDISH` | 2.28.0 | https://bioc-release.r-universe.dev | `scripts/extract_epidish_probes.py`, `scripts/independent_verification.py`, `scripts/methylation_celltypes.py` | called from Python: extract_epidish_probes, independent_verification, methylation_celltypes |
| `ExperimentHub` | 3.2.2 | Bioconductor 3.23 | `R/run_gimicc.R`, `scripts/extract_gimicc_cpgs.py` | gimicc, called from Python: extract_gimicc_cpgs |
| `FARDEEP` | 1.0.1 | CRAN | `R/run_fardeep.R` | fardeep |
| `GIMiCC` | 0.99.1 | Bioconductor | `R/run_gimicc.R` | gimicc |
| `jsonlite` | 2.0.0 | CRAN | `R/common.R`, `R/run_bayesprism_authors.R`, `R/run_gimicc.R`, `R/run_recide.R` | bayesprism_authors, gimicc, recide |
| `LinDeconSeq` | 0.1 | GitHub lihuamei/LinDeconSeq@20f1aec75012 | `R/run_lindeconseq.R` | lindeconseq |
| `linseed` | 0.99.3 | Bioconductor | `R/run_linseed.R` | linseed |
| `Matrix` | 1.7-5 | CRAN | `R/run_bayesprism_authors.R`, `R/run_recide.R` | bayesprism_authors, recide |
| `MIXTURE` | 0.0.1 | local install | `R/run_mixture.R` | mixture |
| `MuSiC` | 1.0.0 | local install | `R/run_music.R` | music |
| `quadprog` | 1.5-8 | CRAN | `R/run_dwls.R` | dwls |
| `quantiseqr` | 1.20.0 | https://bioc-release.r-universe.dev | `R/run_quantiseq.R` | quantiseq |
| `remotes` | 2.5.0 | CRAN | `R/install_deps.R` | support |
| `SCDC` | 0.0.0.9000 | local install | `R/run_scdc.R` | scdc |
| `SingleCellExperiment` | 1.34.0 | https://bioc-release.r-universe.dev | `R/run_music.R` | music |

*'runs': the method name of each `R/run_<method>.R` bridge that loads the package, and the Python scripts that run it through R; 'support' means it is loaded only for installation or data handling.*

## Python packages

| import | distribution | version | loaded by (files) |
|---|---|---|---|
| `anndata` | anndata | 0.13.3.post0 | `ivygap/data/reference.py`, `scripts/atlas_cell_fractions.py`, `scripts/b_profile_tissue_likeness.py`, `scripts/bisque_anchoring.py`, `scripts/gbmap_t_vs_b.py`, `scripts/nk_reannotation.py` and 1 more |
| `ARIC` | ARIC | 1.0.1 | `ivygap/deconv/extension.py` |
| `certifi` | certifi | 2026.6.17 | `ivygap/anatomic/registration.py`, `scripts/build_data_inventory.py`, `scripts/cptac_wgs_purity.py`, `scripts/fetch_cptac_gbm.py` |
| `h5py` | h5py | 3.16.0 | `ivygap/anatomic/run_anatomic.py`, `ivygap/data/reference.py`, `scripts/albiach_constraint_check.py`, `scripts/atlas_cell_fractions.py`, `scripts/b_profile_tissue_likeness.py`, `scripts/bisque_anchoring.py` and 4 more |
| `lifelines` | lifelines | 0.30.3 | `ivygap/survival/run_survival.py` |
| `matplotlib` | matplotlib | 3.10.9 | `ivygap/paper_style.py`, `scripts/build_figure_index.py`, `scripts/build_paper_figures.py`, `scripts/build_stats_figures.py` |
| `numpy` | numpy | 2.4.6 | `ivygap/anatomic/acs.py`, `ivygap/anatomic/agreement.py`, `ivygap/anatomic/control_calibration.py`, `ivygap/anatomic/run_anatomic.py`, `ivygap/bench/calibration.py`, `ivygap/bench/equal_footing.py` and 93 more |
| `pandas` | pandas | 2.3.3 | `ivygap/anatomic/acs.py`, `ivygap/anatomic/agreement.py`, `ivygap/anatomic/control_calibration.py`, `ivygap/anatomic/coverage.py`, `ivygap/anatomic/run_anatomic.py`, `ivygap/bench/calibration.py` and 104 more |
| `pypdf` | pypdf | 6.17.0 | `scripts/klemm_figure1f.py` |
| `requests` | requests | 2.34.2 | `ivygap/data/download_ivygap.py` |
| `rnasieve` | rnasieve | 0.1.4 | `ivygap/deconv/extension.py` |
| `scipy` | scipy | 1.18.0 | `ivygap/anatomic/agreement.py`, `ivygap/anatomic/control_calibration.py`, `ivygap/bench/metrics.py`, `ivygap/data/clinical.py`, `ivygap/deconv/classical.py`, `ivygap/deconv/controls.py` and 44 more |
| `sklearn` | scikit-learn | 1.9.0 | `ivygap/deconv/classical.py`, `ivygap/deconv/elastic_net.py` |

## Vendored code (copied into this repository, with its provenance)

| folder | provenance file |
|---|---|
| `vendor/MIXTURE/` | `vendor/MIXTURE/DESCRIPTION`, `vendor/MIXTURE/PROVENANCE.md`, `vendor/MIXTURE/README.md` |
| `vendor/ReCIDE/` | `vendor/ReCIDE/LICENSE` |
| `vendor/linseed/` | `vendor/linseed/DESCRIPTION`, `vendor/linseed/LICENSE`, `vendor/linseed/PROVENANCE.md`, `vendor/linseed/Readme.md` |

## Data

Every public dataset, with its accession, its source check and the scripts that read it, is in `docs/DATA_INVENTORY.md` (generated by `scripts/build_data_inventory.py`).
