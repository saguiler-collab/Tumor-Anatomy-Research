# References

Numbered in **order of first use in the manuscript**, following the convention of the CJSJ
format reference (Choi, *Columbia Junior Science Journal*): IEEE-style bracketed numerals,
full bibliographic detail, DOI where one exists.

**Two things to know before using this list.**

The manuscript body is still a skeleton — `docs/MANUSCRIPT.md` carries its numbers and its
`> WRITE:` placeholders, not its prose. Only twelve of these sources are named in the draft as
it stands. The order below is therefore the order in which each source is **first needed**,
walking the paper's own section sequence (Introduction → Methods → Results → Discussion →
Limitations). Renumber after the prose is written if the flow moves.

The **Status** column records how each citation was checked. `PDF` means a copy is on disk
and was read. `Crossref` means the record was resolved against the Crossref API on 2026-09-28
— authors, title, journal, volume, pages and year taken from the publisher's own deposited
metadata. `cited` meant quoted from a package's documentation or the source project without being resolved -- the class that hid the [10] error (D25). **On 2026-10-01 every remaining `cited` entry ([5], [7], [11], [18], [20], [21], [23], [28], [29]) was resolved against Crossref; first author, year and title match in all nine, and their status now reads `Crossref`.** Nothing in this list is now unverified: no entry rests on an unresolved citation.

---

## Introduction — the problem this study addresses

| | reference | status |
|---|---|---|
| **[1]** | F. Avila Cobos, J. Alquicira-Hernandez, J. E. Powell, P. Mestdagh, and K. De Preter, "Benchmarking of cell type deconvolution pipelines for transcriptomics data," *Nat. Commun.*, vol. 11, p. 5650, Nov. 2020, doi: 10.1038/s41467-020-19015-1. | PDF |
| **[2]** | G. Sturm, F. Finotello, F. Petitprez, *et al.*, "Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology," *Bioinformatics*, vol. 35, no. 14, pp. i436–i445, Jul. 2019, doi: 10.1093/bioinformatics/btz363. | PDF |
| **[3]** | H. Nguyen, H. Nguyen, D. Tran, S. Draghici, and T. Nguyen, "Fourteen years of cellular deconvolution: methodology, applications, technical evaluation and outstanding challenges," *Nucleic Acids Res.*, vol. 52, no. 9, pp. 4761–4783, 2024, doi: 10.1093/nar/gkae267. | PDF |
| **[4]** | M. Li, Y. Su, Y. Tang, Y. Lee, and W. Tian, "Evaluating deconvolution methods using real bulk RNA-expression data for robust prognostic insights across cancer types," *Genome Biol.*, vol. 27, p. 38, 2026, doi: 10.1186/s13059-026-03942-1. | PDF |
| **[5]** | R. B. Puchalski, N. Shah, J. Miller, *et al.*, "An anatomic transcriptional atlas of human glioblastoma," *Science*, vol. 360, no. 6389, pp. 660–663, May 2018, doi: 10.1126/science.aaf2666. | Crossref |

*[1]–[3] establish that deconvolution benchmarks rest on simulated mixtures or flow cytometry.
[4] is the nearest neighbour to this study and shows that rankings derived from pseudobulk do
not transfer to real bulk. [5] is the Ivy GAP source publication and the origin of the anatomic
structure labels the whole method depends on.*

---

## Methods — data sources

| | reference | status |
|---|---|---|
| **[6]** | M. J. Goldman, B. Craft, M. Hastie, *et al.*, "Visualizing and interpreting cancer genomics data via the Xena platform," *Nat. Biotechnol.*, vol. 38, pp. 675–678, 2020, doi: 10.1038/s41587-020-0546-8. | Crossref |
| **[7]** | S. L. Carter, K. Cibulskis, E. Helman, A. McKenna, *et al.*, "Absolute quantification of somatic DNA alterations in human cancer," *Nat. Biotechnol.*, vol. 30, no. 5, pp. 413–421, May 2012, doi: 10.1038/nbt.2203. | Crossref |
| **[8]** | V. Thorsson, D. L. Gibbs, S. D. Brown, *et al.*, "The Immune Landscape of Cancer," *Immunity*, vol. 48, no. 4, pp. 812–830, Apr. 2018, doi: 10.1016/j.immuni.2018.03.023. | PDF |
| **[9]** | K. Ellrott, M. H. Bailey, G. Saksena, *et al.*, "Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines," *Cell Syst.*, vol. 6, no. 3, pp. 271–281.e7, Mar. 2018, doi: 10.1016/j.cels.2018.03.002. | Crossref |
| **[10]** | C. Ruiz-Moreno, S. M. Salas, E. Samuelsson, M. Minaeva, *et al.*, "Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap," *Neuro-Oncology*, vol. 27, pp. 2281–2295, 2025, doi: 10.1093/neuonc/noaf113. The CELLxGENE dataset used here (collection 999f2a15-3d7e-440b-96ae-2c806799c08c) cites the preprint: "Harmonized single-cell landscape, intercellular crosstalk and tumor architecture of glioblastoma," *bioRxiv*, 2022, doi: 10.1101/2022.08.27.505439. *Corrected 2026-10-01 -- see note below.* | Crossref |
| **[11]** | C. Megill, B. Martin, C. Weaver, *et al.*, "CELLxGENE: a performant, scalable exploration platform for high dimensional sparse matrices," *bioRxiv*, 2021, doi: 10.1101/2021.04.05.438318. | Crossref |
| **[12]** | A. E. Teschendorff, C. E. Breeze, S. C. Zheng, and S. Beck, "A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies" (EpiDISH), *BMC Bioinformatics*, vol. 18, p. 105, 2017, doi: 10.1186/s12859-017-1511-5. | Crossref |

*Correction, 2026-10-01 (`docs/OPEN_DEFECTS.md` D25).* Until today [10] paired GBmap's authors with the title and DOI of a different paper -- Ravi *et al.*, "Spatially resolved multi-omics deciphers bidirectional tumor-host interdependence in glioblastoma", *Cancer Cell* 40:639 (2022) -- under a `cited` status the Crossref pass did not cover. The entry above is resolved against Crossref and against the atlas file's own citation metadata. Every dataset citation is now re-verified on each build of `docs/DATA_INVENTORY.md`.

*[6] serves every TCGA matrix used here. [7] is the DNA copy-number purity that never touches
RNA, which is what makes the accuracy arm independent. [8] supplies the aggregate leukocyte
fraction and the published CIBERSORT comparison. [9] supplies IDH1 status. [10] is the
single-cell atlas every method solves against, obtained through [11]. [12] supplies the
per-cell-type lymphoid truth on which the headline result rests.*

---

## Methods — the deconvolution panel

Eight published tools ran; the other estimators are classical regression baselines written for
this project and have no external source (`docs/METHODS.md`).

| | reference | status |
|---|---|---|
| **[13]** | X. Wang, J. Park, K. Susztak, N. R. Zhang, and M. Li, "Bulk tissue cell type deconvolution with multi-subject single-cell expression reference" (MuSiC), *Nat. Commun.*, vol. 10, p. 380, Jan. 2019, doi: 10.1038/s41467-018-08023-x. | PDF |
| **[14]** | M. Dong, A. Thennavan, E. Urrutia, *et al.*, "SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references," *Brief. Bioinform.*, vol. 22, no. 1, pp. 416–427, 2021, doi: 10.1093/bib/bbz166. | PDF |
| **[15]** | B. Jew, M. Alvarez, E. Rahmani, *et al.*, "Accurate estimation of cell composition in bulk expression through robust integration of single-cell information" (BisqueRNA), *Nat. Commun.*, vol. 11, p. 1971, Apr. 2020, doi: 10.1038/s41467-020-15816-6. | PDF |
| **[16]** | J. Racle and D. Gfeller, "EPIC: A Tool to Estimate the Proportions of Different Cell Types from Bulk Gene Expression Data," in *Bioinformatics for Cancer Immunotherapy*, Methods Mol. Biol., vol. 2120, S. Boegel, Ed., 2020, ch. 17, doi: 10.1007/978-1-0716-0327-7_17. | PDF |
| **[17]** | C. Plattner, F. Finotello, and D. Rieder, "Deconvoluting tumor-infiltrating immune cells from RNA-seq data using quanTIseq," *Methods Enzymol.*, vol. 636, pp. 261–285, 2020, doi: 10.1016/bs.mie.2019.05.056. | PDF |
| **[18]** | D. Tsoucas, R. Dong, H. Chen, Q. Zhu, G. Guo, and G.-C. Yuan, "Accurate estimation of cell-type composition from gene expression data" (DWLS), *Nat. Commun.*, vol. 10, p. 2975, Jul. 2019, doi: 10.1038/s41467-019-10802-z. | Crossref |
| **[19]** | T. Chu, Z. Wang, D. Pe'er, and C. G. Danko, "Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology," *Nat. Cancer*, vol. 3, pp. 505–517, 2022, doi: 10.1038/s43018-022-00356-3. | PDF |
| **[20]** | A. M. Newman, C. B. Steen, C. L. Liu, *et al.*, "Determining cell type abundance and expression from bulk tissues with digital cytometry" (CIBERSORTx), *Nat. Biotechnol.*, vol. 37, pp. 773–782, 2019, doi: 10.1038/s41587-019-0114-2. | Crossref |

---

## Methods — reference-sensitivity arm

| | reference | status |
|---|---|---|
| **[21]** | C. Neftel, J. Laffy, M. G. Filbin, *et al.*, "An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma," *Cell*, vol. 178, no. 4, pp. 835–849.e21, Aug. 2019, doi: 10.1016/j.cell.2019.06.024. | Crossref |
| **[22]** | S. Darmanis, S. A. Sloan, D. Croote, *et al.*, "Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma," *Cell Rep.*, vol. 21, no. 5, pp. 1399–1410, Oct. 2017, doi: 10.1016/j.celrep.2017.10.030. | Crossref |
| **[23]** | A. Mossi Albiach, J. Janusauskas, J. Kjaer, *et al.*, "Futile wound healing drives mesenchymal-like cell phenotypes in human glioblastoma," *bioRxiv*, 2023, doi: 10.1101/2023.09.01.555882. | Crossref |
| **[24]** | K. Siletti, R. Hodge, A. Mossi Albiach, *et al.*, "Transcriptomic diversity of cell types across the adult human brain," *Science*, vol. 382, no. 6667, p. eadd7046, Oct. 2023, doi: 10.1126/science.add7046. | Crossref |

---

## Results — external validation and comparison

| | reference | status |
|---|---|---|
| **[25]** | S. Ajaib, D. Lodha, S. Pollock, *et al.*, "GBMdeconvoluteR accurately infers proportions of neoplastic and immune cell populations from bulk glioblastoma transcriptomic data," *Neuro-Oncology*, vol. 25, no. 7, pp. 1236–1248, 2023, doi: 10.1093/neuonc/noad021. | PDF |
| **[26]** | E. Becht, N. A. Giraldo, L. Lacroix, *et al.*, "Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression" (MCPcounter), *Genome Biol.*, vol. 17, p. 218, 2016, doi: 10.1186/s13059-016-1070-5. | Crossref |
| **[27]** | R. G. W. Verhaak, K. A. Hoadley, E. Purdom, *et al.*, "Integrated Genomic Analysis Identifies Clinically Relevant Subtypes of Glioblastoma Characterized by Abnormalities in PDGFRA, IDH1, EGFR, and NF1," *Cancer Cell*, vol. 17, no. 1, pp. 98–110, Jan. 2010, doi: 10.1016/j.ccr.2009.12.020. | Crossref |

*[25] and [26] supply an external marker-based estimate the anatomic arm is checked against.
[27] defines the mesenchymal subtype used in the failure-factor analysis — and does not exist
for lower-grade glioma, which is why that factor cannot replicate.*

---

## Discussion — tools assessed and excluded

Each was obtained and deliberately not evaluated, with the reason recorded
(`docs/METHOD_LIMITATIONS.md`).

| | reference | status |
|---|---|---|
| **[28]** | K. Kang, Q. Meng, I. Shats, D. M. Umbach, M. Li, Y. Li, X. Li, and L. Li, "CDSeq: A novel complete deconvolution method for dissecting heterogeneous samples using gene expression data," *PLoS Comput. Biol.*, vol. 15, no. 12, p. e1007510, 2019, doi: 10.1371/journal.pcbi.1007510. | Crossref |
| **[29]** | K. Menden, M. Marouf, S. Oller, *et al.*, "Deep learning–based cell composition analysis from tissue expression profiles" (Scaden), *Sci. Adv.*, vol. 6, no. 30, p. eaba2619, 2020, doi: 10.1126/sciadv.aba2619. | Crossref |
| **[30]** | J. Fan, Y. Lyu, Q. Zhang, X. Wang, M. Li, and R. Xiao, "MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data," *Brief. Bioinform.*, vol. 23, no. 6, p. bbac430, 2022, doi: 10.1093/bib/bbac430. | PDF |
| **[31]** | C. B. Steen, B. A. Luca, A. A. Alizadeh, and A. J. Gentles, "Profiling Cellular Ecosystems at Single-Cell Resolution and at Scale with EcoTyper," in *Methods Mol. Biol.*, 2023, pp. 43–71, doi: 10.1007/978-1-0716-2986-4_4. | Crossref |
| **[32]** | G. Sturm, F. Finotello, and M. List, "Immunedeconv: An R Package for Unified Access to Computational Methods for Estimating Immune Cell Fractions from Bulk RNA-Sequencing Data," in *Bioinformatics for Cancer Immunotherapy*, Methods Mol. Biol., vol. 2120, 2020, ch. 16, doi: 10.1007/978-1-0716-0327-7_16. | PDF |

---

## Discussion — positioning and technical background

| | reference | status |
|---|---|---|
| **[33]** | L. C. Gaspard-Boulinc, T. Gautier, S. Guilloux, *et al.*, "Cell-type deconvolution methods for spatial transcriptomics," *Nat. Rev. Genet.*, vol. 26, pp. 828–846, Dec. 2025, doi: 10.1038/s41576-025-00845-y. | PDF |
| **[34]** | F. Liu, H. Qian, and J. Ma, "DNA Methylation-Based Cell Type Deconvolution Reveals the Distinct Cell Composition in Brain Tumor Microenvironment," *bioRxiv*, 2025, doi: 10.1101/2025.01.19.633794. | PDF |
| **[35]** | X. Li, G. Gibson, and P. Qiu, "Gene representation in scRNA-seq is correlated with common motifs at the 3′ end of transcripts," *Front. Bioinform.*, vol. 3, p. 1120290, 2023, doi: 10.3389/fbinf.2023.1120290. | PDF |

---

## Results — post-registration extension panel (added 2026-10-01)

Methods added after registration (`ivygap/deconv/extension.py`; MANUSCRIPT §4.10). Resolved
against Crossref on 2026-10-01 by DOI. **Renumber into first-use order when the prose is
written** — these are appended so the existing numbering does not move.

| | reference | status |
|---|---|---|
| **[36]** | Y. Hao, M. Yan, B. R. Heath, Y. L. Lei, and Y. Xie, "Fast and robust deconvolution of tumor infiltrating lymphocyte from expression profiles using least trimmed squares" (FARDEEP), *PLoS Comput. Biol.*, vol. 15, no. 5, p. e1006976, 2019, doi: 10.1371/journal.pcbi.1006976. | Crossref |
| **[37]** | H. Li, A. Sharma, W. Ming, X. Sun, and H. Liu, "A deconvolution method and its application in analyzing the cellular fractions in acute myeloid leukemia samples" (LinDeconSeq), *BMC Genomics*, vol. 21, p. 652, 2020, doi: 10.1186/s12864-020-06888-1. | Crossref |
| **[38]** | W. Zhang, H. Xu, R. Qiao, B. Zhong, X. Zhang, *et al.*, "ARIC: accurate and robust inference of cell type proportions from bulk gene expression or DNA methylation data," *Brief. Bioinform.*, vol. 23, no. 1, p. bbab362, 2022 (online 2021), doi: 10.1093/bib/bbab362. | Crossref |
| **[39]** | D. D. Erdmann-Pham, J. Fischer, J. Hong, and Y. S. Song, "Likelihood-based deconvolution of bulk gene expression data using single-cell references" (RNA-Sieve), *Genome Res.*, vol. 31, pp. 1794–1806, 2021, doi: 10.1101/gr.272344.120. | Crossref |
| **[40]** | M. Li, Y. Su, Y. Gao, and W. Tian, "ReCIDE: robust estimation of cell type proportions by integrating single-reference-based deconvolutions," *Brief. Bioinform.*, vol. 25, p. bbae422, 2024, doi: 10.1093/bib/bbae422. | Crossref |
| **[41]** | M. P. H. Thomas, S. Ajaib, G. Tanner, A. J. Bulpitt, and L. F. Stead, "GBMPurity: A machine learning tool for estimating glioblastoma tumor purity from bulk RNA-sequencing data," *Neuro-Oncology*, vol. 27, no. 6, pp. 1458–1473, 2025, doi: 10.1093/neuonc/noaf026. | PDF |

*On [40] and [4]:* ReCIDE is the best performer in Li et al. 2026 [4], which comes from the
same laboratory (W. Tian is senior author of both). That is a self-benchmark; an independent
evaluation carries more weight and should be described as such.

---

## Introduction — candidate sources for the glioma immune compartment (added 2026-10-01)

These close the "NEEDS A SOURCE" gap at ¶1.2 of `docs/INTRODUCTION_OUTLINE.md`. They were found
as the primary sources behind a background sentence in Vuyyuru et al. 2026 (*Biology* 15:1624,
on disk), and their metadata is verified — **but none has been read here. Read the passage
you cite before citing it.**

| | reference | status |
|---|---|---|
| **[42]** | N. Abdelfattah, P. Kumar, C. Wang, J.-S. Leu, W. F. Flynn, *et al.*, "Single-cell analysis of human glioma and immune cells identifies S100A4 as an immunotherapy target," *Nat. Commun.*, vol. 13, p. 767, 2022, doi: 10.1038/s41467-022-28372-y. | Crossref |
| **[43]** | N. D. Mathewson, O. Ashenberg, I. Tirosh, S. Gritsch, E. M. Perez, *et al.*, "Inhibitory CD161 receptor identified in glioma-infiltrating T cells by single-cell analysis," *Cell*, vol. 184, pp. 1281–1298.e26, 2021, doi: 10.1016/j.cell.2021.01.022. | Crossref |
| **[44]** | T. Hara, R. Chanoch-Myers, N. D. Mathewson, C. Myskiw, L. Atta, *et al.*, "Interactions between cancer cells and immune cells drive transitions to mesenchymal-like states in glioblastoma," *Cancer Cell*, vol. 39, pp. 779–792.e11, 2021, doi: 10.1016/j.ccell.2021.05.002. | Crossref |

## Extension panel and Extension E2 (added 2026-10-02)

DESeq2 `unmix` entered the post-registration extension panel (MANUSCRIPT §4.10) and carries the
loss-scale stability diagnostic of Extension E2 (§4.12, `docs/EXTENSION_IDENTIFIABILITY.md`).
Resolved against Crossref on 2026-10-02: title, authors, journal, volume, issue, article number
and year match.
[46]–[49] were added the same day for Extension E2's comparison with the 2024–2026 literature. Each was read in full text (Europe PMC open access: PMC12107245, PMC11350143, PMC12837286, PMC11978107), and every passage quoted from it in `docs/EXTENSION_IDENTIFIABILITY.md` §2b was text-searched in that copy.

| | reference | status |
|---|---|---|
| **[45]** | M. I. Love, W. Huber, and S. Anders, "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2," *Genome Biol.*, vol. 15, no. 12, p. 550, 2014, doi: 10.1186/s13059-014-0550-8. | Crossref |
| **[46]** | F. Deng, J. Zou, M. Wang, Y. Gu, J. Wu, L. Gao, *et al.*, "DECEPTICON: a correlation-based strategy for RNA-seq deconvolution inspired by a variation of the Anna Karenina principle," *Brief. Bioinform.*, vol. 26, no. 3, p. bbaf234, 2025, doi: 10.1093/bib/bbaf234. | Crossref |
| **[47]** | B. S. White, A. de Reyniès, A. M. Newman, J. J. Waterfall, A. Lamb, F. Petitprez, *et al.*, "Community assessment of methods to deconvolve cellular composition from bulk gene expression," *Nat. Commun.*, vol. 15, p. 7362, 2024, doi: 10.1038/s41467-024-50618-0. | Crossref |
| **[48]** | A. Dietrich, L. Merotto, K. Pelz, B. Eder, C. Zackl, K. Reinisch, *et al.*, "omnideconv: a unifying framework for using and benchmarking single-cell-informed deconvolution of bulk RNA-seq data," *Genome Biol.*, vol. 27, p. 6, 2026, doi: 10.1186/s13059-026-03955-w. | Crossref |
| **[49]** | L. A. Huuki-Myers, K. D. Montgomery, S. H. Kwon, S. Cinquemani, N. J. Eagles, D. Gonzalez-Padilla, *et al.*, "Benchmark of cellular deconvolution methods using a multi-assay dataset from postmortem human prefrontal cortex," *Genome Biol.*, vol. 26, p. 88, 2025, doi: 10.1186/s13059-025-03552-3. | Crossref |

## Further methods run (added 2026-10-02/03)

MIXTURE completes Nguyen et al.'s [3] seven top-10-in-all-scenarios reference-based methods (extension
panel, MANUSCRIPT §4.10). Linseed is their best reference-free method, run as the reference-free arm's
second instrument beside CDSeq [28]. Both were resolved against Crossref when added; title, authors,
journal, volume and article number match.

| | reference | status |
|---|---|---|
| **[50]** | E. A. Fernández, Y. D. Mahmoud, F. Veigas, D. Rocha, M. Miranda, J. Merlo, *et al.*, "Unveiling the immune infiltrate modulation in cancer and response to immunotherapy by MIXTURE—an enhanced deconvolution method," *Brief. Bioinform.*, vol. 22, no. 4, p. bbaa317, 2021, doi: 10.1093/bib/bbaa317. | Crossref |
| **[51]** | K. Zaitsev, M. Bambouskova, A. Swain, and M. N. Artyomov, "Complete deconvolution of cellular mixtures based on linearity of transcriptional signatures," *Nat. Commun.*, vol. 10, p. 2209, 2019, doi: 10.1038/s41467-019-09990-5. | Crossref |
| **[52]** | S. C. Pike, J. K. Wiencke, Z. Zhang, A. M. Molinaro, H. M. Hansen, D. C. Koestler, *et al.*, "Glioma immune microenvironment composition calculator (GIMiCC): a method of estimating the proportions of eighteen cell types from DNA methylation microarray data," *Acta Neuropathol. Commun.*, vol. 12, no. 1, p. 170, 2024, doi: 10.1186/s40478-024-01874-0. | Crossref |
| **[53]** | F. Klemm, R. R. Maas, R. L. Bowman, M. Kornete, K. Soukup, S. Nassiri, *et al.*, "Interrogation of the microenvironmental landscape in brain tumors reveals disease-specific alterations of immune cells," *Cell*, vol. 181, no. 7, pp. 1643–1660.e17, 2020, doi: 10.1016/j.cell.2020.05.007. | Crossref |


## CPTAC per-sample truth (added 2026-10-06)

The independent cohort and the per-sample DNA truth of the secondary analysis
(`prespecified/cptac_wgs_purity_secondary.md`): the CPTAC glioblastoma study, and the copy-number
method behind GDC's whole-genome purity. All three were resolved against Crossref when added. The
ascatNgs record deposits no page range, so none is given.

| | reference | status |
|---|---|---|
| **[54]** | L.-B. Wang, A. Karpova, M. A. Gritsenko, J. E. Kyle, S. Cao, Y. Li, *et al.*, "Proteogenomic and metabolomic characterization of human glioblastoma," *Cancer Cell*, vol. 39, no. 4, pp. 509–528.e20, 2021, doi: 10.1016/j.ccell.2021.01.006. | Crossref |
| **[55]** | P. Van Loo, S. H. Nordgard, O. C. Lingjærde, H. G. Russnes, I. H. Rye, W. Sun, *et al.*, "Allele-specific copy number analysis of tumors," *Proc. Natl. Acad. Sci. U.S.A.*, vol. 107, no. 39, pp. 16910–16915, 2010, doi: 10.1073/pnas.1009843107. | Crossref |
| **[56]** | K. M. Raine, P. Van Loo, D. C. Wedge, D. Jones, A. Menzies, A. P. Butler, *et al.*, "ascatNgs: identifying somatically acquired copy-number alterations from whole-genome sequencing data," *Curr. Protoc. Bioinformatics*, vol. 56, no. 1, 2016, doi: 10.1002/cpbi.17. | Crossref |
---

## Where each source is load-bearing

Remove any entry in this group and a result disappears.

| reference | what it supplies | result it carries |
|---|---|---|
| **[5]** Ivy GAP | 122 anatomic samples, 10 tumours, 5 structures | every ACS measurement |
| **[7]** ABSOLUTE | DNA copy-number purity, 147 GBM / 496 LGG | the accuracy arm, and its independence |
| **[12]** EpiDISH | per-cell-type T / NK / B truth, 155 + 530 samples | **the headline** — 0 of 12 reproduce T > B |
| **[10]** GBmap | the signature matrix every method solves against | all 15 estimators |
| **[6]** Xena | every TCGA matrix | both validation cohorts |

Everything else is positioning, provenance, or a tool that did not enter the panel.

---

## How many of these does the paper actually need?

The CJSJ format reference (Choi) carries **6 references in a ~2,300-word paper**. This list
has **35**. That gap is worth understanding before cutting, because most of it is forced by
the design rather than by padding.

**The difference is what the two papers are.** Choi tests one technique, so one or two
citations cover the method. This study *benchmarks fifteen estimators*, and every published
tool in the panel must be cited or it is being tested anonymously. That alone is eight
references Choi's design never incurs.

| tier | n | refs | can it be cut? |
|---|---|---|---|
| **Load-bearing data** | 5 | [5] [6] [7] [10] [12] | **No.** Remove one and a result disappears. |
| **Panel methods** | 8 | [13]–[20] | **No.** You cannot benchmark software without citing it. |
| **Positioning** | 4 | [1]–[4] | **Keep.** The Introduction's argument is that benchmarks lean on simulated mixtures; [1]–[3] establish it and [4] independently confirms it. |
| Supporting | 9 | [8] [9] [11] [21] [22] [23] [25] [26] [27] | Cut with the analysis. Each attaches to one arm — reference sensitivity, a failure factor, an external comparison. |
| Optional | 9 | [24] [28]–[35] | Cut freely. Tools assessed and excluded, adjacent fields, technical background. |

**Three defensible targets:**

| target | n | what it costs you |
|---|---|---|
| **13** | load-bearing + panel | the Introduction loses its evidence and reads as assertion |
| **17** | + positioning | **recommended.** Every claim in the paper is still sourced |
| **35** | everything | complete, and long for a venue that expects ~6 |

**At 17 you are still citing nearly three times Choi**, and every one of them is doing work:
five datasets the results rest on, eight pieces of software under test, four benchmarks the
study positions against. That is defensible to any reviewer. Going below 13 is not — it would
mean either an unsourced dataset or an uncited method in the leaderboard.

**What to cut first, if you must.** The Optional tier, in this order: [35] technical
background, [33] [34] adjacent fields, [30] [31] [32] tools assessed but never run, [28] [29]
CDSeq and Scaden — though [29] Scaden is worth keeping if you state the exclusion rationale,
since a deliberate exclusion with a reason reads better than silence.

**What the manuscript currently cites:** [1] [2] [3] [4] [5] [12] [27] [31]. The rest are
available in the list and will attach as the prose is written — a Methods section naming
MuSiC, SCDC, Bisque, EPIC, quanTIseq, DWLS, BayesPrism and CIBERSORTx pulls in [13]–[20] by
itself.

---

## Verification record

**Update 2026-10-02:** [45] (DESeq2, Love et al. 2014) and [46]–[49] (DECEPTICON, the DREAM challenge, omnideconv, the DLPFC benchmark) added and resolved against Crossref; the list now has 49 entries.

**Update 2026-10-03:** [50] (MIXTURE) and [51] (Linseed) added and resolved against Crossref; the list now has 51 entries.

**Update 2026-10-03 (later):** [52] (GIMiCC, Pike et al. 2024) added and resolved against Crossref (all nine authors, journal, volume 12, issue 1, article 170, published 2024-10-28); Europe PMC agrees (PMID 39468647, PMC11514818), and the full text was read there. [53] (Klemm et al. 2020, flow cytometry of brain-tumour immune cells) added and resolved against Crossref; its statement on the lymphocyte compartment was read in the PMC author manuscript (PMC8558904). The list now has 53 entries.

**Update 2026-10-06:** [54] (Wang et al. 2021, CPTAC glioblastoma), [55] (Van Loo et al. 2010, ASCAT) and [56] (Raine et al. 2016, ascatNgs) added and resolved against Crossref (title, first six authors, journal, volume, issue, pages where deposited, year). The list now has 56 entries.

**Update 2026-10-01:** the list now has 44 entries; every one is read from a PDF on disk or
resolved against Crossref (the nine former `cited` entries were resolved today), and [10] was
wrong until today (D25, GBmap). The record below is the 2026-09-28 pass, kept as written.

All 35 entries (as of 2026-09-28) are checked. The eight that were previously reconstructed from memory were
resolved against **Crossref** on 2026-09-28 and are marked `Crossref` in the tables above:

| ref | what changed |
|---|---|
| [6] Goldman, Xena | confirmed exactly as written |
| [9] Ellrott, MC3 | confirmed exactly as written |
| **[12] Teschendorff, EpiDISH** | **confirmed exactly as written** — this is the load-bearing one |
| [22] Darmanis | confirmed exactly as written |
| [24] Siletti | confirmed exactly as written |
| [26] Becht, MCPcounter | confirmed exactly as written |
| [27] Verhaak | **title corrected** — the full title continues *"…Characterized by Abnormalities in PDGFRA, IDH1, EGFR, and NF1"* |
| [31] Steen, EcoTyper | **authors, pages and DOI corrected** — Steen CB, Luca BA, Alizadeh AA, Gentles AJ; pp. 43–71; doi 10.1007/978-1-0716-2986-4_4 |

Six were already exact. Two were wrong in detail and are now right: a truncated title and a
chapter whose author order, page range and DOI had all been reconstructed incorrectly.

`tests/test_references_ordered.py` enforces contiguous numbering, that the stated counts match
the markers, that no panel method carries an unchecked citation, and that every entry has a
year.
