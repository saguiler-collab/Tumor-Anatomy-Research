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
metadata. `cited` means it is quoted from the package's documentation or the source project.
**Nothing in this list is now unverified.**

---

## Introduction — the problem this study addresses

| | reference | status |
|---|---|---|
| **[1]** | F. Avila Cobos, J. Alquicira-Hernandez, J. E. Powell, P. Mestdagh, and K. De Preter, "Benchmarking of cell type deconvolution pipelines for transcriptomics data," *Nat. Commun.*, vol. 11, p. 5650, Nov. 2020, doi: 10.1038/s41467-020-19015-1. | PDF |
| **[2]** | G. Sturm, F. Finotello, F. Petitprez, *et al.*, "Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology," *Bioinformatics*, vol. 35, no. 14, pp. i436–i445, Jul. 2019, doi: 10.1093/bioinformatics/btz363. | PDF |
| **[3]** | H. Nguyen, H. Nguyen, D. Tran, S. Draghici, and T. Nguyen, "Fourteen years of cellular deconvolution: methodology, applications, technical evaluation and outstanding challenges," *Nucleic Acids Res.*, vol. 52, no. 9, pp. 4761–4783, 2024, doi: 10.1093/nar/gkae267. | PDF |
| **[4]** | M. Li, Y. Su, Y. Tang, Y. Lee, and W. Tian, "Evaluating deconvolution methods using real bulk RNA-expression data for robust prognostic insights across cancer types," *Genome Biol.*, vol. 27, p. 38, 2026, doi: 10.1186/s13059-026-03942-1. | PDF |
| **[5]** | R. B. Puchalski, N. Shah, J. Miller, *et al.*, "An anatomic transcriptional atlas of human glioblastoma," *Science*, vol. 360, no. 6389, pp. 660–663, May 2018, doi: 10.1126/science.aaf2666. | cited |

*[1]–[3] establish that deconvolution benchmarks rest on simulated mixtures or flow cytometry.
[4] is the nearest neighbour to this study and shows that rankings derived from pseudobulk do
not transfer to real bulk. [5] is the Ivy GAP source publication and the origin of the anatomic
structure labels the whole method depends on.*

---

## Methods — data sources

| | reference | status |
|---|---|---|
| **[6]** | M. J. Goldman, B. Craft, M. Hastie, *et al.*, "Visualizing and interpreting cancer genomics data via the Xena platform," *Nat. Biotechnol.*, vol. 38, pp. 675–678, 2020, doi: 10.1038/s41587-020-0546-8. | Crossref |
| **[7]** | S. L. Carter, K. Cibulskis, E. Helman, A. McKenna, *et al.*, "Absolute quantification of somatic DNA alterations in human cancer," *Nat. Biotechnol.*, vol. 30, no. 5, pp. 413–421, May 2012, doi: 10.1038/nbt.2203. | cited |
| **[8]** | V. Thorsson, D. L. Gibbs, S. D. Brown, *et al.*, "The Immune Landscape of Cancer," *Immunity*, vol. 48, no. 4, pp. 812–830, Apr. 2018, doi: 10.1016/j.immuni.2018.03.023. | PDF |
| **[9]** | K. Ellrott, M. H. Bailey, G. Saksena, *et al.*, "Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines," *Cell Syst.*, vol. 6, no. 3, pp. 271–281.e7, Mar. 2018, doi: 10.1016/j.cels.2018.03.002. | Crossref |
| **[10]** | C. Ruiz-Moreno, S. F. Salas, E. Samir, *et al.*, "Bidirectional tumor-host interdependence in glioblastoma" (GBmap), *Cancer Cell*, vol. 40, no. 6, p. 639, 2022, doi: 10.1016/j.ccell.2022.05.009. | cited |
| **[11]** | C. Megill, B. Martin, C. Weaver, *et al.*, "CELLxGENE: a performant, scalable exploration platform for high dimensional sparse matrices," *bioRxiv*, 2021, doi: 10.1101/2021.04.05.438318. | cited |
| **[12]** | A. E. Teschendorff, C. E. Breeze, S. C. Zheng, and S. Beck, "A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies" (EpiDISH), *BMC Bioinformatics*, vol. 18, p. 105, 2017, doi: 10.1186/s12859-017-1511-5. | Crossref |

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
| **[18]** | D. Tsoucas, R. Dong, H. Chen, Q. Zhu, G. Guo, and G.-C. Yuan, "Accurate estimation of cell-type composition from gene expression data" (DWLS), *Nat. Commun.*, vol. 10, p. 2975, Jul. 2019, doi: 10.1038/s41467-019-10802-z. | cited |
| **[19]** | T. Chu, Z. Wang, D. Pe'er, and C. G. Danko, "Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology," *Nat. Cancer*, vol. 3, pp. 505–517, 2022, doi: 10.1038/s43018-022-00356-3. | PDF |
| **[20]** | A. M. Newman, C. B. Steen, C. L. Liu, *et al.*, "Determining cell type abundance and expression from bulk tissues with digital cytometry" (CIBERSORTx), *Nat. Biotechnol.*, vol. 37, pp. 773–782, 2019, doi: 10.1038/s41587-019-0114-2. | cited |

---

## Methods — reference-sensitivity arm

| | reference | status |
|---|---|---|
| **[21]** | C. Neftel, J. Laffy, M. G. Filbin, *et al.*, "An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma," *Cell*, vol. 178, no. 4, pp. 835–849.e21, Aug. 2019, doi: 10.1016/j.cell.2019.06.024. | cited |
| **[22]** | S. Darmanis, S. A. Sloan, D. Croote, *et al.*, "Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma," *Cell Rep.*, vol. 21, no. 5, pp. 1399–1410, Oct. 2017, doi: 10.1016/j.celrep.2017.10.030. | Crossref |
| **[23]** | A. Mossi Albiach, J. Janusauskas, J. Kjaer, *et al.*, "Futile wound healing drives mesenchymal-like cell phenotypes in human glioblastoma," *bioRxiv*, 2023, doi: 10.1101/2023.09.01.555882. | cited |
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
| **[28]** | K. Kang, Q. Meng, I. Shats, D. M. Umbach, M. Li, Y. Li, X. Li, and L. Li, "CDSeq: A novel complete deconvolution method for dissecting heterogeneous samples using gene expression data," *PLoS Comput. Biol.*, vol. 15, no. 12, p. e1007510, 2019, doi: 10.1371/journal.pcbi.1007510. | cited |
| **[29]** | K. Menden, M. Marouf, S. Oller, *et al.*, "Deep learning–based cell composition analysis from tissue expression profiles" (Scaden), *Sci. Adv.*, vol. 6, no. 30, p. eaba2619, 2020, doi: 10.1126/sciadv.aba2619. | cited |
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

## Verification record

All 35 entries are checked. The eight that were previously reconstructed from memory were
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
