# Sources

Every reference this project actually depends on. Each DOI below was resolved against
Crossref on 2026-09-10 and the returned author, year, journal and title checked against
the entry — one error was found and corrected that way (an earlier note cited
`10.1186/s13059-014-0550-8` for MAST; that DOI is DESeq2).

Two categories are kept separate on purpose: work whose **result** this project uses,
and work whose **software** it runs. A method's paper and the package that implements it
are not interchangeable, and where this project runs a reimplementation rather than the
published package, that is stated here as it is in the artefacts.

---

## 1 · Data

**Ivy Glioblastoma Atlas Project (Ivy GAP)** — the tissue every anatomic claim is scored on.

> Puchalski RB, Shah N, Miller J, et al. An anatomic transcriptional atlas of human
> glioblastoma. *Science* 360:660–663 (2018). doi:[10.1126/science.aaf2666](https://doi.org/10.1126/science.aaf2666)

- Release used: RNA-Seq of 2014-11-25, `gene_expression_matrix_2014-11-25.zip`.
- Portal: <https://glioblastoma.alleninstitute.org/>
- Clinical and sample metadata: the `/api/v2/gbm/` namespace on that portal
  (`tumor_details.csv`, `rna_seq_samples_details.csv`). See `docs/DATA_SOURCES.md` for how
  those URLs were located after the documented ones went dead.
- Supplies the five anatomic structures — LE, IT, CT, MVP, PAN — assigned by
  H&E-guided laser microdissection.

**GBmap Core** — the single-cell reference. Builds cell-type profiles and the synthetic
mixtures; never scored, and contributes no anatomic claim.

> Ruiz-Moreno C, Salas SM, Samuelsson E, et al. Harmonized single-cell landscape,
> intercellular crosstalk and tumor architecture of glioblastoma. *bioRxiv* (2022).
> doi:[10.1101/2022.08.27.505439](https://doi.org/10.1101/2022.08.27.505439)

- 338,564 cells; obtained via CELLxGENE. 87.1% 10x Chromium by the atlas's own `assay`
  field, which is what makes CIBERSORTx S-mode rather than B-mode the indicated mode.
- One of its constituent studies is cited directly where cross-reference weighting is
  discussed: Neftel C, Laffy J, Filbin MG, et al. An integrative model of cellular
  states, plasticity, and genetics for glioblastoma. *Cell* 178:835–849 (2019).
  doi:[10.1016/j.cell.2019.06.024](https://doi.org/10.1016/j.cell.2019.06.024)

**Allen Human Brain Atlas** — not yet used. Fetched by `scripts/fetch_allen.py` for the
second tissue; the blockers are the constraint file, the roster and the platform, not the
download. <https://human.brain-map.org/static/download>

---

## 2 · Where the constraints come from

The seven pre-registered constraints are not novel biology. Each restates something a
neuropathologist already holds, which is the point — the claim under test is whether
deconvolution methods reproduce them, not whether they are true.

**WHO diagnostic criteria.** Microvascular proliferation and pseudopalisading necrosis
are the two features that define glioblastoma. C3 and C4 (endothelium peaks in MVP) are
close to restatements of a diagnostic criterion.

> Louis DN, Perry A, Wesseling P, et al. The 2021 WHO Classification of Tumors of the
> Central Nervous System: a summary. *Neuro-Oncology* 23:1231–1251 (2021).
> doi:[10.1093/neuonc/noab106](https://doi.org/10.1093/neuonc/noab106)

**Ivy GAP's own sampling definitions.** C1, C2 and C7 follow from how the structures were
selected: LE is defined as infiltrated but preserved brain, CT as dense tumour, IT as
intermediate. Puchalski et al. 2018, above.

**Pseudopalisades are hypoxic** — the premise C5 rests on. This is the paper that
established it, and it is why a myeloid claim about PAN is defensible in advance rather
than a guess.

> Brat DJ, Castellano-Sanchez AA, Hunter SB, et al. Pseudopalisades in glioblastoma are
> hypoxic, express extracellular matrix proteases, and are formed by an actively
> migrating cell population. *Cancer Res* 64:920–927 (2004).
> doi:[10.1158/0008-5472.CAN-03-2073](https://doi.org/10.1158/0008-5472.CAN-03-2073)

**Myeloid cells are recruited into the glioma microenvironment**, including the
perivascular compartment — behind C5 and C6.

> Hambardzumyan D, Gutmann DH, Kettenmann H. The role of microglia and macrophages in
> glioma maintenance and progression. *Nat Neurosci* 19:20–27 (2016).
> doi:[10.1038/nn.4185](https://doi.org/10.1038/nn.4185)

> Hanahan D, Coussens LM. Accessories to the crime: functions of cells recruited to the
> tumor microenvironment. *Cancer Cell* 21:309–322 (2012).
> doi:[10.1016/j.ccr.2012.02.022](https://doi.org/10.1016/j.ccr.2012.02.022)

Per-constraint rationale and evidence are recorded in
`ivygap/anatomic/constraints.py` and travel with every run in
`results/anatomic/constraint_file.json`. Note that those `evidence` fields name the
literature descriptively; this file is where the citations behind them are pinned.

---

## 3 · Deconvolution methods

Ordered as the leaderboard reports them. "Run as" states what actually executed — the
distinction between a published package and a reimplementation is enforced in the
artefacts and repeated here.

| Method | Reference | Run as |
|---|---|---|
| **MuSiC** | Wang X, Park J, Susztak K, Zhang NR, Li M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. *Nat Commun* 10:380 (2019). doi:[10.1038/s41467-018-08023-x](https://doi.org/10.1038/s41467-018-08023-x) | R package `MuSiC` 1.0.0 |
| **DWLS** | Tsoucas D, Dong R, Chen H, Zhu Q, Guo G, Yuan GC. Accurate estimation of cell-type composition from gene expression data. *Nat Commun* 10:2975 (2019). doi:[10.1038/s41467-019-10802-z](https://doi.org/10.1038/s41467-019-10802-z) | R package `DWLS` 0.1.0, with a Python reimplementation as a disclosed wall-clock fallback |
| **Bisque** | Jew B, Alvarez M, Rahmani E, et al. Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. *Nat Commun* 11:1971 (2020). doi:[10.1038/s41467-020-15816-6](https://doi.org/10.1038/s41467-020-15816-6) | R package `BisqueRNA` 1.0.5 — in its no-overlap mode, which is disclosed |
| **SCDC** | Dong M, Thennavan A, Urrutia E, et al. SCDC: bulk gene expression deconvolution by multiple single-cell RNA sequencing references. *Brief Bioinform* 22:416–427 (2021). doi:[10.1093/bib/bbz166](https://doi.org/10.1093/bib/bbz166) | R package `SCDC` 0.0.0.9000. ENSEMBLE with one reference is SCDC, and says so |
| **EPIC** | Racle J, de Jonge K, Baumgaertner P, Speiser DE, Gfeller D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. *eLife* 6:e26476 (2017). doi:[10.7554/eLife.26476](https://doi.org/10.7554/eLife.26476) | R package `EPIC` 1.1.7 |
| **quanTIseq** | Finotello F, Mayer C, Plattner C, et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. *Genome Med* 11:34 (2019). doi:[10.1186/s13073-019-0638-6](https://doi.org/10.1186/s13073-019-0638-6) | R package `quantiseqr` 1.20.0. Ships its own TIL10 signature covering 4 of the 8 roster types, so it is reported but not ranked |
| **BayesPrism** | Chu T, Wang Z, Pe'er D, Danko CG. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. *Nat Cancer* 3:505–517 (2022). doi:[10.1038/s43018-022-00356-3](https://doi.org/10.1038/s43018-022-00356-3) | R package `BayesPrism` 2.2.3, with a disclosed wall-clock fallback |
| **CIBERSORT** (`svr`) | Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. *Nat Methods* 12:453–457 (2015). doi:[10.1038/nmeth.3337](https://doi.org/10.1038/nmeth.3337) | This project's implementation of the published nu-SVR core. The hosted service is licence-gated and is **not** used |
| **CIBERSORTx** (B-mode, S-mode) | Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nat Biotechnol* 37:773–782 (2019). doi:[10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2) | This project's implementation. B-mode from the Methods; S-mode from Supplementary Note 1, p39. Mode selection follows the paper's own Supplementary Table 1d |
| NNLS, Elastic Net, Bayesian, Hierarchical Bayesian | Baselines, not published tools | This project |

**Assessed and deliberately out of scope:** EcoTyper recovers cell *states* and
multicellular *ecotypes*, not the cell-type fractions the constraints order, so it is not
on the leaderboard. Reasoning in `docs/METHODS.md`.

> Luca BA, Steen CB, Matusiak M, et al. Atlas of clinically distinct cell states and
> ecosystems across human solid tumors. *Cell* 184:5482–5496 (2021).
> doi:[10.1016/j.cell.2021.09.014](https://doi.org/10.1016/j.cell.2021.09.014)

---

## 3a · The predecessor project, and what it supplies

The frozen signature matrix in `reference_frozen/` was built by an earlier TCGA-GBM
deconvolution project by the same author and is vendored here unchanged, so the two
projects share a cell-type roster and their results are comparable. That project is
read-only from here.

`Anatomy_Test.md` names a second, orthogonal yardstick that this cohort cannot supply
and TCGA can: **ABSOLUTE** DNA-based tumour purity, which never touches RNA and so is
independent of every method under test. It is listed in `docs/EXTERNAL_ACTIONS.md` as
outstanding rather than used.

> Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA
> alterations in human cancer. *Nat Biotechnol* 30:413–421 (2012).
> doi:[10.1038/nbt.2203](https://doi.org/10.1038/nbt.2203)

> Brennan CW, Verhaak RGW, McKenna A, et al. The somatic genomic landscape of
> glioblastoma. *Cell* 155:462–477 (2013).
> doi:[10.1016/j.cell.2013.09.034](https://doi.org/10.1016/j.cell.2013.09.034)

> Verhaak RGW, Hoadley KA, Purdom E, et al. Integrated genomic analysis identifies
> clinically relevant subtypes of glioblastoma. *Cancer Cell* 17:98–110 (2010).
> doi:[10.1016/j.ccr.2009.12.020](https://doi.org/10.1016/j.ccr.2009.12.020)

---

## 4 · Statistical and computational methods

**ComBat** — empirical-Bayes batch adjustment. Used inside both CIBERSORTx modes.
Implemented here in ~50 lines rather than by adding a Bioconductor dependency, and
labelled as this project's implementation everywhere it appears; `sva` is not installed.

> Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data
> using empirical Bayes methods. *Biostatistics* 8:118–127 (2007).
> doi:[10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

**MAST** — the hurdle model DWLS uses to build its signature, and the step that dominates
its runtime.

> Finak G, McDavid A, Yajima M, et al. MAST: a flexible statistical framework for
> assessing transcriptional changes and characterizing heterogeneity in single-cell RNA
> sequencing data. *Genome Biol* 16:278 (2015).
> doi:[10.1186/s13059-015-0844-5](https://doi.org/10.1186/s13059-015-0844-5)

**LIBSVM** — the nu-SVR solver behind CIBERSORT's core, reached through scikit-learn.

> Chang CC, Lin CJ. LIBSVM: a library for support vector machines. *ACM Trans Intell Syst
> Technol* 2:27 (2011). doi:[10.1145/1961189.1961199](https://doi.org/10.1145/1961189.1961199)

**Not from a paper.** The within-tumour permutation null, the bootstrap over tumours, the
Anatomic Concordance Score and the agreement test are this project's own construction.
They are specified in `Anatomy_Test.md` and implemented in `ivygap/anatomic/`. Anyone
checking them should read the code, not a citation.

---

## 5 · Software

Versions as measured on the analysis machine, 2026-09-10.

| | |
|---|---|
| R | 4.6.0 (2026-04-24) |
| Python | 3.14.5 |
| numpy · pandas · scipy | 2.4.6 · 2.3.3 · 1.18.0 |
| scikit-learn | 1.9.0 |
| anndata · h5py | 0.13.3.post0 · 3.16.0 |
| lifelines | 0.30.3 |
| Bioconductor `Biobase` | 2.72.0 |

**Primary sources consulted directly, beyond the published papers.** These were obtained
and read rather than cited from memory, and two of them changed the implementation:

| Source | What it settled |
|---|---|
| *CIBERSORTx* Supplementary Information (43 pp) | Supplementary Note 1, p39 states the S-mode algorithm in full. S-mode is implemented from it — see `docs/METHODS.md`. |
| *CIBERSORTx* Supplementary Tables 1–4 (`41587_2019_114_MOESM3–6`) | Table 1d records the batch-correction mode used for every deconvolution in the paper, which is the evidence that this configuration calls for S-mode rather than B-mode. |
| `DWLS` CRAN package manual, v0.1.0 (Sistig A, maintainer) | The published defaults and function contract for the R package. |
| `DWLS-master` source (Tsoucas/Yuan lab) | `Deconvolution_functions.R`, the reference implementation. |
| `ecotyper-master` source + Steen CB, Luca BA, Alizadeh AA, Gentles AJ, Newman AM. Profiling cellular ecosystems at single-cell resolution and at scale with EcoTyper. *Methods Mol Biol* (2023), Chapter 4 | Established that EcoTyper recovers cell states and ecotypes rather than cell-type fractions, and is therefore out of scope for this leaderboard. |
| CELLxGENE (Megill C, et al., *bioRxiv* 2021, doi:[10.1101/2021.04.05.438318](https://doi.org/10.1101/2021.04.05.438318)) | The platform GBmap was obtained through; its Ensembl-indexed `var` is why gene-symbol mapping is explicit in `reference.py`. |

R packages are listed with the methods in §3. Every run also writes its own
`method_configs.json` and `implementation_report.json`, which record what actually
executed rather than what was intended — those are authoritative over this table.

---

## 6 · Registration

> Aguilera S. Anatomic concordance: neuropathology as a ground-truth-free criterion for
> selecting a cell-type deconvolution method in glioblastoma. OSF Registries (2026).
> <https://osf.io/nj78m/>

Companion project holding the supplementary table: <https://osf.io/vuh64/>

Constraint freeze hash:
`2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`

---

## A note on what is not cited here

No reference is listed for the central claim — that anatomic concordance can substitute
for ground truth in method selection — because there is no prior work establishing it.
That is the question the study asks. The literature above supplies the tissue, the
reference, the methods under test, and the neuropathology the constraints restate; none
of it supplies the answer.
