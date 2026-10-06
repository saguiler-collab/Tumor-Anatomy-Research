# Data inventory -- every public dataset this study touches

*Generated 2026-10-06 by `scripts/build_data_inventory.py`. Do not edit by hand: every size, hash, consumer and citation check below is measured each run. `docs/supplementary/data_inventory_full.csv` is the same content in one row per dataset; Supplementary Table S7 is its compact form (`scripts/build_supplementary.py`).*

**How to read the checks.** *read by* lists scripts that were opened and found to contain the token that reads the dataset. *source check* compares the local copy with the public server: identical byte size, or, where one side is stored decompressed, the sha256 of the decompressed content on both sides. *citation* is the DOI's Crossref record compared with the first author and year claimed here; the verifier is shown rejecting a known conflation before any citation is trusted.

**Open problems: none.** Every listed file is present, every consumer still reads what it is said to read, every checked copy matches its source, and every citation matches its DOI.

## Summary

| id | dataset | status | accession | used for |
|---|---|---|---|---|
| D01 | Ivy Glioblastoma Atlas Project (Ivy GAP) RNA-seq, normalised FPKM matrix | used | GEO GSE107559; portal well-known file 305873915 | The tissue every anatomic claim is scored on: 270 laser-microdissected samples (122 H&E-selected anatomic samples from 10 tumours across LE, IT, CT, MVP, PAN; 148 ISH-selected cluster samples, never scored as anatomic) |
| D02 | Ivy GAP portal metadata: tumour, RNA-seq sample, sub-block and ISH-gene tables | used | tumor_details.csv, rna_seq_samples_details.csv, sub_block_details.csv, gene_expression_details.csv | Tumour and structure assignment, the donor-tumour join, the clinical table (survival is BLOCKED: no event indicator), and the per-sub-block ISH values used by the no-deconvolution constraint check |
| D03 | Ivy GAP per-sample RSEM gene results (expected counts) | used | 270 x *.genes.results (URLs and sha256 per file in provenance.json) | Fractional expected counts for the count-based arms (CDSeq reference-free ACS, count-model inputs) and closure of the published matrix's provenance chain |
| D04 | TCGA-GBM RNA-seq, STAR gene counts | used | Xena dataset TCGA-GBM.star_counts.tsv | Validation cohort: per-sample Tumor fraction against ABSOLUTE purity (154 samples), the lymphoid ordering against methylation, model fit, extension panel |
| D05 | TCGA-LGG RNA-seq, STAR gene counts | used | Xena dataset TCGA-LGG.star_counts.tsv | Independent second cohort (534 samples; 510 matched to methylation, 496 to ABSOLUTE) for every TCGA analysis |
| D06 | GENCODE v36 comprehensive gene annotation | used | gencode.v36.annotation.gtf.gz | Ensembl id -> gene symbol for both TCGA expression matrices |
| D07 | TCGA-GBM DNA methylation, Illumina HumanMethylation450 | used | Xena dataset TCGA.GBM.sampleMap/HumanMethylation450 | Per-cell-type lymphoid truth (T / NK / B) through EpiDISH -- the methylation side of the headline lymphoid comparison (cohort-level only; see D23) |
| D08 | TCGA-LGG DNA methylation, Illumina HumanMethylation450 | used | Xena dataset TCGA.LGG.sampleMap/HumanMethylation450 | As D07, for LGG |
| D09 | ABSOLUTE tumour purity and ploidy, TCGA PanCanAtlas master calls | used | GDC file 4f277128-f793-4354-a13d-30cc7fe9f6b5 | The DNA ground truth for tumour content, sharing no input with any RNA method (147 GBM / 496 LGG complete cases); ploidy, genome doublings and subclonal fraction as candidate failure factors |
| D10 | MC3 somatic mutation calls, gene level, TCGA-LGG | used | Xena dataset mc3_gene_level/LGG_mc3_gene_level.txt | IDH1 status as a candidate biological failure factor in LGG |
| D11 | Leukocyte fraction from DNA methylation, TCGA PanImmune | used | GDC file 6f75c9d7-5134-4ed1-b8f3-72856c98a4e8 | Aggregate immune truth for the pre-registered immune arm (P1/P2 over-call test) |
| D12 | TCGA-GBM clinical sample annotations (Brennan 2013 study) | used | gbm_tcga_pub2013_clinical_sample.txt (+ _clinical_patient.txt) | Verhaak expression subtype, IDH1 mutation and G-CIMP status as candidate failure factors in GBM |
| D13 | GBmap (core) -- harmonised glioblastoma single-cell atlas | used | collection 999f2a15-3d7e-440b-96ae-2c806799c08c; dataset version 861acfd8-25f0-418b-a445-aa96da232827 (title 'Core GBmap', schema 7.1.0) | The reference every method solves against: the vendored frozen signature (8-type roster, annotation_level_3) and the raw/X-count rebuild that is primary; the source of the synthetic mixtures |
| D14 | Abdelfattah et al. 2022 glioma single-cell RNA-seq (independent atlas) | used | GEO GSE182109 (GSM5518596-GSM5518639; raw supplementary archive GSE182109_RAW); processed matrix + clusters from the Single Cell Portal, whose study accession was not recorded -- content verified against GEO's raw counts | The independent-atlas test of the B-over-T result: a reference not among GBmap's 16 source studies |
| D15 | Neftel et al. 2019 glioblastoma single-cell RNA-seq | used | GEO GSE131928; the Smart-seq2 matrix and cell assignments as distributed on the Broad Single Cell Portal | Reference-dependence arm (a Neftel-only reference) |
| D16 | Darmanis et al. 2017 glioblastoma single-cell RNA-seq (core vs periphery) | used | GEO GSE84465 | Two uses: the only cross-patient core-vs-periphery constraint check (C1, four gates), and a Darmanis-only reference in the reference-dependence arm (a GBmap constituent, so not independent) |
| D17 | Mossi Albiach et al. 2023 spatially sampled glioblastoma single-cell atlas | used | collection 113a558a-e96e-4643-81db-140e95c58578; dataset d45b4ce6-9725-4d79-b97a-70a44158bdbf (donor SL040, schema 7.1.0) | Tests the registered constraints against counted (not inferred) composition per anatomic zone: 3 of 4 testable constraints hold |
| D18 | Ivy GAP in-situ hybridisation (ISH) by anatomic structure | used | gene_expression_details.csv (/api/v2/gbm/ namespace) | A constraint check with no deconvolution anywhere in the chain: CD44 / BIRC5 ISH support C1 and C7 |
| D32 | Flow-cytometry immune composition of brain tumours by IDH status (Klemm et al. 2020, Figure 1F) | used | doi:10.1016/j.cell.2020.05.007 -- Figure 1F (PDF page 3); Table S1 (mmc1.pdf, cohort); Table S2 (mmc2.pdf, gating definitions) | Direct check of the lymphoid truth by IDH status (prespecified/klemm_t_vs_b.md): cohort-mean T, B and NK as % of CD45+ in 17 IDH-mutant and 40 IDH-wildtype gliomas, measured from the figure's vector geometry and validated against two numbers printed in the paper |
| D19 | GBMDeconvoluteR marker sets (Ajaib immune, Ruiz-Moreno level-3 immune, Neftel four-state) and Ajaib's published imaging-mass-cytometry correlations | used | Ajaib_et_al_2022_GBM_Immune_markers.rds, Moreno_et_al_2022_lvl3_immune_markers.rds, Neftel_et_al_2019_four_state_neoplastic_markers.rds | The pre-registered IMC protein test (an external instrument for ACS) and the mesenchymal-program score |
| D20 | MCPcounter default marker genes | used | Signatures/genes.txt | The default-marker arm of the IMC test |
| D21 | EpiDISH blood reference centDHSbloodDMC.m (333 CpGs x 7 blood cell types) | used | EpiDISH::centDHSbloodDMC.m | Turns D07/D08 into per-cell-type lymphoid fractions (RPC mode); a blood reference applied to brain tumour, so used for lymphoid sub-composition only |
| D22 | BayesPrism's bundled gene annotation (gene groups and gene types, GENCODE v22) | used | cleanup.genes() gene groups; select.gene.type() 'protein_coding' | The authors' own gene filtering in the authors'-workflow BayesPrism runs on TCGA (post-registration extension) |
| D31 | GIMiCC reference libraries (glioma-specific hierarchical DNA-methylation deconvolution) | used | EH9483 (GIMiCC_Library.rda, 15 layer libraries); EH9482 (Capper_example_betas.rda, the reproduction control) | The second methylation truth for the lymphoid ordering (post-registration, prespecified/gimicc_truth_confirmation.md): GIMiCC run on the TCGA 450K betas (D07/D08), restricted to its 4,022 library CpGs by scripts/extract_gimicc_cpgs.py |
| D23 | Frozen GBM signature matrix and cell-size factors (from GBmap, by the predecessor project) and its synthetic held-out benchmark | derived | signature_matrix.tsv, cell_size_factors.csv, tcga_benchmark/* | The registered signature; yardstick 1 (500 pseudobulk mixtures, NNLS and SVR only), which shares GBmap with the ACS arm |
| D24 | CIBERSORT relative fractions for TCGA (PanImmune) | on disk, unused | GDC file b3df502e-3594-46ef-9f94-d041a20a0b9a | Obtained for a comparison against a published deconvolution of the same cohort; no script reads it |
| D25 | ABSOLUTE segment tables, TCGA PanCanAtlas | on disk, unused | GDC file 0f4f5701-7b61-41ae-bda9-2805d1ca9781 | Downloaded with D09; nothing reads it |
| D26 | TCGA-LGG GDC clinical and biospecimen tables; GDC whole-exome BAM manifest | on disk, unused | clinical.project-tcga-lgg.2026-09-17, biospecimen.project-tcga-lgg.2026-09-17, gdc_manifest.2026-09-17.191037.txt | Nothing reads the tables |
| D27 | GEO GSE107559 supplementary files (the Ivy GAP portal tables, mirrored) | duplicate | GSE107559_ivygap_rows-genes.csv, GSE107559_ivygap_columns-samples.xlsx | A GEO mirror of D01's tables; nothing reads this copy |
| D28 | Human Brain Cell Atlas v1.0 (Siletti et al.) | considered | collection 283d65eb-dd53-496d-adb7-7570c7caa443 | A normal-brain reference for a second tissue; not obtained |
| D29 | Allen Human Brain Atlas, bulk expression by region | considered | human.brain-map.org downloads | A second tissue for the anatomy test; fetch script exists, the constraint file, roster and platform are the blockers |
| D30 | Spatially resolved multi-omics of glioblastoma (Ravi et al.) | considered |  | Candidate second GBM cohort with spatial structure; not obtained |

## Primary cohort -- the anatomy arm

### D01 · Ivy Glioblastoma Atlas Project (Ivy GAP) RNA-seq, normalised FPKM matrix

- **Status:** used
- **Repository:** Allen Institute for Brain Science, Ivy GAP portal; mirrored at GEO
- **Accession:** GEO GSE107559; portal well-known file 305873915
- **Source:** <https://glioblastoma.alleninstitute.org/static/download.html>
- **Version:** release 2014-11-25 (gene_expression_matrix_2014-11-25.zip)
- **Used for:** The tissue every anatomic claim is scored on: 270 laser-microdissected samples (122 H&E-selected anatomic samples from 10 tumours across LE, IT, CT, MVP, PAN; 148 ISH-selected cluster samples, never scored as anatomic).
- **Citation:** Puchalski RB, et al. Science 360:660-663 (2018), doi:10.1126/science.aaf2666
  - Crossref 10.1126/science.aaf2666: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Ivy GAP', 'anatomic sample'; a mention, not proof of use): II. Methods; 4.9; 4.11; Limitations; Acknowledgements
- **Read by:** `ivygap/data/download_ivygap.py`; `ivygap/data/load_ivygap.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/gene_expression_matrix_2014-11-25.zip` -- the archive as published | 36.0 MB | `fc2ab637b879` | = source size (2026-10-01) |
| `data/raw/ivygap/fpkm_table.csv` -- extracted | 82.2 MB | `27b916d11b73` | -- |
| `data/raw/ivygap/rows-genes.csv` -- extracted | 1.6 MB | `468d0eee470f` | -- |
| `data/raw/ivygap/columns-samples.csv` -- extracted | 49.4 KB | `05bb4dd5cba0` | -- |
| `data/processed/ivygap_bulk_expression.tsv` -- processed: genes x samples | 92.8 MB | `a8ad19045427` | -- |
| `data/processed/ivygap_sample_manifest.tsv` -- processed: one row per sample | 28.7 KB | `eb3d3ff4da64` | -- |

### D02 · Ivy GAP portal metadata: tumour, RNA-seq sample, sub-block and ISH-gene tables

- **Status:** used
- **Repository:** Allen Institute Ivy GAP portal, /api/v2/gbm/ namespace
- **Accession:** tumor_details.csv, rna_seq_samples_details.csv, sub_block_details.csv, gene_expression_details.csv
- **Source:** <https://glioblastoma.alleninstitute.org/api/v2/gbm/tumor_details.csv>
- **Version:** fetched 2026-09-03; sha256 recorded at download (CLINICAL_PROVENANCE.json)
- **Used for:** Tumour and structure assignment, the donor-tumour join, the clinical table (survival is BLOCKED: no event indicator), and the per-sub-block ISH values used by the no-deconvolution constraint check.
- **Citation:** Puchalski RB, et al. Science 360:660-663 (2018), doi:10.1126/science.aaf2666
  - Crossref 10.1126/science.aaf2666: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'survival', 'ISH'; a mention, not proof of use): Title; Abstract; II. Methods; III. Statistical Analysis, Data Analysis and Measurements; 4.1; 4.3; 4.4; 4.7; 4.9; 4.10; 4.11; 4.12; V. Discussion; Limitations; References
- **Read by:** `ivygap/data/portal_metadata.py`; `ivygap/data/clinical.py`; `scripts/fetch_ivygap_ish.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/ivygap/tumor_details.csv` | 2.8 KB | `4accada8f216` | = recorded hash |
| `data/raw/ivygap/rna_seq_samples_details.csv` | 98.5 KB | `6f941bcd4c33` | = recorded hash |
| `data/raw/ivygap/sub_block_details.csv` | 231.4 KB | `1d3ed28782a9` | = recorded hash |
| `data/raw/ivygap/gene_expression_details.csv` | 4.8 MB | `07c5702a2549` | = recorded hash |
| `data/raw/ivygap/CLINICAL_PROVENANCE.json` -- download record | 909 B | `68753389cde8` | -- |

### D03 · Ivy GAP per-sample RSEM gene results (expected counts)

- **Status:** used
- **Repository:** Allen Institute Ivy GAP portal, per-sample well-known files
- **Accession:** 270 x *.genes.results (URLs and sha256 per file in provenance.json)
- **Source:** <https://glioblastoma.alleninstitute.org/static/download.html>
- **Version:** fetched 2026-09-15; GEO GSE107559 states raw data are not provided there
- **Used for:** Fractional expected counts for the count-based arms (CDSeq reference-free ACS, count-model inputs) and closure of the published matrix's provenance chain.
- **Citation:** Puchalski RB, et al. Science 360:660-663 (2018), doi:10.1126/science.aaf2666
  - Crossref 10.1126/science.aaf2666: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'CDSeq'; a mention, not proof of use): 4.11
- **Read by:** `scripts/fetch_ivygap_counts.py`; `scripts/build_ivygap_counts_matrix.py`; `scripts/cdseq_anatomic.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/ivygap_counts` -- 270 per-sample files; each checked against the sha256 recorded when it was fetched | 391.6 MB, 272 files | `not hashed` | 270/270 files = recorded hashes |
| `FPKM/raw_fpkm.csv` -- portal manifest the URLs came from | 72.1 KB | `81037533a5f0` | -- |
| `FPKM/bam_manifest.csv` -- portal BAM manifest -- BAMs deliberately not fetched | 99.8 KB | `fc621e108aca` | -- |
| `data/processed/ivygap_expected_counts.csv.gz` -- assembled matrix | 11.8 MB | `eb4cac2e10b6` | -- |

## Validation cohorts and orthogonal truths -- TCGA

### D04 · TCGA-GBM RNA-seq, STAR gene counts

- **Status:** used
- **Repository:** UCSC Xena, GDC hub (data from the NCI Genomic Data Commons)
- **Accession:** Xena dataset TCGA-GBM.star_counts.tsv
- **Source:** <https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-GBM.star_counts.tsv.gz>
- **Version:** as vendored by the predecessor project (stored decompressed; delivered as log2(x+1) despite the name -- linearised here)
- **Used for:** Validation cohort: per-sample Tumor fraction against ABSOLUTE purity (154 samples), the lymphoid ordering against methylation, model fit, extension panel.
- **Citation:** Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena], doi:10.1038/s41587-020-0546-8
  - Crossref 10.1038/s41587-020-0546-8: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'TCGA-GBM', 'glioblastoma cohort'; a mention, not proof of use): 4.10; 4.12; Limitations; Acknowledgements
- **Read by:** `scripts/build_tcga_bulk.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `/Users/tatopro9130/Downloads/cancer_judging_machine-master 22/public/03_Reference_Free_TME/00_inputs/TCGA-GBM.star_counts.tsv` -- read-only, in the predecessor project | 105.3 MB | `47e40f9e2364` | = source content (2026-10-01) |
| `data/processed/tcga_gbm_bulk_cpm.csv.gz` -- linearised, CPM | 18.0 MB | `8eaa60901344` | -- |
| `data/processed/tcga_gbm_bulk_provenance.json` -- build record | 1.4 KB | `b7465c7033fb` | -- |

### D05 · TCGA-LGG RNA-seq, STAR gene counts

- **Status:** used
- **Repository:** UCSC Xena, GDC hub
- **Accession:** Xena dataset TCGA-LGG.star_counts.tsv
- **Source:** <https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LGG.star_counts.tsv.gz>
- **Version:** downloaded 2026-09-17
- **Used for:** Independent second cohort (534 samples; 510 matched to methylation, 496 to ABSOLUTE) for every TCGA analysis.
- **Citation:** Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena], doi:10.1038/s41587-020-0546-8
  - Crossref 10.1038/s41587-020-0546-8: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'LGG', 'lower-grade'; a mention, not proof of use): Title; Abstract; II. Methods; III. Statistical Analysis, Data Analysis and Measurements; 4.2; 4.3; 4.4; 4.5; 4.6; 4.8; 4.9; 4.10; 4.11; 4.12; V. Discussion; Acknowledgements
- **Read by:** `scripts/build_tcga_bulk.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `TCGA_LGG/TCGA-LGG.star_counts.tsv.gz` | 59.2 MB | `23b5ca39d50a` | = source content (2026-10-01) |
| `data/processed/tcga_lgg_bulk_cpm.csv.gz` -- linearised, CPM | 53.5 MB | `eb6ae17b755e` | -- |
| `data/processed/tcga_lgg_bulk_provenance.json` -- build record. Two text slips in it, both from literals in the generator (fixed 2026-10-01 for future runs; the record is not regenerated mid-analysis): what_this_is says TCGA-GBM, and 'min nonzero 1.000000 = log2(1.5)' should read log2(2). Source path and cohort field say LGG; the inversion applied is unaffected | 1.3 KB | `2fa3dfbb4691` | -- |

### D06 · GENCODE v36 comprehensive gene annotation

- **Status:** used
- **Repository:** GENCODE (EMBL-EBI FTP)
- **Accession:** gencode.v36.annotation.gtf.gz
- **Source:** <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz>
- **Version:** release 36 -- the annotation the Xena GDC hub's STAR counts are keyed to
- **Used for:** Ensembl id -> gene symbol for both TCGA expression matrices.
- **Citation:** Frankish A, et al. Nucleic Acids Res 49:D916-D923 (2021) [GENCODE 2021], doi:10.1093/nar/gkaa1087
  - Crossref 10.1093/nar/gkaa1087: verified (year deposited as 2020 (online-first); claimed 2021)
- **Read by:** `scripts/build_tcga_bulk.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `/Users/tatopro9130/Downloads/cancer_judging_machine-master 22/public/03_Reference_Free_TME/00_inputs/gencode.v36.annotation.gtf.gz` -- read-only, in the predecessor project | 42.4 MB | `81359865bab0` | = source size (2026-10-01) |

### D07 · TCGA-GBM DNA methylation, Illumina HumanMethylation450

- **Status:** used
- **Repository:** UCSC Xena, TCGA hub
- **Accession:** Xena dataset TCGA.GBM.sampleMap/HumanMethylation450
- **Source:** <https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.GBM.sampleMap%2FHumanMethylation450.gz>
- **Version:** downloaded 2026-09-18 (stored decompressed)
- **Used for:** Per-cell-type lymphoid truth (T / NK / B) through EpiDISH -- the methylation side of the headline lymphoid comparison (cohort-level only; see D23).
- **Citation:** Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena], doi:10.1038/s41587-020-0546-8
  - Crossref 10.1038/s41587-020-0546-8: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'methylation', 'EpiDISH'; a mention, not proof of use): Title; Abstract; II. Methods; 4.3; 4.4; 4.11; 4.12; V. Discussion; Limitations; Acknowledgements; References
- **Read by:** `scripts/extract_epidish_probes.py`; `scripts/extract_gimicc_cpgs.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `TCGA_LGG/TCGA.GBM.sampleMap_HumanMethylation450` | 454.4 MB | `f17dc88457a7` | = source content (2026-10-01) |
| `data/processed/gbm_methylation_epidish_probes.csv.gz` -- the EpiDISH CpGs only | 92.1 KB | `442cab01e6bd` | -- |

### D08 · TCGA-LGG DNA methylation, Illumina HumanMethylation450

- **Status:** used
- **Repository:** UCSC Xena, TCGA hub (the same hub as GBM, deliberately)
- **Accession:** Xena dataset TCGA.LGG.sampleMap/HumanMethylation450
- **Source:** <https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LGG.sampleMap%2FHumanMethylation450.gz>
- **Version:** downloaded 2026-09-18
- **Used for:** As D07, for LGG.
- **Citation:** Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena], doi:10.1038/s41587-020-0546-8
  - Crossref 10.1038/s41587-020-0546-8: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'methylation', 'EpiDISH'; a mention, not proof of use): Title; Abstract; II. Methods; 4.3; 4.4; 4.11; 4.12; V. Discussion; Limitations; Acknowledgements; References
- **Read by:** `scripts/extract_epidish_probes.py`; `scripts/extract_gimicc_cpgs.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `TCGA_LGG/TCGA.LGG.sampleMap_HumanMethylation450.gz` | 441.8 MB | `e582907ff45b` | = source size (2026-10-01) |
| `data/processed/lgg_methylation_epidish_probes.csv.gz` -- the EpiDISH CpGs only | 289.7 KB | `110a4ad896d7` | -- |

### D09 · ABSOLUTE tumour purity and ploidy, TCGA PanCanAtlas master calls

- **Status:** used
- **Repository:** NCI GDC, PanCanAtlas publication page
- **Accession:** GDC file 4f277128-f793-4354-a13d-30cc7fe9f6b5
- **Source:** <https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5>
- **Version:** TCGA_mastercalls.abs_tables_JSedit.fixed.txt, downloaded 2026-09-11
- **Used for:** The DNA ground truth for tumour content, sharing no input with any RNA method (147 GBM / 496 LGG complete cases); ploidy, genome doublings and subclonal fraction as candidate failure factors.
- **Citation:** Carter SL, et al. Nat Biotechnol 30:413-421 (2012), doi:10.1038/nbt.2203
  - Crossref 10.1038/nbt.2203: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'ABSOLUTE', 'purity'; a mention, not proof of use): Abstract; II. Methods; 4.1; 4.3; 4.4; 4.5; 4.6; 4.9; 4.10; 4.11; 4.12; Acknowledgements
- **Read by:** `scripts/absolute_purity_yardstick.py`; `scripts/absolute_join_audit.py`; `scripts/gimicc_truth.py`; `scripts/failure_factors.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/tcga/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | 880.7 KB | `f430a975433d` | = recorded hash; = source size (2026-10-01) |
| `pipeline_packages /ABSOLUTE/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` -- the copy as first downloaded | 880.7 KB | `f430a975433d` | -- |

### D10 · MC3 somatic mutation calls, gene level, TCGA-LGG

- **Status:** used
- **Repository:** UCSC Xena, TCGA hub
- **Accession:** Xena dataset mc3_gene_level/LGG_mc3_gene_level.txt
- **Source:** <https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/mc3_gene_level%2FLGG_mc3_gene_level.txt.gz>
- **Version:** downloaded 2026-09-17 (stored decompressed)
- **Used for:** IDH1 status as a candidate biological failure factor in LGG.
- **Citation:** Ellrott K, et al. Cell Syst 6:271-281 (2018), doi:10.1016/j.cels.2018.03.002
  - Crossref 10.1016/j.cels.2018.03.002: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'IDH'; a mention, not proof of use): Title; 4.3; 4.4
- **Read by:** `scripts/failure_factors.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `TCGA_LGG/mc3_gene_level_LGG_mc3_gene_level.txt` | 39.8 MB | `a254bcff3e15` | = source content (2026-10-01) |

### D11 · Leukocyte fraction from DNA methylation, TCGA PanImmune

- **Status:** used
- **Repository:** NCI GDC, PanImmune (Immune Landscape of Cancer) publication page
- **Accession:** GDC file 6f75c9d7-5134-4ed1-b8f3-72856c98a4e8
- **Source:** <https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8>
- **Version:** TCGA_all_leuk_estimate.masked.20170107.tsv
- **Used for:** Aggregate immune truth for the pre-registered immune arm (P1/P2 over-call test).
- **Citation:** Thorsson V, et al. Immunity 48:812-830 (2018), doi:10.1016/j.immuni.2018.03.023
  - Crossref 10.1016/j.immuni.2018.03.023: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'leukocyte', 'Thorsson', 'immune arm'; a mention, not proof of use): II. Methods; 4.4; 4.11; 4.12; V. Discussion; Acknowledgements
- **Read by:** `scripts/immune_arm.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `Immune_Fraction/TCGA_all_leuk_estimate.masked.20170107.tsv` | 547.3 KB | `5a8268caedbf` | = source size (2026-10-01) |

### D12 · TCGA-GBM clinical sample annotations (Brennan 2013 study)

- **Status:** used
- **Repository:** cBioPortal for Cancer Genomics, study gbm_tcga_pub2013
- **Accession:** gbm_tcga_pub2013_clinical_sample.txt (+ _clinical_patient.txt)
- **Source:** <https://www.cbioportal.org/study/summary?id=gbm_tcga_pub2013>
- **Version:** as vendored by the predecessor project
- **Used for:** Verhaak expression subtype, IDH1 mutation and G-CIMP status as candidate failure factors in GBM.
- **Citation:** Brennan CW, et al. Cell 155:462-477 (2013), doi:10.1016/j.cell.2013.09.034
  - Crossref 10.1016/j.cell.2013.09.034: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'subtype', 'G-CIMP', 'Verhaak'; a mention, not proof of use): 4.8; References
- **Read by:** `scripts/failure_factors.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `/Users/tatopro9130/Downloads/cancer_judging_machine-master 22/public/04 Subtypes and Pathways/source/gbm_tcga_pub2013_clinical_sample.txt` -- read-only, in the predecessor project | 69.1 KB | `bae781e3a27e` | -- |
| `/Users/tatopro9130/Downloads/cancer_judging_machine-master 22/public/04 Subtypes and Pathways/source/gbm_tcga_pub2013_clinical_patient.txt` -- read-only, in the predecessor project | 48.9 KB | `82cae4a74413` | -- |

## Single-cell references

### D13 · GBmap (core) -- harmonised glioblastoma single-cell atlas

- **Status:** used
- **Repository:** CZ CELLxGENE Discover
- **Accession:** collection 999f2a15-3d7e-440b-96ae-2c806799c08c; dataset version 861acfd8-25f0-418b-a445-aa96da232827 (title 'Core GBmap', schema 7.1.0)
- **Source:** <https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c>
- **Version:** 338,564 cells x 27,632 genes; identifiers read from the file's own uns/citation
- **Used for:** The reference every method solves against: the vendored frozen signature (8-type roster, annotation_level_3) and the raw/X-count rebuild that is primary; the source of the synthetic mixtures.
- **Citation:** Ruiz-Moreno C, et al. Neuro-Oncology 27:2281-2295 (2025), doi:10.1093/neuonc/noaf113 -- the dataset itself cites the preprint: Ruiz-Moreno C, et al. bioRxiv (2022), doi:10.1101/2022.08.27.505439
  - Crossref 10.1093/neuonc/noaf113: verified (matches)
  - Crossref 10.1101/2022.08.27.505439: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'GBmap'; a mention, not proof of use): II. Methods; 4.4; 4.5; 4.7; 4.11; Limitations; Acknowledgements; References
- **Read by:** `scripts/run_all.py`; `scripts/gbmap_t_vs_b.py`; `ivygap/bench/run_benchmark.py`; `scripts/remeasure_method.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/reference/gbmap_core.h5ad` | 7.6 GB | `459444da84ef` | = source size (2026-10-01) |

### D14 · Abdelfattah et al. 2022 glioma single-cell RNA-seq (independent atlas)

- **Status:** used
- **Repository:** GEO (raw and sample identifiers); processed matrix and clusters from the Broad Single Cell Portal
- **Accession:** GEO GSE182109 (GSM5518596-GSM5518639; raw supplementary archive GSE182109_RAW); processed matrix + clusters from the Single Cell Portal, whose study accession was not recorded -- content verified against GEO's raw counts
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109>
- **Version:** Processed_matrix.mtx.gz (Seurat log-normalised, inverted to counts here), genes.tsv, barcodes.tsv, Cluster_GBM.txt
- **Used for:** The independent-atlas test of the B-over-T result: a reference not among GBmap's 16 source studies. All 201,986 barcodes carry one of the series' 44 GSM prefixes; none comes from Neftel's GSE131928. Provenance closed against GEO (scripts/verify_gse182109_raw.py): 201,893 of the 201,986 processed cells are in their own sample's Cell Ranger barcodes (the 93 others, in two samples, look like a more permissive cell call), and on the samples checked every recovered count equals GEO's raw count and every library equals the raw UMI total.
- **Citation:** Abdelfattah N, et al. Nat Commun 13:767 (2022), doi:10.1038/s41467-022-28372-y
  - Crossref 10.1038/s41467-022-28372-y: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Abdelfattah', 'independent atlas', 'GSE182109'; a mention, not proof of use): 4.4; 4.5; 4.11; Acknowledgements
- **Read by:** `scripts/build_abdelfattah_reference.py`; `scripts/atlas_t_vs_b.py`; `scripts/independent_atlas_test.py`; `scripts/verify_gse182109_raw.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `GSE182109/Processed_matrix.mtx.gz.download/Processed_matrix.mtx.gz` | 1.9 GB | `not hashed` | -- |
| `GSE182109/genes.tsv` | 635.2 KB | `99cad0d3699b` | -- |
| `GSE182109/barcodes.tsv` | 5.8 MB | `4b412821a02f` | -- |
| `GSE182109/Cluster_GBM.txt` | 18.9 MB | `1b53ec791284` | -- |
| `GSE182109/gsm_to_sample.tsv` -- GSM -> sample name, from GEO | 926 B | `54c8723ef988` | -- |
| `GSE182109/GSE182109_RAW` -- GEO's raw archive: 44 Cell Ranger outputs (5.0.0 x19, 5.0.1 x25), 132 files; sha256 of each in results/gse182109_raw_verification.json | 2.3 GB, 132 files | `not hashed` | -- |
| `data/reference/abdelfattah_2022/profile.csv` -- built reference | 2.9 MB | `aa0ff9912939` | -- |
| `data/reference/abdelfattah_2022/provenance.json` -- build record | 4.5 KB | `183d3d416604` | -- |

### D15 · Neftel et al. 2019 glioblastoma single-cell RNA-seq

- **Status:** used
- **Repository:** GEO
- **Accession:** GEO GSE131928; the Smart-seq2 matrix and cell assignments as distributed on the Broad Single Cell Portal
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928>
- **Version:** processed TPM tables and authors' metadata
- **Used for:** Reference-dependence arm (a Neftel-only reference). One of GBmap's constituent studies, so NOT an independent reference.
- **Citation:** Neftel C, et al. Cell 178:835-849 (2019), doi:10.1016/j.cell.2019.06.024
  - Crossref 10.1016/j.cell.2019.06.024: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Neftel'; a mention, not proof of use): 4.11; Acknowledgements
- **Read by:** `scripts/build_neftel_reference.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages / repos/SCDC/neftel_data/IDHwtGBM.processed.SS2.logTPM.txt` -- Smart-seq2 log-TPM, Broad Single Cell Portal copy -- the matrix the build reads (its GEO-named twin GSM3828672 was byte-identical and deleted, docs/DELETED_FILES.md) | 938.9 MB | `not hashed` | -- |
| `pipeline_packages / repos/SCDC/neftel_data/GBM_Metadata` -- authors' cell assignments (SCP) | 1.4 MB | `b119ac2d7c28` | -- |
| `pipeline_packages / repos/SCDC/neftel_data/GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv` -- 10x subset -- read by nothing | 1.5 GB | `12efdba4d2e6` | -- |
| `pipeline_packages / repos/SCDC/neftel_data/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx` -- GEO sample sheet -- read by nothing | 828.2 KB | `4bcee2ff8b09` | -- |

### D16 · Darmanis et al. 2017 glioblastoma single-cell RNA-seq (core vs periphery)

- **Status:** used
- **Repository:** GEO
- **Accession:** GEO GSE84465
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465>
- **Version:** GSE84465_GBM_All_data.csv.gz and the series matrix, fetched 2026-09-14
- **Used for:** Two uses: the only cross-patient core-vs-periphery constraint check (C1, four gates), and a Darmanis-only reference in the reference-dependence arm (a GBmap constituent, so not independent).
- **Citation:** Darmanis S, et al. Cell Rep 21:1399-1410 (2017), doi:10.1016/j.celrep.2017.10.030
  - Crossref 10.1016/j.celrep.2017.10.030: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Darmanis'; a mention, not proof of use): 4.11; Acknowledgements
- **Read by:** `scripts/build_darmanis_reference.py`; `scripts/darmanis_constraint_check.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/darmanis_2017/GSE84465_GBM_All_data.csv.gz` | 19.6 MB | `6431e67742c1` | = source size (2026-10-01) |
| `data/raw/darmanis_2017/GSE84465_series_matrix.txt.gz` | 142.1 KB | `440b07f02aee` | = source size (2026-10-01) |
| `data/raw/darmanis_2017/cell_metadata.csv` -- parsed from the series matrix | 360.3 KB | `4958456b561e` | -- |

## Constraint checks that use no deconvolution

### D17 · Mossi Albiach et al. 2023 spatially sampled glioblastoma single-cell atlas

- **Status:** used
- **Repository:** CZ CELLxGENE Discover
- **Accession:** collection 113a558a-e96e-4643-81db-140e95c58578; dataset d45b4ce6-9725-4d79-b97a-70a44158bdbf (donor SL040, schema 7.1.0)
- **Source:** <https://cellxgene.cziscience.com/collections/113a558a-e96e-4643-81db-140e95c58578>
- **Version:** 135,482 cells x 58,234 genes, one donor, 27 samples in 12 locations
- **Used for:** Tests the registered constraints against counted (not inferred) composition per anatomic zone: 3 of 4 testable constraints hold.
- **Citation:** Mossi Albiach A, et al. bioRxiv (2023), doi:10.1101/2023.09.01.555882
  - Crossref 10.1101/2023.09.01.555882: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Albiach'; a mention, not proof of use): 4.11; Acknowledgements
- **Read by:** `scripts/albiach_constraint_check.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages / repos/SCDC/albiach_data/d45b4ce6-9725-4d79-b97a-70a44158bdbf.h5ad` | 1.4 GB | `e68dc89ef023` | = source size (2026-10-01) |

### D18 · Ivy GAP in-situ hybridisation (ISH) by anatomic structure

- **Status:** used
- **Repository:** Allen Institute Ivy GAP portal: the per-(gene x sub-block) ISH table in D02
- **Accession:** gene_expression_details.csv (/api/v2/gbm/ namespace)
- **Source:** <https://glioblastoma.alleninstitute.org/api/v2/gbm/gene_expression_details.csv>
- **Version:** the portal table, sha256 recorded at download; the API-based route (scripts/fetch_ivygap_ish.py) is SUPERSEDED and its output is not used
- **Used for:** A constraint check with no deconvolution anywhere in the chain: CD44 / BIRC5 ISH support C1 and C7.
- **Citation:** Puchalski RB, et al. Science 360:660-663 (2018), doi:10.1126/science.aaf2666
  - Crossref 10.1126/science.aaf2666: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'ISH', 'in-situ'; a mention, not proof of use): Title; Abstract; II. Methods; III. Statistical Analysis, Data Analysis and Measurements; 4.1; 4.3; 4.4; 4.7; 4.9; 4.10; 4.11; 4.12; V. Discussion; Limitations; References
- **Read by:** `scripts/ish_constraint_check.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/ivygap/gene_expression_details.csv` -- listed under D02 | 4.8 MB | `07c5702a2549` | -- |

### D32 · Flow-cytometry immune composition of brain tumours by IDH status (Klemm et al. 2020, Figure 1F)

- **Status:** used
- **Repository:** Cell (Elsevier); article PDF and supplements placed locally by the user (gitignored)
- **Accession:** doi:10.1016/j.cell.2020.05.007 -- Figure 1F (PDF page 3); Table S1 (mmc1.pdf, cohort); Table S2 (mmc2.pdf, gating definitions)
- **Source:** <https://doi.org/10.1016/j.cell.2020.05.007>
- **Version:** article PDF sha256 df824704317a2f84ec9a49e9d8ebff648789cd59b95ce6800f93f58f612b6503
- **Used for:** Direct check of the lymphoid truth by IDH status (prespecified/klemm_t_vs_b.md): cohort-mean T, B and NK as % of CD45+ in 17 IDH-mutant and 40 IDH-wildtype gliomas, measured from the figure's vector geometry and validated against two numbers printed in the paper.
- **Citation:** Klemm F, et al. Cell 181:1643-1660 (2020), doi:10.1016/j.cell.2020.05.007
  - Crossref 10.1016/j.cell.2020.05.007: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'Klemm', 'flow cytometry'; a mention, not proof of use): Title; I. Introduction; 4.4; 4.12
- **Read by:** `scripts/klemm_figure1f.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `celldecov_reference_papers/Klemm_et_al_2020/PIIS0092867420305699.pdf` -- the article PDF (publisher copyright; not redistributed) | 20.3 MB | `df824704317a` | -- |
| `celldecov_reference_papers/Klemm_et_al_2020/mmc1.pdf` -- Table S1, the flow-cytometry cohort | 62.9 KB | `85bc295ab237` | -- |
| `celldecov_reference_papers/Klemm_et_al_2020/mmc2.pdf` -- Table S2, gating definitions | 30.5 KB | `d636cae87a4e` | -- |

## Marker sets and reference data bundled with software

### D19 · GBMDeconvoluteR marker sets (Ajaib immune, Ruiz-Moreno level-3 immune, Neftel four-state) and Ajaib's published imaging-mass-cytometry correlations

- **Status:** used
- **Repository:** GBMDeconvoluteR source repository (vendored); Ajaib et al. 2023 paper
- **Accession:** Ajaib_et_al_2022_GBM_Immune_markers.rds, Moreno_et_al_2022_lvl3_immune_markers.rds, Neftel_et_al_2019_four_state_neoplastic_markers.rds
- **Source:** <https://github.com/GliomaGenomics/GBMDeconvoluteR>
- **Version:** vendored copy in pipeline_packages / repos/GBMDeconvoluteR
- **Used for:** The pre-registered IMC protein test (an external instrument for ACS) and the mesenchymal-program score.
- **Citation:** Ajaib S, et al. Neuro-Oncology 25:1236-1248 (2023), doi:10.1093/neuonc/noad021
  - Crossref 10.1093/neuonc/noad021: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'IMC', 'MCP', 'mesenchymal'; a mention, not proof of use): 4.8; 4.11; V. Discussion; Limitations
- **Read by:** `scripts/imc_anchored_test.py`; `scripts/mes_score.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages / repos/GBMDeconvoluteR` -- vendored repository | 52.0 MB, 76 files | `not hashed` | -- |

### D20 · MCPcounter default marker genes

- **Status:** used
- **Repository:** MCPcounter source repository
- **Accession:** Signatures/genes.txt
- **Source:** <https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt>
- **Version:** fetched 2026-09-15
- **Used for:** The default-marker arm of the IMC test.
- **Citation:** Becht E, et al. Genome Biol 17:218 (2016), doi:10.1186/s13059-016-1070-5
  - Crossref 10.1186/s13059-016-1070-5: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'MCP'; a mention, not proof of use): 4.11
- **Read by:** `scripts/imc_anchored_test.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `data/raw/mcp_default_genes.txt` | 5.2 KB | `408f6c5d02c8` | = source size (2026-10-01) |

### D21 · EpiDISH blood reference centDHSbloodDMC.m (333 CpGs x 7 blood cell types)

- **Status:** used
- **Repository:** Bioconductor package EpiDISH (data object)
- **Accession:** EpiDISH::centDHSbloodDMC.m
- **Source:** <https://bioconductor.org/packages/EpiDISH>
- **Version:** EpiDISH 2.28.0 as installed (packageVersion, read 2026-10-01; the methylation artefacts do not record it)
- **Used for:** Turns D07/D08 into per-cell-type lymphoid fractions (RPC mode); a blood reference applied to brain tumour, so used for lymphoid sub-composition only.
- **Citation:** Teschendorff AE, et al. BMC Bioinformatics 18:105 (2017), doi:10.1186/s12859-017-1511-5
  - Crossref 10.1186/s12859-017-1511-5: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'EpiDISH'; a mention, not proof of use): Title; II. Methods; 4.4; 4.12; V. Discussion; References
- **Read by:** `scripts/methylation_celltypes.py`; `scripts/extract_epidish_probes.py`

### D22 · BayesPrism's bundled gene annotation (gene groups and gene types, GENCODE v22)

- **Status:** used
- **Repository:** BayesPrism source repository (vendored, Danko lab)
- **Accession:** cleanup.genes() gene groups; select.gene.type() 'protein_coding'
- **Source:** <https://github.com/Danko-Lab/BayesPrism>
- **Version:** BayesPrism 2.2.3 as installed (packageVersion, read 2026-10-01); source vendored in pipeline_packages / repos/BayesPrism-main
- **Used for:** The authors' own gene filtering in the authors'-workflow BayesPrism runs on TCGA (post-registration extension).
- **Citation:** Chu T, et al. Nat Cancer 3:505-517 (2022), doi:10.1038/s43018-022-00356-3
  - Crossref 10.1038/s43018-022-00356-3: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'BayesPrism'; a mention, not proof of use): 4.5; 4.7; 4.10; 4.12; V. Discussion
- **Read by:** `R/run_bayesprism_authors.R`; `R/run_bayesprism_authors.R`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages / repos/BayesPrism-main` -- vendored repository | 77.4 MB, 43 files | `not hashed` | -- |

### D31 · GIMiCC reference libraries (glioma-specific hierarchical DNA-methylation deconvolution)

- **Status:** used
- **Repository:** Bioconductor ExperimentHub; github.com/SalasLab/GIMiCC
- **Accession:** EH9483 (GIMiCC_Library.rda, 15 layer libraries); EH9482 (Capper_example_betas.rda, the reproduction control)
- **Source:** <https://github.com/SalasLab/GIMiCC>
- **Version:** GIMiCC 0.99.1 @26cb8a15 (the user's copy, byte-identical); ExperimentHub cache sha256 EH9483 bdfbf2849ed6513c3afa84560f24b4b49d74245e20cd2f2a15ea59395a0bdb14, EH9482 47da7b0f14925d16fe17ac2e3fc975083e663142e1026ec0d00ab7dcc030e140 (both added 2024-04-02)
- **Used for:** The second methylation truth for the lymphoid ordering (post-registration, prespecified/gimicc_truth_confirmation.md): GIMiCC run on the TCGA 450K betas (D07/D08), restricted to its 4,022 library CpGs by scripts/extract_gimicc_cpgs.py.
- **Citation:** Pike SC, et al. Acta Neuropathol Commun 12:170 (2024), doi:10.1186/s40478-024-01874-0
  - Crossref 10.1186/s40478-024-01874-0: verified (matches)
- **Manuscript sections mentioning it** (keyword scan for 'GIMiCC'; a mention, not proof of use): Title; 4.4; 4.12; V. Discussion; Limitations
- **Read by:** `R/run_gimicc.R`; `scripts/extract_gimicc_cpgs.py`; `scripts/gimicc_truth.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages / repos/GIMiCC-main` -- the package source (the user's copy) | 62.0 KB, 17 files | `not hashed` | -- |

## Derived by this project or its predecessor (not public data)

### D23 · Frozen GBM signature matrix and cell-size factors (from GBmap, by the predecessor project) and its synthetic held-out benchmark

- **Status:** derived
- **Repository:** this repository, reference_frozen/ (vendored unchanged)
- **Accession:** signature_matrix.tsv, cell_size_factors.csv, tcga_benchmark/*
- **Version:** frozen 2026-07-12T03:36:26Z; sha256 recorded in PROVENANCE.json
- **Used for:** The registered signature; yardstick 1 (500 pseudobulk mixtures, NNLS and SVR only), which shares GBmap with the ACS arm.
- **Manuscript sections mentioning it** (keyword scan for 'frozen signature'; a mention, not proof of use): II. Methods; 4.4; 4.10; V. Discussion
- **Read by:** `ivygap/config.py`

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `reference_frozen/signature_matrix.tsv` | 2.4 MB | `57d0fe572e11` | = recorded hash |
| `reference_frozen/cell_size_factors.csv` | 1.2 KB | `8d1621fd5b7c` | = recorded hash |
| `reference_frozen/tcga_benchmark/benchmark_summary.csv` | 2.7 KB | `961b3577aba7` | = recorded hash |

## Present on disk but read by nothing

### D24 · CIBERSORT relative fractions for TCGA (PanImmune)

- **Status:** on disk, unused
- **Repository:** NCI GDC, PanImmune publication page
- **Accession:** GDC file b3df502e-3594-46ef-9f94-d041a20a0b9a
- **Source:** <https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a>
- **Version:** TCGA.Kallisto.fullIDs.cibersort.relative.tsv
- **Used for:** Obtained for a comparison against a published deconvolution of the same cohort; no script reads it. DATA_SOURCES.md A2 lists it as a supply -- that comparison was never built.
- **Citation:** Thorsson V, et al. Immunity 48:812-830 (2018), doi:10.1016/j.immuni.2018.03.023
  - Crossref 10.1016/j.immuni.2018.03.023: verified (matches)

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `Immune_Fraction/TCGA.Kallisto.fullIDs.cibersort.relative.tsv` | 3.6 MB | `8bf3d2cb511d` | = source size (2026-10-01) |

### D25 · ABSOLUTE segment tables, TCGA PanCanAtlas

- **Status:** on disk, unused
- **Repository:** NCI GDC, PanCanAtlas publication page
- **Accession:** GDC file 0f4f5701-7b61-41ae-bda9-2805d1ca9781
- **Source:** <https://api.gdc.cancer.gov/data/0f4f5701-7b61-41ae-bda9-2805d1ca9781>
- **Version:** TCGA_mastercalls.abs_segtabs.fixed.txt
- **Used for:** Downloaded with D09; nothing reads it.
- **Citation:** Carter SL, et al. Nat Biotechnol 30:413-421 (2012), doi:10.1038/nbt.2203
  - Crossref 10.1038/nbt.2203: verified (matches)

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages /ABSOLUTE/TCGA_mastercalls.abs_segtabs.fixed.txt` | 241.3 MB | `04080d770547` | = source size (2026-10-01) |

### D26 · TCGA-LGG GDC clinical and biospecimen tables; GDC whole-exome BAM manifest

- **Status:** on disk, unused
- **Repository:** NCI GDC Data Portal
- **Accession:** clinical.project-tcga-lgg.2026-09-17, biospecimen.project-tcga-lgg.2026-09-17, gdc_manifest.2026-09-17.191037.txt
- **Source:** <https://portal.gdc.cancer.gov/projects/TCGA-LGG>
- **Version:** exported 2026-09-17
- **Used for:** Nothing reads the tables. The manifest lists controlled-access BAMs that were never downloaded and are not needed.
- **Citation:** Grossman RL, et al. N Engl J Med 375:1109-1112 (2016) [GDC], doi:10.1056/NEJMp1607591
  - Crossref 10.1056/NEJMp1607591: verified (matches)

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `TCGA_LGG/project/clinical.project-tcga-lgg.2026-09-17/clinical.tsv` | 4.2 MB | `a2e0028fd5af` | -- |
| `TCGA_LGG/project/gdc_manifest.2026-09-17.191037.txt` | 5.0 MB | `6dab6669d4ec` | -- |

### D27 · GEO GSE107559 supplementary files (the Ivy GAP portal tables, mirrored)

- **Status:** duplicate
- **Repository:** GEO
- **Accession:** GSE107559_ivygap_rows-genes.csv, GSE107559_ivygap_columns-samples.xlsx
- **Source:** <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107559>
- **Used for:** A GEO mirror of D01's tables; nothing reads this copy.
- **Citation:** Puchalski RB, et al. Science 360:660-663 (2018), doi:10.1126/science.aaf2666
  - Crossref 10.1126/science.aaf2666: verified (matches)

| local file | size | sha256 (first 12) | check |
|---|---|---|---|
| `pipeline_packages /GSE107559/GSE107559_ivygap_rows-genes.csv` | 1.6 MB | `468d0eee470f` | -- |
| `pipeline_packages /GSE107559/GSE107559_ivygap_columns-samples.xlsx` | 53.1 KB | `e384afa62521` | -- |

## Identified and assessed, not used

### D28 · Human Brain Cell Atlas v1.0 (Siletti et al.)

- **Status:** considered
- **Repository:** CZ CELLxGENE Discover
- **Accession:** collection 283d65eb-dd53-496d-adb7-7570c7caa443
- **Source:** <https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443>
- **Used for:** A normal-brain reference for a second tissue; not obtained.
- **Citation:** Siletti K, et al. Science 382:eadd7046 (2023), doi:10.1126/science.add7046
  - Crossref 10.1126/science.add7046: verified (matches)

### D29 · Allen Human Brain Atlas, bulk expression by region

- **Status:** considered
- **Repository:** Allen Institute for Brain Science
- **Accession:** human.brain-map.org downloads
- **Source:** <https://human.brain-map.org/static/download>
- **Used for:** A second tissue for the anatomy test; fetch script exists, the constraint file, roster and platform are the blockers.
- **Citation:** Hawrylycz MJ, et al. Nature 489:391-399 (2012), doi:10.1038/nature11405
  - Crossref 10.1038/nature11405: verified (matches)

### D30 · Spatially resolved multi-omics of glioblastoma (Ravi et al.)

- **Status:** considered
- **Repository:** publication data
- **Used for:** Candidate second GBM cohort with spatial structure; not obtained.
- **Citation:** Ravi VM, et al. Cancer Cell 40:639-655 (2022), doi:10.1016/j.ccell.2022.05.009
  - Crossref 10.1016/j.ccell.2022.05.009: verified (matches)

## Negative control for the citation check

| DOI | claimed first author | rejected? | why |
|---|---|---|---|
| 10.1016/j.ccell.2022.05.009 | Ruiz-Moreno | yes | first author is 'Ravi', not 'Ruiz-Moreno' |

## Data availability statement (draft for the manuscript)

All data analysed in this study are publicly available. Bulk RNA-seq of laser-microdissected glioblastoma structures is from the Ivy Glioblastoma Atlas Project (Allen Institute for Brain Science; release 2014-11-25; GEO GSE107559). TCGA-GBM and TCGA-LGG RNA-seq (STAR counts) and HumanMethylation450 data were obtained from the UCSC Xena GDC and TCGA hubs; ABSOLUTE purity calls, leukocyte fractions and MC3 mutation calls from the TCGA PanCanAtlas resources at the NCI Genomic Data Commons and UCSC Xena; TCGA-GBM clinical annotations from cBioPortal (gbm_tcga_pub2013). Single-cell references are GBmap (CZ CELLxGENE Discover collection 999f2a15-3d7e-440b-96ae-2c806799c08c, dataset 861acfd8-25f0-418b-a445-aa96da232827), Abdelfattah et al. (GEO GSE182109), Neftel et al. (GEO GSE131928), Darmanis et al. (GEO GSE84465) and Mossi Albiach et al. (CZ CELLxGENE Discover collection 113a558a-e96e-4643-81db-140e95c58578). Exact files, versions, sizes and checksums are listed in Supplementary Table S7.
