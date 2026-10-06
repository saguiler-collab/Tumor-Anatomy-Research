#!/usr/bin/env python3
"""
build_data_inventory.py -- every public dataset this study uses, measured, not recalled.

The reader of a paper is owed an exact answer to "what data, from where, which version,
used for what". This script produces it from the files themselves:

  * every local file: present?, bytes, sha256 (cached by path + size + mtime);
  * every claim that a script reads a dataset: the script is opened and the token that
    reads it must be found -- a stale path or a renamed folder shows up here, not in a
    reviewer's attempt to reproduce;
  * every citation: resolved against Crossref (`--verify-dois`, cached) and its first
    author and year compared with what this file claims, so a citation that pairs one
    paper's authors with another paper's DOI fails;
  * every public source: re-checked against the server (`--verify-remote`, cached):
    byte size, or for a file stored here decompressed, the sha256 of the decompressed
    stream on both sides;
  * where each dataset's results appear in docs/MANUSCRIPT.md, found by scanning it.

Status vocabulary, used exactly:
  used                -- read by code whose results are reported
  derived             -- built by this project or its predecessor from a dataset above
  on disk, unused     -- present, read by nothing; listed so nobody assumes otherwise
  duplicate           -- a second copy of a file listed elsewhere
  considered          -- identified and assessed, not obtained or not used

    python3 scripts/build_data_inventory.py                    # offline, uses the cache
    python3 scripts/build_data_inventory.py --verify-dois      # + Crossref
    python3 scripts/build_data_inventory.py --verify-remote    # + source servers
    python3 scripts/build_data_inventory.py --hash-large       # + sha256 of files > 600 MB

Writes docs/DATA_INVENTORY.md, docs/supplementary/data_inventory_full.csv and the cache
docs/supplementary/data_inventory_cache.json.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import re
import sys
import unicodedata
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

CACHE = ROOT / "docs" / "supplementary" / "data_inventory_cache.json"
OUT_MD = ROOT / "docs" / "DATA_INVENTORY.md"
OUT_CSV = ROOT / "docs" / "supplementary" / "data_inventory_full.csv"   # S7 is built from it
MANUSCRIPT = ROOT / "docs" / "MANUSCRIPT.md"
LARGE = 600 * 1024 ** 2
UA = "anatomy-test-data-inventory/1.0 (research provenance check)"

try:   # the python.org macOS build ships no CA bundle of its own; certifi's is the fix
    import certifi
    import ssl
    _CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:  # pragma: no cover
    _CTX = None

R, P, REF = config.RAW_DIR, config.PROCESSED_DIR, config.REFERENCE_DIR
TCGA, IMM, GEO182 = config.TCGA_DOWNLOADS_DIR, config.IMMUNE_FRACTION_DIR, config.GSE182109_DIR
VEN, REPOS, PRED = config.VENDORED_DIR, config.VENDORED_REPOS_DIR, config.PREDECESSOR_ROOT
XENA_GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
XENA_TCGA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
GDC_DATA = "https://api.gdc.cancer.gov/data/"


def F(path, note="", remote=None, compare="size", recorded_sha=None, status=None):
    """One file. `compare`: size | gunzip_remote (remote is .gz of this file) |
    gunzip_both (both gzip; compare decompressed) | none."""
    return {"path": Path(path), "note": note, "remote": remote, "compare": compare,
            "recorded_sha": recorded_sha, "status": status}


# =============================================================================
# THE REGISTRY. Order is the order a reader needs: the cohort the claims rest on,
# the truths it is checked against, the references every method solves against,
# the checks that use no deconvolution, then what is bundled, derived, unused.
# =============================================================================

GROUPS = [
    "Primary cohort -- the anatomy arm",
    "Validation cohorts and orthogonal truths -- TCGA",
    "Single-cell references",
    "Constraint checks that use no deconvolution",
    "Marker sets and reference data bundled with software",
    "Derived by this project or its predecessor (not public data)",
    "Present on disk but read by nothing",
    "Identified and assessed, not used",
]

REGISTRY: list[dict] = [
    # ---------------------------------------------------------------- primary
    dict(id="D01", group=0, status="used",
         name="Ivy Glioblastoma Atlas Project (Ivy GAP) RNA-seq, normalised FPKM matrix",
         repository="Allen Institute for Brain Science, Ivy GAP portal; mirrored at GEO",
         accession="GEO GSE107559; portal well-known file 305873915",
         url="https://glioblastoma.alleninstitute.org/static/download.html",
         version="release 2014-11-25 (gene_expression_matrix_2014-11-25.zip)",
         role="The tissue every anatomic claim is scored on: 270 laser-microdissected samples "
              "(122 H&E-selected anatomic samples from 10 tumours across LE, IT, CT, MVP, PAN; "
              "148 ISH-selected cluster samples, never scored as anatomic).",
         citation="Puchalski RB, et al. Science 360:660-663 (2018)",
         doi="10.1126/science.aaf2666", first_author="Puchalski", year=2018,
         files=[F(config.IVYGAP_ZIP_PATH, "the archive as published",
                  remote="http://api.brain-map.org/api/v2/well_known_file_download/305873915"),
                F(config.IVYGAP_FPKM_PATH, "extracted"),
                F(config.IVYGAP_GENES_PATH, "extracted"),
                F(config.IVYGAP_SAMPLES_PATH, "extracted"),
                F(config.BULK_EXPRESSION_PATH, "processed: genes x samples", status="derived"),
                F(config.SAMPLE_MANIFEST_PATH, "processed: one row per sample", status="derived")],
         consumers=[("ivygap/data/download_ivygap.py", "305873915"),
                    ("ivygap/data/load_ivygap.py", "IVYGAP_FPKM_PATH")],
         keywords=["Ivy GAP", "anatomic sample"]),
    dict(id="D02", group=0, status="used",
         name="Ivy GAP portal metadata: tumour, RNA-seq sample, sub-block and ISH-gene tables",
         repository="Allen Institute Ivy GAP portal, /api/v2/gbm/ namespace",
         accession="tumor_details.csv, rna_seq_samples_details.csv, sub_block_details.csv, "
                   "gene_expression_details.csv",
         url="https://glioblastoma.alleninstitute.org/api/v2/gbm/tumor_details.csv",
         version="fetched 2026-09-03; sha256 recorded at download (CLINICAL_PROVENANCE.json)",
         role="Tumour and structure assignment, the donor-tumour join, the clinical table "
              "(survival is BLOCKED: no event indicator), and the per-sub-block ISH values "
              "used by the no-deconvolution constraint check.",
         citation="Puchalski RB, et al. Science 360:660-663 (2018)",
         doi="10.1126/science.aaf2666", first_author="Puchalski", year=2018,
         files=[F(R / "ivygap" / n, recorded_sha=("json", R / "ivygap" / "CLINICAL_PROVENANCE.json",
                                                   [n, "sha256"]))
                for n in ("tumor_details.csv", "rna_seq_samples_details.csv",
                          "sub_block_details.csv", "gene_expression_details.csv")]
               + [F(R / "ivygap" / "CLINICAL_PROVENANCE.json", "download record")],
         consumers=[("ivygap/data/portal_metadata.py", "rna_seq_samples_details"),
                    ("ivygap/data/clinical.py", "tumor_details"),
                    ("scripts/fetch_ivygap_ish.py", "gene_expression_details")],
         keywords=["survival", "ISH"]),
    dict(id="D03", group=0, status="used",
         name="Ivy GAP per-sample RSEM gene results (expected counts)",
         repository="Allen Institute Ivy GAP portal, per-sample well-known files",
         accession="270 x *.genes.results (URLs and sha256 per file in provenance.json)",
         url="https://glioblastoma.alleninstitute.org/static/download.html",
         version="fetched 2026-09-15; GEO GSE107559 states raw data are not provided there",
         role="Fractional expected counts for the count-based arms (CDSeq reference-free ACS, "
              "count-model inputs) and closure of the published matrix's provenance chain.",
         citation="Puchalski RB, et al. Science 360:660-663 (2018)",
         doi="10.1126/science.aaf2666", first_author="Puchalski", year=2018,
         files=[F(R / "ivygap_counts", "270 per-sample files; each checked against the "
                  "sha256 recorded when it was fetched", compare="per_file_provenance"),
                F(config.PROJECT_ROOT / "FPKM" / "raw_fpkm.csv", "portal manifest the URLs came from"),
                F(config.PROJECT_ROOT / "FPKM" / "bam_manifest.csv",
                  "portal BAM manifest -- BAMs deliberately not fetched"),
                F(P / "ivygap_expected_counts.csv.gz", "assembled matrix", status="derived")],
         consumers=[("scripts/fetch_ivygap_counts.py", "genes.results"),
                    ("scripts/build_ivygap_counts_matrix.py", "genes.results"),
                    ("scripts/cdseq_anatomic.py", "ivygap_expected_counts")],
         keywords=["CDSeq"]),
    # ---------------------------------------------------------------- TCGA
    dict(id="D04", group=1, status="used",
         name="TCGA-GBM RNA-seq, STAR gene counts",
         repository="UCSC Xena, GDC hub (data from the NCI Genomic Data Commons)",
         accession="Xena dataset TCGA-GBM.star_counts.tsv",
         url=XENA_GDC + "TCGA-GBM.star_counts.tsv.gz",
         version="as vendored by the predecessor project (stored decompressed; delivered as "
                 "log2(x+1) despite the name -- linearised here)",
         role="Validation cohort: per-sample Tumor fraction against ABSOLUTE purity (154 "
              "samples), the lymphoid ordering against methylation, model fit, extension panel.",
         citation="Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena]",
         doi="10.1038/s41587-020-0546-8", first_author="Goldman", year=2020,
         files=[F(PRED / "public/03_Reference_Free_TME/00_inputs/TCGA-GBM.star_counts.tsv",
                  "read-only, in the predecessor project",
                  remote=XENA_GDC + "TCGA-GBM.star_counts.tsv.gz", compare="gunzip_remote"),
                F(P / "tcga_gbm_bulk_cpm.csv.gz", "linearised, CPM", status="derived"),
                F(P / "tcga_gbm_bulk_provenance.json", "build record", status="derived")],
         consumers=[("scripts/build_tcga_bulk.py", "star_counts")],
         keywords=["TCGA-GBM", "glioblastoma cohort"]),
    dict(id="D05", group=1, status="used",
         name="TCGA-LGG RNA-seq, STAR gene counts",
         repository="UCSC Xena, GDC hub",
         accession="Xena dataset TCGA-LGG.star_counts.tsv",
         url=XENA_GDC + "TCGA-LGG.star_counts.tsv.gz",
         version="downloaded 2026-09-17",
         role="Independent second cohort (534 samples; 510 matched to methylation, 496 to "
              "ABSOLUTE) for every TCGA analysis.",
         citation="Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena]",
         doi="10.1038/s41587-020-0546-8", first_author="Goldman", year=2020,
         files=[F(TCGA / "TCGA-LGG.star_counts.tsv.gz", remote=XENA_GDC + "TCGA-LGG.star_counts.tsv.gz",
                  compare="gunzip_both"),
                F(P / "tcga_lgg_bulk_cpm.csv.gz", "linearised, CPM", status="derived"),
                F(P / "tcga_lgg_bulk_provenance.json", "build record. Two text slips in it, "
                  "both from literals in the generator (fixed 2026-10-01 for future runs; the "
                  "record is not regenerated mid-analysis): what_this_is says TCGA-GBM, and "
                  "'min nonzero 1.000000 = log2(1.5)' should read log2(2). Source path and "
                  "cohort field say LGG; the inversion applied is unaffected",
                  status="derived")],
         consumers=[("scripts/build_tcga_bulk.py", "star_counts")],
         keywords=["LGG", "lower-grade"]),
    dict(id="D06", group=1, status="used",
         name="GENCODE v36 comprehensive gene annotation",
         repository="GENCODE (EMBL-EBI FTP)",
         accession="gencode.v36.annotation.gtf.gz",
         url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/"
             "gencode.v36.annotation.gtf.gz",
         version="release 36 -- the annotation the Xena GDC hub's STAR counts are keyed to",
         role="Ensembl id -> gene symbol for both TCGA expression matrices.",
         citation="Frankish A, et al. Nucleic Acids Res 49:D916-D923 (2021) [GENCODE 2021]",
         doi="10.1093/nar/gkaa1087", first_author="Frankish", year=2021,
         files=[F(PRED / "public/03_Reference_Free_TME/00_inputs/gencode.v36.annotation.gtf.gz",
                  "read-only, in the predecessor project",
                  remote="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/"
                         "gencode.v36.annotation.gtf.gz")],
         consumers=[("scripts/build_tcga_bulk.py", "gencode.v36")],
         keywords=[]),
    dict(id="D07", group=1, status="used",
         name="TCGA-GBM DNA methylation, Illumina HumanMethylation450",
         repository="UCSC Xena, TCGA hub",
         accession="Xena dataset TCGA.GBM.sampleMap/HumanMethylation450",
         url=XENA_TCGA + "TCGA.GBM.sampleMap%2FHumanMethylation450.gz",
         version="downloaded 2026-09-18 (stored decompressed)",
         role="Per-cell-type lymphoid truth (T / NK / B) through EpiDISH -- the methylation "
              "side of the headline lymphoid comparison (cohort-level only; see D23).",
         citation="Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena]",
         doi="10.1038/s41587-020-0546-8", first_author="Goldman", year=2020,
         files=[F(TCGA / "TCGA.GBM.sampleMap_HumanMethylation450",
                  remote=XENA_TCGA + "TCGA.GBM.sampleMap%2FHumanMethylation450.gz",
                  compare="gunzip_remote"),
                F(P / "gbm_methylation_epidish_probes.csv.gz", "the EpiDISH CpGs only",
                  status="derived")],
         consumers=[("scripts/extract_epidish_probes.py", "HumanMethylation450"),
                    ("scripts/extract_gimicc_cpgs.py", "TCGA_HM450")],
         keywords=["methylation", "EpiDISH"]),
    dict(id="D08", group=1, status="used",
         name="TCGA-LGG DNA methylation, Illumina HumanMethylation450",
         repository="UCSC Xena, TCGA hub (the same hub as GBM, deliberately)",
         accession="Xena dataset TCGA.LGG.sampleMap/HumanMethylation450",
         url=XENA_TCGA + "TCGA.LGG.sampleMap%2FHumanMethylation450.gz",
         version="downloaded 2026-09-18",
         role="As D07, for LGG.",
         citation="Goldman MJ, et al. Nat Biotechnol 38:675-678 (2020) [Xena]",
         doi="10.1038/s41587-020-0546-8", first_author="Goldman", year=2020,
         files=[F(TCGA / "TCGA.LGG.sampleMap_HumanMethylation450.gz",
                  remote=XENA_TCGA + "TCGA.LGG.sampleMap%2FHumanMethylation450.gz"),
                F(P / "lgg_methylation_epidish_probes.csv.gz", "the EpiDISH CpGs only",
                  status="derived")],
         consumers=[("scripts/extract_epidish_probes.py", "HumanMethylation450"),
                    ("scripts/extract_gimicc_cpgs.py", "TCGA_HM450")],
         keywords=["methylation", "EpiDISH"]),
    dict(id="D09", group=1, status="used",
         name="ABSOLUTE tumour purity and ploidy, TCGA PanCanAtlas master calls",
         repository="NCI GDC, PanCanAtlas publication page",
         accession="GDC file 4f277128-f793-4354-a13d-30cc7fe9f6b5",
         url=GDC_DATA + "4f277128-f793-4354-a13d-30cc7fe9f6b5",
         version="TCGA_mastercalls.abs_tables_JSedit.fixed.txt, downloaded 2026-09-11",
         role="The DNA ground truth for tumour content, sharing no input with any RNA method "
              "(147 GBM / 496 LGG complete cases); ploidy, genome doublings and subclonal "
              "fraction as candidate failure factors.",
         citation="Carter SL, et al. Nat Biotechnol 30:413-421 (2012)",
         doi="10.1038/nbt.2203", first_author="Carter", year=2012,
         files=[F(R / "tcga" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                  remote=GDC_DATA + "4f277128-f793-4354-a13d-30cc7fe9f6b5",
                  recorded_sha=("json", config.FROZEN_REFERENCE_DIR / "tcga_benchmark" /
                                "ABSOLUTE_PROVENANCE.json", ["sha256"])),
                F(VEN / "ABSOLUTE" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                  "the copy as first downloaded", status="duplicate")],
         consumers=[("scripts/absolute_purity_yardstick.py", "mastercalls"),
                    ("scripts/absolute_join_audit.py", "ABS_T"), ("scripts/gimicc_truth.py", "ABS_T"),
                    ("scripts/failure_factors.py", "mastercalls")],
         keywords=["ABSOLUTE", "purity"]),
    dict(id="D10", group=1, status="used",
         name="MC3 somatic mutation calls, gene level, TCGA-LGG",
         repository="UCSC Xena, TCGA hub",
         accession="Xena dataset mc3_gene_level/LGG_mc3_gene_level.txt",
         url=XENA_TCGA + "mc3_gene_level%2FLGG_mc3_gene_level.txt.gz",
         version="downloaded 2026-09-17 (stored decompressed)",
         role="IDH1 status as a candidate biological failure factor in LGG.",
         citation="Ellrott K, et al. Cell Syst 6:271-281 (2018)",
         doi="10.1016/j.cels.2018.03.002", first_author="Ellrott", year=2018,
         files=[F(TCGA / "mc3_gene_level_LGG_mc3_gene_level.txt",
                  remote=XENA_TCGA + "mc3_gene_level%2FLGG_mc3_gene_level.txt.gz",
                  compare="gunzip_remote")],
         consumers=[("scripts/failure_factors.py", "mc3_gene_level")],
         keywords=["IDH"]),
    dict(id="D11", group=1, status="used",
         name="Leukocyte fraction from DNA methylation, TCGA PanImmune",
         repository="NCI GDC, PanImmune (Immune Landscape of Cancer) publication page",
         accession="GDC file 6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
         url=GDC_DATA + "6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
         version="TCGA_all_leuk_estimate.masked.20170107.tsv",
         role="Aggregate immune truth for the pre-registered immune arm (P1/P2 over-call test).",
         citation="Thorsson V, et al. Immunity 48:812-830 (2018)",
         doi="10.1016/j.immuni.2018.03.023", first_author="Thorsson", year=2018,
         files=[F(IMM / "TCGA_all_leuk_estimate.masked.20170107.tsv",
                  remote=GDC_DATA + "6f75c9d7-5134-4ed1-b8f3-72856c98a4e8")],
         consumers=[("scripts/immune_arm.py", "TCGA_all_leuk_estimate")],
         keywords=["leukocyte", "Thorsson", "immune arm"]),
    dict(id="D12", group=1, status="used",
         name="TCGA-GBM clinical sample annotations (Brennan 2013 study)",
         repository="cBioPortal for Cancer Genomics, study gbm_tcga_pub2013",
         accession="gbm_tcga_pub2013_clinical_sample.txt (+ _clinical_patient.txt)",
         url="https://www.cbioportal.org/study/summary?id=gbm_tcga_pub2013",
         version="as vendored by the predecessor project",
         role="Verhaak expression subtype, IDH1 mutation and G-CIMP status as candidate "
              "failure factors in GBM.",
         citation="Brennan CW, et al. Cell 155:462-477 (2013)",
         doi="10.1016/j.cell.2013.09.034", first_author="Brennan", year=2013,
         files=[F(PRED / "public/04 Subtypes and Pathways/source/gbm_tcga_pub2013_clinical_sample.txt",
                  "read-only, in the predecessor project", compare="none"),
                F(PRED / "public/04 Subtypes and Pathways/source/gbm_tcga_pub2013_clinical_patient.txt",
                  "read-only, in the predecessor project", compare="none")],
         consumers=[("scripts/failure_factors.py", "gbm_tcga_pub2013")],
         keywords=["subtype", "G-CIMP", "Verhaak"]),
    # ---------------------------------------------------------------- single-cell
    dict(id="D13", group=2, status="used",
         name="GBmap (core) -- harmonised glioblastoma single-cell atlas",
         repository="CZ CELLxGENE Discover",
         accession="collection 999f2a15-3d7e-440b-96ae-2c806799c08c; dataset version "
                   "861acfd8-25f0-418b-a445-aa96da232827 (title 'Core GBmap', schema 7.1.0)",
         url="https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c",
         version="338,564 cells x 27,632 genes; identifiers read from the file's own uns/citation",
         role="The reference every method solves against: the vendored frozen signature "
              "(8-type roster, annotation_level_3) and the raw/X-count rebuild that is "
              "primary; the source of the synthetic mixtures.",
         citation="Ruiz-Moreno C, et al. Neuro-Oncology 27:2281-2295 (2025)",
         cite_note="the dataset itself cites the preprint: Ruiz-Moreno C, et al. bioRxiv (2022), "
                   "doi:10.1101/2022.08.27.505439",
         doi="10.1093/neuonc/noaf113", first_author="Ruiz-Moreno", year=2025,
         extra_dois=[("10.1101/2022.08.27.505439", "Ruiz-Moreno", 2022)],
         files=[F(REF / "gbmap_core.h5ad",
                  remote="https://datasets.cellxgene.cziscience.com/"
                         "861acfd8-25f0-418b-a445-aa96da232827.h5ad")],
         consumers=[("scripts/run_all.py", "gbmap_core"),
                    ("scripts/gbmap_t_vs_b.py", "gbmap_core"),
                    ("ivygap/bench/run_benchmark.py", "gbmap_core"),
                    ("scripts/remeasure_method.py", "gbmap_core")],
         keywords=["GBmap"]),
    dict(id="D14", group=2, status="used",
         name="Abdelfattah et al. 2022 glioma single-cell RNA-seq (independent atlas)",
         repository="GEO (raw and sample identifiers); processed matrix and clusters from the "
                    "Broad Single Cell Portal",
         accession="GEO GSE182109 (GSM5518596-GSM5518639; raw supplementary archive "
                   "GSE182109_RAW); processed matrix + clusters from the Single Cell Portal, whose "
                   "study accession was not recorded -- content verified against GEO's raw counts",
         url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109",
         version="Processed_matrix.mtx.gz (Seurat log-normalised, inverted to counts here), "
                 "genes.tsv, barcodes.tsv, Cluster_GBM.txt",
         role="The independent-atlas test of the B-over-T result: a reference not among "
              "GBmap's 16 source studies. All 201,986 barcodes carry one of the series' 44 "
              "GSM prefixes; none comes from Neftel's GSE131928. Provenance closed against GEO "
              "(scripts/verify_gse182109_raw.py): 201,893 of the 201,986 processed cells are in "
              "their own sample's Cell Ranger barcodes (the 93 others, in two samples, look like a "
              "more permissive cell call), and on the samples checked every recovered count equals "
              "GEO's raw count and every library equals the raw UMI total.",
         citation="Abdelfattah N, et al. Nat Commun 13:767 (2022)",
         doi="10.1038/s41467-022-28372-y", first_author="Abdelfattah", year=2022,
         files=[F(GEO182 / "Processed_matrix.mtx.gz.download" / "Processed_matrix.mtx.gz",
                  compare="none"),
                F(GEO182 / "genes.tsv", compare="none"),
                F(GEO182 / "barcodes.tsv", compare="none"),
                F(GEO182 / "Cluster_GBM.txt", compare="none"),
                F(GEO182 / "gsm_to_sample.tsv", "GSM -> sample name, from GEO", compare="none"),
                F(GEO182 / "GSE182109_RAW", "GEO's raw archive: 44 Cell Ranger outputs (5.0.0 x19, "
                  "5.0.1 x25), 132 files; sha256 of each in results/gse182109_raw_verification.json",
                  compare="none"),

                F(REF / "abdelfattah_2022" / "profile.csv", "built reference", status="derived"),
                F(REF / "abdelfattah_2022" / "provenance.json", "build record", status="derived")],
         consumers=[("scripts/build_abdelfattah_reference.py", "Cluster_GBM"),
                    ("scripts/atlas_t_vs_b.py", "Cluster_GBM"),
                    ("scripts/independent_atlas_test.py", "abdelfattah"),
                    ("scripts/verify_gse182109_raw.py", "GSE182109_RAW")],
         keywords=["Abdelfattah", "independent atlas", "GSE182109"]),
    dict(id="D15", group=2, status="used",
         name="Neftel et al. 2019 glioblastoma single-cell RNA-seq",
         repository="GEO",
         accession="GEO GSE131928; the Smart-seq2 matrix and cell assignments as distributed on "
                   "the Broad Single Cell Portal",
         url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928",
         version="processed TPM tables and authors' metadata",
         role="Reference-dependence arm (a Neftel-only reference). One of GBmap's constituent "
              "studies, so NOT an independent reference.",
         citation="Neftel C, et al. Cell 178:835-849 (2019)",
         doi="10.1016/j.cell.2019.06.024", first_author="Neftel", year=2019,
         files=[F(REPOS / "SCDC" / "neftel_data" / "IDHwtGBM.processed.SS2.logTPM.txt",
                  "Smart-seq2 log-TPM, Broad Single Cell Portal copy -- the matrix the build reads "
                  "(its GEO-named twin GSM3828672 was byte-identical and deleted, docs/DELETED_FILES.md)",
                  compare="none"),
                F(REPOS / "SCDC" / "neftel_data" / "GBM_Metadata", "authors' cell assignments (SCP)",
                  compare="none"),
                F(REPOS / "SCDC" / "neftel_data" / "GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv",
                  "10x subset -- read by nothing", compare="none"),
                F(REPOS / "SCDC" / "neftel_data" / "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx",
                  "GEO sample sheet -- read by nothing", compare="none")],
         consumers=[("scripts/build_neftel_reference.py", "IDHwtGBM.processed.SS2.logTPM.txt")],
         keywords=["Neftel"]),
    dict(id="D16", group=2, status="used",
         name="Darmanis et al. 2017 glioblastoma single-cell RNA-seq (core vs periphery)",
         repository="GEO",
         accession="GEO GSE84465",
         url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465",
         version="GSE84465_GBM_All_data.csv.gz and the series matrix, fetched 2026-09-14",
         role="Two uses: the only cross-patient core-vs-periphery constraint check (C1, four "
              "gates), and a Darmanis-only reference in the reference-dependence arm (a GBmap "
              "constituent, so not independent).",
         citation="Darmanis S, et al. Cell Rep 21:1399-1410 (2017)",
         doi="10.1016/j.celrep.2017.10.030", first_author="Darmanis", year=2017,
         files=[F(R / "darmanis_2017" / "GSE84465_GBM_All_data.csv.gz",
                  remote="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/suppl/"
                         "GSE84465_GBM_All_data.csv.gz"),
                F(R / "darmanis_2017" / "GSE84465_series_matrix.txt.gz",
                  remote="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/matrix/"
                         "GSE84465_series_matrix.txt.gz"),
                F(R / "darmanis_2017" / "cell_metadata.csv", "parsed from the series matrix",
                  status="derived")],
         consumers=[("scripts/build_darmanis_reference.py", "GSE84465_GBM_All_data"),
                    ("scripts/darmanis_constraint_check.py", "GSE84465")],
         keywords=["Darmanis"]),
    # ---------------------------------------------------------------- no-deconvolution checks
    dict(id="D17", group=3, status="used",
         name="Mossi Albiach et al. 2023 spatially sampled glioblastoma single-cell atlas",
         repository="CZ CELLxGENE Discover",
         accession="collection 113a558a-e96e-4643-81db-140e95c58578; dataset "
                   "d45b4ce6-9725-4d79-b97a-70a44158bdbf (donor SL040, schema 7.1.0)",
         url="https://cellxgene.cziscience.com/collections/113a558a-e96e-4643-81db-140e95c58578",
         version="135,482 cells x 58,234 genes, one donor, 27 samples in 12 locations",
         role="Tests the registered constraints against counted (not inferred) composition "
              "per anatomic zone: 3 of 4 testable constraints hold.",
         citation="Mossi Albiach A, et al. bioRxiv (2023)",
         doi="10.1101/2023.09.01.555882", first_author="Mossi Albiach", year=2023,
         files=[F(REPOS / "SCDC" / "albiach_data" / "d45b4ce6-9725-4d79-b97a-70a44158bdbf.h5ad",
                  remote="https://datasets.cellxgene.cziscience.com/"
                         "d45b4ce6-9725-4d79-b97a-70a44158bdbf.h5ad")],
         consumers=[("scripts/albiach_constraint_check.py", "albiach_data")],
         keywords=["Albiach"]),
    dict(id="D18", group=3, status="used",
         name="Ivy GAP in-situ hybridisation (ISH) by anatomic structure",
         repository="Allen Institute Ivy GAP portal: the per-(gene x sub-block) ISH table in D02",
         accession="gene_expression_details.csv (/api/v2/gbm/ namespace)",
         url="https://glioblastoma.alleninstitute.org/api/v2/gbm/gene_expression_details.csv",
         version="the portal table, sha256 recorded at download; the API-based route "
                 "(scripts/fetch_ivygap_ish.py) is SUPERSEDED and its output is not used",
         role="A constraint check with no deconvolution anywhere in the chain: CD44 / BIRC5 "
              "ISH support C1 and C7.",
         citation="Puchalski RB, et al. Science 360:660-663 (2018)",
         doi="10.1126/science.aaf2666", first_author="Puchalski", year=2018,
         files=[F(R / "ivygap" / "gene_expression_details.csv", "listed under D02",
                  compare="none")],
         consumers=[("scripts/ish_constraint_check.py", "gene_expression_details")],
         keywords=["ISH", "in-situ"]),
    # ---------------------------------------------------------------- bundled
    dict(id="D32", group=3, status="used",
         name="Flow-cytometry immune composition of brain tumours by IDH status (Klemm et al. 2020, Figure 1F)",
         repository="Cell (Elsevier); article PDF and supplements placed locally by the user (gitignored)",
         accession="doi:10.1016/j.cell.2020.05.007 -- Figure 1F (PDF page 3); Table S1 (mmc1.pdf, cohort); "
                   "Table S2 (mmc2.pdf, gating definitions)",
         url="https://doi.org/10.1016/j.cell.2020.05.007",
         version="article PDF sha256 df824704317a2f84ec9a49e9d8ebff648789cd59b95ce6800f93f58f612b6503",
         role="Direct check of the lymphoid truth by IDH status (prespecified/klemm_t_vs_b.md): cohort-mean "
              "T, B and NK as % of CD45+ in 17 IDH-mutant and 40 IDH-wildtype gliomas, measured from the "
              "figure's vector geometry and validated against two numbers printed in the paper.",
         citation="Klemm F, et al. Cell 181:1643-1660 (2020)",
         doi="10.1016/j.cell.2020.05.007", first_author="Klemm", year=2020,
         files=[F(config.KLEMM_2020_PDF, "the article PDF (publisher copyright; not redistributed)", compare="none"),
                F(config.KLEMM_2020_PDF.parent / "mmc1.pdf", "Table S1, the flow-cytometry cohort", compare="none"),
                F(config.KLEMM_2020_PDF.parent / "mmc2.pdf", "Table S2, gating definitions", compare="none")],
         consumers=[("scripts/klemm_figure1f.py", "KLEMM_2020_PDF")],
         keywords=["Klemm", "flow cytometry"]),
    dict(id="D19", group=4, status="used",
         name="GBMDeconvoluteR marker sets (Ajaib immune, Ruiz-Moreno level-3 immune, Neftel "
              "four-state) and Ajaib's published imaging-mass-cytometry correlations",
         repository="GBMDeconvoluteR source repository (vendored); Ajaib et al. 2023 paper",
         accession="Ajaib_et_al_2022_GBM_Immune_markers.rds, "
                   "Moreno_et_al_2022_lvl3_immune_markers.rds, "
                   "Neftel_et_al_2019_four_state_neoplastic_markers.rds",
         url="https://github.com/GliomaGenomics/GBMDeconvoluteR",
         version="vendored copy in pipeline_packages / repos/GBMDeconvoluteR",
         role="The pre-registered IMC protein test (an external instrument for ACS) and the "
              "mesenchymal-program score.",
         citation="Ajaib S, et al. Neuro-Oncology 25:1236-1248 (2023)",
         doi="10.1093/neuonc/noad021", first_author="Ajaib", year=2023,
         files=[F(REPOS / "GBMDeconvoluteR", "vendored repository", compare="none")],
         consumers=[("scripts/imc_anchored_test.py", "Ajaib_et_al_2022_GBM_Immune_markers"),
                    ("scripts/mes_score.py", "Neftel_et_al_2019_four_state_neoplastic_markers")],
         keywords=["IMC", "MCP", "mesenchymal"]),
    dict(id="D20", group=4, status="used",
         name="MCPcounter default marker genes",
         repository="MCPcounter source repository",
         accession="Signatures/genes.txt",
         url="https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt",
         version="fetched 2026-09-15",
         role="The default-marker arm of the IMC test.",
         citation="Becht E, et al. Genome Biol 17:218 (2016)",
         doi="10.1186/s13059-016-1070-5", first_author="Becht", year=2016,
         files=[F(R / "mcp_default_genes.txt",
                  remote="https://raw.githubusercontent.com/ebecht/MCPcounter/master/"
                         "Signatures/genes.txt")],
         consumers=[("scripts/imc_anchored_test.py", "mcp_default_genes")],
         keywords=["MCP"]),
    dict(id="D21", group=4, status="used",
         name="EpiDISH blood reference centDHSbloodDMC.m (333 CpGs x 7 blood cell types)",
         repository="Bioconductor package EpiDISH (data object)",
         accession="EpiDISH::centDHSbloodDMC.m",
         url="https://bioconductor.org/packages/EpiDISH",
         version="EpiDISH 2.28.0 as installed (packageVersion, read 2026-10-01; the "
                 "methylation artefacts do not record it)",
         role="Turns D07/D08 into per-cell-type lymphoid fractions (RPC mode); a blood "
              "reference applied to brain tumour, so used for lymphoid sub-composition only.",
         citation="Teschendorff AE, et al. BMC Bioinformatics 18:105 (2017)",
         doi="10.1186/s12859-017-1511-5", first_author="Teschendorff", year=2017,
         files=[],
         consumers=[("scripts/methylation_celltypes.py", "centDHSbloodDMC"),
                    ("scripts/extract_epidish_probes.py", "centDHSbloodDMC")],
         keywords=["EpiDISH"]),
    dict(id="D22", group=4, status="used",
         name="BayesPrism's bundled gene annotation (gene groups and gene types, GENCODE v22)",
         repository="BayesPrism source repository (vendored, Danko lab)",
         accession="cleanup.genes() gene groups; select.gene.type() 'protein_coding'",
         url="https://github.com/Danko-Lab/BayesPrism",
         version="BayesPrism 2.2.3 as installed (packageVersion, read 2026-10-01); source "
                 "vendored in pipeline_packages / repos/BayesPrism-main",
         role="The authors' own gene filtering in the authors'-workflow BayesPrism runs on "
              "TCGA (post-registration extension).",
         citation="Chu T, et al. Nat Cancer 3:505-517 (2022)",
         doi="10.1038/s43018-022-00356-3", first_author="Chu", year=2022,
         files=[F(REPOS / "BayesPrism-main", "vendored repository", compare="none")],
         consumers=[("R/run_bayesprism_authors.R", "cleanup.genes"),
                    ("R/run_bayesprism_authors.R", "select.gene.type")],
         keywords=["BayesPrism"]),
    dict(id="D31", group=4, status="used",
         name="GIMiCC reference libraries (glioma-specific hierarchical DNA-methylation deconvolution)",
         repository="Bioconductor ExperimentHub; github.com/SalasLab/GIMiCC",
         accession="EH9483 (GIMiCC_Library.rda, 15 layer libraries); EH9482 (Capper_example_betas.rda, "
                   "the reproduction control)",
         url="https://github.com/SalasLab/GIMiCC",
         version="GIMiCC 0.99.1 @26cb8a15 (the user's copy, byte-identical); ExperimentHub cache sha256 "
                 "EH9483 bdfbf2849ed6513c3afa84560f24b4b49d74245e20cd2f2a15ea59395a0bdb14, "
                 "EH9482 47da7b0f14925d16fe17ac2e3fc975083e663142e1026ec0d00ab7dcc030e140 (both added 2024-04-02)",
         role="The second methylation truth for the lymphoid ordering (post-registration, "
              "prespecified/gimicc_truth_confirmation.md): GIMiCC run on the TCGA 450K betas (D07/D08), "
              "restricted to its 4,022 library CpGs by scripts/extract_gimicc_cpgs.py.",
         citation="Pike SC, et al. Acta Neuropathol Commun 12:170 (2024)",
         doi="10.1186/s40478-024-01874-0", first_author="Pike", year=2024,
         files=[F(REPOS / "GIMiCC-main", "the package source (the user's copy)", compare="none")],
         consumers=[("R/run_gimicc.R", "GIMiCC_Deconvo"),
                    ("scripts/extract_gimicc_cpgs.py", "EH9483"),
                    ("scripts/gimicc_truth.py", "run_gimicc.R")],
         keywords=["GIMiCC"]),
    # ---------------------------------------------------------------- derived
    dict(id="D23", group=5, status="derived",
         name="Frozen GBM signature matrix and cell-size factors (from GBmap, by the "
              "predecessor project) and its synthetic held-out benchmark",
         repository="this repository, reference_frozen/ (vendored unchanged)",
         accession="signature_matrix.tsv, cell_size_factors.csv, tcga_benchmark/*",
         url="",
         version="frozen 2026-07-12T03:36:26Z; sha256 recorded in PROVENANCE.json",
         role="The registered signature; yardstick 1 (500 pseudobulk mixtures, NNLS and SVR "
              "only), which shares GBmap with the ACS arm.",
         citation="", doi="", first_author="", year=None,
         files=[F(config.FROZEN_SIGNATURE_PATH,
                  recorded_sha=("json", config.FROZEN_REFERENCE_DIR / "PROVENANCE.json",
                                ["files", "signature_matrix.tsv", "sha256"])),
                F(config.FROZEN_CELL_SIZE_PATH,
                  recorded_sha=("json", config.FROZEN_REFERENCE_DIR / "PROVENANCE.json",
                                ["files", "cell_size_factors.csv", "sha256"])),
                F(config.FROZEN_REFERENCE_DIR / "tcga_benchmark" / "benchmark_summary.csv",
                  recorded_sha=("json", config.FROZEN_REFERENCE_DIR / "tcga_benchmark" /
                                "PROVENANCE.json", ["files", "benchmark_summary.csv", "sha256"]))],
         consumers=[("ivygap/config.py", "FROZEN_SIGNATURE_PATH")],
         keywords=["frozen signature"]),
    # ---------------------------------------------------------------- unused
    dict(id="D24", group=6, status="on disk, unused",
         name="CIBERSORT relative fractions for TCGA (PanImmune)",
         repository="NCI GDC, PanImmune publication page",
         accession="GDC file b3df502e-3594-46ef-9f94-d041a20a0b9a",
         url=GDC_DATA + "b3df502e-3594-46ef-9f94-d041a20a0b9a",
         version="TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
         role="Obtained for a comparison against a published deconvolution of the same "
              "cohort; no script reads it. DATA_SOURCES.md A2 lists it as a supply -- that "
              "comparison was never built.",
         citation="Thorsson V, et al. Immunity 48:812-830 (2018)",
         doi="10.1016/j.immuni.2018.03.023", first_author="Thorsson", year=2018,
         files=[F(IMM / "TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
                  remote=GDC_DATA + "b3df502e-3594-46ef-9f94-d041a20a0b9a")],
         consumers=[], keywords=[]),
    dict(id="D25", group=6, status="on disk, unused",
         name="ABSOLUTE segment tables, TCGA PanCanAtlas",
         repository="NCI GDC, PanCanAtlas publication page",
         accession="GDC file 0f4f5701-7b61-41ae-bda9-2805d1ca9781",
         url=GDC_DATA + "0f4f5701-7b61-41ae-bda9-2805d1ca9781",
         version="TCGA_mastercalls.abs_segtabs.fixed.txt",
         role="Downloaded with D09; nothing reads it.",
         citation="Carter SL, et al. Nat Biotechnol 30:413-421 (2012)",
         doi="10.1038/nbt.2203", first_author="Carter", year=2012,
         files=[F(VEN / "ABSOLUTE" / "TCGA_mastercalls.abs_segtabs.fixed.txt",
                  remote=GDC_DATA + "0f4f5701-7b61-41ae-bda9-2805d1ca9781")],
         consumers=[], keywords=[]),
    dict(id="D26", group=6, status="on disk, unused",
         name="TCGA-LGG GDC clinical and biospecimen tables; GDC whole-exome BAM manifest",
         repository="NCI GDC Data Portal",
         accession="clinical.project-tcga-lgg.2026-09-17, biospecimen.project-tcga-lgg."
                   "2026-09-17, gdc_manifest.2026-09-17.191037.txt",
         url="https://portal.gdc.cancer.gov/projects/TCGA-LGG",
         version="exported 2026-09-17",
         role="Nothing reads the tables. The manifest lists controlled-access BAMs that were "
              "never downloaded and are not needed.",
         citation="Grossman RL, et al. N Engl J Med 375:1109-1112 (2016) [GDC]",
         doi="10.1056/NEJMp1607591", first_author="Grossman", year=2016,
         files=[F(TCGA / "project" / "clinical.project-tcga-lgg.2026-09-17" / "clinical.tsv",
                  compare="none"),
                F(TCGA / "project" / "gdc_manifest.2026-09-17.191037.txt", compare="none")],
         consumers=[], keywords=[]),
    dict(id="D27", group=6, status="duplicate",
         name="GEO GSE107559 supplementary files (the Ivy GAP portal tables, mirrored)",
         repository="GEO",
         accession="GSE107559_ivygap_rows-genes.csv, GSE107559_ivygap_columns-samples.xlsx",
         url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107559",
         version="", role="A GEO mirror of D01's tables; nothing reads this copy.",
         citation="Puchalski RB, et al. Science 360:660-663 (2018)",
         doi="10.1126/science.aaf2666", first_author="Puchalski", year=2018,
         files=[F(VEN / "GSE107559" / "GSE107559_ivygap_rows-genes.csv", compare="none"),
                F(VEN / "GSE107559" / "GSE107559_ivygap_columns-samples.xlsx", compare="none")],
         consumers=[], keywords=[]),
    # ---------------------------------------------------------------- considered
    dict(id="D28", group=7, status="considered",
         name="Human Brain Cell Atlas v1.0 (Siletti et al.)",
         repository="CZ CELLxGENE Discover",
         accession="collection 283d65eb-dd53-496d-adb7-7570c7caa443",
         url="https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443",
         version="", role="A normal-brain reference for a second tissue; not obtained.",
         citation="Siletti K, et al. Science 382:eadd7046 (2023)",
         doi="10.1126/science.add7046", first_author="Siletti", year=2023,
         files=[], consumers=[], keywords=[]),
    dict(id="D29", group=7, status="considered",
         name="Allen Human Brain Atlas, bulk expression by region",
         repository="Allen Institute for Brain Science",
         accession="human.brain-map.org downloads",
         url="https://human.brain-map.org/static/download",
         version="", role="A second tissue for the anatomy test; fetch script exists, the "
                          "constraint file, roster and platform are the blockers.",
         citation="Hawrylycz MJ, et al. Nature 489:391-399 (2012)",
         doi="10.1038/nature11405", first_author="Hawrylycz", year=2012,
         files=[], consumers=[], keywords=[]),
    dict(id="D30", group=7, status="considered",
         name="Spatially resolved multi-omics of glioblastoma (Ravi et al.)",
         repository="publication data", accession="", url="",
         version="", role="Candidate second GBM cohort with spatial structure; not obtained.",
         citation="Ravi VM, et al. Cancer Cell 40:639-655 (2022)",
         doi="10.1016/j.ccell.2022.05.009", first_author="Ravi", year=2022,
         files=[], consumers=[], keywords=[]),
]

#: Citations that MUST fail the author check -- the negative control for the citation
#: verifier. This is the exact conflation found in docs/REFERENCES.md [10] on 2026-10-01:
#: GBmap's authors attached to the DOI of Ravi et al.'s Cancer Cell paper.
CONTROL_CONFLATIONS = [("10.1016/j.ccell.2022.05.009", "Ruiz-Moreno", 2022)]


# =============================================================================
# measurement
# =============================================================================

def _open_nocache(path: Path):
    fh = open(path, "rb")
    try:                                  # macOS: bypass the buffer cache, so hashing an
        fcntl.fcntl(fh.fileno(), 48, 1)   # 8 GB atlas does not evict a running job's pages
    except OSError:
        pass
    return fh


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with _open_nocache(path) as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_gunzip_file(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_gunzip_url(url: str) -> str:
    h = hashlib.sha256()
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120, context=_CTX) as r, gzip.GzipFile(fileobj=r) as g:
        for chunk in iter(lambda: g.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def remote_size(url: str) -> int | None:
    """Byte size the server reports, via a one-byte range request (GDC rejects HEAD)."""
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Range": "bytes=0-0"})
    try:
        with urllib.request.urlopen(req, timeout=60, context=_CTX) as r:
            cr = r.headers.get("Content-Range")
            if cr and "/" in cr:
                return int(cr.rsplit("/", 1)[1])
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def _get(d: dict, keys: list):
    for k in keys:
        d = d[k]
    return d


def measure_file(f: dict, cache: dict, hash_large: bool) -> dict:
    p: Path = f["path"]
    rel = str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
    out = {"path": rel, "exists": p.exists(), "note": f["note"]}
    if not p.exists():
        return out
    if p.is_dir():
        files = [q for q in p.rglob("*") if q.is_file() and not q.name.startswith(".")]
        out.update(kind="dir", n_files=len(files), bytes=sum(q.stat().st_size for q in files))
        if f["compare"] == "per_file_provenance":
            prov = json.loads((p / "provenance.json").read_text())
            want = {s["rna_well"]: s["sha256"] for s in prov["samples"]}
            bad, missing = [], []
            for well, sha in want.items():
                q = p / f"{well}.genes.results"
                if not q.exists():
                    missing.append(well)
                    continue
                key = f"{q}|{q.stat().st_size}|{int(q.stat().st_mtime)}"
                got = cache["sha256"].get(key) or sha256_file(q)
                cache["sha256"][key] = got
                if got != sha:
                    bad.append(well)
            out["per_file_check"] = {"n_recorded": len(want), "n_missing": len(missing),
                                     "n_mismatch": len(bad)}
        return out
    st = p.stat()
    out.update(kind="file", bytes=st.st_size)
    key = f"{p}|{st.st_size}|{int(st.st_mtime)}"
    if key in cache["sha256"]:
        out["sha256"] = cache["sha256"][key]
    elif st.st_size <= LARGE or (hash_large and f["status"] != "duplicate"):
        out["sha256"] = cache["sha256"][key] = sha256_file(p)
    if f["recorded_sha"]:
        _, src, keys = f["recorded_sha"]
        rec = _get(json.loads(Path(src).read_text()), keys)
        out["recorded_sha256"] = rec
        out["matches_recorded"] = (out.get("sha256") == rec) if out.get("sha256") else None
    return out


def verify_remote(f: dict, cache: dict) -> dict | None:
    url, how, p = f["remote"], f["compare"], f["path"]
    if not url or how == "none" or not p.exists():
        return None
    st = p.stat()
    if how == "size":
        rs = remote_size(url)
        return {"how": "byte size", "remote_bytes": rs, "local_bytes": st.st_size,
                "identical_size": rs == st.st_size}
    local_key = f"gunzip|{p}|{st.st_size}|{int(st.st_mtime)}"
    if how == "gunzip_remote":
        local = cache["sha256"].get(f"{p}|{st.st_size}|{int(st.st_mtime)}") or sha256_file(p)
    else:
        local = cache["sha256"].get(local_key) or sha256_gunzip_file(p)
        cache["sha256"][local_key] = local
    remote = sha256_gunzip_url(url)
    return {"how": "sha256 of the decompressed content, both sides",
            "remote_sha256": remote, "local_sha256": local, "identical_content": remote == local}


def crossref(doi: str) -> dict | None:
    req = urllib.request.Request(f"https://api.crossref.org/works/{doi}",
                                 headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60, context=_CTX) as r:
            m = json.load(r)["message"]
    except Exception as e:
        return {"error": str(e)}
    au = m.get("author") or []
    first = (au[0].get("family") or au[0].get("name") or "") if au else ""
    return {"title": (m.get("title") or [""])[0], "first_author": first,
            "year": ((m.get("issued") or {}).get("date-parts") or [[None]])[0][0],
            "venue": (m.get("container-title") or [""])[0],
            "checked": dt.date.today().isoformat()}


def _fold(s: str) -> str:
    s = unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z]", "", s)


def citation_ok(rec: dict | None, first_author: str, year: int | None) -> tuple[bool, str]:
    """Does the Crossref record for a DOI match the first author and year claimed?"""
    if not rec or "error" in rec:
        return False, "no Crossref record"
    fa, claimed = _fold(rec["first_author"]), _fold(first_author)
    # "Mossi Albiach" is deposited as family "Mossi Albiach" or split; accept containment.
    if not fa or not (fa == claimed or fa in claimed or claimed in fa):
        return False, f"first author is {rec['first_author']!r}, not {first_author!r}"
    if year is not None and rec["year"] not in (year, year - 1, year + 1):
        return False, f"year is {rec['year']}, not {year}"
    if year is not None and rec["year"] != year:
        return True, f"year deposited as {rec['year']} (online-first); claimed {year}"
    return True, "matches"


def consumer_check(consumers: list[tuple[str, str]]) -> list[dict]:
    out = []
    for rel, token in consumers:
        p = ROOT / rel
        found = p.exists() and token in p.read_text(errors="ignore")
        out.append({"file": rel, "token": token, "found": bool(found)})
    return out


def manuscript_sections(keywords: list[str]) -> list[str]:
    if not keywords or not MANUSCRIPT.exists():
        return []
    secs, cur, buf = [], "front matter", []
    for line in MANUSCRIPT.read_text().splitlines():
        if line.startswith("## ") or line.startswith("### "):
            secs.append((cur, "\n".join(buf)))
            cur, buf = line.lstrip("# ").strip(), []
        else:
            buf.append(line)
    secs.append((cur, "\n".join(buf)))
    hits = []
    for name, body in secs:
        if any(k.lower() in body.lower() for k in keywords):
            short = name.split(" · ")[0] if " · " in name else name
            if short not in hits:
                hits.append(short)
    return hits


# =============================================================================
# main
# =============================================================================

def human(n: int | None) -> str:
    if n is None:
        return "?"
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return str(n)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify-dois", action="store_true")
    ap.add_argument("--verify-remote", action="store_true")
    ap.add_argument("--hash-large", action="store_true")
    args = ap.parse_args()

    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    for k in ("sha256", "remote", "crossref"):
        cache.setdefault(k, {})

    if args.verify_dois:
        dois = {e["doi"] for e in REGISTRY if e["doi"]}
        dois |= {d for e in REGISTRY for d, _, _ in e.get("extra_dois", [])}
        dois |= {d for d, _, _ in CONTROL_CONFLATIONS}
        import time  # noqa: PLC0415
        for d in sorted(dois):
            # Crossref rate-limits bursts (HTTP 429). A failed fetch must never overwrite a good cached
            # record -- that is how seven verified citations became "no Crossref record" on
            # 2026-10-03. Back off and retry; keep the old record if the fetch still fails.
            rec = crossref(d)
            for wait in (5, 15, 45):
                if not (rec and "429" in str(rec.get("error", ""))):
                    break
                time.sleep(wait)
                rec = crossref(d)
            old_rec = cache["crossref"].get(d)
            if rec and "error" in rec and old_rec and "error" not in old_rec:
                print(f"crossref {d}: fetch failed ({rec['error']}); keeping the cached record", flush=True)
            else:
                cache["crossref"][d] = rec
            print(f"crossref {d}: {cache['crossref'][d].get('first_author', '?')} "
                  f"{cache['crossref'][d].get('year', '?')}", flush=True)
            time.sleep(1.0)                 # polite pacing

    rows, problems = [], []
    for e in REGISTRY:
        measured = [measure_file(f, cache, args.hash_large) for f in e["files"]]
        if args.verify_remote:
            for f, m in zip(e["files"], measured):
                r = verify_remote(f, cache)
                if r is not None:
                    r["checked"] = dt.date.today().isoformat()
                    cache["remote"][m["path"]] = r
                    print(f"remote {m['path']}: {r}", flush=True)
        for m in measured:
            m["remote_check"] = cache["remote"].get(m["path"])
        cons = consumer_check(e["consumers"])
        cites = []
        if e["doi"]:
            cites.append((e["doi"], *citation_ok(cache["crossref"].get(e["doi"]),
                                                 e["first_author"], e["year"])))
        for d, fa, yr in e.get("extra_dois", []):
            cites.append((d, *citation_ok(cache["crossref"].get(d), fa, yr)))
        secs = manuscript_sections(e["keywords"])
        rows.append({**{k: e[k] for k in ("id", "status", "name", "repository", "accession",
                                          "url", "version", "role", "citation", "doi")},
                     "cite_note": e.get("cite_note", ""),
                     "group": GROUPS[e["group"]], "files": measured, "consumers": cons,
                     "citation_checks": cites, "manuscript_sections": secs,
                     "keywords": e["keywords"]})
        # what counts as a problem
        for m in measured:
            fstat = next(f["status"] for f in e["files"] if
                         (str(f["path"].relative_to(ROOT)) if f["path"].is_relative_to(ROOT)
                          else str(f["path"])) == m["path"])
            if (not m["exists"] and e["status"] in ("used", "derived")
                    and fstat in (None, "derived")):
                problems.append(f"{e['id']}: missing {m['path']}")
            if m.get("matches_recorded") is False:
                problems.append(f"{e['id']}: {m['path']} sha256 differs from the recorded hash")
            pfc = m.get("per_file_check")
            if pfc and (pfc["n_missing"] or pfc["n_mismatch"]):
                problems.append(f"{e['id']}: {pfc}")
            rc = m.get("remote_check") or {}
            if rc.get("identical_size") is False or rc.get("identical_content") is False:
                problems.append(f"{e['id']}: {m['path']} differs from the public source ({rc})")
        for c in cons:
            if not c["found"]:
                problems.append(f"{e['id']}: {c['file']} no longer contains {c['token']!r}")
        for d, ok, why in cites:
            if not ok and cache["crossref"].get(d):
                problems.append(f"{e['id']}: citation {d}: {why}")

    control_results = []
    for d, fa, yr in CONTROL_CONFLATIONS:
        ok, why = citation_ok(cache["crossref"].get(d), fa, yr)
        control_results.append({"doi": d, "claimed_first_author": fa, "rejected": not ok,
                                "why": why})
        if ok and cache["crossref"].get(d):
            problems.append(f"NEGATIVE CONTROL FAILED: conflated citation {d} was accepted")

    CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    write_csv(rows)
    write_md(rows, problems, control_results)
    print(f"{len(rows)} datasets; {len(problems)} problem(s)")
    for p in problems:
        print("  PROBLEM:", p)
    return 0


def write_csv(rows: list[dict]) -> None:
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "group", "status", "dataset", "repository", "accession", "url",
                    "version", "used_for", "citation", "doi", "local_files", "total_bytes",
                    "read_by", "manuscript_sections_keyword_scan"])
        for r in rows:
            present = [m for m in r["files"] if m["exists"]]
            w.writerow([r["id"], r["group"], r["status"], r["name"], r["repository"],
                        r["accession"], r["url"], r["version"], r["role"], r["citation"],
                        r["doi"], "; ".join(m["path"] for m in r["files"]),
                        sum(m.get("bytes", 0) for m in present),
                        "; ".join(c["file"] for c in r["consumers"] if c["found"]),
                        "; ".join(r["manuscript_sections"])])


def write_md(rows: list[dict], problems: list[str], controls: list[dict]) -> None:
    today = dt.date.today().isoformat()
    L = ["# Data inventory -- every public dataset this study touches",
         "",
         f"*Generated {today} by `scripts/build_data_inventory.py`. Do not edit by hand: "
         "every size, hash, consumer and citation check below is measured each run. "
         "`docs/supplementary/data_inventory_full.csv` is the same content in one row per "
         "dataset; Supplementary Table S7 is its compact form (`scripts/build_supplementary.py`).*",
         "",
         "**How to read the checks.** *read by* lists scripts that were opened and found to "
         "contain the token that reads the dataset. *source check* compares the local copy "
         "with the public server: identical byte size, or, where one side is stored "
         "decompressed, the sha256 of the decompressed content on both sides. *citation* is "
         "the DOI's Crossref record compared with the first author and year claimed here; "
         "the verifier is shown rejecting a known conflation before any citation is trusted.",
         ""]
    if problems:
        L += ["## Open problems", ""] + [f"- {p}" for p in problems] + [""]
    else:
        L += ["**Open problems: none.** Every listed file is present, every consumer still "
              "reads what it is said to read, every checked copy matches its source, and "
              "every citation matches its DOI.", ""]
    L += ["## Summary", "",
          "| id | dataset | status | accession | used for |",
          "|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['id']} | {r['name']} | {r['status']} | {r['accession']} | "
                 f"{r['role'].split('. ')[0].rstrip('.')} |")
    L.append("")
    for g in GROUPS:
        group_rows = [r for r in rows if r["group"] == g]
        if not group_rows:
            continue
        L += [f"## {g}", ""]
        for r in group_rows:
            L += [f"### {r['id']} · {r['name']}", ""]
            L += [x for x in (
                f"- **Status:** {r['status']}",
                f"- **Repository:** {r['repository']}",
                f"- **Accession:** {r['accession']}" if r["accession"] else "",
                f"- **Source:** <{r['url']}>" if r["url"] else "",
                f"- **Version:** {r['version']}" if r["version"] else "",
                f"- **Used for:** {r['role']}",
                (f"- **Citation:** {r['citation']}" + (f", doi:{r['doi']}" if r["doi"] else "")
                 + (f" -- {r['cite_note']}" if r.get("cite_note") else ""))
                if r["citation"] else "") if x]
            for d, ok, why in r["citation_checks"]:
                L.append(f"  - Crossref {d}: {'verified' if ok else 'NOT VERIFIED'} ({why})")
            if r["manuscript_sections"]:
                L.append("- **Manuscript sections mentioning it** (keyword scan for "
                         + ", ".join(f"'{k}'" for k in r["keywords"]) + "; a mention, "
                         "not proof of use): " + "; ".join(r["manuscript_sections"]))
            if r["consumers"]:
                L.append("- **Read by:** " + "; ".join(
                    f"`{c['file']}`" + ("" if c["found"] else " **(token missing)**")
                    for c in r["consumers"]))
            if r["files"]:
                L += ["", "| local file | size | sha256 (first 12) | check |", "|---|---|---|---|"]
                for m in r["files"]:
                    if not m["exists"]:
                        L.append(f"| `{m['path']}` | -- | -- | **absent** |")
                        continue
                    checks = []
                    if m.get("matches_recorded") is True:
                        checks.append("= recorded hash")
                    elif m.get("matches_recorded") is False:
                        checks.append("**differs from recorded hash**")
                    pfc = m.get("per_file_check")
                    if pfc:
                        checks.append(f"{pfc['n_recorded'] - pfc['n_missing'] - pfc['n_mismatch']}"
                                      f"/{pfc['n_recorded']} files = recorded hashes")
                    rc = m.get("remote_check")
                    if rc:
                        if "identical_size" in rc:
                            checks.append("= source size" if rc["identical_size"] else
                                          f"**source is {rc['remote_bytes']} bytes**")
                        else:
                            checks.append("= source content" if rc["identical_content"] else
                                          "**source content differs**")
                        checks[-1] += f" ({rc.get('checked', '')})"
                    size = human(m.get("bytes"))
                    if m.get("kind") == "dir":
                        size += f", {m['n_files']} files"
                    sha = (m.get("sha256") or "")[:12] or "not hashed"
                    note = f" -- {m['note']}" if m["note"] else ""
                    L.append(f"| `{m['path']}`{note} | {size} | `{sha}` | "
                             f"{'; '.join(checks) or '--'} |")
            L.append("")
    L += ["## Negative control for the citation check", "",
          "| DOI | claimed first author | rejected? | why |", "|---|---|---|---|"]
    for c in controls:
        L.append(f"| {c['doi']} | {c['claimed_first_author']} | "
                 f"{'yes' if c['rejected'] else '**NO**'} | {c['why']} |")
    L += ["", "## Data availability statement (draft for the manuscript)", "",
          "All data analysed in this study are publicly available. Bulk RNA-seq of "
          "laser-microdissected glioblastoma structures is from the Ivy Glioblastoma Atlas "
          "Project (Allen Institute for Brain Science; release 2014-11-25; GEO GSE107559). "
          "TCGA-GBM and TCGA-LGG RNA-seq (STAR counts) and HumanMethylation450 data were "
          "obtained from the UCSC Xena GDC and TCGA hubs; ABSOLUTE purity calls, leukocyte "
          "fractions and MC3 mutation calls from the TCGA PanCanAtlas resources at the NCI "
          "Genomic Data Commons and UCSC Xena; TCGA-GBM clinical annotations from cBioPortal "
          "(gbm_tcga_pub2013). Single-cell references are GBmap (CZ CELLxGENE Discover "
          "collection 999f2a15-3d7e-440b-96ae-2c806799c08c, dataset 861acfd8-25f0-418b-a445-"
          "aa96da232827), Abdelfattah et al. (GEO GSE182109), Neftel et al. (GEO GSE131928), "
          "Darmanis et al. (GEO GSE84465) and Mossi Albiach et al. (CZ CELLxGENE Discover "
          "collection 113a558a-e96e-4643-81db-140e95c58578). Exact files, versions, sizes and "
          "checksums are listed in Supplementary Table S7.", ""]
    OUT_MD.write_text(re.sub(r"\n{3,}", "\n\n", "\n".join(L)))


if __name__ == "__main__":
    raise SystemExit(main())
