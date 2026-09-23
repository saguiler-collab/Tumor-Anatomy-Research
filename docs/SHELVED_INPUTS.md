# Shelved raw inputs — large files kept on disk but out of git

Two raw input files are **present on this machine, gzipped, and deliberately not tracked by
git**. Neither is lost. Neither is needed to reproduce any number in the paper, because both
have already been consumed into derived artefacts that *are* kept.

They were removed from git history rather than merely untracked, because GitHub rejects any
push containing a blob over **100 MB** and both had already been committed into the branch's
unpushed history. Untracking them at the tip would not have helped: the blobs live in the
commits, not only in the working tree.

## What was shelved

| file | raw size | now | sha256 (of the uncompressed file) |
|---|---|---|---|
| `TCGA_LGG/TCGA-LGG.star_counts.tsv` | 318 MB | `…star_counts.tsv.gz` | `b12726cf554d158fe252cba19d185efbf19bffc9ed1fa07c52635da1e26fe80d` |
| `ivygap/ivygap_downloads/gene_expression_matrix_2014-11-25/fpkm_table.csv` | 82 MB | `…fpkm_table.csv.gz` | `27b916d11b73e9195bbbb5026592e129dfd3cfd50b46ecef25be496f0bf05750` |

Restore either with `gunzip -k <file>.gz` — `-k` keeps the archive.

## Why shelving them costs nothing

**The Ivy GAP FPKM table is a duplicate.** The pipeline reads
`data/raw/ivygap/fpkm_table.csv`, set in `ivygap/config.py` as `IVYGAP_FPKM_PATH`, and that
copy is untouched and already gitignored. The shelved copy is the pristine download bundle's
version. The two were verified **byte-identical** — same sha256, above — so nothing reads the
shelved one and no code path changes.

**The TCGA-LGG STAR table has already been consumed.** It is the raw input to
`scripts/build_tcga_bulk.py`, whose output `data/processed/tcga_lgg_bulk_cpm.csv.gz`
(53.5 MB, 534 samples, 59,427 gene symbols after collapse) is what every downstream stage
actually reads. `data/processed/tcga_lgg_bulk_provenance.json` records the source path, the
`log2(x + 1)` scale correction, the gene mapping and the duplicate-symbol rule, so the
derived matrix is auditable without the raw table.

`build_tcga_bulk.py` resolves the `.gz` transparently (`_resolve`), so even re-running the
build from scratch works against the shelved file. If neither form is present it names the
file and points here rather than failing on a bare `FileNotFoundError`.

## Where to re-download, if this machine is ever lost

- **TCGA-LGG `star_counts`** — UCSC Xena GDC hub, TCGA-LGG cohort, `STAR - Counts` dataset.
  Despite the name the delivered values are `log2(x + 1)`, not counts; the build measures
  this rather than assuming it (100% of small inverted values sit on the 0.5 grid).
- **Ivy GAP `fpkm_table.csv`** — the `gene_expression_matrix_2014-11-25` bundle from the Ivy
  Glioblastoma Atlas Project, as recorded in `docs/DATA_SOURCES.md`.

## What is NOT shelved

Everything the results depend on. The frozen archives under `results_archive/` carry every
scored artefact with a hash per file, the derived bulk matrices are on disk, and
`docs/DATA_SOURCES.md` remains the record of what the archive contains.
