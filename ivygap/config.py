"""
config.py — the single control panel for the Ivy GAP GBM anatomic deconvolution project.

WHAT THIS FILE IS
-----------------
Every other module reads its paths, its cell-type roster, its anatomic-structure
roster and its analysis knobs from here. Nothing in this file downloads data, does
math, or produces a result — it only *declares* shared truth. Change something once
here and the whole pipeline changes with it.

This is the Ivy GAP successor to the TCGA project's `Step 0 Setup/config.py`. The
cell-type roster is deliberately carried over unchanged so results remain comparable
across the two cohorts; what changes is the *unit of analysis*.

    TCGA-GBM      : 1 bulk sample per patient  -> patient-level composition
    Ivy GAP       : many laser-microdissected anatomic samples per tumor
                    -> (patient x anatomic structure) composition, nested design

That nesting is the whole point of this rebuild, and it is also the main statistical
hazard: tumors contribute unequal numbers of blocks. Every scoring path in this
project therefore aggregates *within patient first*, then across patients with equal
weight per patient. See `ivygap/bench/equal_footing.py`, which enforces it.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

# =============================================================================
# SECTION 1 — WHERE THINGS LIVE
# =============================================================================

PACKAGE_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_ROOT.parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"                  # untouched downloads
PROCESSED_DIR = DATA_DIR / "processed"      # aligned/standardised intermediates
REFERENCE_DIR = DATA_DIR / "reference"      # single-cell references + signature matrices
PSEUDOBULK_DIR = DATA_DIR / "pseudobulk"    # synthetic mixtures with known truth

RESULTS_DIR = PROJECT_ROOT / "results"
ESTIMATES_DIR = RESULTS_DIR / "estimates"   # one composition table per method
BENCH_DIR = RESULTS_DIR / "benchmark"
ANATOMIC_DIR = RESULTS_DIR / "anatomic"
SURVIVAL_DIR = RESULTS_DIR / "survival"
RELEASE_DIR = PROJECT_ROOT / "release"      # the stable, public-facing snapshot
FIGURES_DIR = RESULTS_DIR / "figures"

ALL_OUTPUT_DIRS = [
    DATA_DIR, RAW_DIR, PROCESSED_DIR, REFERENCE_DIR, PSEUDOBULK_DIR,
    RESULTS_DIR, ESTIMATES_DIR, BENCH_DIR, ANATOMIC_DIR, SURVIVAL_DIR,
    RELEASE_DIR, FIGURES_DIR,
]

# --- Ivy GAP raw inputs -------------------------------------------------------
IVYGAP_ZIP_PATH = RAW_DIR / "gene_expression_matrix_2014-11-25.zip"
IVYGAP_FPKM_PATH = RAW_DIR / "ivygap" / "fpkm_table.csv"
IVYGAP_GENES_PATH = RAW_DIR / "ivygap" / "rows-genes.csv"
IVYGAP_SAMPLES_PATH = RAW_DIR / "ivygap" / "columns-samples.csv"
IVYGAP_TUMOR_DETAILS_PATH = RAW_DIR / "ivygap" / "tumor_details.csv"
# The portal's per-RNA-seq-sample metadata. Authoritative for protocol step 3's
# reconciliation gate, and the only file carrying the donor_id <-> tumor_id join.
IVYGAP_RNA_SEQ_DETAILS_PATH = RAW_DIR / "ivygap" / "rna_seq_samples_details.csv"

# --- processed intermediates --------------------------------------------------
BULK_EXPRESSION_PATH = PROCESSED_DIR / "ivygap_bulk_expression.tsv"   # genes x samples
SAMPLE_MANIFEST_PATH = PROCESSED_DIR / "ivygap_sample_manifest.tsv"   # one row per RNA-seq sample
CLINICAL_PATH = PROCESSED_DIR / "ivygap_clinical.tsv"                 # one row per tumor/patient
ALIGNED_BULK_PATH = PROCESSED_DIR / "bulk_aligned.tsv"
ALIGNED_SIGNATURE_PATH = PROCESSED_DIR / "signature_aligned.tsv"

# --- the vendored frozen reference (protocol step 1) -------------------------
# Committed to the repo with recorded hashes so the pipeline runs without rebuilding a
# reference from a multi-gigabyte atlas. See reference_frozen/PROVENANCE.json for what
# it is, where it came from, and what it cannot support.
FROZEN_REFERENCE_DIR = PROJECT_ROOT / "reference_frozen"
FROZEN_SIGNATURE_PATH = FROZEN_REFERENCE_DIR / "signature_matrix.tsv"
FROZEN_CELL_SIZE_PATH = FROZEN_REFERENCE_DIR / "cell_size_factors.csv"

# --- signature matrices (one per single-cell reference) -----------------------
SIGNATURE_MATRIX_PATHS = {
    "gbmap": REFERENCE_DIR / "signature_gbmap.tsv",
    "neftel": REFERENCE_DIR / "signature_neftel.tsv",
}
PRIMARY_REFERENCE = "gbmap"   # the frozen reference; the others exist so SCDC ENSEMBLE
                              # has real cross-dataset weighting to do, and so that
                              # "robust to reference choice" is a *measurable* property.

# =============================================================================
# SECTION 2 — THE CELL-TYPE ROSTER
# -----------------------------------------------------------------------------
# Carried over unchanged from the TCGA project so the two cohorts stay comparable.
# Fixed order: columns of the signature matrix, columns of every output table, and
# the slice order in every stacked-bar figure.
#
# Astrocyte is the eighth, DIAGNOSTIC-ONLY type. In the TCGA project GFAP alone
# accounted for 38.8% of squared reconstruction error and Astrocyte proved
# unidentifiable against Tumor under NNLS/SVR, so it was demoted to a sidecar.
# It is retained here *as a sidecar* for two reasons:
#   1. Ivy GAP's Leading Edge is genuinely astrocyte-rich normal brain, so the
#      identifiability question is materially different from TCGA's cellular-tumor
#      samples and deserves to be re-asked rather than assumed.
#   2. Tumor/Astrocyte collinearity is exactly the failure mode that Elastic Net
#      (L2 shrinkage across collinear columns) and the Bayesian model (an explicit
#      prior) are supposed to improve on. Dropping it would delete the test case.
# It is never summed into the seven primary types.
# =============================================================================

PRIMARY_CELL_TYPES = [
    "Tumor",                 # malignant glioma cells (all molecular states pooled)
    "Macrophage_Microglia",  # tumour-associated myeloid cells
    "T_cell",
    "NK_cell",
    "B_cell",
    "Endothelial",
    "Oligodendrocyte",
]
SIDECAR_CELL_TYPES = ["Astrocyte"]
CELL_TYPES = PRIMARY_CELL_TYPES + SIDECAR_CELL_TYPES

MARKER_GENES = {
    "Tumor":                ["EGFR", "SOX2", "PDGFRA", "CDK4", "OLIG2"],
    "Macrophage_Microglia": ["CD68", "AIF1", "ITGAM", "CX3CR1", "P2RY12", "TMEM119", "CD14"],
    "T_cell":               ["CD3D", "CD3E", "CD2", "CD8A", "IL7R"],
    "NK_cell":              ["NCAM1", "NKG7", "KLRD1", "GNLY", "KLRF1"],
    "B_cell":               ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Endothelial":          ["PECAM1", "VWF", "CLDN5", "FLT1", "CDH5"],
    "Oligodendrocyte":      ["MBP", "PLP1", "MOG", "MAG", "SOX10"],
    "Astrocyte":            ["GFAP", "AQP4", "S100B", "SLC1A3", "ALDH1L1"],
}

# Mapping from GBmap `annotation_level_3` labels onto the roster above. Carried over
# from the TCGA project verbatim so the reference is built identically.
GBMAP_CELL_TYPE_MAP = {
    "AC-like": "Tumor", "MES-like": "Tumor", "NPC-like": "Tumor", "OPC-like": "Tumor",
    "TAM-BDM": "Macrophage_Microglia", "TAM-MG": "Macrophage_Microglia",
    "CD4/CD8": "T_cell",
    "NK": "NK_cell",
    "B cell": "B_cell",
    "Endothelial": "Endothelial",
    "Oligodendrocyte": "Oligodendrocyte",
    "Astrocyte": "Astrocyte",
}

# =============================================================================
# SECTION 3 — THE ANATOMIC-STRUCTURE ROSTER  (this is what Ivy GAP adds)
# -----------------------------------------------------------------------------
# Ivy GAP laser-microdissected each tumour into histologically defined structures
# and RNA-sequenced them separately. `structure_abbreviation` in columns-samples.csv
# carries the label. These five are the ones the Ivy GAP reference set is built on
# (Puchalski et al., Science 2018) and the ones this project analyses.
#
# Ordered from tumour periphery to tumour core, because several analyses treat the
# order as an ordinal axis.
# =============================================================================

PRIMARY_STRUCTURES = [
    "LE",   # Leading Edge — outermost; mostly normal brain infiltrated by few tumour cells
    "IT",   # Infiltrating Tumor — transitional
    "CT",   # Cellular Tumor — the dense malignant core
    "MVP",  # Microvascular Proliferation
    "PAN",  # Peri-necrotic / pseudopalisading cells around necrosis
]

# -----------------------------------------------------------------------------
# HOW IVY GAP ACTUALLY LABELS ITS SAMPLES  (verified against the 2014-11-25 release)
# -----------------------------------------------------------------------------
# The archive's 270 RNA-seq samples are TWO different studies sharing one file, and
# telling them apart is the single most important thing this loader does.
#
#   122 samples / 10 tumours   anatomic structures, selected by H&E histology.
#                              `structure_abbreviation` ends in "-reference-histology"
#                              e.g. "CT-reference-histology", "CTmvp-reference-histology"
#   148 samples / 34 tumours   putative cancer-stem-cell clusters, selected by ISH with
#                              an 18-probe reference set. Labels carry the probe(s)
#                              instead: "CT-CD44", "CTpnz-PI3", "CT-control-TGFBR2".
#
# Only the first study is the Anatomy Test's instrument. Its structure labels were
# assigned by a neuropathologist reading H&E, independent of any RNA — which is exactly
# what makes them usable as a test rather than as more data. The ISH-cluster labels were
# assigned using expression, so scoring anatomic constraints on them would be circular.
#
# A naive match on the acronym prefix would silently pull all 148 cluster samples into
# the anatomic cohort — "CT-CD44" starts with "CT" — roughly tripling the apparent
# sample size with samples that do not belong. Hence the explicit suffix requirement.
ANATOMIC_SUFFIX = "-reference-histology"

# Within the anatomic study the archive uses CTmvp/CTpan; the constraint file and the
# protocol use MVP/PAN. Same structures, one vocabulary on disk.
STRUCTURE_ALIASES = {
    "CT": "CT", "IT": "IT", "LE": "LE",
    "CTMVP": "MVP", "MVP": "MVP",
    "CTPAN": "PAN", "PAN": "PAN",
}


def is_anatomic_sample(raw: str) -> bool:
    """True for the 122 H&E-selected anatomic-structure samples, false for ISH clusters."""
    return isinstance(raw, str) and raw.strip().endswith(ANATOMIC_SUFFIX)


def canonical_structure(raw: str) -> str:
    """
    Fold a raw Ivy GAP structure label onto the roster above.

    Anatomic samples become LE / IT / CT / MVP / PAN. Cancer-stem-cell cluster samples
    keep their raw label, which will not match PRIMARY_STRUCTURES, so they are excluded
    from the anatomic analysis by construction rather than by a filter someone has to
    remember to apply.
    """
    if not isinstance(raw, str):
        return "UNKNOWN"
    raw = raw.strip()
    if not is_anatomic_sample(raw):
        return raw
    stem = raw[: -len(ANATOMIC_SUFFIX)]
    key = stem.replace("-", "").replace("_", "").replace(" ", "").upper()
    return STRUCTURE_ALIASES.get(key, stem)
STRUCTURE_LONG_NAMES = {
    "LE": "Leading Edge",
    "IT": "Infiltrating Tumor",
    "CT": "Cellular Tumor",
    "MVP": "Microvascular Proliferation",
    "PAN": "Pseudopalisading Cells Around Necrosis",
    "CTpnz": "Perinecrotic Zone",
    "CThbv": "Hyperplastic Blood Vessels",
    "CTne": "Necrosis",
}
# CTpnz / CThbv appear in the archive ONLY with an ISH probe suffix (e.g. "CTpnz-PI3"),
# never as reference histology, so there are no pure anatomic samples of them. They are
# kept here so the roster documents that they exist and are deliberately out of scope,
# not because any sample will ever be assigned to them.
SECONDARY_STRUCTURES = ["CTpnz", "CThbv", "CTne"]

# The reference / baseline structure that paired contrasts are taken against.
REFERENCE_STRUCTURE = "CT"

# =============================================================================
# SECTION 4 — ANALYSIS SETTINGS
# =============================================================================

# Linear CPM/TPM-like scale, never log. S*w = b is a LINEAR additive mixing model;
# log(a+b) != log(a)+log(b), so log-transforming before solving does not correspond
# to the model at all. (Verified the hard way on TCGA: log2 inputs drove T/NK/B to
# ~0% cohort-wide and put myeloid below endothelial, contradicting GBM biology.)
# Ivy GAP ships FPKM, which is already a linear per-sample-normalised scale; we
# rescale each column to sum to 1e6 so it is comparable to CPM.
NORMALIZATION = "cpm"
FPKM_INPUT = True          # Ivy GAP distributes FPKM, not integer counts
RESCALE_COLUMNS_TO = 1e6

RANDOM_SEED = 0

# Gene-space construction
MIN_GENES_SHARED = 500              # abort if bulk and signature overlap less than this
MAX_MISSING_FRACTION_PER_GENE = 0.2
SIGNATURE_GENES_PER_TYPE = 200      # top-N differentially expressed genes per cell type

# Deconvolution is run on the informative subset of the signature, not on every gene the
# reference and the bulk happen to share. Two reasons, in order of importance:
#
#   1. Statistics. Least squares is dominated by the largest residuals, which belong to
#      the most highly expressed genes — mostly housekeeping genes that are similar in
#      every cell type and carry no information about composition. CIBERSORT's LM22 uses
#      547 genes for exactly this reason.
#   2. Tractability. nu-SVR is superlinear in the number of rows; on the full 16,007-gene
#      intersection a single cohort takes hours, and the fit is not better for it.
#
# The equal-footing contract is unaffected either way: every method receives the same
# gene set, and its hash is in the certificate. What matters is that the choice is made
# once, here, and not per method.
USE_SIGNATURE_GENE_SUBSET = True

# Rows of every composition table must sum to 100 within this tolerance (percent).
PERCENT_SUM_TOLERANCE = 1.0

# Pseudobulk benchmark (known-truth mixtures used to SELECT the method)
N_TRAIN_PSEUDOBULK = 2000
N_TEST_PSEUDOBULK = 500
PSEUDOBULK_TEST_DONOR_FRAC = 0.2    # held out by DONOR, never by cell
PSEUDOBULK_TOTAL_CELLS = (300, 1000)
PSEUDOBULK_MAX_ZEROED = 3           # let some types be genuinely absent, so structural
                                    # zeros and limit-of-detection are measurable
PSEUDOBULK_DIRICHLET_ALPHA = 1.0

# Cross-validation. ALWAYS grouped by patient — never by sample. A tumour's blocks
# are not independent observations of each other.
N_CV_FOLDS = 5
N_CV_REPEATS = 10

# =============================================================================
# SECTION 5 — HELPERS
# =============================================================================


#: Which kind of run is writing. Set once, at startup, before anything else runs.
RUN_LABEL = "real"


def use_synthetic_paths() -> None:
    """
    Repoint every output directory at a synthetic-only tree.

    WHY THIS EXISTS
    ---------------
    A fixture run and a real run originally wrote to the same paths, so
    `run_all.py --synthetic` silently destroyed the results of a real run that had taken
    half an hour. The release manifest still said SYNTHETIC FIXTURE, so nothing was
    *mislabelled* — but the real numbers were gone, and the only reason they survived was
    that they had already been transcribed into RESULTS.md by hand.

    Separating the trees means a fixture run can never cost real work. Call this before
    `ensure_dirs()` and before importing anything that resolves a path at import time
    (nothing here does; every writer reads `config.<DIR>` at call time, which is what
    makes this rebinding work).
    """
    global RUN_LABEL, RESULTS_DIR, ESTIMATES_DIR, BENCH_DIR, ANATOMIC_DIR
    global SURVIVAL_DIR, RELEASE_DIR, FIGURES_DIR, ALL_OUTPUT_DIRS

    RUN_LABEL = "synthetic"
    RESULTS_DIR = PROJECT_ROOT / "results_synthetic"
    ESTIMATES_DIR = RESULTS_DIR / "estimates"
    BENCH_DIR = RESULTS_DIR / "benchmark"
    ANATOMIC_DIR = RESULTS_DIR / "anatomic"
    SURVIVAL_DIR = RESULTS_DIR / "survival"
    FIGURES_DIR = RESULTS_DIR / "figures"
    RELEASE_DIR = PROJECT_ROOT / "release_synthetic"

    ALL_OUTPUT_DIRS = [
        DATA_DIR, RAW_DIR, PROCESSED_DIR, REFERENCE_DIR, PSEUDOBULK_DIR,
        RESULTS_DIR, ESTIMATES_DIR, BENCH_DIR, ANATOMIC_DIR, SURVIVAL_DIR,
        RELEASE_DIR, FIGURES_DIR,
    ]


class ResultsTreeBusy(RuntimeError):
    """Another run holds the lock on this results tree."""


def acquire_run_lock(force: bool = False):
    """
    Take an exclusive lock on the results tree, and return a release callable.

    WHY THIS EXISTS
    ---------------
    This project has now been bitten twice by two runs sharing one output tree. The
    first time, a `--synthetic` fixture overwrote a completed 30-minute real run, and
    the numbers survived only because someone had copied them into RESULTS.md by hand;
    `use_synthetic_paths()` was the fix for that one. The second time, a background run
    reported as killed had not actually died — only its wrapper had — and it was still
    deconvolving into the same `results/` a fresh run was writing to.

    Labelling is not isolation, and neither is assuming a process is dead because
    something said so. The lock records the pid and start time, and a stale lock whose
    pid is gone is reclaimed automatically, so a crashed run does not block the next one.
    """
    import json as _json
    import os
    import time

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    lock = RESULTS_DIR / ".run.lock"

    if lock.exists() and not force:
        try:
            info = _json.loads(lock.read_text())
            pid = int(info.get("pid", -1))
        except Exception:                                  # noqa: BLE001
            info, pid = {}, -1
        alive = False
        if pid > 0:
            try:
                os.kill(pid, 0)
                alive = True
            except OSError:
                alive = False
        if alive:
            raise ResultsTreeBusy(
                f"{lock} is held by pid {pid}, started {info.get('started')}, "
                f"writing to {RESULTS_DIR}. Two runs sharing one results tree is how "
                f"this project lost a completed run once already. Stop that process, or "
                f"pass --force-lock if you are certain it is gone."
            )
        lock.unlink(missing_ok=True)

    lock.write_text(_json.dumps({
        "pid": os.getpid(),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "results_dir": str(RESULTS_DIR),
    }, indent=2))

    def release() -> None:
        lock.unlink(missing_ok=True)

    return release


def ensure_dirs() -> None:
    """Create every output directory. Safe to call repeatedly."""
    for d in ALL_OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    """Content hash of a file, used throughout for provenance."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_frame(df) -> str:
    """
    Stable content hash of a DataFrame, used by the equal-footing guard to prove that
    two methods really did receive byte-identical inputs. Index and columns are
    included, and float formatting is fixed, so the hash does not drift with pandas'
    display settings.
    """
    import numpy as np
    h = hashlib.sha256()
    h.update("\t".join(map(str, df.columns)).encode())
    h.update(b"\n")
    h.update("\t".join(map(str, df.index)).encode())
    h.update(b"\n")
    arr = np.ascontiguousarray(df.to_numpy(dtype="float64"))
    h.update(np.round(arr, 10).tobytes())
    return h.hexdigest()


def validate_config() -> None:
    """Fail loudly at import time if the declarations above contradict each other."""
    missing = [c for c in CELL_TYPES if c not in MARKER_GENES]
    if missing:
        raise ValueError(f"MARKER_GENES is missing entries for: {missing}")
    extra = [c for c in MARKER_GENES if c not in CELL_TYPES]
    if extra:
        raise ValueError(f"MARKER_GENES has entries not in CELL_TYPES: {extra}")
    if len(set(CELL_TYPES)) != len(CELL_TYPES):
        raise ValueError("CELL_TYPES contains duplicates")
    if REFERENCE_STRUCTURE not in PRIMARY_STRUCTURES:
        raise ValueError(f"REFERENCE_STRUCTURE {REFERENCE_STRUCTURE!r} is not a primary structure")
    if set(PRIMARY_STRUCTURES) & set(SECONDARY_STRUCTURES):
        raise ValueError("A structure cannot be both primary and secondary")
    unmapped = set(GBMAP_CELL_TYPE_MAP.values()) - set(CELL_TYPES)
    if unmapped:
        raise ValueError(f"GBMAP_CELL_TYPE_MAP maps to unknown cell types: {unmapped}")
    if PRIMARY_REFERENCE not in SIGNATURE_MATRIX_PATHS:
        raise ValueError(f"PRIMARY_REFERENCE {PRIMARY_REFERENCE!r} has no signature path")


validate_config()
