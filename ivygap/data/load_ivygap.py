"""
load_ivygap.py — turn the raw Allen Institute CSVs into the two tables the rest of
the pipeline consumes.

    bulk       genes (HGNC symbol) x samples (rna_well_id), linear CPM-like scale
    manifest   one row per RNA-seq sample: sample_id, patient_id, structure, block

A NOTE ON TOLERANT COLUMN RESOLUTION
------------------------------------
The Allen Institute has shipped `columns-samples.csv` with several different column
namings across releases (`structure_acronym` vs `structure_abbreviation`, `tumor_id`
vs `donor_id`, and so on). Rather than pin one spelling and break on the next
release, each field below is resolved from a list of accepted names and fails with an
explicit message naming what was actually present. Guessing silently — or worse,
falling back to positional indexing — is how a pipeline ends up analysing the wrong
column without anyone noticing.

CONTRACT VALIDATION
-------------------
Counts (samples, patients, structures) are measured from the data and recorded, never
asserted against hard-coded literals. The published Ivy GAP RNA-seq release is ~270
samples from ~37 tumours, but this loader reports what it actually found; a test that
hard-codes 270 to make itself pass would hide a truncated download.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from ivygap import config

# Accepted spellings for each field we need, in priority order.
_SAMPLE_ID_COLS = ["rna_well_id", "rna_well", "well_id", "sample_id"]
_STRUCTURE_COLS = ["structure_acronym", "structure_abbreviation", "structure_abbrev",
                   "structure_acronym_abbrev"]
_STRUCTURE_NAME_COLS = ["structure_name", "structure_full_name"]
_PATIENT_COLS = ["tumor_id", "donor_id", "tumor_name", "specimen_name", "donor_name"]
_BLOCK_COLS = ["block_name", "sub_block_name", "specimen_name"]
_GENE_ID_COLS = ["gene_id", "gene_idx", "id"]
_GENE_SYMBOL_COLS = ["gene_symbol", "symbol", "gene_name"]


@dataclass
class LoadReport:
    """Measured properties of what was actually loaded. Written to disk as evidence."""
    n_genes_raw: int
    n_genes_after_symbol_collapse: int
    n_samples: int
    n_patients: int
    structures_found: dict          # structure -> sample count
    samples_per_patient: dict       # min / median / max
    n_primary_structure_samples: int
    n_anatomic_study_samples: int
    n_anatomic_study_tumors: int
    n_cluster_study_samples: int
    n_cluster_study_tumors: int
    patients_with_reference_structure: int
    structures_per_tumor: dict
    bulk_sha256: str
    manifest_sha256: str

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def _resolve(df: pd.DataFrame, candidates: list[str], what: str) -> str:
    """Pick the first accepted column name that is present, or fail informatively."""
    lower = {c.lower().strip(): c for c in df.columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    raise KeyError(
        f"could not find a column for {what}. Accepted names: {candidates}. "
        f"The file actually has: {list(df.columns)}"
    )


def _normalise_structure(value: str) -> str:
    """
    Map a raw structure label onto the roster in config.

    Allen ships acronyms in several decorations across releases ('CT-mvp', 'CTmvp',
    'CT_mvp'), and the protocol calls those same structures MVP and PAN. One shared
    canonicaliser lives in config so the loader and the constraint file cannot drift
    into different vocabularies and split one structure's blocks across two labels.
    """
    return config.canonical_structure(value)


def load_raw() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the three raw CSVs exactly as distributed, with no transformation."""
    for path in (config.IVYGAP_FPKM_PATH, config.IVYGAP_GENES_PATH,
                 config.IVYGAP_SAMPLES_PATH):
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not found. Run: python -m ivygap.data.download_ivygap"
            )
    fpkm = pd.read_csv(config.IVYGAP_FPKM_PATH, index_col=0)
    genes = pd.read_csv(config.IVYGAP_GENES_PATH)
    samples = pd.read_csv(config.IVYGAP_SAMPLES_PATH)
    return fpkm, genes, samples


def build_manifest(samples: pd.DataFrame) -> pd.DataFrame:
    """
    One tidy row per RNA-seq sample.

    `patient_id` is the grouping variable for every fold, every contrast and every
    aggregation in this project. Getting it from the wrong column would silently
    destroy the equal-footing guarantee, so it is resolved explicitly and checked for
    plausibility (a patient id that is unique per sample means we picked a
    sample-level column by mistake).
    """
    sid = _resolve(samples, _SAMPLE_ID_COLS, "the RNA-seq sample id")
    struct = _resolve(samples, _STRUCTURE_COLS, "the anatomic structure acronym")
    pid = _resolve(samples, _PATIENT_COLS, "the patient/tumour id")

    manifest = pd.DataFrame({
        "sample_id": samples[sid].astype(str),
        "patient_id": samples[pid].astype(str),
        "structure": samples[struct].map(_normalise_structure),
    })

    try:
        block = _resolve(samples, _BLOCK_COLS, "the tissue block name")
        manifest["block"] = samples[block].astype(str)
    except KeyError:
        manifest["block"] = pd.NA

    try:
        sname = _resolve(samples, _STRUCTURE_NAME_COLS, "the structure long name")
        manifest["structure_name"] = samples[sname].astype(str)
    except KeyError:
        manifest["structure_name"] = manifest["structure"].map(
            config.STRUCTURE_LONG_NAMES).fillna(manifest["structure"])

    # The archive holds two studies. `is_anatomic_sample` distinguishes them from the
    # raw label, before canonicalisation collapses it.
    manifest["is_anatomic_study"] = samples[struct].map(config.is_anatomic_sample).to_numpy()
    manifest["is_primary_structure"] = (
        manifest["is_anatomic_study"]
        & manifest["structure"].isin(config.PRIMARY_STRUCTURES)
    )

    if manifest["sample_id"].duplicated().any():
        dupes = manifest.loc[manifest["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"duplicate sample ids in columns-samples.csv: {dupes[:10]}")

    # If every sample has its own "patient", we resolved a sample-level column and the
    # nested design — the entire premise of this project — would be lost.
    if manifest["patient_id"].nunique() == len(manifest):
        raise ValueError(
            f"column {pid!r} gives one distinct patient per sample "
            f"({len(manifest)} of {len(manifest)}). That is a sample-level identifier, "
            f"not a patient identifier — the nested design would be lost. "
            f"Available columns: {list(samples.columns)}"
        )
    return manifest.set_index("sample_id")


def build_bulk(fpkm: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """
    Map Allen gene ids to HGNC symbols, collapse duplicate symbols, and put every
    sample on a common linear scale.

    Duplicate symbols are collapsed by SUM, not mean: FPKM is an additive abundance
    measure, and two probes for the same gene represent parts of the same transcript
    pool. Averaging them would halve that gene's contribution to the mixing model.
    """
    gid = _resolve(genes, _GENE_ID_COLS, "the gene id")
    gsym = _resolve(genes, _GENE_SYMBOL_COLS, "the gene symbol")

    id_to_symbol = dict(zip(genes[gid].astype(str), genes[gsym].astype(str)))
    bulk = fpkm.copy()
    bulk.index = [id_to_symbol.get(str(i), None) for i in bulk.index]

    n_raw = len(bulk)
    bulk = bulk[[i is not None and isinstance(i, str) and i.strip() not in ("", "nan")
                 for i in bulk.index]]
    bulk.index.name = "gene_symbol"
    bulk = bulk.groupby(level=0).sum()

    bulk.columns = [str(c).strip() for c in bulk.columns]
    bulk = bulk.astype("float64")

    # Ivy GAP ships FPKM: already within-sample length- and depth-normalised, but the
    # column sums are not identical across samples. Rescaling each column to a common
    # total puts them on a CPM-comparable footing, which is what the linear mixing
    # model S*w = b assumes. Columns that are entirely zero are left alone rather than
    # producing NaN.
    if config.FPKM_INPUT:
        totals = bulk.sum(axis=0)
        safe = totals.replace(0.0, np.nan)
        bulk = bulk.div(safe, axis=1) * config.RESCALE_COLUMNS_TO
        bulk = bulk.fillna(0.0)

    return bulk, n_raw


def load(write: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, LoadReport]:
    """
    Full load: raw CSVs -> (bulk, manifest, report).

    Bulk columns and manifest rows are intersected and put in identical order, so
    every downstream method receives the same samples in the same order by
    construction rather than by convention.
    """
    config.ensure_dirs()
    fpkm, genes, samples = load_raw()

    manifest = build_manifest(samples)
    bulk, n_genes_raw = build_bulk(fpkm, genes)

    shared = [s for s in bulk.columns if s in manifest.index]
    if not shared:
        raise ValueError(
            "no overlap between fpkm_table.csv columns and columns-samples.csv ids. "
            f"Example expression column: {list(bulk.columns)[:3]}; "
            f"example manifest id: {list(manifest.index)[:3]}"
        )
    bulk = bulk[shared]
    manifest = manifest.loc[shared]

    per_patient = manifest.groupby("patient_id").size()
    ref = config.REFERENCE_STRUCTURE
    patients_with_ref = manifest.loc[
        manifest["is_anatomic_study"] & (manifest["structure"] == ref), "patient_id"
    ].nunique()

    anat = manifest[manifest["is_anatomic_study"]]
    cluster = manifest[~manifest["is_anatomic_study"]]

    report = LoadReport(
        n_genes_raw=int(n_genes_raw),
        n_genes_after_symbol_collapse=int(bulk.shape[0]),
        n_samples=int(bulk.shape[1]),
        n_patients=int(manifest["patient_id"].nunique()),
        structures_found=anat["structure"].value_counts().to_dict(),
        samples_per_patient={
            "min": int(per_patient.min()),
            "median": float(per_patient.median()),
            "max": int(per_patient.max()),
        },
        n_primary_structure_samples=int(manifest["is_primary_structure"].sum()),
        n_anatomic_study_samples=int(len(anat)),
        n_anatomic_study_tumors=int(anat["patient_id"].nunique()),
        n_cluster_study_samples=int(len(cluster)),
        n_cluster_study_tumors=int(cluster["patient_id"].nunique()),
        patients_with_reference_structure=int(patients_with_ref),
        structures_per_tumor={
            str(t): sorted(g["structure"].unique())
            for t, g in anat.groupby("patient_id")
        },
        bulk_sha256=config.sha256_frame(bulk),
        manifest_sha256=config.sha256_frame(
            manifest[["is_primary_structure"]].astype(float)
        ),
    )

    if write:
        bulk.to_csv(config.BULK_EXPRESSION_PATH, sep="\t")
        manifest.to_csv(config.SAMPLE_MANIFEST_PATH, sep="\t")
        (config.PROCESSED_DIR / "load_report.json").write_text(report.to_json())

    return bulk, manifest, report


def load_cached() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the processed tables written by `load()`, without re-parsing the raws."""
    if not config.BULK_EXPRESSION_PATH.exists():
        return load(write=True)[:2]
    bulk = pd.read_csv(config.BULK_EXPRESSION_PATH, sep="\t", index_col=0)
    manifest = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", index_col=0)
    manifest.index = manifest.index.astype(str)
    bulk.columns = [str(c) for c in bulk.columns]
    return bulk, manifest


def main() -> int:
    bulk, manifest, report = load(write=True)
    print(report.to_json())
    print(f"\nwrote {config.BULK_EXPRESSION_PATH}")
    print(f"wrote {config.SAMPLE_MANIFEST_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
