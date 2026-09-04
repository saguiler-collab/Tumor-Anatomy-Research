"""
reference.py — build signature matrices and reference bundles from single-cell data.

WHAT A REFERENCE HAS TO CARRY
-----------------------------
A signature matrix alone is not enough for this project. MuSiC needs the cross-donor
variance of each gene within each cell type; Bisque and SCDC need donor-level
profiles; the cell-size correction needs mean transcript content per cell type. All of
those are properties of the *cells*, and are destroyed by averaging into a signature
matrix. So the reference is built once, in this module, and everything derived from it
comes out of the same pass — which is also what stops two methods being handed subtly
different references.

TWO WAYS IN
-----------
`build_from_h5ad` is the real path: a GBM single-cell atlas (GBmap by default),
annotated with donors and cell types, collapsed onto the project's eight-type roster.

`build_synthetic` generates a reference with known ground truth for tests and for
end-to-end validation before any multi-gigabyte download has happened. It is clearly
labelled as such at every point it touches disk; a synthetic reference must never be
mistaken for a real one in a results table, so its name carries the `synthetic_`
prefix all the way into the manifest.

DONOR BALANCE
-------------
Cells are capped per donor per cell type before averaging. Without a cap, one deeply
sequenced donor supplies most of the cells of a common type and the "reference" becomes
that donor's profile — which then inflates the apparent cross-donor consistency that
MuSiC is trying to measure, exactly backwards.
"""

from __future__ import annotations

import json

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from ivygap import config
from ivygap.deconv.base import ReferenceBundle


def _collapse_to_roster(labels: pd.Series, mapping: dict[str, str]) -> pd.Series:
    """Map fine-grained annotations onto the project roster; unmapped cells become NaN."""
    return labels.map(mapping)


def build_reference(expression: pd.DataFrame, cell_meta: pd.DataFrame, name: str,
                    cells_per_donor_per_type: int = 500,
                    seed: int = config.RANDOM_SEED,
                    cell_totals: pd.Series | None = None) -> ReferenceBundle:
    """
    Core builder: cells -> ReferenceBundle.

    Parameters
    ----------
    expression
        genes x cells, linear scale (CPM-like), non-negative.
    cell_meta
        indexed by cell id, with columns `donor` and `cell_type`. `cell_type` must
        already be collapsed onto `config.CELL_TYPES`.
    """
    rng = np.random.default_rng(seed)

    missing_cols = {"donor", "cell_type"} - set(cell_meta.columns)
    if missing_cols:
        raise ValueError(f"cell_meta is missing column(s): {sorted(missing_cols)}")

    shared = [c for c in expression.columns if c in cell_meta.index]
    if not shared:
        raise ValueError("no cell ids shared between expression and cell_meta")
    expression = expression[shared]
    cell_meta = cell_meta.loc[shared]

    known = cell_meta["cell_type"].isin(config.CELL_TYPES)
    expression = expression.loc[:, known.to_numpy()]
    cell_meta = cell_meta.loc[known]
    if expression.shape[1] == 0:
        raise ValueError(
            f"no cells carry a label in the roster {config.CELL_TYPES}. "
            f"Labels present: {sorted(cell_meta['cell_type'].unique())[:10]}"
        )

    # --- donor-balanced subsampling -----------------------------------------
    keep: list[str] = []
    for (donor, ctype), idx in cell_meta.groupby(["donor", "cell_type"]).groups.items():
        idx = list(idx)
        if len(idx) > cells_per_donor_per_type:
            idx = list(rng.choice(idx, size=cells_per_donor_per_type, replace=False))
        keep.extend(idx)
    expression = expression[keep]
    cell_meta = cell_meta.loc[keep]

    genes = list(expression.index)
    types = list(config.CELL_TYPES)
    donors = sorted(cell_meta["donor"].astype(str).unique())

    # --- per-donor, per-type mean profiles -----------------------------------
    donor_profiles: dict[str, pd.DataFrame] = {}
    for donor in donors:
        cells = cell_meta.index[cell_meta["donor"].astype(str) == donor]
        block = pd.DataFrame(0.0, index=genes, columns=types)
        for ctype in types:
            sel = [c for c in cells if cell_meta.at[c, "cell_type"] == ctype]
            if sel:
                block[ctype] = expression[sel].mean(axis=1)
        donor_profiles[donor] = block

    # --- the signature: mean ACROSS DONORS, not across cells -----------------
    # Averaging across cells would weight donors by how many cells they contributed;
    # averaging the donor means gives each donor one vote, which is the quantity whose
    # variance MuSiC then uses.
    stack = np.stack([donor_profiles[d].to_numpy(dtype="float64") for d in donors])

    # A donor that contributed no cells of a type has a structural zero there, not a
    # measured zero. Masking those out keeps a missing donor from dragging the mean
    # toward zero and from inflating the variance.
    present = np.stack([
        np.array([
            (cell_meta.loc[cell_meta["donor"].astype(str) == d, "cell_type"] == t).any()
            for t in types
        ]) for d in donors
    ])
    mask = present[:, None, :]
    n_present = mask.sum(axis=0)

    with np.errstate(invalid="ignore", divide="ignore"):
        profile = np.where(n_present > 0, (stack * mask).sum(axis=0) / np.maximum(n_present, 1), 0.0)
        centred = (stack - profile[None]) * mask
        sigma = np.where(n_present > 1,
                         (centred ** 2).sum(axis=0) / np.maximum(n_present - 1, 1),
                         0.0)

    profile_df = pd.DataFrame(profile, index=genes, columns=types)
    sigma_df = pd.DataFrame(sigma, index=genes, columns=types)

    # --- cell size: mean total expression per cell of each type --------------
    # This is what converts an RNA proportion into a cell proportion. It must be
    # strictly positive; a type with no cells at all gets the roster median so the
    # division stays defined and the type is not silently amplified to infinity.
    # `cell_totals` is each cell's total counts BEFORE per-cell normalisation. It has to
    # be supplied, because normalising every cell to a common library size makes every
    # type's mean total identical by construction — the cell-size factors then come out
    # all equal and `to_cell_fractions` silently becomes the identity, so RNA
    # proportions are reported as cell proportions. That is exactly the correction this
    # project's own config says it needs: a tumour cell and a lymphocyte do not carry
    # comparable amounts of mRNA.
    totals = (cell_totals.reindex(expression.columns) if cell_totals is not None
              else expression.sum(axis=0))
    sizes = {}
    for ctype in types:
        sel = cell_meta.index[cell_meta["cell_type"] == ctype]
        sizes[ctype] = float(totals[sel].mean()) if len(sel) else np.nan
    size_series = pd.Series(sizes, index=types, dtype="float64")
    if size_series.isna().all():
        raise ValueError("no cell type has any cells; cannot compute cell sizes")
    size_series = size_series.fillna(size_series.median())
    size_series = size_series.replace(0.0, size_series[size_series > 0].median())

    return ReferenceBundle(
        name=name,
        profile=profile_df,
        sigma=sigma_df,
        cell_size=size_series,
        n_donors=len(donors),
        donor_profiles=donor_profiles,
    )


#: Cap on cells taken from any one (donor, cell type) group.
#:
#: A real atlas cannot be densified. Core GBmap is 338,564 cells x ~30,000 genes; as a
#: dense float32 matrix that is roughly 40 TB, and the machine this was written on has
#: 9 GB of RAM. `X.toarray()` on the full matrix is not a slow path, it is an impossible
#: one, and it was what the previous implementation did unconditionally.
#:
#: Subsampling is also what the predecessor project did — `reference_frozen/
#: PROVENANCE.json` records the frozen signature as "100 donors, donor-balanced" — so
#: this follows an established choice rather than inventing one. Balancing by
#: (donor, cell type) matters for more than memory: without it the reference's cell-type
#: means are dominated by whichever donors were sequenced most deeply, and MuSiC's
#: cross-donor variance is estimated from an unbalanced panel.
MAX_CELLS_PER_DONOR_TYPE = 50

#: Overall target, applied after the per-group cap. 20,000 cells against ~16,000 shared
#: genes is about 1.3 GB as float32, which fits alongside the rest of a run.
#:
#: A target rather than a cap: the trim is proportional within cell type and keeps a
#: floor of 25 cells for any type, so a rare population cannot be thinned out of
#: existence by an abundant one. That floor can push the total a little over. The
#: sampling record reports the overshoot rather than hiding it.
MAX_TOTAL_CELLS = 20_000


def _balanced_cell_sample(obs: pd.DataFrame, max_cells_per_donor_type: int,
                          max_total_cells: int,
                          seed: int) -> tuple[list, dict]:
    """
    Choose cells donor-balanced within cell type, deterministically.

    Returns (index labels to keep, a provenance dict). The provenance is written to
    disk because "which cells built the reference" is part of every downstream result,
    and a subsample nobody recorded is a subsample nobody can reproduce.
    """
    rng = np.random.default_rng(seed)

    picked: list = []
    for (_, _), block in obs.groupby(["donor", "cell_type"], observed=True):
        if len(block) <= max_cells_per_donor_type:
            picked.extend(block.index.tolist())
        else:
            picked.extend(rng.choice(block.index.to_numpy(),
                                     size=max_cells_per_donor_type,
                                     replace=False).tolist())

    # The overall cap is applied proportionally WITHIN cell type rather than by a flat
    # random draw, so a rare type is not thinned out of existence by an abundant one.
    # NOTE this is a TARGET, not a hard cap: the per-type floor below can push the
    # total slightly over it, and protecting a rare cell type is worth more than hitting
    # a round number exactly. The provenance records both so the overshoot is visible.
    if max_total_cells and len(picked) > max_total_cells:
        sub = obs.loc[picked]
        frac = max_total_cells / len(picked)
        trimmed: list = []
        for _, block in sub.groupby("cell_type", observed=True):
            # Never drop a type below a floor it can still be estimated from.
            k = max(int(round(len(block) * frac)), min(len(block), 25))
            trimmed.extend(rng.choice(block.index.to_numpy(), size=k,
                                      replace=False).tolist()
                           if k < len(block) else block.index.tolist())
        picked = trimmed

    kept = obs.loc[picked]
    return picked, {
        "strategy": "donor-balanced within cell type, seeded",
        "seed": int(seed),
        "max_cells_per_donor_type": int(max_cells_per_donor_type),
        "target_total_cells": int(max_total_cells) if max_total_cells else None,
        "target_exceeded_by": (max(0, int(len(kept)) - int(max_total_cells))
                               if max_total_cells else 0),
        "per_type_floor": 25,
        "n_cells_kept": int(len(kept)),
        "n_donors_kept": int(kept["donor"].nunique()),
        "cells_per_type": {str(k): int(v)
                           for k, v in kept["cell_type"].value_counts().items()},
        "donors_per_type": {
            str(t): int(b["donor"].nunique())
            for t, b in kept.groupby("cell_type", observed=True)},
        "why": "a real atlas cannot be densified on this hardware, and an unbalanced "
               "panel would let deeply-sequenced donors dominate every cell-type mean",
    }


def build_from_h5ad(path: Path, name: str = "gbmap",
                    annotation_column: str = "annotation_level_3",
                    donor_column: str = "donor_id",
                    mapping: dict[str, str] | None = None,
                    max_cells_per_donor_type: int = MAX_CELLS_PER_DONOR_TYPE,
                    max_total_cells: int = MAX_TOTAL_CELLS,
                    gene_name_column: str | None = "feature_name",
                    restrict_to_genes=None,
                    export: bool = True,
                    seed: int = config.RANDOM_SEED,
                    **kwargs) -> tuple[ReferenceBundle, pd.DataFrame, pd.DataFrame]:
    """
    Build a reference from a single-cell atlas in .h5ad form.

    Returns (bundle, expression, meta) — the same triple as `build_synthetic` — because
    the benchmark needs the CELLS, not just the collapsed bundle: it re-splits donors
    and rebuilds the reference from training donors only, and it pools held-out donors'
    cells into test mixtures. Returning only the bundle would make the donor-held-out
    design impossible to construct downstream.

    Requires `anndata`, which is an optional dependency: the synthetic path and the
    whole benchmark run without it, and there is no reason to force a scanpy stack on
    someone reproducing the fixture results.
    """
    try:
        import h5py
        from anndata.io import read_elem, sparse_dataset
    except ImportError as exc:                            # pragma: no cover
        raise ImportError(
            "reading .h5ad needs `anndata` (pip install anndata). The synthetic "
            "reference path and the full benchmark do not require it."
        ) from exc

    # Read obs, var and X — and NOTHING else.
    #
    # `anndata.read_h5ad(backed="r")` leaves only X on disk and eagerly loads every
    # other element. Core GBmap carries a `layers['scaled']` matrix and a `raw/X`
    # matrix, each the size of X itself, so "backed" mode still tries to pull roughly
    # 16 GB into 9 GB of RAM and simply hangs. Reading the three elements that matter
    # through the low-level API skips all of it.
    with h5py.File(path, "r") as f:
        obs_raw = read_elem(f["obs"])
        var_raw = read_elem(f["var"])

        for col in (annotation_column, donor_column):
            if col not in obs_raw.columns:
                raise KeyError(
                    f"{path.name} has no obs column {col!r}. Available: "
                    f"{list(obs_raw.columns)[:20]}"
                )

        obs = pd.DataFrame({
            "donor": obs_raw[donor_column].astype(str).to_numpy(),
            "cell_type": _collapse_to_roster(
                obs_raw[annotation_column].astype(str),
                mapping or config.GBMAP_CELL_TYPE_MAP,
            ).to_numpy(),
        }, index=obs_raw.index.astype(str))
        obs["_row"] = np.arange(len(obs))
        n_all = len(obs)
        obs = obs.dropna(subset=["cell_type"])
        n_mapped = len(obs)

        keep_rows, sampling = _balanced_cell_sample(
            obs, max_cells_per_donor_type=max_cells_per_donor_type,
            max_total_cells=max_total_cells, seed=seed)
        obs = obs.loc[keep_rows]

        # CELLxGENE atlases are indexed by Ensembl id (ENSG...), while Ivy GAP's bulk
        # matrix — and this project's whole gene space — uses HGNC symbols. Matching the
        # two directly yields an empty intersection, so the symbols are taken from the
        # `feature_name` var column when it exists.
        #
        # This is not a cosmetic rename. Without it `restrict_to_genes` selects nothing
        # and the reference is built on whatever survives, which is why the disjoint
        # case raises below rather than quietly proceeding.
        var_names = var_raw.index.astype(str)
        gene_id_space = "index"
        if gene_name_column and gene_name_column in var_raw.columns:
            symbols = var_raw[gene_name_column].astype(str)
            if symbols.notna().any() and symbols.nunique() > 1:
                var_names = pd.Index(symbols.to_numpy(), dtype=object)
                gene_id_space = gene_name_column

        # A symbol can appear on several Ensembl ids. Keep the first occurrence rather
        # than summing, because summing two rows that a reference chose to keep separate
        # invents an expression value neither of them reported.
        dup_mask = ~pd.Index(var_names).duplicated(keep="first")
        n_dup = int((~dup_mask).sum())

        # Gene restriction BEFORE materialising anything. Handing in the bulk's gene
        # index turns a ~36,000-column read into a ~16,000-column one, and the genes
        # that are dropped could not have contributed to a deconvolution anyway.
        if restrict_to_genes is not None:
            wanted = set(map(str, restrict_to_genes))
            gene_mask = np.array([g in wanted for g in var_names]) & dup_mask
            if not gene_mask.any():
                raise ValueError(
                    f"restrict_to_genes shares no gene with the atlas. The atlas is "
                    f"indexed in the {gene_id_space!r} space, e.g. "
                    f"{list(var_names[:3])}; the requested genes look like "
                    f"{list(map(str, restrict_to_genes))[:3]}. Check symbols vs "
                    f"Ensembl ids.")
        else:
            gene_mask = dup_mask.copy()

        # ALIGNMENT. `sparse_dataset[...]` returns rows in ascending file order, so the
        # metadata has to be put in that same order before its index is used to label
        # the columns. `_balanced_cell_sample` returns cells grouped by (donor, type),
        # not sorted, so using its order here attaches every cell's expression to a
        # DIFFERENT cell's donor and cell type — which scrambles every per-type mean
        # while leaving the column ORDER and every shape check intact.
        #
        # Measured consequence when this was wrong: 23 of 24 marker genes landed in the
        # wrong cell-type column, and the resulting leaderboard put every real method
        # below a random control.
        obs = obs.sort_values("_row")
        rows = obs["_row"].to_numpy()
        assert np.all(np.diff(rows) > 0), "row indices must be strictly increasing"

        X = sparse_dataset(f["X"])[rows]                  # only the chosen cells
        X = X[:, gene_mask]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        X = np.asarray(X, dtype="float32")

    expression = pd.DataFrame(X.T, index=var_names[gene_mask],
                              columns=obs.index.to_numpy())
    del X
    # Each cell's library size, captured BEFORE normalisation. Normalising first and
    # measuring after would make every type's mean total 1e6 and turn the cell-size
    # correction into a no-op.
    raw_totals = expression.sum(axis=0)

    # Normalise each cell to a common total so a deeply sequenced cell does not
    # dominate its type's mean profile.
    totals = raw_totals.replace(0.0, np.nan)
    expression = (expression.div(totals, axis=1) * 1e6).fillna(0.0)

    meta = obs.drop(columns=["_row"])
    sampling.update({
        "source_file": str(path),
        "sha256": config.sha256_file(path),
        "n_cells_in_file": int(n_all),
        "n_cells_after_roster_mapping": int(n_mapped),
        "n_genes_in_file": int(len(var_names)),
        "n_genes_kept": int(gene_mask.sum()),
        "gene_id_space": gene_id_space,
        "n_duplicate_symbols_dropped": n_dup,
        "annotation_column": annotation_column,
        "donor_column": donor_column,
    })

    ref = build_reference(expression, meta, name=name, cell_totals=raw_totals, **kwargs)
    # ReferenceBundle is frozen by design, so the sampling record lives on disk rather
    # than on the object. That is the right place for it anyway: "which cells built this
    # reference" has to survive the process that built them.
    (config.REFERENCE_DIR / f"reference_sampling_{name}.json").write_text(
        json.dumps(sampling, indent=2, default=str))
    if export:
        export_cell_level(expression, meta, name)
    return ref, expression, meta


def export_cell_level(expression: pd.DataFrame, meta: pd.DataFrame, name: str) -> None:
    """
    Write the cell-level export the real R packages need.

    MuSiC, Bisque and SCDC compute their own cross-donor statistics and cannot work
    from a collapsed signature matrix, so without this export `r_bridge` correctly
    reports itself unavailable and the run falls back to Python.
    """
    config.ensure_dirs()
    counts_path = config.REFERENCE_DIR / f"sc_counts_{name}.csv.gz"
    meta_path = config.REFERENCE_DIR / f"sc_meta_{name}.csv"
    # Write-then-rename: an interrupted export otherwise leaves a truncated file that
    # `r_bridge.check()` will read as an available reference.
    tmp_path = counts_path.with_suffix(counts_path.suffix + ".part")
    with gzip.open(tmp_path, "wt") as fh:
        expression.to_csv(fh)
    tmp_path.replace(counts_path)
    meta.to_csv(meta_path)
    print(f"  wrote cell-level export: {counts_path.name}, {meta_path.name}")


def select_signature_genes(ref: ReferenceBundle,
                           n_per_type: int = config.SIGNATURE_GENES_PER_TYPE) -> list[str]:
    """
    Pick the genes that actually distinguish cell types.

    Score = a type's mean expression divided by the mean across all other types, so a
    gene that is high everywhere scores low. Marker genes from the config are always
    included: they are the biological sanity check that makes a signature matrix
    inspectable, and losing GFAP or PECAM1 to a scoring quirk would make the reference
    much harder to argue about.
    """
    profile = ref.profile
    scores = {}
    for ctype in profile.columns:
        others = profile.drop(columns=[ctype]).mean(axis=1)
        scores[ctype] = profile[ctype] / (others + 1.0)

    chosen: list[str] = []
    for ctype, s in scores.items():
        chosen.extend(s.sort_values(ascending=False).head(n_per_type).index.tolist())
    for markers in config.MARKER_GENES.values():
        chosen.extend(g for g in markers if g in profile.index)

    seen, out = set(), []
    for g in chosen:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def frozen_gene_space(frozen: ReferenceBundle, bulk_index) -> tuple[list[str], dict]:
    """
    The gene space a run uses when the benchmark stage is skipped.

    WHY THIS IS A FUNCTION AND NOT THREE LINES IN run_all.py
    --------------------------------------------------------
    It used to be three lines in run_all.py, and they ignored
    `config.USE_SIGNATURE_GENE_SUBSET`. The switch was applied only inside
    `run_benchmark`, which does not run without a cell-level single-cell reference — so
    on exactly the configuration this project ships in, the flag was declared, documented
    as load-bearing, and never read. Every method received the full 16,007-gene
    intersection.

    config.py already states both halves of why that is wrong:

      * least squares is dominated by the largest residuals, which belong to the most
        highly expressed genes — mostly housekeeping genes similar in every cell type,
        carrying no information about composition;
      * nu-SVR is superlinear in the number of rows, and on the full intersection a
        single cohort takes hours.

    Both were measured on this machine rather than assumed: one NuSVR fit on 16,007
    genes takes 21 s at nu=0.25 and 44 s at nu=0.75, all three hitting the 200,000
    iteration cap, against 366 fits per cohort.

    Returned alongside the genes is a provenance dict, because "which genes did every
    method see" is part of the equal-footing claim and belongs in the artefacts rather
    than in a print statement.
    """
    # Hoisted: `g in set(bulk_index)` inside the comprehension rebuilds the set on
    # every one of ~18,000 iterations.
    bulk_genes = set(bulk_index)
    shared = [g for g in frozen.profile.index if g in bulk_genes]
    if not config.USE_SIGNATURE_GENE_SUBSET:
        return shared, {
            "strategy": "full intersection",
            "n_shared": len(shared),
            "n_selected": len(shared),
            "why": "config.USE_SIGNATURE_GENE_SUBSET is False",
        }

    # Rank within the SHARED space, not the whole frozen profile. Ranking first and
    # intersecting afterwards silently discards top-scoring genes the bulk does not
    # carry, so the subset ends up smaller than requested for a reason that has nothing
    # to do with how informative the genes are.
    selected = select_signature_genes(frozen.subset_genes(shared),
                                      n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    shared_set = set(shared)
    genes = [g for g in selected if g in shared_set]
    return genes, {
        "strategy": f"top {config.SIGNATURE_GENES_PER_TYPE} per cell type, "
                    f"plus every marker gene, ranked within the shared space",
        "n_shared": len(shared),
        "n_selected": len(genes),
        "why": "least squares is dominated by highly expressed housekeeping genes that "
               "carry no compositional information; and nu-SVR is superlinear in rows",
    }


def build_synthetic(n_genes: int = 1200, n_donors: int = 8,
                    cells_per_type_per_donor: int = 40,
                    seed: int = config.RANDOM_SEED,
                    name: str = "synthetic_a") -> tuple[ReferenceBundle, pd.DataFrame, pd.DataFrame]:
    """
    A synthetic single-cell reference with a known generative structure.

    Used by the tests and by the fixture end-to-end run. Every cell type gets a block
    of genes it expresses highly, plus a shared housekeeping background, plus
    donor-level offsets — so cross-donor variance is real and MuSiC's weighting has
    something genuine to find rather than fitting noise.

    Returns (bundle, expression, meta) so callers can also build matched pseudobulk.
    """
    rng = np.random.default_rng(seed)
    types = list(config.CELL_TYPES)
    genes = [f"GENE{i:05d}" for i in range(n_genes)]

    # Give the real marker symbols to the first genes of each type's block, so the
    # fixture reference is inspectable with the same biology checks as a real one.
    block = n_genes // (len(types) + 1)
    base = np.full((n_genes, len(types)), 5.0)
    for k, ctype in enumerate(types):
        lo, hi = k * block, (k + 1) * block
        base[lo:hi, k] += rng.gamma(shape=3.0, scale=40.0, size=hi - lo)
        for j, marker in enumerate(config.MARKER_GENES[ctype]):
            if lo + j < n_genes:
                genes[lo + j] = marker
    housekeeping = rng.gamma(shape=2.0, scale=10.0, size=n_genes)
    base += housekeeping[:, None]

    # Cell types differ in transcript content; that is what makes the cell-size
    # correction non-trivial rather than a no-op.
    size_factor = pd.Series(
        rng.uniform(0.6, 1.8, size=len(types)), index=types).to_dict()

    rows, meta_rows = {}, {}
    for d in range(n_donors):
        donor = f"D{d:02d}"
        donor_shift = rng.lognormal(mean=0.0, sigma=0.25, size=(n_genes, len(types)))
        for k, ctype in enumerate(types):
            mu = base[:, k] * donor_shift[:, k] * size_factor[ctype]
            for c in range(cells_per_type_per_donor):
                cid = f"{donor}_{ctype}_{c}"
                rows[cid] = rng.poisson(np.clip(mu, 0, None)).astype("float64")
                meta_rows[cid] = {"donor": donor, "cell_type": ctype}

    expression = pd.DataFrame(rows, index=genes)
    meta = pd.DataFrame.from_dict(meta_rows, orient="index")
    bundle = build_reference(expression, meta, name=name, seed=seed)
    return bundle, expression, meta


def load_frozen_reference(name: str = "gbmap_frozen") -> ReferenceBundle:
    """
    Load the signature matrix vendored in `reference_frozen/`.

    This is protocol step 1 — "restore the frozen inputs" — and it is what lets the
    pipeline run on real Ivy GAP tissue without first rebuilding a reference from a
    multi-gigabyte single-cell atlas.

    WHAT THIS REFERENCE CANNOT SUPPORT, AND WHY THAT IS RECORDED RATHER THAN PATCHED
    -------------------------------------------------------------------------------
    A signature matrix is a collapsed object: it is the mean profile per cell type, and
    the per-donor spread that produced it has been averaged away. Two methods in the
    roster need exactly that spread.

    MuSiC weights each gene by its cross-donor variance within a cell type — that is the
    whole method. Bisque estimates its assay transform from donor-level profiles. With a
    collapsed reference, `sigma` is unavailable.

    Filling sigma with zeros and saying nothing would be the tempting move, and it would
    silently turn "MuSiC" into "residual-weighted NNLS" while still printing the name
    MuSiC on a leaderboard. Instead sigma is zeroed AND `has_cross_donor_variance` is set
    false, which the methods check and report. A degenerate method is reported as
    degenerate.
    """
    if not config.FROZEN_SIGNATURE_PATH.exists():
        raise FileNotFoundError(
            f"{config.FROZEN_SIGNATURE_PATH} not found. It should be committed to the "
            "repository; see reference_frozen/PROVENANCE.json."
        )

    profile = pd.read_csv(config.FROZEN_SIGNATURE_PATH, sep="\t", index_col=0)
    missing = [c for c in config.CELL_TYPES if c not in profile.columns]
    if missing:
        raise ValueError(
            f"the frozen signature is missing cell type(s) {missing}; it has "
            f"{list(profile.columns)}"
        )
    profile = profile[list(config.CELL_TYPES)].astype("float64")
    profile = profile[~profile.index.duplicated(keep="first")]

    sizes = pd.read_csv(config.FROZEN_CELL_SIZE_PATH, index_col=0)
    # L_median matches the source project's frozen CELL_SIZE_DIVISOR = "median".
    # Picking a different column here would silently change the estimand.
    col = "L_median" if "L_median" in sizes.columns else sizes.columns[0]
    cell_size = sizes[col].reindex(list(config.CELL_TYPES)).astype("float64")
    if cell_size.isna().any():
        raise ValueError(
            f"cell size factors are missing for "
            f"{cell_size[cell_size.isna()].index.tolist()}"
        )

    n_donors = int(sizes["n_donors"].max()) if "n_donors" in sizes.columns else 0

    bundle = ReferenceBundle(
        name=name,
        profile=profile,
        sigma=pd.DataFrame(0.0, index=profile.index, columns=profile.columns),
        cell_size=cell_size,
        n_donors=n_donors,
        donor_profiles=None,
    )
    object.__setattr__(bundle, "has_cross_donor_variance", False)
    return bundle
