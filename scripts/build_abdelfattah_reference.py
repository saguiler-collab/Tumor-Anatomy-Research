"""Build an INDEPENDENT reference from Abdelfattah et al. 2022 (GSE182109), streaming -- prespecified/abdelfattah_cluster_mapping.md.

Abdelfattah 2022 is not one of GBmap's 16 source studies (GBmap `author` field, checked 2026-10-01)
and is lymphocyte-rich, so it can test whether the lymphoid inversion is a property of GBmap's
lymphoid profiles. This script only BUILDS the reference; scripts/independent_atlas_test.py refits.

INPUT, as obtained (2026-10-01): the Single Cell Portal files -- Processed_matrix.mtx.gz (38,224
genes x 201,986 cells, 542.7 M nonzeros), genes.tsv, barcodes.tsv, Cluster_GBM.txt (authors'
cluster ids C1-C12) -- plus gsm_to_sample.tsv from GEO's file list.

THE MATRIX IS LOG-SCALE, and is NOT used as such (the D16 lesson). Its values are Seurat's
log1p(count / library x 1e4). Inverting is exact: within a cell the smallest nonzero expm1 value is
one count, so counts = expm1(v) / min(expm1(v)) and library = 1e4 / min(expm1(v)). Verified on
1,220 cells before this script was written: max relative deviation from integers 1.6e-15. Every
cell is checked again here, and the recovered counts must also sum to the recovered library.

The matrix is too large to load (8 GB RAM), so it is streamed twice, in chunks, relying on the
file's column-major order (checked): pass 1 recovers library sizes and the marker-panel counts
needed to name clusters; pass 2 collects only the selected cells.

Cluster naming follows the PRE-SPECIFIED precedence: the authors' labels if a metadata file is
present, else the fixed marker rule (margin 0.5 on mean log1p-CPM). Construction follows this
project's builder: CPM per cell, at most 50 cells per (patient, type), seeded.
Writes data/reference/abdelfattah_2022/{profile,sigma,cell_size}.csv + provenance.json; with
--seeds, other seeds write abdelfattah_2022_seed<k>/ from the same streaming pass.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                         # noqa: E402

SRC = config.GSE182109_DIR      # never a literal: the folder has been renamed once (D25 audit)
OUT = config.REFERENCE_DIR / "abdelfattah_2022"
PANELS = {   # fixed in prespecified/abdelfattah_cluster_mapping.md before any expression was seen
    "Tumor": ["SOX2", "PTPRZ1", "EGFR", "NES", "PDGFRA"],
    "Macrophage_Microglia": ["CD14", "CD68", "AIF1", "C1QA", "CSF1R"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "FLT1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "MAG"],
    "Astrocyte": ["AQP4", "GJA1", "ALDH1L1"],
    "_excluded": ["RGS5", "PDGFRB", "JCHAIN", "MZB1", "TPSAB1", "CPA3", "CLEC9A", "FCER1A"],
}
MARGIN, CAP, CHUNK = 0.5, 50, 8_000_000


def matrix_path() -> Path:
    for p in (SRC / "Processed_matrix.mtx.gz", SRC / "Processed_matrix.mtx.gz.download" / "Processed_matrix.mtx.gz",
              SRC / "Processed_matrix.mtx", SRC / "matrix.mtx"):
        if p.exists() and p.stat().st_size > 0:
            return p
    raise SystemExit(f"BLOCKED: no Processed_matrix in {SRC}/")


def stream(path: Path):
    """Yield (gene_idx0, cell_idx0, value) arrays per chunk; header parsed first."""
    fh = gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    dims = None
    for line in fh:
        if not line.startswith("%"):
            dims = tuple(int(x) for x in line.split()); break
    yield dims
    for ch in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=CHUNK,
                          dtype={"g": np.int32, "c": np.int32, "v": np.float64}):
        yield ch.g.to_numpy() - 1, ch.c.to_numpy() - 1, ch.v.to_numpy()


def per_cell_units(g, c, v, carry):
    """Join the previous chunk's partial last cell, return complete-cell arrays and the new carry."""
    if carry is not None:
        g, c, v = np.concatenate([carry[0], g]), np.concatenate([carry[1], c]), np.concatenate([carry[2], v])
    last = c[-1]
    done = c != last
    return (g[done], c[done], v[done]), (g[~done], c[~done], v[~done])


def out_dir(seed: int) -> Path:
    """The registered build keeps its path; a seed-sensitivity build gets its own."""
    return OUT if seed == config.RANDOM_SEED else OUT.parent / f"{OUT.name}_seed{seed}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seeds", type=int, nargs="+", default=[config.RANDOM_SEED],
                    help="cell-sampling seeds; all share one streaming pass. The registered "
                         "build is config.RANDOM_SEED; others write to abdelfattah_2022_seed<k>/")
    args = ap.parse_args()
    genes = [l.rstrip("\n").split("\t")[-1] for l in open(SRC / "genes.tsv")]
    cells = [l.strip() for l in open(SRC / "barcodes.tsv")]
    cl = pd.read_csv(SRC / "Cluster_GBM.txt", sep=",", skiprows=[1], header=None,
                     names=["NAME", "N2", "X", "Y", "Category"], low_memory=False)
    cl = cl[cl.NAME != "NAME"]
    assert list(cl.NAME) == cells, "barcodes.tsv and Cluster_GBM.txt disagree on cell order"
    clusters = cl.Category.to_numpy()
    gsm2s = dict(l.rstrip("\n").split("\t") for l in open(SRC / "gsm_to_sample.tsv"))
    patient = np.array(["-".join(gsm2s[n.split("_")[0]].split("-")[:2]) for n in cells])
    meta_files = [p for p in SRC.iterdir() if re.search(r"meta", p.name, re.I) and p.is_file()]
    if meta_files:
        raise SystemExit(f"authors' metadata present ({meta_files[0].name}): the pre-specified precedence "
                         f"puts their labels first -- inspect and wire it before running the marker rule")
    panel_genes = sorted({g for p in PANELS.values() for g in p if g in set(genes)})
    pidx = {genes.index(g): k for k, g in enumerate(panel_genes)}
    is_panel = np.zeros(len(genes), bool); is_panel[list(pidx)] = True
    n_cells = len(cells)
    lib = np.full(n_cells, np.nan); total = np.zeros(n_cells); maxdev = np.zeros(n_cells)
    panel_counts = np.zeros((n_cells, len(panel_genes)), np.float32)

    # ---------------- pass 1 ----------------
    path = matrix_path(); print(f"streaming {path.name} (pass 1)", flush=True)
    it = stream(path); dims = next(it)
    assert dims[0] == len(genes) and dims[1] == n_cells, f"matrix dims {dims} vs files"
    carry, seen = None, 0
    for g, c, v in it:
        seen += len(v)
        (g, c, v), carry = per_cell_units(g, c, v, carry)
        if len(c) == 0:
            continue
        e = np.expm1(v)
        starts = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
        ids = c[starts]; units = np.minimum.reduceat(e, starts)
        reps = np.repeat(units, np.diff(np.r_[starts, len(c)]))
        cnt = e / reps
        lib[ids] = 1e4 / units
        total[ids] = np.add.reduceat(cnt, starts)
        maxdev[ids] = np.maximum.reduceat(np.abs(cnt - np.round(cnt)) / np.maximum(cnt, 1), starts)
        m = is_panel[g]
        panel_counts[c[m], [pidx[x] for x in g[m]]] = np.round(cnt[m])
        print(f"  {seen / dims[2]:6.1%}", flush=True) if (seen // CHUNK) % 10 == 0 else None
    if carry is not None and len(carry[1]):                      # the final cell
        g, c, v = carry; e = np.expm1(v); u = e.min(); cid = c[0]; cnt = e / u
        lib[cid] = 1e4 / u; total[cid] = cnt.sum()
        maxdev[cid] = float((np.abs(cnt - np.round(cnt)) / np.maximum(cnt, 1)).max())
        m = is_panel[g]; panel_counts[cid, [pidx[x] for x in g[m]]] = np.round(cnt[m])
    assert seen == dims[2], f"read {seen:,} entries, header says {dims[2]:,}"
    ok_int = bool(np.nanmax(maxdev) < 1e-6)
    lib_err = float(np.nanmax(np.abs(total - lib) / lib))
    print(f"pass 1: {np.isfinite(lib).sum():,} cells; integral counts {ok_int} (max rel dev {np.nanmax(maxdev):.1e}); "
          f"recovered totals vs library max rel error {lib_err:.1e}", flush=True)
    if not ok_int or lib_err > 1e-6:
        raise SystemExit("BLOCKED: the log-normalisation identity does not hold for every cell")

    # ---------------- name the clusters (pre-declared rule) ----------------
    lcpm = np.log1p(panel_counts / lib[:, None] * 1e6)
    col = {g: k for k, g in enumerate(panel_genes)}
    mapping, scores = {}, {}
    for cid in sorted(set(clusters), key=lambda s: int(s[1:])):
        rows = clusters == cid
        sc = {t: float(lcpm[rows][:, [col[g] for g in p if g in col]].mean()) for t, p in PANELS.items()}
        best = sorted(sc.items(), key=lambda kv: -kv[1])
        scores[cid] = {t: round(x, 3) for t, x in sc.items()}
        if best[0][0] in ("T_cell", "NK_cell"):
            # OPEN_DEFECTS D28: a T/NK cluster is judged against the best panel OUTSIDE {T, NK}; with
            # the plain runner-up margin the declared T/NK split below could never trigger. No cluster
            # of this atlas changes (C3: |T - NK| = 3.64; no other cluster is T- or NK-best).
            other = max(v for k, v in sc.items() if k not in ("T_cell", "NK_cell"))
            good = max(sc["T_cell"], sc["NK_cell"]) - other >= MARGIN
        else:
            good = best[0][0] != "_excluded" and best[0][1] - best[1][1] >= MARGIN
        mapping[cid] = best[0][0] if good else None
        print(f"  {cid:<4} n={rows.sum():>6}  best {best[0][0]:<20} {best[0][1]:.2f}  runner-up {best[1][0]:<20} "
              f"{best[1][1]:.2f}  -> {mapping[cid]}", flush=True)
    label = np.array([mapping[x] for x in clusters], dtype=object)
    tcols = [col[g] for g in PANELS["T_cell"] if g in col]; ncols = [col[g] for g in PANELS["NK_cell"] if g in col]
    tnk = [cid for cid, t in mapping.items() if t in ("T_cell", "NK_cell")
           and abs(scores[cid]["T_cell"] - scores[cid]["NK_cell"]) < MARGIN]
    for cid in tnk:
        rows = np.flatnonzero(clusters == cid)
        d = lcpm[rows][:, tcols].mean(1) - lcpm[rows][:, ncols].mean(1)
        label[rows] = np.where(d >= 0, "T_cell", "NK_cell")

    # ---------------- select cells: cap per (patient, type), seeded ----------------
    df = pd.DataFrame({"type": label, "patient": patient}); df = df[df.type.notna()]
    takes = {}
    for seed in args.seeds:
        rng = np.random.default_rng(seed)
        take = []
        for _, grp in df.groupby(["patient", "type"]):
            ix = grp.index.to_numpy(); take += list(rng.choice(ix, CAP, replace=False) if len(ix) > CAP else ix)
        takes[seed] = np.sort(np.array(take))
        print(f"seed {seed}: selected {len(takes[seed]):,} cells across "
              f"{df.loc[takes[seed]].patient.nunique()} patients", flush=True)
    want = np.zeros(n_cells, bool)
    for take in takes.values():
        want[take] = True

    # ---------------- pass 2: collect selected cells ----------------
    print("streaming (pass 2)", flush=True)
    it = stream(path); next(it)
    rows_g, rows_c, vals = [], [], []
    for g, c, v in it:
        m = want[c]
        if m.any():
            rows_g.append(g[m]); rows_c.append(c[m]); vals.append(v[m])
    g = np.concatenate(rows_g); c = np.concatenate(rows_c)
    cnt = np.round(np.expm1(np.concatenate(vals)) * lib[c] / 1e4)
    cpm = cnt / lib[c] * 1e6
    for seed, take in takes.items():
        write_build(seed, take, df, g, c, cpm, lib, genes, patient, path, dims, maxdev, lib_err,
                    mapping, scores, tnk)
    return 0


def write_build(seed, take, df, g, c, cpm, lib, genes, patient, path, dims, maxdev, lib_err,
                mapping, scores, tnk) -> None:
    OUTS = out_dir(seed)
    types = [t for t in config.CELL_TYPES if t in set(df.loc[take, "type"])]
    prof, sig, size, ncell, npat = {}, {}, {}, {}, {}
    for t in types:
        tc = take[df.loc[take, "type"].to_numpy() == t]; n = len(tc)
        m = np.isin(c, tc)
        s1 = np.bincount(g[m], weights=cpm[m], minlength=len(genes))
        s2 = np.bincount(g[m], weights=cpm[m] ** 2, minlength=len(genes))
        mu = s1 / n
        prof[t], sig[t] = mu, np.clip(s2 / n - mu ** 2, 0, None)
        size[t], ncell[t] = float(lib[tc].mean()), int(n)
        npat[t] = int(pd.Series(patient[tc]).nunique())
    OUTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(prof, index=genes).groupby(level=0).sum().to_csv(OUTS / "profile.csv")
    pd.DataFrame(sig, index=genes).groupby(level=0).sum().to_csv(OUTS / "sigma.csv")
    pd.Series(size).rename("cell_size").to_csv(OUTS / "cell_size.csv")
    prov = {"source": "Abdelfattah et al. 2022 Nat Commun 13:767; GEO GSE182109; Single Cell Portal processed files",
            "independent_of_gbmap": "not among GBmap's 16 source studies (author field), checked 2026-10-01",
            "matrix": {"file": path.name, "dims": dims, "scale": "Seurat log1p(count/library x 1e4)",
                       "inverted_to_counts": True, "max_rel_dev_from_integers": float(np.nanmax(maxdev)),
                       "max_rel_error_totals_vs_library": lib_err},
            "labelling": "pre-declared marker rule (prespecified/abdelfattah_cluster_mapping.md); no authors' metadata file",
            "cluster_mapping": mapping, "cluster_panel_scores": scores, "t_nk_split_clusters": tnk,
            "types": types, "cells_per_type": ncell, "patients_per_type": npat,
            "cap_per_patient_type": CAP, "seed": seed,
            "n_patients": int(len(set(patient)))}
    (OUTS / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"seed {seed} -> {OUTS}")
    print(json.dumps({k: prov[k] for k in ("cluster_mapping", "cells_per_type", "patients_per_type")}, indent=1))


if __name__ == "__main__":
    raise SystemExit(main())
