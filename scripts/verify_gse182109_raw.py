"""
verify_gse182109_raw.py -- does the processed Abdelfattah matrix derive from GEO's raw counts?

The independent atlas (scripts/build_abdelfattah_reference.py) was built from processed files
obtained from the Single Cell Portal, whose study accession was never recorded. GEO GSE182109's
own supplementary archive -- one Cell Ranger output per sample (matrix.mtx.gz, features.tsv.gz,
barcodes.tsv.gz) -- was added on 2026-10-01 as GSE182109/GSE182109_RAW/. This script examines all
of it and tests the processed matrix against it, so the reference's provenance rests on a GEO
accession rather than on an unrecorded portal page.

  A. Every one of the 44 samples: the file triplet is complete; the sample name in the file name
     matches GEO's GSM -> sample table; features, barcodes and the matrix header agree in
     dimension; Cell Ranger version; sha256 of every file.
  B. Every processed cell (all 201,986): its barcode exists among its own sample's raw barcodes.
  C. The processed gene list against the raw feature lists (Seurat makes duplicated symbols
     unique with .1, .2 suffixes; the same rule is applied here before matching).
  D. Exactness, on the first --n-samples samples (the processed matrix is column-major in sample
     order, so their cells are streamed without reading the rest): invert Seurat's
     log1p(count / library x 1e4) exactly as the builder does and require every recovered
     count to EQUAL the raw count for the same cell and gene, and the recovered library to equal
     the raw UMI total over the genes the processed object kept.

    python3 scripts/verify_gse182109_raw.py      # writes results/gse182109_raw_verification.json
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402

SRC = config.GSE182109_DIR
RAW = SRC / "GSE182109_RAW"
OUT = config.RESULTS_DIR / "gse182109_raw_verification.json"
CHUNK = 8_000_000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def mtx_header(path: Path) -> tuple[tuple[int, int, int], str | None, str]:
    meta, field = None, ""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("%%MatrixMarket"):
                field = line.split()[3] if len(line.split()) > 3 else ""
            elif line.startswith("%metadata_json"):
                try:
                    meta = json.loads(line.split(":", 1)[1]).get("software_version")
                except Exception:                                   # noqa: BLE001
                    meta = None
            elif not line.startswith("%"):
                return tuple(int(x) for x in line.split()), meta, field
    raise ValueError(f"no header in {path}")


def make_unique(names: list[str]) -> list[str]:
    """Seurat/R make.unique: second and later copies of a name get .1, .2, ..."""
    seen: dict[str, int] = {}
    out = []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}.{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--n-samples", type=int, default=3, help="samples checked for exactness (part D)")
    a = ap.parse_args()
    if not RAW.is_dir():
        print(f"BLOCKED: {RAW} not found"); return 2

    gsm2s = dict(l.rstrip("\n").split("\t") for l in open(SRC / "gsm_to_sample.tsv"))
    files = sorted(RAW.iterdir())
    by_gsm: dict[str, dict[str, Path]] = defaultdict(dict)
    name_in_file: dict[str, str] = {}
    unexpected = []
    for f in files:
        m = re.match(r"^(GSM\d+)_(.+)_(matrix\.mtx|features\.tsv|barcodes\.tsv)\.gz$", f.name)
        if not m:
            unexpected.append(f.name); continue
        by_gsm[m.group(1)][m.group(3).split(".")[0]] = f
        name_in_file[m.group(1)] = m.group(2)

    # ---------------- A: every sample ----------------
    samples, ref_sets = {}, Counter()
    feature_lists: dict[str, list[str]] = {}
    for gsm in sorted(gsm2s):
        trip = by_gsm.get(gsm, {})
        rec = {"sample": gsm2s[gsm], "name_in_file": name_in_file.get(gsm),
               "complete_triplet": set(trip) == {"matrix", "features", "barcodes"}}
        rec["name_matches_geo_table"] = rec["name_in_file"] == gsm2s[gsm]
        if rec["complete_triplet"]:
            feats = pd.read_csv(trip["features"], sep="\t", header=None, compression="gzip")
            n_bc = sum(1 for _ in gzip.open(trip["barcodes"], "rt"))
            (nr, nc, nnz), sw, field = mtx_header(trip["matrix"])
            fid = hashlib.sha256("\n".join(feats[0].astype(str)).encode()).hexdigest()[:16]
            ref_sets[fid] += 1
            feature_lists.setdefault(fid, list(feats[1].astype(str)))
            rec.update({"n_features": int(len(feats)), "feature_list_id": fid,
                        "n_barcodes": int(n_bc), "matrix_dims": [nr, nc], "nnz": nnz,
                        "value_field": field, "cellranger": sw,
                        "dims_agree": nr == len(feats) and nc == n_bc,
                        "bytes": {k: v.stat().st_size for k, v in trip.items()},
                        "sha256": {k: sha256(v) for k, v in trip.items()}})
        samples[gsm] = rec
        print(f"A {gsm} {rec['sample']:<12} triplet {rec['complete_triplet']} "
              f"name {rec['name_matches_geo_table']} {rec.get('n_barcodes', '?'):>6} barcodes "
              f"{rec.get('n_features', '?')} features {rec.get('cellranger', '?')}", flush=True)

    # ---------------- B: every processed cell ----------------
    proc_bc = [l.strip() for l in open(SRC / "barcodes.tsv")]
    per_gsm_proc: dict[str, list[str]] = defaultdict(list)
    for b in proc_bc:
        g, bc = b.split("_", 1)
        per_gsm_proc[g].append(bc)
    cell_check = {}
    for gsm, bcs in per_gsm_proc.items():
        raw_bc = set(l.strip() for l in gzip.open(by_gsm[gsm]["barcodes"], "rt"))
        found = sum(bc in raw_bc for bc in bcs)
        cell_check[gsm] = {"n_processed": len(bcs), "n_found_in_raw": found,
                           "n_raw": len(raw_bc), "retained_fraction": round(len(bcs) / len(raw_bc), 4)}
    n_found = sum(v["n_found_in_raw"] for v in cell_check.values())
    print(f"B processed cells found among their sample's raw barcodes: {n_found:,} / {len(proc_bc):,}")

    # ---------------- C: genes ----------------
    proc_genes = [l.rstrip("\n").split("\t")[0] for l in open(SRC / "genes.tsv")]
    union, raw_unique = [], {}
    for fid, syms in feature_lists.items():
        raw_unique[fid] = make_unique(syms)
    union_set = set().union(*map(set, raw_unique.values()))
    gene_check = {"n_processed_genes": len(proc_genes), "n_distinct_feature_lists": len(ref_sets),
                  "feature_list_usage": dict(ref_sets),
                  "n_raw_union_after_make_unique": len(union_set),
                  "processed_genes_in_raw_union": sum(g in union_set for g in proc_genes),
                  "processed_genes_not_in_raw": [g for g in proc_genes if g not in union_set][:25]}
    print(f"C genes: processed {len(proc_genes):,}; raw feature lists {len(ref_sets)} "
          f"(usage {dict(ref_sets)}); processed genes found in raw union "
          f"{gene_check['processed_genes_in_raw_union']:,}", flush=True)

    # ---------------- D: exactness on the first n samples ----------------
    order = list(dict.fromkeys(b.split("_", 1)[0] for b in proc_bc))
    first = order[: a.n_samples]
    n_cells_first = sum(len(per_gsm_proc[g]) for g in first)
    path = SRC / "Processed_matrix.mtx.gz.download" / "Processed_matrix.mtx.gz"
    if not path.exists():
        path = SRC / "Processed_matrix.mtx.gz"
    rows_g, rows_c, vals = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for ch in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=CHUNK,
                              dtype={"g": np.int32, "c": np.int32, "v": np.float64}):
            m = ch.c.to_numpy() <= n_cells_first
            rows_g.append(ch.g.to_numpy()[m] - 1); rows_c.append(ch.c.to_numpy()[m] - 1)
            vals.append(ch.v.to_numpy()[m])
            if not m[-1]:                     # column-major: past the last wanted cell
                break
    g_ = np.concatenate(rows_g); c_ = np.concatenate(rows_c); v_ = np.expm1(np.concatenate(vals))
    exact = {}
    offset = 0
    for gsm in first:
        bcs = per_gsm_proc[gsm]
        sel = (c_ >= offset) & (c_ < offset + len(bcs))
        cg, cc, cv = g_[sel], c_[sel] - offset, v_[sel]
        mins = np.full(len(bcs), np.inf)
        np.minimum.at(mins, cc, cv)
        rec_counts = cv / mins[cc]                       # the builder's inversion
        lib = 1e4 / mins                                 # Seurat's library size, recovered
        trip = by_gsm[gsm]
        X = sio.mmread(trip["matrix"]).tocsc()           # features x barcodes, integer counts
        fid = samples[gsm]["feature_list_id"]
        raw_names = raw_unique[fid]
        raw_bc = [l.strip() for l in gzip.open(trip["barcodes"], "rt")]
        col_of = {b: i for i, b in enumerate(raw_bc)}
        row_of = {n: i for i, n in enumerate(raw_names)}
        p2r = np.array([row_of.get(n, -1) for n in proc_genes])
        kept = np.flatnonzero(p2r >= 0)
        sub = X[:, [col_of[b] for b in bcs]].tocsc()
        raw_vals = np.asarray(sub[p2r[cg].clip(min=0), cc]).ravel()
        mapped = p2r[cg] >= 0
        diff = np.abs(rec_counts[mapped] - raw_vals[mapped])
        totals_all = np.asarray(sub.sum(axis=0)).ravel()
        totals_kept = np.asarray(sub[p2r[kept], :].sum(axis=0)).ravel()
        n_raw_nz_kept = int(sub[p2r[kept], :].count_nonzero())
        exact[gsm] = {
            "sample": gsm2s[gsm], "n_cells": len(bcs), "n_processed_nonzeros": int(len(cv)),
            "n_raw_nonzeros_on_kept_genes": n_raw_nz_kept,
            "nonzeros_match_count": int(len(cv)) == n_raw_nz_kept,
            "max_abs_count_difference": float(diff.max()) if diff.size else None,
            "fraction_entries_exact": float(np.mean(diff < 1e-6)) if diff.size else None,
            "max_rel_library_vs_raw_total_kept_genes": float(np.max(np.abs(lib - totals_kept) / totals_kept)),
            "max_rel_library_vs_raw_total_all_genes": float(np.max(np.abs(lib - totals_all) / totals_all)),
        }
        print(f"D {gsm} {gsm2s[gsm]}: {len(bcs):,} cells, {len(cv):,} nonzeros; max |count diff| "
              f"{exact[gsm]['max_abs_count_difference']:.3g}; exact {exact[gsm]['fraction_entries_exact']:.6f}; "
              f"library vs raw total (kept genes) max rel {exact[gsm]['max_rel_library_vs_raw_total_kept_genes']:.2g}",
              flush=True)
        offset += len(bcs)

    report = {
        "what_this_is": __doc__.split("\n")[1],
        "raw_dir": str(RAW.relative_to(config.PROJECT_ROOT)),
        "n_files": len(files), "unexpected_files": unexpected,
        "n_gsm_in_geo_table": len(gsm2s), "n_gsm_with_complete_triplet":
            sum(r["complete_triplet"] for r in samples.values()),
        "all_names_match_geo_table": all(r["name_matches_geo_table"] for r in samples.values()),
        "all_dims_agree": all(r.get("dims_agree") for r in samples.values()),
        "cellranger_versions": dict(Counter(r.get("cellranger") for r in samples.values())),
        "total_raw_barcodes": sum(r.get("n_barcodes", 0) for r in samples.values()),
        "processed_cells": len(proc_bc), "processed_cells_found_in_raw": n_found,
        "samples": samples, "cells_per_sample": cell_check, "genes": gene_check,
        "exactness_samples": exact,
    }
    OUT.write_text(json.dumps(report, indent=2))
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
