#!/usr/bin/env python3
"""
imc_anchored_test.py — does ACS reproduce a ranking established by PROTEIN ground truth?

THE POINT
---------
Until now this project has had no external ground truth. The accuracy arm is synthetic
pseudobulk built from the same atlas the ACS arm deconvolves against — the shared dependence
recorded as C4 — so "ACS agrees with accuracy" has always been two GBmap-based rankings
agreeing with each other.

Ajaib et al., Neuro-Oncology 2023;25(7):1236-1248, supply an external one: imaging mass
cytometry, 33 antibodies, single-cell resolution, on ten IDHwt GBM samples with matched bulk
RNA-seq on the same tissue. They scored deconvolution approaches against it, and the core of
their result is a near-controlled contrast — ONE algorithm, MCPcounter, differing only in which
marker set it is handed:

    MCP_GBM      GBM-specific markers    r = 0.37 immune, 0.43 neoplastic
    MCP_default  MCPcounter's own        r = 0.27 immune
    MCP_GBmap    GBmap-derived markers   r = 0.06 immune, 0.22 neoplastic

GBmap is the atlas this entire leaderboard rests on.

This runs the same contrast through this project's pipeline on Ivy GAP and asks whether ACS
recovers the same ordering. The prediction, the arms, the preprocessing and the constraint
subsets were all fixed in `prespecified/imc_anchored_prediction.md` and committed BEFORE any
ACS here was computed.

WHAT IS HELD FIXED
------------------
One algorithm (genuine MCPcounter 1.2.0), one bulk matrix, one log transform applied identically
to every arm, one frozen constraint file. Only the marker set varies.

Every pairwise comparison is scored on the constraints BOTH arms can express, by setting
non-shared cell types to NaN so the scorer drops exactly those pairs — the same partial-coverage
path quanTIseq already goes through. `constraints.py` is not touched and its hash is unchanged.

    python scripts/imc_anchored_test.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402
from ivygap.anatomic.acs import score as acs_score                 # noqa: E402
from ivygap.data.load_ivygap import load_cached                    # noqa: E402

#: NOTE THE SPACES. The vendored tree is literally named "pipeline_packages " with a trailing
#: space, containing " repos" with a LEADING space. Writing "repos/..." as one segment silently
#: drops the leading space and the path does not exist.
GBMD = (ROOT / "pipeline_packages " / " repos" / "GBMDeconvoluteR"
        / "GBMDeconvoluteR-main" / "data")
WORK = config.RESULTS_DIR / "imc_anchored"
RS = ROOT / "scripts" / "mcp_marker_source.R"

#: Marker populations -> this project's roster. Populations with no roster equivalent are
#: DROPPED, never folded into a neighbouring type: forcing monocytes, DC or mast cells into an
#: eight-type roster would invent a correspondence the panel does not have.
AJAIB = {"TAM": "Macrophage_Microglia", "Microglia": "Macrophage_Microglia",
         "T cells": "T_cell", "NK cells": "NK_cell", "B cells": "B_cell"}
MORENO = {**AJAIB, "Endothelial": "Endothelial", "Oligodendrocyte": "Oligodendrocyte",
          "Astrocyte": "Astrocyte"}
#: "Monocytic lineage" -> Macrophage_Microglia is an APPROXIMATION and is declared as one:
#: MCPcounter's default panel has no tumour-associated macrophage or microglia population,
#: because it was not built for brain tissue. That is the point of the contrast, not a flaw
#: in how it is being run.
DEFAULT = {"Endothelial cells": "Endothelial", "Monocytic lineage": "Macrophage_Microglia",
           "T cells": "T_cell", "NK cells": "NK_cell", "B lineage": "B_cell"}
NEFTEL_TO_TUMOR = ("AC", "MES", "NPC", "OPC")


def read_rds(name: str) -> pd.DataFrame:
    """Read one of GBMdeconvoluteR's shipped .rds marker objects via R."""
    out = WORK / f"_{name}.tsv"
    subprocess.run(["Rscript", "-e", f'''
        x <- readRDS("{GBMD / name}")
        if (is.list(x) && !is.data.frame(x)) {{
          x <- data.frame("HUGO symbols" = unlist(x, use.names = FALSE),
                          "Cell population" = rep(names(x), lengths(x)),
                          check.names = FALSE)
        }}
        write.table(x, "{out}", sep = "\\t", row.names = FALSE, quote = FALSE)
    '''], check=True, capture_output=True)
    return pd.read_csv(out, sep="\t", dtype=str)


def to_roster(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    """Rewrite a marker table into the roster's vocabulary, dropping unmapped populations."""
    d = df.copy()
    d.columns = ["HUGO symbols", "Cell population"][:len(d.columns)]
    d["Cell population"] = d["Cell population"].map(mapping)
    return d.dropna(subset=["Cell population"]).drop_duplicates()


def gbmap_markers(n_per_type: int = config.SIGNATURE_GENES_PER_TYPE) -> pd.DataFrame:
    """
    GBmap's own markers, PER TYPE.

    The same score `select_signature_genes` uses — a type's mean divided by the mean across
    the others — kept per type instead of flattened into one list, which is the only change.
    """
    from ivygap.data.reference import build_reference          # noqa: PLC0415
    prof = pd.read_csv(ROOT / "data/reference/gbmap_linear/profile.csv", index_col=0)
    rows = []
    for ctype in prof.columns:
        others = prof.drop(columns=[ctype]).mean(axis=1)
        s = (prof[ctype] / (others + 1.0)).sort_values(ascending=False)
        for g in s.head(n_per_type).index:
            rows.append((g, ctype))
    return pd.DataFrame(rows, columns=["HUGO symbols", "Cell population"])


def run_mcp(tag: str, markers: pd.DataFrame, bulk_csv: Path) -> pd.DataFrame:
    mk = WORK / f"markers_{tag}.tsv"
    markers.to_csv(mk, sep="\t", index=False)
    out = WORK / f"scores_{tag}.csv"
    r = subprocess.run(["Rscript", str(RS), str(mk), str(bulk_csv), str(out)],
                       capture_output=True, text=True)
    print(f"  [{tag}] {r.stdout.strip()}")
    if r.returncode != 0:
        print(f"  [{tag}] FAILED: {r.stderr.strip()[-400:]}")
        return pd.DataFrame()
    s = pd.read_csv(out, index_col=0)                # populations x samples
    est = s.T                                        # samples x populations
    for t in config.CELL_TYPES:                      # every roster column, NaN where absent
        if t not in est.columns:
            est[t] = np.nan
    est = est[list(config.CELL_TYPES)]
    est.index.name = "sample_id"
    return est


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    expr, meta = load_cached()
    anat = [s for s in expr.columns
            if bool(meta.loc[s, "is_anatomic_study"])
            and meta.loc[s, "structure"] in config.PRIMARY_STRUCTURES]
    print(f"Ivy GAP anatomic cohort: {len(anat)} samples, {expr.shape[0]:,} genes")
    bulk_csv = WORK / "bulk_anatomic.csv"
    expr.loc[:, anat].to_csv(bulk_csv)
    man = meta.loc[anat]

    neftel = read_rds("Neftel_et_al_2019_four_state_neoplastic_markers.rds")
    tumor = to_roster(neftel, {k: "Tumor" for k in NEFTEL_TO_TUMOR})
    ajaib = pd.concat([to_roster(read_rds("Ajaib_et_al_2022_GBM_Immune_markers.rds"), AJAIB),
                       tumor])
    moreno = pd.concat([to_roster(read_rds("Moreno_et_al_2022_lvl3_immune_markers.rds"), MORENO),
                        tumor])
    dflt = to_roster(pd.read_csv(WORK.parent / "_mcp_default_genes.tsv", sep="\t", dtype=str)
                     if (WORK.parent / "_mcp_default_genes.tsv").exists()
                     else pd.read_csv(ROOT / "data/raw/mcp_default_genes.txt", sep="\t",
                                      dtype=str)[["HUGO symbols", "Cell population"]], DEFAULT)
    gbmap = gbmap_markers()

    arms = {"MCP_GBM": ajaib, "MCP_GBM_moreno": moreno,
            "MCP_default": dflt, "MCP_GBmap": gbmap}
    print("\nmarker sets, in the roster's vocabulary:")
    for tag, mk in arms.items():
        cov = sorted(mk["Cell population"].unique())
        print(f"  {tag:16s} {len(mk):5d} genes  covers {len(cov)}: {', '.join(cov)}")

    print("\nrunning genuine MCPcounter 1.2.0, one marker set at a time:")
    est = {t: run_mcp(t, mk, bulk_csv) for t, mk in arms.items()}
    est = {t: e for t, e in est.items() if not e.empty}
    for t, e in est.items():
        e.to_csv(WORK / f"estimates_{t}.csv")

    # --- the comparisons, each on the constraints BOTH arms express -------------------
    CONSTRAINT_TYPES = {"Tumor", "Oligodendrocyte", "Endothelial", "Macrophage_Microglia"}
    have = {t: {c for c in CONSTRAINT_TYPES if e[c].notna().any()} for t, e in est.items()}
    print("\nconstraint-relevant coverage per arm:")
    for t, h in have.items():
        print(f"  {t:16s} {sorted(h)}")

    def acs_on(tag: str, types: set[str], label: str) -> dict:
        e = est[tag].copy()
        for c in config.CELL_TYPES:                 # blank everything outside the subset
            if c not in types:
                e[c] = np.nan
        res = acs_score(e, man, method=f"{tag}__{label}",
                        n_permutations=10000, n_boot=2000)
        s = res.summary()
        return {"acs": round(float(s["acs"]), 4),
                "ci": [round(float(s["ci_low"]), 4), round(float(s["ci_high"]), 4)],
                "null_p": float(s["null_p"]),
                "n_pairs": int(s["n_constraint_tumor_pairs"]),
                "n_tumors": int(s["n_tumors"])}

    report = {"arms": {t: sorted(h) for t, h in have.items()}, "comparisons": {}}
    pairs = [("MCP_GBM", "MCP_GBmap", "PRIMARY — the pre-specified test"),
             ("MCP_GBM_moreno", "MCP_GBmap", "secondary — full-coverage GBM markers"),
             ("MCP_default", "MCP_GBmap", "secondary — MCPcounter's own markers"),
             ("MCP_GBM", "MCP_default", "secondary — GBM-specific vs generic")]
    for a, b in [(a, b) for a, b, _ in pairs]:
        if a not in est or b not in est:
            continue
    for a, b, why in pairs:
        if a not in est or b not in est:
            print(f"\n[{a} vs {b}] SKIPPED: an arm failed")
            continue
        shared = have[a] & have[b]
        if not shared:
            print(f"\n[{a} vs {b}] SKIPPED: no shared constraint type")
            continue
        label = "+".join(sorted(shared))
        print(f"\n=== {a} vs {b} — {why} ===")
        print(f"    shared cell types: {sorted(shared)}")
        ra, rb = acs_on(a, shared, label), acs_on(b, shared, label)
        print(f"    {a:16s} ACS {ra['acs']:.4f}  CI [{ra['ci'][0]:.4f}, {ra['ci'][1]:.4f}]  "
              f"p={ra['null_p']:.4f}  pairs={ra['n_pairs']}")
        print(f"    {b:16s} ACS {rb['acs']:.4f}  CI [{rb['ci'][0]:.4f}, {rb['ci'][1]:.4f}]  "
              f"p={rb['null_p']:.4f}  pairs={rb['n_pairs']}")
        print(f"    -> {a} {'>' if ra['acs'] > rb['acs'] else '<=' } {b}")
        report["comparisons"][f"{a}_vs_{b}"] = {
            "why": why, "shared_cell_types": sorted(shared),
            a: ra, b: rb, "winner": a if ra["acs"] > rb["acs"] else
            (b if rb["acs"] > ra["acs"] else "tie")}

    report["external_ground_truth"] = {
        "source": "Ajaib et al., Neuro-Oncology 2023;25(7):1236-1248",
        "modality": "imaging mass cytometry, 33 antibodies, single-cell resolution",
        "cohort": "10 IDHwt GBM (5 paired primary/recurrent) with matched bulk RNA-seq",
        "pearson_r_vs_IMC": {"MCP_GBM": {"immune": 0.37, "neoplastic": 0.43},
                             "MCP_default": {"immune": 0.27},
                             "MCP_GBmap": {"immune": 0.06, "neoplastic": 0.22},
                             "CIBERSORTx": {"immune": 0.05, "neoplastic": 0.02}},
        "prespecified": "prespecified/imc_anchored_prediction.md",
        "prediction": "ACS(MCP_GBM) > ACS(MCP_GBmap) on the shared constraints",
    }
    (config.RESULTS_DIR / "imc_anchored_test.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {(config.RESULTS_DIR / 'imc_anchored_test.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
