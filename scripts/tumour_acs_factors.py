"""
tumour_acs_factors.py -- which tumours do deconvolution methods reproduce anatomically, and why?

Implements prespecified/tumour_acs_factors.md exactly (read it first: ACS is agreement with anatomy,
not accuracy, so nothing here is a driver of deconvolution accuracy).

  outcome   per-tumour ACS minus that tumour's own chance ACS (2,000 within-tumour permutations),
            averaged over the comparable registered methods
  factors   F1 anatomy strength in the tumour's own RNA (marker-panel contrasts, deconvolution-free)
            F1b marker ACS (the constraints scored on marker panels, deconvolution-free)
            F2 sampling design, F3 sequencing quality, F4 clinical/molecular (descriptive)
  tests     Spearman across the 9 tumours, exact two-sided p by full enumeration of 9! orderings,
            Holm across the 9 continuous factors; the broken-input controls as negative control

Control first: each tumour's ACS recomputed here must equal the registered per-tumour table.

    python3 scripts/tumour_acs_factors.py      # writes results/tumour_acs_factors.json
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402
from ivygap.anatomic import constraints as K  # noqa: E402
from ivygap.anatomic.acs import (_acs_from_table, _score_table, permutation_null,  # noqa: E402
                                 tumor_structure_means)

PANELS = {   # prespecified/abdelfattah_cluster_mapping.md, fixed before any data
    "Tumor": ["SOX2", "PTPRZ1", "EGFR", "NES", "PDGFRA"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "MAG"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "FLT1"],
    "Macrophage_Microglia": ["CD14", "CD68", "AIF1", "C1QA", "CSF1R"],
}
N_NULL = 2000
OUT = config.RESULTS_DIR / "tumour_acs_factors.json"


def load_manifest() -> pd.DataFrame:
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t", dtype={"sample_id": str, "patient_id": str})
    man = man.set_index("sample_id")
    return man[man["structure"].isin(config.PRIMARY_STRUCTURES) & man["is_anatomic_study"].astype(bool)]


def panel_matrix(samples: list[str]) -> pd.DataFrame:
    """FPKM of the panel genes, genes x samples, streamed from the processed bulk table."""
    want = {g for p in PANELS.values() for g in p}
    rows = {}
    with open(config.BULK_EXPRESSION_PATH) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = [header.index(s) for s in samples]
        for line in fh:
            g, rest = line.split("\t", 1)
            if g in want:
                vals = line.rstrip("\n").split("\t")
                rows[g] = [float(vals[i]) for i in idx]
    missing = sorted(want - set(rows))
    if missing:
        print("panel genes absent from the bulk table:", missing)
    return pd.DataFrame(rows, index=samples).T


def structure_panel_means(fpkm: pd.DataFrame, man: pd.DataFrame, tumour: str) -> dict:
    """{cell_type: {structure: [mean FPKM per panel gene]}} for one tumour."""
    s_of = man["structure"]
    out = {}
    for ct, genes in PANELS.items():
        g = [x for x in genes if x in fpkm.index]
        cols = man.index[man["patient_id"] == tumour]
        out[ct] = {s: fpkm.loc[g, [c for c in cols if s_of[c] == s]].mean(axis=1).to_numpy()
                   for s in sorted(set(s_of[cols]))}
    return out


def _num(v) -> float:
    """First number in a metadata field ("17 yrs" -> 17.0); NaN when there is none."""
    import re                                                       # noqa: PLC0415
    m = re.search(r"\d+(?:\.\d+)?", str(v)) if v is not None else None
    return float(m.group()) if m else float("nan")


def contrast(hi: np.ndarray, lo: np.ndarray) -> float:
    return float(np.mean(np.log2(hi + 1) - np.log2(lo + 1)))


def f1_tumour(means: dict) -> tuple[float, float, int]:
    """(weighted mean contrast, marker ACS, n evaluable) over the constraints ACS evaluates."""
    num = den = sat = 0.0
    n = 0
    for c in K.CONSTRAINTS:
        m = means[c.cell_type]
        if c.kind == "pairwise":
            if not set(c.structures) <= set(m):
                continue
            hi, lo = (m[s] for s in c.structures)
            con = contrast(hi, lo)
            ok = con > 0
        elif c.kind == "maximum":
            t = c.structures[0]
            others = [s for s in c.among if s != t and s in m]
            if t not in m or len(others) < 2:
                continue
            con = contrast(m[t], np.mean([m[s] for s in others], axis=0))
            score = {s: np.mean(np.log2(m[s] + 1)) for s in [t] + others}
            ok = score[t] > max(score[s] for s in others)
        else:                                                  # monotone chain
            if not set(c.structures) <= set(m):
                continue
            steps = list(zip(c.structures, c.structures[1:]))
            con = float(np.mean([contrast(m[b], m[a]) for a, b in steps]))
            score = [np.mean(np.log2(m[s] + 1)) for s in c.structures]
            ok = all(x < y for x, y in zip(score, score[1:]))
        num += c.weight * con
        den += c.weight
        sat += c.weight * float(ok)
        n += 1
    return (num / den, sat / den, n) if den else (float("nan"), float("nan"), 0)


def exact_spearman_p(x: np.ndarray, y: np.ndarray, perms: np.ndarray) -> tuple[float, float]:
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    if den == 0:
        return float("nan"), float("nan")
    obs = float((rx * ry).sum() / den)
    null = (ry[perms] * rx).sum(axis=1) / den
    return obs, float(np.mean(np.abs(null) >= abs(obs) - 1e-12))


def holm(ps: dict) -> dict:
    items = sorted((p, k) for k, p in ps.items() if np.isfinite(p))
    m, out, running = len(items), {}, 0.0
    for i, (p, k) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        out[k] = running
    return out


def main() -> int:
    man = load_manifest()
    board = pd.read_csv(config.ANATOMIC_DIR / "acs_leaderboard.csv").set_index("method")
    per_t = pd.read_csv(config.ANATOMIC_DIR / "acs_per_tumor.csv", dtype={"tumor_id": str})
    comparable = [m for m in board.index if not board.loc[m, "is_control"]
                  and board.loc[m, "comparable"] and m != "scdc_ensemble"]
    controls = [m for m in board.index if board.loc[m, "is_control"]]
    tumours = sorted(per_t["tumor_id"].unique())

    # ---- per-tumour ACS (control: must equal the registered table) and its chance level
    excess, raw, ctrl_excess, checks = {}, {}, {}, []
    for m in comparable + controls:
        est = pd.read_csv(config.ESTIMATES_DIR / f"ivygap_{m}.csv", index_col=0)
        est.index = est.index.astype(str)
        for t in tumours:
            ids = [s for s in man.index[man["patient_id"] == t] if s in est.index]
            sub, msub = est.loc[ids], man.loc[ids]
            acs_t = _acs_from_table(_score_table(tumor_structure_means(sub, msub)))
            reg = per_t[(per_t["method"] == m) & (per_t["tumor_id"] == t)]["acs"]
            if len(reg) and not np.isclose(acs_t, float(reg.iloc[0]), atol=1e-9):
                checks.append(f"{m}/{t}: recomputed {acs_t} != registered {float(reg.iloc[0])}")
            null = permutation_null(sub, msub, n_permutations=N_NULL, seed=config.RANDOM_SEED)
            chance = float(np.nanmean(null))
            target = ctrl_excess if m in controls else excess
            target.setdefault(t, []).append(acs_t - chance)
            if m not in controls:
                raw.setdefault(t, []).append(acs_t)
        print(f"  {m}: per-tumour ACS recomputed and chance-corrected", flush=True)
    if checks:
        print("ABORT -- per-tumour ACS did not reproduce:\n  " + "\n  ".join(checks[:10]))
        return 2

    # ---- factors
    fpkm = panel_matrix(list(man.index))
    prov = json.loads((config.RAW_DIR / "ivygap_counts" / "provenance.json").read_text())
    q = {s["rna_well"]: s for s in prov["samples"]}
    clin = pd.read_csv(config.IVYGAP_RNA_SEQ_DETAILS_PATH, dtype={"tumor_id": str})
    clin = clin.drop_duplicates("tumor_id").set_index("tumor_id")
    rows = {}
    for t in tumours:
        ids = list(man.index[man["patient_id"] == t])
        f1, f1b, n_ev = f1_tumour(structure_panel_means(fpkm, man, t))
        reads = [q[s]["rnaseq_total_reads"] for s in ids if s in q]
        aligned = [q[s]["percent_aligned_mrna"] for s in ids if s in q]
        c = clin.loc[t] if t in clin.index else {}
        rows[t] = {
            "acs_excess": float(np.mean(excess[t])), "acs_raw": float(np.mean(raw[t])),
            "control_acs_excess": float(np.mean(ctrl_excess[t])),
            "F1_rna_anatomy_strength": f1, "F1b_marker_acs": f1b, "n_constraints_evaluable": n_ev,
            "F2_n_samples": len(ids), "F2_n_structures": int(man.loc[ids, "structure"].nunique()),
            "F2_weight_evaluated": float(per_t[(per_t.tumor_id == t) & (per_t.method == comparable[0])]
                                          ["weight_evaluated"].iloc[0]),
            "F3_mean_total_reads": float(np.mean(reads)) if reads else float("nan"),
            "F3_mean_pct_aligned_mrna": float(np.mean(aligned)) if aligned else float("nan"),
            "F4_age": _num(c.get("age_in_years")) if len(c) else float("nan"),   # stored as "17 yrs"
            "F4_initial_kps": _num(c.get("initial_kps")) if len(c) else float("nan"),
            "F4_mgmt": str(c.get("mgmt_methylation", "")) if len(c) else "",
            "F4_egfr_amp": str(c.get("egfr_amplification", "")) if len(c) else "",
            "F4_subtype": str(c.get("molecular_subtype", "")) if len(c) else "",
            "F4_resection": str(c.get("extent_of_resection", "")) if len(c) else "",
            "F4_surgery": str(c.get("surgery", "")) if len(c) else "",
        }
    df = pd.DataFrame(rows).T
    print(df[["acs_raw", "acs_excess", "control_acs_excess", "F1_rna_anatomy_strength", "F1b_marker_acs",
              "F2_n_samples", "F2_n_structures"]].to_string())

    # ---- tests
    perms = np.array(list(itertools.permutations(range(len(tumours)))), dtype=np.int8)
    cont = ["F1_rna_anatomy_strength", "F1b_marker_acs", "F2_n_samples", "F2_n_structures",
            "F2_weight_evaluated", "F3_mean_total_reads", "F3_mean_pct_aligned_mrna", "F4_age",
            "F4_initial_kps"]
    y = df["acs_excess"].astype(float).to_numpy()
    yc = df["control_acs_excess"].astype(float).to_numpy()
    res, ps = {}, {}
    for f in cont:
        x = pd.to_numeric(df[f], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < len(x):
            res[f] = {"note": f"{(~ok).sum()} tumour(s) lack this factor; not tested"}
            continue
        rho, p = exact_spearman_p(x, y, perms)
        rho_c, p_c = exact_spearman_p(x, yc, perms)
        res[f] = {"rho": rho, "p_exact": p, "control_rho": rho_c, "control_p_exact": p_c}
        ps[f] = p
    for f, padj in holm(ps).items():
        res[f]["p_holm"] = padj
        res[f]["associated"] = bool(padj < 0.05 and res[f]["control_p_exact"] >= 0.05)
    cats = {}
    for f in ("F4_mgmt", "F4_egfr_amp", "F4_subtype", "F4_resection", "F4_surgery"):
        g = df.groupby(df[f].replace("", "missing").fillna("missing"))["acs_excess"]
        cats[f] = {k: {"n": int(len(v)), "mean_acs_excess": round(float(v.mean()), 4)} for k, v in g}
    out = {"rule": "prespecified/tumour_acs_factors.md", "n_tumours": len(tumours),
           "comparable_methods": comparable, "controls": controls, "n_null": N_NULL,
           "per_tumour_acs_reproduced": True, "tumours": rows, "tests": res,
           "categorical_descriptive": cats}
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print("\nfactor                        rho     p_exact  p_holm  | control rho  p   | associated")
    for f in cont:
        r = res[f]
        if "rho" in r:
            print(f"{f:28s} {r['rho']:+.3f}  {r['p_exact']:.4f}  {r.get('p_holm', float('nan')):.4f}  | "
                  f"{r['control_rho']:+.3f}  {r['control_p_exact']:.3f} | {r.get('associated')}")
        else:
            print(f"{f:28s} {r['note']}")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
