"""
cptac_per_sample.py -- prespecified/cptac_per_sample_truth.md: per-sample accuracy against single-nucleus
counts from the SAME cryopulverised tissue (CPTAC glioblastoma, Wang et al. 2021; 17 tumours).

Stages (each writes under results/cptac/, so the job resumes):
  bulk     GDC STAR counts ('unstranded') -> CPM on GENCODE v36 symbols, the TCGA bulk's construction
  truth    per-sample nucleus composition: GDC's Seurat clusters, named by the marker rule declared on
           2026-10-01 for the Abdelfattah atlas (PANELS, margin 0.5), unchanged
  methods  the registered panel under the frozen signature and the donor-level raw/X reference, run
           exactly as the TCGA arms run it (absolute_purity_yardstick.py), all samples required finite
  analyse  controls first (truth vs DNA purity; mapping coverage), then C1-C4

    python3 scripts/cptac_per_sample.py --stage bulk|truth|methods|analyse|all
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402

SRC = config.CPTAC_GBM_DIR
OUT = config.RESULTS_DIR / "cptac"
BULK = config.PROCESSED_DIR / "cptac_gbm_bulk_cpm.csv.gz"
MIN_MAPPED = 0.70          # prespecified: a sample with < 70% of nuclei mapped is excluded
N_PERM = 10_000


def manifest() -> dict:
    return json.loads((SRC / "MANIFEST.json").read_text())


def files(kind: str) -> dict[str, Path]:
    """case -> local file, for one kind; a case with two files of a kind is an error, never a choice."""
    out: dict[str, Path] = {}
    for r in manifest()["files"]:
        if r["kind"] == kind:
            if r["case"] in out:
                raise SystemExit(f"BLOCKED: case {r['case']} has two {kind} files; the rule does not choose")
            out[r["case"]] = config.PROJECT_ROOT / r["local"]
    return out


# ------------------------------------------------------------------------------------------- bulk
def stage_bulk() -> None:
    cols = {}
    for case, p in sorted(files("bulk").items()):
        df = pd.read_csv(p, sep="\t", comment="#")
        df = df[~df["gene_id"].astype(str).str.startswith("N_")]
        s = df.groupby("gene_name")["unstranded"].mean()          # duplicates collapsed by mean
        cols[case] = s
    expr = pd.DataFrame(cols).fillna(0.0)
    cpm = expr.div(expr.sum(axis=0), axis=1) * 1e6
    BULK.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(BULK, "wt") as fh:
        cpm.to_csv(fh, float_format="%.4f")
    print(f"bulk: {cpm.shape[0]:,} genes x {cpm.shape[1]} samples -> {BULK.relative_to(config.PROJECT_ROOT)}")


# ------------------------------------------------------------------------------------------- truth
def read_mex(tgz: Path):
    """CellRanger filtered matrix from a .tar.gz: returns (sparse genes x nuclei counts, gene symbols,
    barcodes). Member names are matched by suffix, so the archive's folder layout does not matter."""
    import scipy.io  # noqa: PLC0415
    with tarfile.open(tgz, "r:gz") as tf:
        mem = {m.name.split("/")[-1]: m for m in tf.getmembers() if m.isfile()}
        def get(stem):
            for k in (stem, stem + ".gz"):
                if k in mem:
                    raw = tf.extractfile(mem[k]).read()
                    return gzip.decompress(raw) if k.endswith(".gz") else raw
            raise SystemExit(f"BLOCKED: {tgz.name} has no {stem}[.gz]: {sorted(mem)}")
        mtx = scipy.io.mmread(io.BytesIO(get("matrix.mtx"))).tocsc()
        feats = pd.read_csv(io.BytesIO(get("features.tsv")), sep="\t", header=None)
        bcs = pd.read_csv(io.BytesIO(get("barcodes.tsv")), sep="\t", header=None)[0].astype(str).tolist()
    sym = feats[1].astype(str).tolist() if feats.shape[1] > 1 else feats[0].astype(str).tolist()
    return mtx, sym, bcs


def gdc_clusters(tsv: Path) -> pd.Series:
    """GDC's Seurat cluster per barcode. Exactly one cluster column must be identifiable; anything
    else is BLOCKED (prespecified: never re-cluster, never pick among clusterings after looking)."""
    df = pd.read_csv(tsv, sep="\t")
    cands = [c for c in df.columns if "cluster" in c.lower()]
    if len(cands) != 1:
        raise SystemExit(f"BLOCKED: {tsv.name}: cluster columns {cands} (need exactly one) -- columns {list(df.columns)}")
    bc = next((c for c in df.columns if c.lower() in ("barcode", "barcodes", "cell", "cell_barcode", "cellbarcode")),
              df.columns[0])
    return pd.Series(df[cands[0]].astype(str).to_numpy(), index=df[bc].astype(str).to_numpy(), name=cands[0])


def label_sample(mtx, sym: list[str], bcs: list[str], clusters: pd.Series, panels: dict | None = None
                 ) -> tuple[pd.Series, dict]:
    """The Abdelfattah rule, unchanged: per-cluster mean log1p-CPM of each PANEL; best panel must beat the
    runner-up by MARGIN; excluded or ambiguous clusters unmapped; a T/NK cluster within the margin is
    split per nucleus by the sign of (T panel - NK panel)."""
    from build_abdelfattah_reference import MARGIN, PANELS  # noqa: PLC0415
    PANELS = panels or PANELS  # noqa: N806
    idx = {g: i for i, g in enumerate(sym)}
    lib = np.asarray(mtx.sum(axis=0)).ravel()
    keep = [i for i, b in enumerate(bcs) if b in clusters.index]
    cl = clusters.reindex([bcs[i] for i in keep]).to_numpy()
    def panel_lcpm(genes):
        rows = [idx[g] for g in genes if g in idx]
        if not rows:
            return np.full(len(keep), np.nan)
        sub = mtx[rows][:, keep].toarray().astype(float)
        return np.log1p(sub / lib[keep] * 1e6).mean(axis=0)
    scores = {t: panel_lcpm(g) for t, g in PANELS.items()}
    label = np.array([None] * len(keep), dtype=object)
    per_cluster = {}
    for c in sorted(set(cl)):
        r = cl == c
        sc = {t: float(np.nanmean(v[r])) for t, v in scores.items()}
        best = sorted(sc.items(), key=lambda kv: -kv[1])
        top = best[0][0]
        if top in ("T_cell", "NK_cell"):
            # Addendum 1 (the declared rule, implemented as declared): the margin is between the
            # lymphoid T/NK pair and the best OTHER panel; within the pair, |T - NK| < MARGIN splits
            # the cluster per nucleus by the sign of (T panel - NK panel)
            other = max(v for k, v in sc.items() if k not in ("T_cell", "NK_cell"))
            if max(sc["T_cell"], sc["NK_cell"]) - other < MARGIN:
                t = None
                label[r] = None
            elif abs(sc["T_cell"] - sc["NK_cell"]) >= MARGIN:
                t = top
                label[r] = t
            else:
                d = scores["T_cell"][r] - scores["NK_cell"][r]
                label[r] = np.where(d >= 0, "T_cell", "NK_cell")
                t = "T/NK split"
        else:
            ok = top != "_excluded" and best[0][1] - best[1][1] >= MARGIN
            t = top if ok else None
            label[r] = t
        per_cluster[c] = {"n": int(r.sum()), "label": t, "scores": {k: round(v, 3) for k, v in sc.items()}}
    lab = pd.Series(label, index=[bcs[i] for i in keep])
    # Diagnostic M1 (Addendum 2): nuclei inside tumour-labelled clusters whose OWN myeloid panel beats
    # their tumour panel by >= MARGIN -- myeloid nuclei the clustering probably absorbed
    in_tumour = np.array([x == "Tumor" for x in label])
    m1 = (scores["Macrophage_Microglia"] - scores["Tumor"] >= MARGIN) & in_tumour
    per_cluster["_M1"] = {"n_in_tumour_clusters": int(in_tumour.sum()), "n_myeloid_like": int(m1.sum()),
                          "share": round(float(m1.sum() / max(in_tumour.sum(), 1)), 4)}
    return lab, per_cluster


def stage_truth(variant: str = "registered") -> None:
    """variant 'registered': the declared panels. 's1' (Addendum 2): the NK panel without NCAM1."""
    from build_abdelfattah_reference import PANELS  # noqa: PLC0415
    panels = dict(PANELS)
    if variant == "s1":
        panels["NK_cell"] = [g for g in PANELS["NK_cell"] if g != "NCAM1"]
    tag = "" if variant == "registered" else f"_{variant}"
    mats, tabs = files("snrna"), files("snrna_gdc")
    rows, detail = {}, {}
    for case in sorted(set(mats) & set(tabs)):
        mtx, sym, bcs = read_mex(mats[case])
        clusters = gdc_clusters(tabs[case])
        lab, per_cluster = label_sample(mtx, sym, bcs, clusters, panels)
        m1 = per_cluster.pop("_M1")
        n = len(lab)
        mapped = lab.dropna()
        frac = mapped.value_counts().reindex(list(config.CELL_TYPES), fill_value=0) / max(len(mapped), 1)
        rows[case] = frac
        detail[case] = {"n_nuclei_filtered": len(bcs), "n_with_gdc_cluster": n, "n_mapped": int(len(mapped)),
                        "share_mapped": round(len(mapped) / n, 4) if n else 0.0,
                        "counts": {k: int(v) for k, v in mapped.value_counts().items()}, "M1": m1,
                        "clusters": per_cluster}
        print(f"  {case}: {len(bcs)} nuclei, {n} clustered, {len(mapped)} mapped ({detail[case]['share_mapped']:.0%})", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).T.rename_axis("case").to_csv(OUT / f"snrna_truth{tag}.csv")
    (OUT / f"snrna_truth{tag}_detail.json").write_text(json.dumps(detail, indent=2))
    print(f"wrote {(OUT / f'snrna_truth{tag}.csv').relative_to(config.PROJECT_ROOT)}")


# ------------------------------------------------------------------------------------------- methods
def stage_methods(reference: str) -> None:
    from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: PLC0415
    from ivygap.deconv.base import DeconvolutionInput  # noqa: PLC0415
    from ivygap.deconv.registry import build_methods  # noqa: PLC0415
    with gzip.open(BULK, "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    if reference == "h5ad":
        from ivygap.data.reference import build_from_h5ad  # noqa: PLC0415
        from ivygap.deconv import r_bridge  # noqa: PLC0415
        ref, sc_expr, sc_meta = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad", matrix="raw/X",
                                                restrict_to_genes=bulk.index, export=False)
        r_bridge.set_cell_source(config.PRIMARY_REFERENCE, sc_expr, sc_meta)
    else:
        ref = load_frozen_reference()
    shared = [g for g in ref.profile.index if g in bulk.index]
    picked = select_signature_genes(ref.subset_genes(shared), n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    genes = [g for g in picked if g in bulk.index]
    cols = list(bulk.columns)
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6
    manifest_df = pd.DataFrame({"patient_id": cols}, index=cols)      # one sample per patient here
    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),), manifest=manifest_df,
                              cell_types=tuple(ref.cell_types), bulk_full=bulk.loc[:, cols])
    print(f"{reference}: {len(shared):,} shared genes; {len(genes)} markers from the reference alone; {len(cols)} samples")
    rows, report = [], {}
    for m in build_methods(prefer_r=True):
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                       # noqa: BLE001
            report[m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:160]}"}
            print(f"  {m.name:24s} FAILED {type(exc).__name__}"); continue
        full = est.reindex(index=cols, columns=list(config.CELL_TYPES))
        finite = int(np.isfinite(full.to_numpy(float)).all(axis=1).sum())
        if finite < len(cols):                                         # prespecified small-cohort rule
            report[m.name] = {"skipped": f"{finite} of {len(cols)} samples finite"}
            print(f"  {m.name:24s} SKIPPED ({finite}/{len(cols)} finite)"); continue
        report[m.name] = {"implementation": getattr(m, "implementation_", "python"),
                          "degenerate": bool(getattr(m, "degenerate_", False)),
                          "is_control": m.name.startswith("control_")}
        rows.append(full.assign(method=m.name).rename_axis("sample").reset_index())
        print(f"  {m.name:24s} ok ({report[m.name]['implementation']})", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat(rows).to_csv(OUT / f"estimates_{reference}.csv", index=False)
    (OUT / f"methods_{reference}.json").write_text(json.dumps({"reference": reference, "n_markers": len(genes),
                                                               "methods": report}, indent=2))
    print(f"wrote {(OUT / f'estimates_{reference}.csv').relative_to(config.PROJECT_ROOT)}")


# ------------------------------------------------------------------------------------------- methylation (control 1)
def stage_methyl() -> None:
    """GIMiCC's Layer-0 purity (InfiniumPurify) per case from the CPTAC methylation: the DNA side of control 1.
    Library CpGs only (results/gimicc/gimicc_library_cpgs.txt); incomplete CpGs dropped by run_gimicc.R."""
    import subprocess  # noqa: PLC0415
    keep = set((config.RESULTS_DIR / "gimicc" / "gimicc_library_cpgs.txt").read_text().split())
    cols = {}
    for case, f in sorted(files("methylation").items()):
        d = pd.read_csv(f, sep="\t", header=None, index_col=0)
        d = d[d.index.isin(keep)]
        cols[case] = pd.to_numeric(d.iloc[:, 0], errors="coerce")
    m = pd.DataFrame(cols)
    OUT.mkdir(parents=True, exist_ok=True)
    src = OUT / "cptac_gimicc_cpgs.tsv"
    m.rename_axis("sample").to_csv(src, sep="\t", na_rep="NA")
    args = {"in": str(src), "tumor_type": "GBM", "h": 4, "out": str(OUT / "gimicc_h4.csv")}
    (OUT / "gimicc_args.json").write_text(json.dumps(args))
    r = subprocess.run(["Rscript", str(ROOT / "R" / "run_gimicc.R"), str(OUT / "gimicc_args.json")],
                       capture_output=True, text=True)
    (OUT / "gimicc.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0:
        raise SystemExit(f"BLOCKED: GIMiCC failed on CPTAC methylation; see {OUT / 'gimicc.log'}")
    print(f"GIMiCC on {m.shape[1]} cases, {m.shape[0]:,} library CpGs present -> {args['out']}")


# ------------------------------------------------------------------------------------------- analyse
def perm_p(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    rho = float(stats.spearmanr(x, y).statistic)
    if not np.isfinite(rho):
        return float("nan"), float("nan")
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    # the same permutations, drawn in the same order from the same generator as a per-draw loop;
    # only the correlations are computed at once (Pearson of ranks = mean product of z-scores)
    P = np.array([rng.permutation(ry) for _ in range(N_PERM)])
    zx = (rx - rx.mean()) / rx.std()
    zp = (P - P.mean(axis=1, keepdims=True)) / P.std(axis=1, keepdims=True)
    null = zp @ zx / len(rx)
    return rho, float((1 + np.sum(np.abs(null) >= abs(rho) - 1e-12)) / (N_PERM + 1))


COMPARTMENTS = {t: [t] for t in config.CELL_TYPES}
COMPARTMENTS.update({"Leukocytes": ["Macrophage_Microglia", "T_cell", "NK_cell", "B_cell"],
                     "Lymphoid": ["T_cell", "NK_cell", "B_cell"]})


def analyse_truth(tag: str, rng: np.random.Generator) -> dict:
    """Controls and C1-C4 against one version of the truth ('' registered; '_s1' Addendum 2)."""
    truth = pd.read_csv(OUT / f"snrna_truth{tag}.csv", index_col=0)
    det = json.loads((OUT / f"snrna_truth{tag}_detail.json").read_text())
    share = pd.Series({k: v["share_mapped"] for k, v in det.items()})
    excluded = sorted(share[share < MIN_MAPPED].index)
    cases = [c for c in truth.index if c not in excluded]
    out: dict = {"n_cases_with_truth": int(len(truth)), "excluded_low_mapping": excluded, "n_cases_scored": len(cases)}
    out["control_mapping"] = {"median_share_mapped": round(float(share.median()), 4), "pass": bool(share.median() >= MIN_MAPPED)}
    g = pd.read_csv(OUT / "gimicc_h4.csv", index_col=0)
    g.index = [str(i) for i in g.index]
    common = [c for c in cases if c in g.index and np.isfinite(g.loc[c, "Tumor"])]
    r1, p1 = perm_p(truth.loc[common, "Tumor"].to_numpy(float), g.loc[common, "Tumor"].to_numpy(float), rng)
    out["control_truth_vs_dna_purity"] = {"spearman": round(r1, 4), "perm_p": p1, "n": len(common), "pass": bool(r1 >= 0.40)}
    out["controls_pass"] = out["control_mapping"]["pass"] and out["control_truth_vs_dna_purity"]["pass"]
    m1 = {c: det[c]["M1"]["share"] for c in cases}
    out["M1_myeloid_absorbed_share"] = m1
    out["M1_flag_myeloid_unreliable"] = bool(max(m1.values()) > 0.10) if m1 else None
    lb = pd.read_csv(config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv")
    acs = dict(zip(lb["method"], lb["acs"]))
    reg12 = list(json.loads((config.RESULTS_DIR / "yardstick_agreement.json").read_text())["per_method"])
    for ref in ("frozen", "h5ad"):
        f = OUT / f"estimates_{ref}.csv"
        if not f.exists():
            continue
        est = pd.read_csv(f)
        R: dict = {"methods": {}}
        for m, gm in est.groupby("method"):
            gm = gm.set_index("sample").reindex(cases)
            rec = {}
            for comp, types in COMPARTMENTS.items():
                x = gm[types].sum(axis=1).to_numpy(float)
                y = truth.loc[cases, types].sum(axis=1).to_numpy(float)
                if np.nanstd(x) == 0 or np.nanstd(y) == 0:
                    rec[comp] = {"rho": None, "why": "constant"}
                    continue
                r, p = perm_p(x, y, rng)
                rec[comp] = {"rho": round(r, 4), "perm_p": p, "mean_est": round(float(np.mean(x)), 4),
                             "mean_truth": round(float(np.mean(y)), 4)}
            L = gm[["T_cell", "NK_cell", "B_cell"]]
            Ls = L[L.sum(axis=1) > 0]
            mu = (Ls.div(Ls.sum(axis=1), axis=0)).mean() if len(Ls) else pd.Series(dtype=float)
            rec["_within_lymphoid_mean"] = {k: round(float(v), 4) for k, v in mu.items()}
            rec["_is_control"] = m.startswith("control_")
            R["methods"][m] = rec
        real = {m: v for m, v in R["methods"].items() if not v["_is_control"]}
        R["compartment_median_rho"] = {
            comp: (round(float(np.median([v[comp]["rho"] for v in real.values() if v[comp].get("rho") is not None])), 4)
                   if any(v[comp].get("rho") is not None for v in real.values()) else None)
            for comp in COMPARTMENTS}
        R["C1_reading"] = {comp: ("tracks" if (v is not None and v >= 0.40) else "not recovered")
                           for comp, v in R["compartment_median_rho"].items()}
        ctrl = {m: {c: v[c].get("rho") for c in ("Tumor", "Leukocytes")} for m, v in R["methods"].items() if v["_is_control"]}
        R["negative_control_methods"] = ctrl
        ms = [m for m in reg12 if m in real and real[m]["Tumor"].get("rho") is not None and m in acs]
        if len(ms) >= 6:
            r2 = stats.spearmanr([acs[m] for m in ms], [real[m]["Tumor"]["rho"] for m in ms])
            R["C2_acs_vs_tumour_accuracy"] = {"n_methods": len(ms), "spearman": round(float(r2.statistic), 4),
                                              "p": round(float(r2.pvalue), 4),
                                              "reading": "anatomy ranks methods by accuracy here" if r2.statistic >= 0.60
                                              else "anatomy does not rank methods by accuracy here"}
        out[ref] = R
    if "frozen" in out and "h5ad" in out:
        out["C3_reference_effect_tumour_rho"] = {
            m: {"frozen": out["frozen"]["methods"][m]["Tumor"].get("rho"), "h5ad": out["h5ad"]["methods"][m]["Tumor"].get("rho")}
            for m in out["frozen"]["methods"] if m in out["h5ad"]["methods"]}
    tot = {"T_cell": 0, "NK_cell": 0, "B_cell": 0}
    n_with_b = 0
    for c in cases:
        cnt = det[c]["counts"]
        for k in tot:
            tot[k] += cnt.get(k, 0)
        n_with_b += cnt.get("B_cell", 0) > 0
    s_ = sum(tot.values())
    out["C4_snrna_lymphoid"] = {"pooled_nuclei": tot,
                                "pooled_within_lymphoid_share": {k: round(v / s_, 4) for k, v in tot.items()} if s_ else None,
                                "n_cases_with_any_B": int(n_with_b),
                                "per_sample_B_comparison": ("UNRESOLVED (fewer than 5 cases with B nuclei)" if n_with_b < 5
                                                            else "evaluable")}
    return out


def stage_analyse() -> int:
    rng = np.random.default_rng(config.RANDOM_SEED)
    out = {"rule": "prespecified/cptac_per_sample_truth.md (with Addenda 1-2)",
           "registered_truth": analyse_truth("", rng), "s1_truth_nk_without_ncam1": analyse_truth("_s1", rng)}
    (config.RESULTS_DIR / "cptac_per_sample_truth.json").write_text(json.dumps(out, indent=2, default=float))
    for name in ("registered_truth", "s1_truth_nk_without_ncam1"):
        o = out[name]
        print(f"=== {name}: {o['n_cases_scored']} cases scored; excluded {o['excluded_low_mapping']}")
        print("  controls:", {k: o[k] for k in ("control_mapping", "control_truth_vs_dna_purity", "controls_pass")})
        print("  M1 flag (myeloid unreliable):", o["M1_flag_myeloid_unreliable"], "| max M1", max(o["M1_myeloid_absorbed_share"].values()))
        for ref in ("frozen", "h5ad"):
            if ref in o:
                print(f"  {ref} median rho by compartment:", o[ref]["compartment_median_rho"])
                print(f"  {ref} C2:", o[ref].get("C2_acs_vs_tumour_accuracy"))
                print(f"  {ref} negative controls:", o[ref]["negative_control_methods"])
        print("  C4:", o["C4_snrna_lymphoid"])
    print("wrote results/cptac_per_sample_truth.json")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stage", choices=["bulk", "truth", "methyl", "methods", "analyse", "all"], default="all")
    ap.add_argument("--reference", choices=["frozen", "h5ad", "both"], default="both")
    a = ap.parse_args()
    if a.stage in ("bulk", "all"):
        stage_bulk()
    if a.stage in ("truth", "all"):
        stage_truth("registered")
        stage_truth("s1")
    if a.stage in ("methyl", "all"):
        stage_methyl()
    if a.stage in ("methods", "all"):
        for r in (["frozen", "h5ad"] if a.reference == "both" else [a.reference]):
            stage_methods(r)
    if a.stage in ("analyse", "all"):
        return stage_analyse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
