"""
repaired_methods.py -- re-run ONLY the repaired methods (docs/METHOD_REPAIRS.md) on one TCGA arm, under
IVYGAP_REPAIRED=1, and score them against every truth the evaluation matrix uses.

Writes results/repaired/tcga_<cohort>_<reference>.json and results/repaired/estimates_<cohort>_<reference>.csv.
It never writes a registered artefact.

The inputs are prepared exactly as scripts/absolute_purity_yardstick.py prepares them (same bulk, ABSOLUTE join,
reference build, marker rule and manifest). The run aborts if its sample or gene count differs from the registered
artefact's, so a repaired number is always comparable with the registered one beside it.

    python3 scripts/repaired_methods.py --cohort gbm --reference frozen
    python3 scripts/repaired_methods.py --cohort lgg --reference h5ad --methods dwls quantiseq
"""
from __future__ import annotations

import os

os.environ["IVYGAP_REPAIRED"] = "1"            # before ivygap.config is imported: this run IS the repaired one

import argparse  # noqa: E402
import gzip  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import identifiability_diagnostics as idd  # noqa: E402
from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: E402

assert config.REPAIRED_METHODS, "IVYGAP_REPAIRED was not seen by ivygap.config"
OUT = config.RESULTS_DIR / "repaired"
DEFAULT_METHODS = ["dwls", "quantiseq", "bayesian_hierarchical"]
COMPS = [("Tumor", "ABSOLUTE"), ("Leukocytes", "LF"), ("Lymphoid", "EpiDISH"), ("T_cell", "EpiDISH"),
         ("B_cell", "EpiDISH"), ("NK_cell", "EpiDISH"), ("Myeloid", "none")]


def prepare(cohort: str, reference: str):
    """absolute_purity_yardstick.main's preparation, line for line; returns the DeconvolutionInput and checks."""
    from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: PLC0415
    from ivygap.deconv.base import DeconvolutionInput  # noqa: PLC0415
    bulk_path = config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz"
    with gzip.open(bulk_path, "rt") as fh:
        bulk = pd.read_csv(fh, index_col=0)
    absol = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    absol["k"] = absol["sample"].map(key4)
    absol = absol.drop_duplicates("k").set_index("k")
    bmap = {key4(c): c for c in bulk.columns}
    shared = sorted(set(absol.index) & set(bmap))
    cols = [bmap[k] for k in shared]
    if reference == "h5ad":
        from ivygap.data.reference import build_from_h5ad  # noqa: PLC0415
        from ivygap.deconv import r_bridge  # noqa: PLC0415
        ref, sc_expr, sc_meta = build_from_h5ad(config.REFERENCE_DIR / "gbmap_core.h5ad", matrix="raw/X",
                                                restrict_to_genes=bulk.index, export=False)
        r_bridge.set_cell_source(config.PRIMARY_REFERENCE, sc_expr, sc_meta)
    else:
        ref = load_frozen_reference()
    shared_genes = [g for g in ref.profile.index if g in bulk.index]
    picked = select_signature_genes(ref.subset_genes(shared_genes), n_per_type=config.SIGNATURE_GENES_PER_TYPE)
    genes = [g for g in picked if g in bulk.index]
    sub = bulk.loc[genes, cols]
    sub = sub / sub.sum(axis=0) * 1e6
    manifest = pd.DataFrame({"patient_id": ["-".join(str(c).split("-")[:3]) for c in sub.columns]},
                            index=sub.columns)
    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),), manifest=manifest,
                              cell_types=tuple(ref.cell_types), bulk_full=bulk.loc[:, cols])
    return data, len(cols), len(genes)


def prepare_cptac(reference: str):
    """cptac_per_sample.stage_methods's preparation, line for line (CPTAC bulk; the same reference build and
    marker rule; one sample per patient)."""
    from ivygap.data.reference import load_frozen_reference, select_signature_genes  # noqa: PLC0415
    from ivygap.deconv.base import DeconvolutionInput  # noqa: PLC0415
    import cptac_per_sample as cp  # noqa: PLC0415
    with gzip.open(cp.BULK, "rt") as fh:
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
    manifest = pd.DataFrame({"patient_id": cols}, index=cols)
    data = DeconvolutionInput(bulk=sub, references=(ref.subset_genes(genes),), manifest=manifest,
                              cell_types=tuple(ref.cell_types), bulk_full=bulk.loc[:, cols])
    return data, len(cols), len(genes)


def run_cptac(reference: str, methods: list[str]) -> int:
    """CPTAC (18 glioblastomas): the repaired methods against whole-genome AscatNGS purity (secondary analysis S2's
    truth), beside each method's registered value there."""
    from ivygap.deconv import r_bridge  # noqa: PLC0415
    from ivygap.deconv.registry import build_methods  # noqa: PLC0415
    import cptac_per_sample as cp  # noqa: PLC0415
    r_bridge.R_METHOD_TIMEOUTS.update({"dwls": 10**9, "quantiseq": 10**9})
    reg_meta = json.loads((cp.OUT / f"methods_{reference}.json").read_text())
    reg_rho = json.loads((config.RESULTS_DIR / "cptac_wgs_purity.json").read_text())["all_cases"][reference]["methods"]
    data, n_samples, n_genes = prepare_cptac(reference)
    if (n_samples, n_genes) != (18, reg_meta["n_markers"]):
        print(f"BLOCKED: prepared {n_samples} x {n_genes}; the registered CPTAC arm has 18 x {reg_meta['n_markers']}")
        return 2
    purity = pd.read_csv(cp.OUT / "wgs_purity_per_sample.csv", index_col=0)["wgs_purity"]
    rng = np.random.default_rng(config.RANDOM_SEED)
    report: dict = {"what_this_is": "repaired methods on CPTAC against whole-genome purity", "cohort": "cptac",
                    "reference": reference, "n_samples": n_samples, "n_genes": n_genes, "repaired": True, "methods": {}}
    frames = []
    for m in build_methods(prefer_r=True):
        if m.name not in methods:
            continue
        t0 = time.time()
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                       # noqa: BLE001
            report["methods"][m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:200]}"}
            continue
        rec = {"implementation": getattr(m, "implementation_", "python"),
               "fallback_reason": getattr(m, "fallback_reason_", None), "seconds": round(time.time() - t0, 1)}
        t = est["Tumor"].reindex(purity.index).to_numpy(float)
        ok = np.isfinite(t) & np.isfinite(purity.to_numpy(float))
        if ok.sum() >= 4 and np.ptp(t[ok]) > 0:
            r, pval = cp.perm_p(t[ok], purity.to_numpy(float)[ok], rng)
            rec["Tumor|WGS"] = {"rho": round(r, 4), "perm_p": pval, "n": int(ok.sum())}
        else:
            rec["Tumor|WGS"] = {"rho": None, "n": int(ok.sum()), "why": "not modelled or constant"}
        rec["registered"] = {"implementation": (reg_rho.get(m.name) or {}).get("implementation"),
                             "tumour_rho": (reg_rho.get(m.name) or {}).get("rho"),
                             "status": ("ok" if "implementation" in reg_meta["methods"].get(m.name, {})
                                        else next(iter(reg_meta["methods"].get(m.name, {"absent": ""}).keys())))}
        report["methods"][m.name] = rec
        frames.append(est.assign(method=m.name).rename_axis("sample_id").reset_index())
        print(f"  {m.name:24s} {rec['implementation']}: tumour vs WGS {rec['Tumor|WGS'].get('rho')} "
              f"(registered {rec['registered']['tumour_rho']}) [{rec['seconds']} s]", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"cptac_{reference}.json").write_text(json.dumps(report, indent=2, default=float))
    if frames:
        pd.concat(frames).to_csv(OUT / f"estimates_cptac_{reference}.csv", index=False)
    print(f"wrote results/repaired/cptac_{reference}.json")
    return 0


def lymphoid_direction(est: pd.DataFrame, tr: dict) -> dict:
    """Cohort-mean within-lymphoid T, NK, B over the methylation-matched samples, as lymphoid_ordering reports it."""
    e = est.copy()
    e.index = [idd.k4s(i) for i in e.index]
    e = e[~e.index.duplicated()]
    L = e.loc[e.index.intersection(tr["epi_T"].index), ["T_cell", "NK_cell", "B_cell"]]
    if L.isna().all().all():
        return {"why": "the method returns no lymphoid estimate"}
    tot = L.sum(axis=1)
    sh = L[tot > 0].div(tot[tot > 0], axis=0).mean()
    return {"n": int(len(L)), "T": round(float(sh["T_cell"]), 4), "NK": round(float(sh["NK_cell"]), 4),
            "B": round(float(sh["B_cell"]), 4), "T_exceeds_B": bool(sh["T_cell"] > sh["B_cell"]),
            "samples_with_no_lymphoid_signal": int((tot <= 0).sum())}


def main() -> int:
    from ivygap.deconv import r_bridge  # noqa: PLC0415
    from ivygap.deconv.registry import build_methods  # noqa: PLC0415
    # NO BUDGET for the genuine packages in a repaired run (the user's 2026-10-01 directive). The registered
    # pipeline's 40-minute DWLS budget would cut a 510-sample run off and silently substitute the
    # reimplementation; this changes the limit in this process only.
    r_bridge.R_METHOD_TIMEOUTS.update({"dwls": 10**9, "quantiseq": 10**9})
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", choices=["gbm", "lgg", "cptac"], required=True)
    ap.add_argument("--reference", choices=["frozen", "h5ad"], required=True)
    ap.add_argument("--methods", nargs="+", default=DEFAULT_METHODS)
    ap.add_argument("--python-only", action="store_true",
                    help="run this project's Python versions (the repaired DWLS reimplements the package to 5e-11 "
                         "given its dampening constant); labelled as reimplementations, never as the package")
    a = ap.parse_args()
    if a.cohort == "cptac":
        return run_cptac(a.reference, a.methods)
    tag = ("" if a.cohort == "gbm" else "_lgg") + ("_h5ad" if a.reference == "h5ad" else "")
    registered = json.loads((config.RESULTS_DIR / f"absolute_purity_yardstick{tag}.json").read_text())
    data, n_samples, n_genes = prepare(a.cohort, a.reference)
    if (n_samples, n_genes) != (registered["n_samples"], registered["n_genes"]):
        print(f"BLOCKED: prepared {n_samples} samples x {n_genes} genes; the registered arm has "
              f"{registered['n_samples']} x {registered['n_genes']}. Not comparable, so nothing is written.")
        return 2
    print(f"{a.cohort} {a.reference}: {n_samples} samples x {n_genes} genes (= the registered arm)")
    tr = idd.truths(a.cohort)
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {"what_this_is": __doc__.split("\n")[1].strip(), "cohort": a.cohort, "reference": a.reference,
                    "n_samples": n_samples, "n_genes": n_genes, "repaired": True, "methods": {}}
    frames = []
    published = {"music", "dwls", "bisque", "scdc", "scdc_ensemble", "epic", "quantiseq", "bayesprism"}
    for m in build_methods(prefer_r=not a.python_only):
        if m.name not in a.methods:
            continue
        t0 = time.time()
        try:
            est = m.fit_predict(data)
        except Exception as exc:                                       # noqa: BLE001
            report["methods"][m.name] = {"failed": f"{type(exc).__name__}: {str(exc)[:200]}"}
            print(f"  {m.name:24s} FAILED {type(exc).__name__}: {str(exc)[:120]}", flush=True)
            continue
        impl = getattr(m, "implementation_", None)
        if impl is None:   # a Python class run directly (--python-only): never labelled as the package
            impl = ("python-reimplementation" + getattr(m, "variant_suffix_", "")) if m.name in published else "python"
        rec = {"implementation": impl,
               "fallback_reason": getattr(m, "fallback_reason_", None),
               "degenerate": bool(getattr(m, "degenerate_", False)),
               "seconds": round(time.time() - t0, 1),
               "finite_samples_by_type": {c: int(np.isfinite(est[c].to_numpy(float)).sum()) for c in est.columns}}
        for comp, truth in COMPS:
            if truth == "none":
                continue
            cols = idd.COMPARTMENTS[comp]
            x = est[cols]
            if x.isna().all().all():
                rec[f"{comp}|{truth}"] = {"rho": None, "why": "not modelled by this method"}
                continue
            r, n = idd.truth_rho(x.sum(axis=1), comp, tr)
            rec[f"{comp}|{truth}"] = {"rho": None if not np.isfinite(r) else round(r, 4), "n": n}
        rec["T>B"] = lymphoid_direction(est, tr)
        reg = registered["methods"].get(m.name, {})
        rec["registered"] = {"implementation": reg.get("implementation"), "tumour_rho": reg.get("spearman_vs_purity"),
                             "status": reg.get("skipped") or reg.get("failed") or "ok"}
        report["methods"][m.name] = rec
        frames.append(est.assign(method=m.name).rename_axis("sample_id").reset_index())
        print(f"  {m.name:24s} {rec['implementation']}: tumour {rec['Tumor|ABSOLUTE'].get('rho')} "
              f"(registered {rec['registered']['tumour_rho']}), leukocytes {rec['Leukocytes|LF'].get('rho')}, "
              f"T>B {rec['T>B'].get('T_exceeds_B')} [{rec['seconds']} s]", flush=True)
    stem = f"tcga_{a.cohort}_{a.reference}"
    (OUT / f"{stem}.json").write_text(json.dumps(report, indent=2, default=float))
    if frames:
        pd.concat(frames).to_csv(OUT / f"estimates_{a.cohort}_{a.reference}.csv", index=False)
    print(f"wrote results/repaired/{stem}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
