"""
cptac_wgs_purity.py -- prespecified/cptac_wgs_purity_secondary.md (S2, SECONDARY: registered after the primary
CPTAC DNA control failed). Whole-genome AscatNGS tumour purity as the per-sample truth for CPTAC tumour content.

Stages:
  fetch    GDC file records: tumor_purity / tumor_ploidy of each case's one open AscatNGS segment file, and the
           pathology percent_tumor_nuclei of the slides of each case's analysed tumour sample(s)
  gimicc   S2.2(b): GIMiCC re-run on the samples covering >= 90% of the library CpGs (complete cases among them)
  analyse  S2.1-S2.6 -> results/cptac_wgs_purity.json

    python3 scripts/cptac_wgs_purity.py --stage fetch|gimicc|analyse|all
"""
from __future__ import annotations

import argparse
import json
import ssl
import subprocess
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import cptac_per_sample as cp  # noqa: E402

SRC = config.CPTAC_GBM_DIR
OUT = cp.OUT
PURITY = SRC / "wgs_ascat_purity.tsv"
PATHOLOGY = SRC / "pathology_tumor_nuclei.tsv"
GDC = "https://api.gdc.cancer.gov"
UA = {"User-Agent": "TumorAnatomyResearch/1.0 (open-data research use)"}
MIN_COVERAGE = 0.90        # S2.2(b), registered: any threshold from 68% to 95% selects the same samples
EDGE = 0.99                # purity values at or above this may be fits to a genome without aberrations
BAR_TRUTH, BAR_AGREE = 0.40, 0.60


def _get(endpoint: str, filt: dict, fields: list[str], size: int = 500) -> list[dict]:
    import certifi  # noqa: PLC0415
    q = urllib.parse.urlencode({"filters": json.dumps(filt), "size": str(size), "format": "json",
                                "fields": ",".join(fields)})
    req = urllib.request.Request(f"{GDC}/{endpoint}?{q}", headers=UA)
    ctx = ssl.create_default_context(cafile=certifi.where())
    return json.load(urllib.request.urlopen(req, context=ctx, timeout=120))["data"]["hits"]


def _tumour_samples(hit: dict) -> list[str]:
    return sorted({s["submitter_id"] for c in hit["cases"] for s in c.get("samples", [])
                   if "Tumor" in (s.get("sample_type") or "")})


# ------------------------------------------------------------------------------------------- fetch
def stage_fetch() -> None:
    man = cp.manifest()
    cases = sorted(man["cases"]) if isinstance(man["cases"][0], str) else sorted(c["case"] for c in man["cases"])
    samp = ["cases.submitter_id", "cases.samples.submitter_id", "cases.samples.sample_type"]
    hits = _get("files", {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.submitter_id", "value": cases}},
        {"op": "=", "content": {"field": "data_type", "value": "Allele-specific Copy Number Segment"}},
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": "AscatNGS"}},
        {"op": "=", "content": {"field": "access", "value": "open"}}]},
        ["file_id", "file_name", "md5sum", "data_release", "tumor_purity", "tumor_ploidy",
         "analysis.workflow_type", "analysis.workflow_version"] + samp)
    rows = {}
    for h in hits:
        case = h["cases"][0]["submitter_id"]
        if case in rows:
            raise SystemExit(f"BLOCKED: {case} has two AscatNGS records; the rule does not choose")
        rows[case] = {"case": case, "tumor_purity": h.get("tumor_purity"), "tumor_ploidy": h.get("tumor_ploidy"),
                      "tumour_samples": ";".join(_tumour_samples(h)), "file_id": h["file_id"],
                      "file_name": h["file_name"], "md5sum": h.get("md5sum"), "data_release": h.get("data_release"),
                      "workflow": f"{h['analysis']['workflow_type']} {h['analysis'].get('workflow_version', '')}".strip()}
    missing = [c for c in cases if c not in rows]
    pur = pd.DataFrame([rows[c] for c in cases if c in rows])
    pur.to_csv(PURITY, sep="\t", index=False)
    # the analysed tumour sample(s) per case are those of its bulk RNA file (identical to the WGS samples, verified)
    bulk_ids = {r["file_id"]: r["case"] for r in man["files"] if r["kind"] == "bulk"}
    analysed = {bulk_ids[h["file_id"]]: _tumour_samples(h)
                for h in _get("files", {"op": "in", "content": {"field": "file_id", "value": list(bulk_ids)}},
                              ["file_id"] + samp)}
    slides = _get("cases", {"op": "in", "content": {"field": "submitter_id", "value": cases}},
                  ["submitter_id", "samples.submitter_id", "samples.portions.slides.submitter_id",
                   "samples.portions.slides.percent_tumor_nuclei", "samples.portions.slides.section_location"])
    prow = []
    for c in slides:
        for s in c.get("samples", []):
            if s["submitter_id"] not in analysed.get(c["submitter_id"], []):
                continue
            for p in s.get("portions", []) or []:
                for sl in p.get("slides", []) or []:
                    prow.append({"case": c["submitter_id"], "sample": s["submitter_id"],
                                 "slide": sl.get("submitter_id"), "section_location": sl.get("section_location"),
                                 "percent_tumor_nuclei": sl.get("percent_tumor_nuclei")})
    pd.DataFrame(prow).sort_values(["case", "slide"]).to_csv(PATHOLOGY, sep="\t", index=False)
    prov = {"source": "NCI GDC API (open access), project CPTAC-3", "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "purity_field": "files.tumor_purity on data_type 'Allele-specific Copy Number Segment', workflow AscatNGS",
            "pathology_field": "cases.samples.portions.slides.percent_tumor_nuclei, slides of the bulk file's tumour samples",
            "cases_requested": len(cases), "cases_with_record": len(rows), "cases_without_record": missing,
            "data_releases": sorted({str(r["data_release"]) for r in rows.values()}),
            "rule": "prespecified/cptac_wgs_purity_secondary.md"}
    (SRC / "wgs_ascat_purity.provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"AscatNGS records: {len(rows)} of {len(cases)} cases; without: {missing or 'none'}")
    print(f"pathology slides of analysed samples: {len(prow)} rows -> {PATHOLOGY.relative_to(config.PROJECT_ROOT)}")


# ------------------------------------------------------------------------------------------- gimicc (S2.2 b)
def covered(m: pd.DataFrame, threshold: float = MIN_COVERAGE) -> list[str]:
    """Samples (columns) whose share of non-missing library CpGs is >= threshold. Nothing is imputed."""
    cov = m.notna().mean()
    return sorted(cov[cov >= threshold].index)


def reading_s22(ra: float | None, rb: float | None) -> str:
    """S2.2's registered reading: attributed to coverage only if the complete-sample run reaches the bar and the
    all-sample run does not."""
    if rb is not None and rb >= BAR_TRUTH and (ra is None or ra < BAR_TRUTH):
        return "failed primary control attributed to methylation coverage"
    return "not attributed to methylation coverage by this rule"


def stage_gimicc() -> None:
    m = pd.read_csv(OUT / "cptac_gimicc_cpgs.tsv", sep="\t", index_col=0, na_values=["NA"])
    keep = covered(m)
    print(f"library coverage >= {MIN_COVERAGE:.0%}: {len(keep)} of {m.shape[1]} samples; "
          f"left out {sorted(set(m.columns) - set(keep))}")
    src = OUT / "cptac_gimicc_cpgs_cov90.tsv"
    m[keep].rename_axis("sample").to_csv(src, sep="\t", na_rep="NA")
    args = {"in": str(src), "tumor_type": "GBM", "h": 4, "out": str(OUT / "gimicc_h4_cov90.csv")}
    (OUT / "gimicc_args_cov90.json").write_text(json.dumps(args))
    r = subprocess.run(["Rscript", str(ROOT / "R" / "run_gimicc.R"), str(OUT / "gimicc_args_cov90.json")],
                       capture_output=True, text=True)
    (OUT / "gimicc_cov90.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0:
        raise SystemExit(f"BLOCKED: GIMiCC failed; see {OUT / 'gimicc_cov90.log'}")
    print([ln for ln in r.stdout.splitlines() if "gimicc" in ln.lower()])


# ------------------------------------------------------------------------------------------- analyse
def _rho(x, y, rng) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    # np.ptp, not np.std: the std of identical floats can be ~1e-16 rather than 0
    if ok.sum() < 4 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return {"rho": None, "n": int(ok.sum()), "why": "too few or constant"}
    r, p = cp.perm_p(x[ok], y[ok], rng)
    if not np.isfinite(r):
        return {"rho": None, "n": int(ok.sum()), "why": "undefined"}
    return {"rho": round(r, 4), "perm_p": p, "n": int(ok.sum())}


def _gimicc(name: str) -> pd.Series:
    g = pd.read_csv(OUT / name, index_col=0)
    g.index = [str(i) for i in g.index]
    return g["Tumor"] / 100.0                      # GIMiCC reports percent


def _analyse(purity: pd.Series, rng, label: str) -> dict:
    """S2.1-S2.5 on one set of cases (all, or without edge values)."""
    out: dict = {"cases": label, "n_cases_with_purity": int(purity.notna().sum())}
    pooled = set()
    for c, s in pd.read_csv(PURITY, sep="\t", index_col=0)["tumour_samples"].items():
        if isinstance(s, str) and ";" in s:
            pooled.add(c)
    # S2.1 the nuclei truth against WGS purity
    s21 = {}
    for tag, name in (("", "registered"), ("_s1", "s1_nk_without_ncam1")):
        truth = pd.read_csv(OUT / f"snrna_truth{tag}.csv", index_col=0)
        det = json.loads((OUT / f"snrna_truth{tag}_detail.json").read_text())
        scored = [c for c in truth.index if det[c]["share_mapped"] >= cp.MIN_MAPPED and c in purity.dropna().index]
        single = [c for c in scored if c not in pooled]
        s21[name] = {"all_scored": _rho(truth.loc[scored, "Tumor"], purity[scored], rng),
                     "single_piece_cases (descriptive)": _rho(truth.loc[single, "Tumor"], purity[single], rng),
                     "scored_cases": scored}
    r = s21["registered"]["all_scored"]["rho"]
    s21["pass"] = bool(r is not None and r >= BAR_TRUTH)
    out["S2_1_nuclei_truth_vs_wgs"] = s21
    # S2.2 the methylation instrument
    ga, gb = _gimicc("gimicc_h4.csv"), (_gimicc("gimicc_h4_cov90.csv") if (OUT / "gimicc_h4_cov90.csv").exists() else None)
    ia = [c for c in ga.index if c in purity.dropna().index]
    s22 = {"a_complete_case_769_cpgs": _rho(ga[ia], purity[ia], rng)}
    if gb is not None:
        ib = [c for c in gb.index if c in purity.dropna().index]
        s22["b_coverage_ge_90pct"] = _rho(gb[ib], purity[ib], rng)
        s22["reading"] = reading_s22(s22["a_complete_case_769_cpgs"]["rho"], s22["b_coverage_ge_90pct"]["rho"])
        # the registered claim that any threshold from 68% to 95% selects the same samples, checked on the input
        m = pd.read_csv(OUT / "cptac_gimicc_cpgs.tsv", sep="\t", index_col=0, na_values=["NA"])
        cov = m.notna().mean()
        sets = {f"{t:.2f}": covered(m, t) for t in (0.68, MIN_COVERAGE, 0.95)}
        s22["threshold_check"] = {"coverage_by_sample": {k: round(float(v), 4) for k, v in cov.sort_values().items()},
                                  "same_samples_from_0.68_to_0.95": len({tuple(v) for v in sets.values()}) == 1}
        # post hoc, descriptive: the nuclei truth against the working methylation purity
        truth = pd.read_csv(OUT / "snrna_truth.csv", index_col=0)
        ic = [c for c in s21["registered"]["scored_cases"] if c in gb.index]
        s22["nuclei_vs_b (post hoc, descriptive)"] = _rho(truth.loc[ic, "Tumor"], gb[ic], rng)
    out["S2_2_methylation_instrument"] = s22
    # S2.3-S2.5 bulk methods against WGS purity
    lb = pd.read_csv(config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv")
    acs = dict(zip(lb["method"], lb["acs"]))
    reg12 = list(json.loads((config.RESULTS_DIR / "yardstick_agreement.json").read_text())["per_method"])
    for ref, yard in (("frozen", "absolute_purity_yardstick.json"), ("h5ad", "absolute_purity_yardstick_h5ad.json")):
        est = pd.read_csv(OUT / f"estimates_{ref}.csv")
        meta = json.loads((OUT / f"methods_{ref}.json").read_text())["methods"]
        yard_methods = json.loads((config.RESULTS_DIR / yard).read_text())["methods"]
        tcga = {m: v["spearman_vs_purity"] for m, v in yard_methods.items() if v.get("spearman_vs_purity") is not None}
        R: dict = {"methods": {}}
        for m, gm in est.groupby("method"):
            gm = gm.set_index("sample")
            cs = [c for c in purity.dropna().index if c in gm.index]
            x, y = gm.loc[cs, "Tumor"].to_numpy(float), purity[cs].to_numpy(float)
            rec = _rho(x, y, rng)
            rec.update({"recovery_slope": round(float(np.polyfit(y, x, 1)[0]), 4),
                        "mean_bias": round(float(np.mean(x - y)), 4),
                        "degenerate": meta[m].get("degenerate"), "implementation": meta[m].get("implementation"),
                        "tcga_gbm_rho_vs_absolute": tcga.get(m)})
            R["methods"][m] = rec
        rhos = {m: v["rho"] for m, v in R["methods"].items() if v["rho"] is not None}
        both = [m for m in rhos if m in tcga]
        med = float(np.median(list(rhos.values())))
        R["S2_3"] = {"median_rho": round(med, 4), "n_methods": len(rhos),
                     "reading": "tracks" if med >= BAR_TRUTH else "not recovered",
                     "tcga_gbm_median_same_methods": round(float(np.median([tcga[m] for m in both])), 4),
                     "median_recovery_slope": round(float(np.median([v["recovery_slope"] for v in R["methods"].values()])), 4)}
        if len(both) >= 6:
            t = stats.spearmanr([rhos[m] for m in both], [tcga[m] for m in both])
            R["S2_4_ranking_transfer"] = {"n_methods": len(both), "spearman": round(float(t.statistic), 4),
                                          "p": round(float(t.pvalue), 4),
                                          "reading": "transfers" if t.statistic >= BAR_AGREE else "does not transfer"}
            # post hoc, descriptive: only methods run by the same implementation in both cohorts
            impl = {m: v.get("implementation") for m, v in yard_methods.items()}
            same = [m for m in both if impl.get(m) == R["methods"][m]["implementation"]]
            if len(same) >= 6 and len(same) < len(both):
                t = stats.spearmanr([rhos[m] for m in same], [tcga[m] for m in same])
                R["S2_4_matched_implementation (post hoc, descriptive)"] = {
                    "n_methods": len(same), "left_out": sorted(set(both) - set(same)),
                    "spearman": round(float(t.statistic), 4), "p": round(float(t.pvalue), 4)}
        ms = [m for m in reg12 if m in rhos and m in acs]
        if len(ms) >= 6:
            t = stats.spearmanr([acs[m] for m in ms], [rhos[m] for m in ms])
            R["S2_5_acs_vs_wgs_accuracy"] = {"n_methods": len(ms), "methods": ms, "spearman": round(float(t.statistic), 4),
                                             "p": round(float(t.pvalue), 4),
                                             "reading": ("anatomy ranks methods by accuracy here" if t.statistic >= BAR_AGREE
                                                         else "anatomy does not rank methods by accuracy here")}
        out[ref] = R
    return out


def stage_analyse() -> int:
    rng = np.random.default_rng(config.RANDOM_SEED)
    pur = pd.read_csv(PURITY, sep="\t", index_col=0)
    purity = pd.to_numeric(pur["tumor_purity"], errors="coerce")
    edge = sorted(purity[purity >= EDGE].index)
    res: dict = {"rule": "prespecified/cptac_wgs_purity_secondary.md (SECONDARY; registered after the primary was read)",
                 "truth": {"n_cases": int(len(purity)), "n_with_value": int(purity.notna().sum()),
                           "excluded_no_value": sorted(purity[purity.isna()].index),
                           "range": [round(float(purity.min()), 4), round(float(purity.max()), 4)],
                           "median": round(float(purity.median()), 4), "cases_at_or_above_0.99": edge,
                           "ploidy_range": [float(pur["tumor_ploidy"].min()), float(pur["tumor_ploidy"].max())]},
                 "all_cases": _analyse(purity, rng, "all cases with a value")}
    if edge:
        res["without_edge_values"] = _analyse(purity.drop(edge), rng, f"without {edge}")
    # S2.6 pathology, descriptive
    if PATHOLOGY.exists():
        pa = pd.read_csv(PATHOLOGY, sep="\t")
        pn = pd.to_numeric(pa["percent_tumor_nuclei"], errors="coerce").groupby(pa["case"]).mean() / 100.0
        truth = pd.read_csv(OUT / "snrna_truth.csv", index_col=0)
        sc = res["all_cases"]["S2_1_nuclei_truth_vs_wgs"]["registered"]["scored_cases"]
        cw = [c for c in pn.dropna().index if c in purity.dropna().index]
        cs = [c for c in sc if c in pn.dropna().index]
        res["S2_6_pathology_descriptive"] = {
            "n_cases": int(pn.notna().sum()), "range": [round(float(pn.min()), 3), round(float(pn.max()), 3)],
            "vs_wgs_purity": _rho(pn[cw], purity[cw], rng), "vs_snrna_tumour_share": _rho(pn[cs], truth.loc[cs, "Tumor"], rng)}
    # per-sample table for the figure
    tab = pd.DataFrame({"wgs_purity": purity})
    for tag in ("", "_s1"):
        tab[f"snrna_tumour{tag or '_registered'}"] = pd.read_csv(OUT / f"snrna_truth{tag}.csv", index_col=0)["Tumor"]
    tab["gimicc_769cpg"] = _gimicc("gimicc_h4.csv")
    if (OUT / "gimicc_h4_cov90.csv").exists():
        tab["gimicc_cov90"] = _gimicc("gimicc_h4_cov90.csv")
    if PATHOLOGY.exists():
        tab["pathology_tumor_nuclei"] = pn
    tab["pooled_pieces"] = pur["tumour_samples"].astype(str).str.contains(";")
    m = pd.read_csv(OUT / "cptac_gimicc_cpgs.tsv", sep="\t", index_col=0, na_values=["NA"])
    tab["methylation_library_coverage"] = m.notna().mean()
    tab.rename_axis("case").to_csv(OUT / "wgs_purity_per_sample.csv")
    (config.RESULTS_DIR / "cptac_wgs_purity.json").write_text(json.dumps(res, indent=2, default=float))
    a = res["all_cases"]
    print("truth:", res["truth"])
    print("S2.1:", {k: v["all_scored"] for k, v in a["S2_1_nuclei_truth_vs_wgs"].items() if isinstance(v, dict)},
          "pass:", a["S2_1_nuclei_truth_vs_wgs"]["pass"])
    print("S2.2:", a["S2_2_methylation_instrument"])
    for ref in ("frozen", "h5ad"):
        print(f"{ref} S2.3:", a[ref]["S2_3"])
        print(f"{ref} S2.4:", a[ref].get("S2_4_ranking_transfer"))
        print(f"{ref} S2.5:", {k: v for k, v in a[ref].get("S2_5_acs_vs_wgs_accuracy", {}).items() if k != "methods"})
    if "S2_6_pathology_descriptive" in res:
        print("S2.6:", res["S2_6_pathology_descriptive"])
    print("wrote results/cptac_wgs_purity.json and results/cptac/wgs_purity_per_sample.csv")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stage", choices=["fetch", "gimicc", "analyse", "all"], default="all")
    a = ap.parse_args()
    if a.stage in ("fetch", "all"):
        stage_fetch()
    if a.stage in ("gimicc", "all"):
        stage_gimicc()
    if a.stage in ("analyse", "all"):
        return stage_analyse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
