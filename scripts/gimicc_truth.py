"""
gimicc_truth.py -- prespecified/gimicc_truth_confirmation.md (with Addenda 1 and 2).

Does a glioma-specific methylation method (GIMiCC) confirm the lymphoid truth the study's headline
rests on -- methylation's T > B ordering, until now measured with EpiDISH and a BLOOD reference?

Runs (genuine GIMiCC 0.99.1 through R/run_gimicc.R; every output kept, so the job resumes):
  primary      each cohort at its default tumour type ("GBM" / "AST"), levels 4 and 5
  level 3      the default tumour type at level 3: T, NK and B straight from Layer 3A, which does not
               depend on the CD4/CD8 projection (Addendum 2)
  sensitivity  the other three tumour types at level 4 (structural invariance of lymphoid shares)
  negative     CpG labels permuted (must break GIMiCC Tumor vs ABSOLUTE)
  robustness   20 repeats dropping a further random 20% of the available library CpGs

Analysis, controls before the primary reading:
  structural invariance; tumour layer vs ABSOLUTE; Immune + Microglia vs the Thorsson leukocyte
  fraction; the negative control; probe-drop robustness; then Q1 (cohort-level T vs B), with EpiDISH
  per-sample agreement and the RNA-matched subsets.

    python3 scripts/gimicc_truth.py [--stage runs|analyse|all]
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402

G = config.RESULTS_DIR / "gimicc"
RUNS = G / "runs"
OUT = config.RESULTS_DIR / "gimicc_truth_confirmation.json"
COHORTS = ("gbm", "lgg")
DEFAULT_TYPE = {"gbm": "GBM", "lgg": "AST"}
TYPES = ("GBM", "AST", "AST-HG", "OLG")
N_REPEATS, DROP = 20, 0.20
LYMPH = ["T", "NK", "B"]
IMMUNE = ["CD4Tcell", "CD8Tcell", "Bcell", "NK", "Mono", "Neu"]
IMMUNE_H3 = ["Tcell", "Bcell", "NK", "Mono", "Neu"]


def immune_total(df: pd.DataFrame, cols: list[str]) -> tuple[pd.Series, int]:
    """Immune + Microglia per sample; a sample with a NaN in any of the columns is EXCLUDED (counted),
    never summed as if the NaN were zero (Addendum 2)."""
    sub = df[cols + ["Microglia"]]
    ok = sub.notna().all(axis=1)
    return sub[ok].sum(axis=1), int((~ok).sum())


def run(cohort: str, tumor_type: str, h: int, tag: str, src: Path | None = None, shuffle_seed: int | None = None) -> Path:
    out = RUNS / f"{cohort}_{tag}.csv"
    if out.exists():
        return out
    args = {"in": str(src or G / f"tcga_{cohort}_gimicc_cpgs.tsv"), "tumor_type": tumor_type, "h": h, "out": str(out)}
    if shuffle_seed is not None:
        args["shuffle_seed"] = shuffle_seed
    a = RUNS / f"{cohort}_{tag}.args.json"
    a.write_text(json.dumps(args))
    with open(RUNS / f"{cohort}_{tag}.log", "w") as log:
        r = subprocess.run(["Rscript", str(ROOT / "R" / "run_gimicc.R"), str(a)], stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0 or not out.exists():
        raise RuntimeError(f"GIMiCC run {cohort}/{tag} failed; see {RUNS / f'{cohort}_{tag}.log'}")
    return out


def robustness_inputs(cohort: str) -> list[Path]:
    """Repeat k keeps a random 80% of the complete-case library CpGs (seed RANDOM_SEED + k)."""
    src = G / f"tcga_{cohort}_gimicc_cpgs.tsv"
    full = pd.read_csv(src, sep="\t", index_col=0, na_values=["", "NA", "NaN"])
    cc = full.dropna(axis=0, how="any")
    paths = []
    for k in range(N_REPEATS):
        p = RUNS / f"{cohort}_robust{k:02d}.tsv"
        if not p.exists():
            rng = np.random.default_rng(config.RANDOM_SEED + k)
            keep = rng.choice(cc.index.to_numpy(), size=int(round(len(cc) * (1 - DROP))), replace=False)
            cc.loc[sorted(keep)].to_csv(p, sep="\t")
        paths.append(p)
    return paths


def stage_runs() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    for c in COHORTS:
        d = DEFAULT_TYPE[c]
        run(c, d, 4, f"{d}_h4")
        run(c, d, 5, f"{d}_h5")
        run(c, d, 3, f"{d}_h3")
        for t in TYPES:
            run(c, t, 4, f"{t}_h4")
        run(c, d, 4, f"{d}_h4_shuffled", shuffle_seed=config.RANDOM_SEED)
        print(f"{c}: primary, sensitivity and negative-control runs done", flush=True)
    for c in COHORTS:
        for k, p in enumerate(robustness_inputs(c)):
            run(c, DEFAULT_TYPE[c], 4, f"robust{k:02d}_h4", src=p)
        print(f"{c}: {N_REPEATS} robustness repeats done", flush=True)


# ---------------------------------------------------------------- analysis
def load(cohort: str, tag: str) -> pd.DataFrame:
    df = pd.read_csv(RUNS / f"{cohort}_{tag}.csv", index_col=0) / 100.0      # GIMiCC reports percent
    df.index = df.index.astype(str)
    return df


def lymph_shares(df: pd.DataFrame) -> pd.DataFrame:
    """Within-lymphoid composition per sample: T = CD4T + CD8T at level 4 (registered), or Layer 3A's
    Tcell at level 3 (Addendum 2). Samples with no lymphoid signal, or a NaN (noisy) row, are dropped --
    counted by the caller, never zero-filled."""
    t = df["Tcell"] if "Tcell" in df.columns else df["CD4Tcell"] + df["CD8Tcell"]
    L = pd.DataFrame({"T": t, "NK": df["NK"], "B": df["Bcell"]})
    tot = L.sum(axis=1)
    return L[tot > 0].div(tot[tot > 0], axis=0).dropna()


def ordering(mu: pd.Series) -> str:
    return ">".join(mu.sort_values(ascending=False).index)


def spearman(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    idx = a.dropna().index.intersection(b.dropna().index)
    if len(idx) < 10:
        return float("nan"), len(idx)
    return float(stats.spearmanr(a[idx], b[idx]).statistic), len(idx)


def truths(cohort: str) -> dict:
    from absolute_purity_yardstick import ABS_T, PURITY_COL, key4  # noqa: PLC0415
    from lymphoid_ordering import k4s  # noqa: PLC0415
    base = "" if cohort == "gbm" else "_lgg"
    a = pd.read_csv(ABS_T, sep="\t").dropna(subset=[PURITY_COL])
    a["k"] = a["sample"].map(key4)
    purity = a.drop_duplicates("k").set_index("k")[PURITY_COL].astype(float)   # key4: RNA barcodes (01A)
    # KEY FIX (2026-10-03, after the first analysis run matched 0 samples): ABSOLUTE's sample IDs carry
    # the vial letter (TCGA-xx-xxxx-01A-...), the Xena methylation IDs do not (TCGA-xx-xxxx-01), so
    # key4 can never match them. Methylation-side comparisons key both by k4s, as every other
    # methylation comparison in the project does; the first vial is kept where two exist.
    # D27 (found the same day by tests/test_gimicc_independent.py): 795 ABSOLUTE rows carry a Broad-
    # internal ID in `sample` ("GBM-TCGA-06-5416-Tumor-SM-1QETM") that no key function can parse; their
    # `array` column is the TCGA barcode in exactly the methylation IDs' form (TCGA-xx-xxxx-01).
    purity_meth = a.drop_duplicates("array").set_index("array")[PURITY_COL].astype(float)
    lf = pd.read_csv(config.IMMUNE_FRACTION_DIR / "TCGA_all_leuk_estimate.masked.20170107.tsv", sep="\t",
                     header=None, names=["study", "barcode", "lf"])
    lf["k"] = lf["barcode"].map(k4s)
    lf = pd.to_numeric(lf.drop_duplicates("k").set_index("k")["lf"], errors="coerce").dropna()
    m = pd.read_csv(config.RESULTS_DIR / f"methylation_celltypes{base}.csv", index_col=0)
    m.index = [k4s(i) for i in m.index]
    m = m[~m.index.duplicated()]
    epi = pd.DataFrame({"T": m["CD4T"] + m["CD8T"], "NK": m["NK"], "B": m["B"]})
    epi = epi.div(epi.sum(axis=1), axis=0)
    # RNA-matched = the registered comparisons' set: RNA samples with ABSOLUTE purity (as
    # extension_tcga selects them); intersected with GIMiCC's samples by the caller.
    with gzip.open(config.PROCESSED_DIR / f"tcga_{cohort}_bulk_cpm.csv.gz", "rt") as fh:
        rna = {k4s(c) for c in pd.read_csv(fh, index_col=0, nrows=1).columns if key4(c) in purity.index}
    return {"purity": purity, "purity_meth": purity_meth, "lf": lf, "epi": epi, "rna_samples": rna,
            "key4": key4, "k4s": k4s}


def rekey(s: pd.Series | pd.DataFrame, f) -> pd.Series | pd.DataFrame:
    s = s.copy()
    s.index = [f(i) for i in s.index]
    return s[~s.index.duplicated()]


def q1(df: pd.DataFrame, subset: set | None = None, k4s=None) -> dict:
    L = lymph_shares(df)
    if subset is not None:
        L = L[[k4s(i) in subset for i in L.index]]
    mu = L.mean()
    return {"n_samples_with_lymphoid": int(len(L)), "mean_within_lymphoid_share": mu.round(4).to_dict(),
            "ordering": ordering(mu), "T_exceeds_B": bool(mu["T"] > mu["B"]),
            "share_samples_T_over_B": round(float((L["T"] > L["B"]).mean()), 4) if len(L) else float("nan")}


def reading_of(ok_controls: bool, t_over_b: list[bool], robust: bool) -> tuple[str, str]:
    """The registered reading (the document's text, Addendum 2 item 3) and the symmetric one.
    A failing control -> INCONCLUSIVE; B >= T in either cohort -> NOT CONFIRMED (fragility stated with
    it); T > B in both and robust -> CONFIRMED; T > B in both but fragile -> INCONCLUSIVE."""
    if not ok_controls:
        reading = "INCONCLUSIVE"
    elif not all(t_over_b):
        reading = "NOT CONFIRMED"
    else:
        reading = "CONFIRMED" if robust else "INCONCLUSIVE"
    symmetric = "INCONCLUSIVE" if (not ok_controls or not robust) else ("CONFIRMED" if all(t_over_b) else "NOT CONFIRMED")
    return reading, symmetric


def stage_analyse() -> int:
    out: dict = {"rule": "prespecified/gimicc_truth_confirmation.md (with Addendum 1)",
                 "implementation": "R:GIMiCC 0.99.1 (github SalasLab/GIMiCC @26cb8a15; library EH9483)",
                 "extraction": json.loads((G / "extraction_summary.json").read_text()), "cohorts": {}}
    ok_controls = True
    for c in COHORTS:
        d = DEFAULT_TYPE[c]
        tr = truths(c)
        prim = load(c, f"{d}_h4")
        n_nan = int(prim.isna().any(axis=1).sum())
        res: dict = {"default_tumor_type": d, "n_samples": int(len(prim)), "n_noisy_nan_rows": n_nan}
        # 1 structural invariance across tumour types
        Ls = {t: lymph_shares(load(c, f"{t}_h4")) for t in TYPES}
        common = set.intersection(*(set(v.index) for v in Ls.values()))
        diffs = [float((Ls[t].loc[sorted(common)] - Ls[d].loc[sorted(common)]).abs().to_numpy().max()) for t in TYPES]
        res["structural_invariance"] = {"max_abs_diff_within_lymphoid_across_types": max(diffs),
                                        "n_samples_compared": len(common), "pass": max(diffs) < 1e-9}
        res["tissue_fractions_by_type"] = {t: {"Tumor_mean": round(float(load(c, f"{t}_h4")["Tumor"].mean()), 4)} for t in TYPES}
        # 2 tumour layer vs ABSOLUTE
        r2, n2 = spearman(rekey(prim["Tumor"], tr["k4s"]), tr["purity_meth"])
        res["control_tumor_vs_absolute"] = {"spearman": round(r2, 4), "n": n2, "pass": r2 >= 0.40}
        # 3 immune total vs Thorsson LF (NaN rows excluded and counted -- Addendum 2)
        h3 = load(c, f"{d}_h3")
        imm_mg, n_ex = immune_total(prim, IMMUNE)
        r3, n3 = spearman(rekey(imm_mg, tr["k4s"]), tr["lf"])
        imm_only = imm_mg - prim.loc[imm_mg.index, "Microglia"]
        r3b, _ = spearman(rekey(imm_only, tr["k4s"]), tr["lf"])
        imm3, n_ex3 = immune_total(h3, IMMUNE_H3)
        r3c, n3c = spearman(rekey(imm3, tr["k4s"]), tr["lf"])
        res["control_immune_vs_leukocyte_fraction"] = {"spearman_immune_plus_microglia": round(r3, 4),
                                                       "spearman_immune_only": round(r3b, 4), "n": n3,
                                                       "n_excluded_nan": n_ex, "pass": r3 >= 0.40,
                                                       "level3_spearman_immune_plus_microglia": round(r3c, 4),
                                                       "level3_n": n3c, "level3_n_excluded_nan": n_ex3}
        # 4 negative control
        sh = load(c, f"{d}_h4_shuffled")
        r4, n4 = spearman(rekey(sh["Tumor"], tr["k4s"]), tr["purity_meth"])
        res["negative_control_shuffled_cpgs"] = {"tumor_vs_absolute_spearman": round(r4, 4), "n": n4, "pass": abs(r4) < 0.20}
        for name, n_ in (("tumour layer vs ABSOLUTE", n2), ("immune vs leukocyte fraction", n3), ("shuffled CpGs", n4)):
            if n_ < 100:
                raise SystemExit(f"BLOCKED: control '{name}' matched only {n_} samples in {c} -- a sample-key "
                                 f"defect, not a result; fix the matching before any reading")
        # robustness (Addendum 1)
        reps = [q1(load(c, f"robust{k:02d}_h4"))["T_exceeds_B"] for k in range(N_REPEATS)
                if (RUNS / f"{c}_robust{k:02d}_h4.csv").exists()]
        res["probe_drop_robustness"] = {"n_repeats": len(reps), "n_T_over_B": int(sum(reps))}
        ok_controls &= res["control_tumor_vs_absolute"]["pass"] and res["control_immune_vs_leukocyte_fraction"]["pass"] \
            and res["negative_control_shuffled_cpgs"]["pass"] and res["structural_invariance"]["pass"]
        # Q1
        res["Q1_all_methylation_samples"] = q1(prim)
        res["Q1_rna_matched_samples"] = q1(prim, tr["rna_samples"], tr["k4s"])
        # Addendum 2: level 3 (Layer 3A's T, independent of the CD4/CD8 projection), and the samples
        # whose level-4 T is lost (NaN) or silently zero while Layer 3A's T is positive
        t4 = prim["CD4Tcell"] + prim["CD8Tcell"]
        t3 = h3["Tcell"].reindex(prim.index)
        res["Q1_level3_sensitivity"] = q1(h3)
        res["Q1_level3_rna_matched"] = q1(h3, tr["rna_samples"], tr["k4s"])
        res["level4_T_lost_while_level3_T_positive"] = {
            "n_nan_at_level4": int((t4.isna() & (t3 > 0)).sum()),
            "n_zero_at_level4": int(((t4 == 0) & (t3 > 0)).sum()),
            "n_level3_T_positive": int((t3 > 0).sum()),
            "max_abs_diff_T_where_both_defined": float((t4 - t3)[t4.notna() & t3.notna() & (t4 > 0)].abs().max())
            if bool((t4.notna() & t3.notna() & (t4 > 0)).any()) else float("nan")}
        res["level3_direction_agrees_with_level4"] = (res["Q1_level3_sensitivity"]["T_exceeds_B"]
                                                      == res["Q1_all_methylation_samples"]["T_exceeds_B"])
        L = rekey(lymph_shares(prim), tr["k4s"])
        epi = tr["epi"]
        rT, nT = spearman(L["T"], epi["T"])
        rB, nB = spearman(L["B"], epi["B"])
        mu_epi = epi.loc[epi.index.intersection(L.index)].mean()
        res["agreement_with_epidish"] = {"spearman_T_share": round(rT, 4), "spearman_B_share": round(rB, 4), "n": nT,
                                         "epidish_ordering_same_samples": ordering(mu_epi),
                                         "epidish_T_exceeds_B_same_samples": bool(mu_epi["T"] > mu_epi["B"])}
        h5 = load(c, f"{d}_h5")
        res["CD4_vs_CD8_mean_tissue_fraction"] = {"CD4": round(float(prim["CD4Tcell"].mean()), 6),
                                                  "CD8": round(float(prim["CD8Tcell"].mean()), 6)}
        res["level5_mean_tissue_fraction"] = {k: round(float(h5[k].mean()), 6) for k in
                                              ("CD4nv", "CD4mem", "Treg", "CD8nv", "CD8mem", "Bnv", "Bmem", "NK")}
        res["lymphoid_tissue_fraction_mean"] = round(float(prim[["CD4Tcell", "CD8Tcell", "Bcell", "NK"]].sum(axis=1).mean()), 6)
        out["cohorts"][c] = res
    both = [out["cohorts"][c]["Q1_all_methylation_samples"]["T_exceeds_B"] for c in COHORTS]
    for c in COHORTS:
        pr = out["cohorts"][c]["probe_drop_robustness"]
        held = pr["n_T_over_B"] if out["cohorts"][c]["Q1_all_methylation_samples"]["T_exceeds_B"] else N_REPEATS - pr["n_T_over_B"]
        pr["n_direction_held"] = int(held)
        pr["robust"] = bool(pr["n_repeats"] == N_REPEATS and held >= 18)
    robust = all(out["cohorts"][c]["probe_drop_robustness"]["robust"] for c in COHORTS)
    reading, symmetric = reading_of(ok_controls, both, robust)
    out["controls_pass"], out["robust"], out["Q1_reading"] = ok_controls, robust, reading
    out["Q1_reading_if_robustness_required_in_both_directions"] = symmetric
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({c: {k: v for k, v in out["cohorts"][c].items() if k.startswith(("control", "negative", "structural", "probe"))}
                      for c in COHORTS}, indent=1, default=float))
    print("controls pass:", ok_controls, "| robust:", robust)
    print("Q1:", {c: out["cohorts"][c]["Q1_all_methylation_samples"] for c in COHORTS})
    print("Q1 level 3 (Addendum 2):", {c: out["cohorts"][c]["Q1_level3_sensitivity"] for c in COHORTS})
    print("level-4 T lost while level-3 T positive:", {c: out["cohorts"][c]["level4_T_lost_while_level3_T_positive"] for c in COHORTS})
    print("Q1 reading:", reading, "| symmetric-robustness reading:", symmetric)
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stage", choices=["runs", "analyse", "all"], default="all")
    a = ap.parse_args()
    if a.stage in ("runs", "all"):
        stage_runs()
    if a.stage in ("analyse", "all"):
        return stage_analyse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
