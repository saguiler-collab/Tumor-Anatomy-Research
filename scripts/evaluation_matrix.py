"""
evaluation_matrix.py -- every method against every truth, in both cohorts and under both reference builds,
with the implementation that ran, the sample count, and the reason for every empty cell.

It reads saved estimates and truths only and computes no new estimate. The sample joins are the published
ones (`identifiability_diagnostics.truths` / `truth_rho`), and every tumour-content cell is checked against
the registered yardstick (`absolute_purity_yardstick*.json`), so a disagreement in the joins cannot pass
silently.

Truths (columns):
  ACS          anatomic concordance on Ivy GAP -- glioblastoma only, donor-level reference (registered leaderboard)
  ABSOLUTE     TCGA DNA copy-number purity, against the tumour estimate
  LF           Thorsson et al.'s methylation leukocyte fraction, against all leukocytes
  EpiDISH      methylation T, B, NK (the registered lymphoid truth), against lymphoid / T / B / NK
  T>B          does the method's cohort-mean T exceed its B? (flow cytometry and two atlases: it should)
  CPTAC WGS    AscatNGS whole-genome purity, 18 glioblastomas (secondary analysis S2)

    python3 scripts/evaluation_matrix.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import identifiability_diagnostics as idd  # noqa: E402

R = config.RESULTS_DIR
OUT = R / "evaluation_matrix.json"
DOC = config.PROJECT_ROOT / "docs" / "EVALUATION_MATRIX.md"
ARMS = [("gbm", "frozen", ""), ("gbm", "h5ad", "_h5ad"), ("lgg", "frozen", "_lgg"), ("lgg", "h5ad", "_lgg_h5ad")]
COMPS = [("Tumor", "ABSOLUTE"), ("Leukocytes", "LF"), ("Lymphoid", "EpiDISH"), ("T_cell", "EpiDISH"),
         ("B_cell", "EpiDISH"), ("NK_cell", "EpiDISH")]
TOL = 1e-3          # the recomputed tumour rho must equal the registered yardstick's to this


GIM_RUN = {"gbm": "GBM_h4", "lgg": "AST_h4"}          # gimicc_truth.DEFAULT_TYPE at level 4: the primary runs
GIM_COMPS = ["Tumor", "Leukocytes", "Lymphoid", "T_cell", "B_cell", "NK_cell"]


def gimicc_truths(cohort: str) -> dict:
    """GIMiCC's per-sample fractions (gimicc_truth.load; its immune total excludes, never zero-fills, a NaN row),
    keyed by k4s as every methylation comparison in the project is."""
    import gimicc_truth as gt  # noqa: PLC0415
    df = gt.load(cohort, GIM_RUN[cohort])
    imm, _ = gt.immune_total(df, gt.IMMUNE)
    t = df["CD4Tcell"] + df["CD8Tcell"]
    out = {"Tumor": df["Tumor"], "Leukocytes": imm, "Lymphoid": t + df["Bcell"] + df["NK"], "T_cell": t,
           "B_cell": df["Bcell"], "NK_cell": df["NK"]}
    rek = {}
    for k, v in out.items():
        v = v.dropna().copy()
        v.index = [idd.k4s(i) for i in v.index]
        rek[k] = v[~v.index.duplicated()]
    return rek


def gim_rho(est_c: pd.Series, truth: pd.Series) -> tuple[float, int]:
    e = est_c.copy()
    e.index = [idd.k4s(i) for i in e.index]
    e = e[~e.index.duplicated()].dropna()
    shared = e.index.intersection(truth.index)
    if len(shared) < 10 or e[shared].nunique() < 2:
        return float("nan"), len(shared)
    from scipy import stats  # noqa: PLC0415
    return float(stats.spearmanr(e[shared], truth[shared]).statistic), len(shared)


def J(name: str) -> dict:
    p = R / name
    return json.loads(p.read_text()) if p.exists() else {}


def build() -> dict:
    out: dict = {"what_this_is": __doc__.split("\n")[1].strip(), "arms": {}, "gaps": [], "join_checks": []}
    truth_cache: dict = {}
    gim_cache = {c: gimicc_truths(c) for c in ("gbm", "lgg")}
    out["gimicc_runs"] = GIM_RUN
    for cohort, ref, tag in ARMS:
        arm = f"{cohort}|{ref}"
        yard = J(f"absolute_purity_yardstick{tag}.json").get("methods", {})
        lymph = J(f"lymphoid_ordering{tag}.json").get("methods", {})
        est = pd.read_csv(R / f"estimates_full{tag}.csv")
        est = est.rename(columns={"sample": "sample_id"})          # the donor-level arms name the column 'sample'
        if cohort not in truth_cache:
            truth_cache[cohort] = idd.truths(cohort)
        tr = truth_cache[cohort]
        cells: dict = {}
        for m in sorted(set(yard) | set(est["method"].unique())):
            y = yard.get(m, {})
            c: dict = {"implementation": y.get("implementation"), "degenerate": bool(y.get("degenerate", False))}
            gm = est[est["method"] == m].set_index("sample_id")
            if gm.empty:
                why = y.get("skipped") or y.get("failed") or y.get("why") or "no estimates written"
                c["status"] = f"no estimates: {why}"
                out["gaps"].append({"arm": arm, "method": m, "what": "every truth", "why": why})
                cells[m] = c
                continue
            c["status"] = "ok"
            X = gm[list(config.CELL_TYPES)].to_numpy(float)
            Xf = X[~np.isnan(X).any(axis=1)]
            c["per_sample_solves"] = {          # OPEN_DEFECTS D7: degenerate solves that would otherwise pass unseen
                "n": int(len(X)), "nan_rows": int(len(X) - len(Xf)),
                "all_in_one_type": int(np.isclose(Xf.max(axis=1), 1.0, atol=1e-6).sum()) if len(Xf) else 0,
                "exactly_uniform": int(np.all(np.abs(Xf - 1 / Xf.shape[1]) < 1e-6, axis=1).sum()) if len(Xf) else 0,
                "duplicate_rows": int(len(Xf) - len(np.unique(np.round(Xf, 8), axis=0))) if len(Xf) else 0}
            for comp, truth in COMPS:
                r, n = idd.truth_rho(idd.comp(gm, comp), comp, tr)
                c[f"{comp}|{truth}"] = {"rho": None if not np.isfinite(r) else round(r, 4), "n": n}
                if not np.isfinite(r):
                    out["gaps"].append({"arm": arm, "method": m, "what": f"{comp} vs {truth}",
                                        "why": f"undefined (n = {n}; constant estimate or < 10 shared samples)"})
            for comp in GIM_COMPS:
                r, n = gim_rho(idd.comp(gm, comp), gim_cache[cohort][comp])
                c[f"{comp}|GIMiCC"] = {"rho": None if not np.isfinite(r) else round(r, 4), "n": n}
            reg = y.get("spearman_vs_purity")
            mine = c["Tumor|ABSOLUTE"]["rho"]
            if reg is not None and mine is not None:
                ok = abs(reg - mine) <= TOL
                out["join_checks"].append({"arm": arm, "method": m, "registered": reg, "recomputed": mine,
                                           "n_registered": y.get("n"), "n_recomputed": c["Tumor|ABSOLUTE"]["n"],
                                           "agree": ok})
            L = lymph.get(m, {})
            if L:
                c["T>B"] = {"T_exceeds_B": L.get("T_exceeds_B"), "ordering": L.get("ordering"),
                            "samples_with_no_lymphoid_signal": L.get("n_no_lymphoid_signal"), "n": L.get("n")}
            cells[m] = c
        out["arms"][arm] = cells
    # anatomy: glioblastoma, donor-level reference only
    lb = pd.read_csv(R / "anatomic" / "acs_leaderboard.csv")
    out["acs"] = {r.method: {"acs": round(float(r.acs), 4), "ci": [round(float(r.ci_low), 4), round(float(r.ci_high), 4)],
                             "comparable": bool(r.comparable), "is_control": bool(r.is_control),
                             "degenerate": bool(r.degenerate), "implementation": r.implementation,
                             "tie_group": None if str(r.acs_tie_group) in ("-", "nan") else r.acs_tie_group}
                  for r in lb.itertuples()}
    # CPTAC whole-genome purity (secondary)
    cw = J("cptac_wgs_purity.json").get("all_cases", {})
    out["cptac_wgs"] = {ref: {m: {"rho": v["rho"], "n": v["n"], "implementation": v.get("implementation")}
                              for m, v in cw.get(ref, {}).get("methods", {}).items()} for ref in ("frozen", "h5ad")}
    # implementation consistency: does a method run the same code in the anatomy arm and the TCGA donor-level arm?
    out["implementation_mismatches"] = [
        {"method": m, "acs_arm": a["implementation"], "tcga_gbm_donor_level": out["arms"]["gbm|h5ad"].get(m, {}).get("implementation"),
         "cptac_donor_level": out["cptac_wgs"]["h5ad"].get(m, {}).get("implementation")}
        for m, a in out["acs"].items() if not a["is_control"]
        and len({a["implementation"], out["arms"]["gbm|h5ad"].get(m, {}).get("implementation"),
                 out["cptac_wgs"]["h5ad"].get(m, {}).get("implementation")} - {None}) > 1]
    out["instrument_agreement"] = instrument_agreement(out)
    out["repaired"] = repaired_layer(out, truth_cache)
    out["repaired_agreement"] = repaired_agreement(out)
    out["n_join_checks"] = len(out["join_checks"])
    out["n_join_disagreements"] = sum(not j["agree"] for j in out["join_checks"])
    return out


def instrument_agreement(o: dict) -> dict:
    """Post hoc, descriptive: do two independent DNA instruments rank the methods alike? Spearman, across methods,
    between accuracy against ABSOLUTE and against GIMiCC's purity (tumour), and against Thorsson's LF and GIMiCC's
    immune total (leukocytes). scdc_ensemble is left out: it duplicates scdc."""
    from scipy import stats  # noqa: PLC0415
    out = {}
    for arm, cells in o["arms"].items():
        rec = {}
        for label, a_key, b_key in (("tumour", "Tumor|ABSOLUTE", "Tumor|GIMiCC"),
                                    ("leukocytes", "Leukocytes|LF", "Leukocytes|GIMiCC")):
            ms = [m for m, c in cells.items() if c["status"] == "ok" and m != "scdc_ensemble"
                  and c[a_key]["rho"] is not None and c[b_key]["rho"] is not None]
            if len(ms) >= 6:
                r = stats.spearmanr([cells[m][a_key]["rho"] for m in ms], [cells[m][b_key]["rho"] for m in ms])
                rec[label] = {"n_methods": len(ms), "spearman": round(float(r.statistic), 4),
                              "p": round(float(r.pvalue), 4)}
        out[arm] = rec
    return out


def _score(est: pd.DataFrame, tr: dict) -> dict:
    rec = {}
    for comp, truth in COMPS:
        x = est[idd.COMPARTMENTS[comp]]
        if x.isna().all().all():
            rec[f"{comp}|{truth}"] = {"rho": None, "n": 0, "why": "not modelled"}
            continue
        r, n = idd.truth_rho(x.sum(axis=1), comp, tr)
        rec[f"{comp}|{truth}"] = {"rho": None if not np.isfinite(r) else round(r, 4), "n": n}
    return rec


def repaired_layer(o: dict, truth_cache: dict) -> dict:
    """Registered value beside the repaired (results/repaired/) or genuine-package value, per method and arm."""
    rep: dict = {"tcga": {}, "genuine_bayesprism_tcga": {}, "acs": {}, "cptac": {}}
    for ref in ("frozen", "h5ad"):
        f = R / "repaired" / f"cptac_{ref}.json"
        if f.exists():
            rep["cptac"][ref] = json.loads(f.read_text())["methods"]
    for cohort, ref, _ in ARMS:
        f = R / "repaired" / f"tcga_{cohort}_{ref}.json"
        if f.exists():
            rep["tcga"][f"{cohort}|{ref}"] = json.loads(f.read_text())["methods"]
    # BayesPrism: the genuine package, run unbudgeted (results/extension, 2026-10-01/02), donor-level reference
    for cohort in ("gbm", "lgg"):
        f = R / "extension" / f"bayesprism_tcga_{cohort}_pipeline.csv"
        if f.exists():
            est = pd.read_csv(f, index_col=0)
            rep["genuine_bayesprism_tcga"][f"{cohort}|h5ad"] = {
                "implementation": "R:BayesPrism (pipeline configuration, unbudgeted)",
                "source": str(f.relative_to(config.PROJECT_ROOT)), **_score(est, truth_cache[cohort])}
    for name, f in (("bayesprism", R / "bayesprism_remeasured.json"),
                    ("dwls", R / "repaired" / "ivygap_dwls_genuine_rescaled.json")):
        if f.exists():
            d = json.loads(f.read_text())
            rep["acs"][name] = {"acs": d.get("acs"), "ci": d.get("ci"), "n_pairs": d.get("n_pairs"),
                                "n_samples_unestimated": d.get("n_samples_unestimated"),
                                "source": str(f.relative_to(config.PROJECT_ROOT))}
    return rep


def repaired_agreement(o: dict) -> dict:
    """POST HOC SENSITIVITY of the registered R2 test (ACS against tumour accuracy on TCGA-GBM; bar 0.60) to the
    implementation repairs: the same methods, with each broken implementation's value replaced by its repaired or
    genuine-package value where one exists. Every substitution is listed. The registered reading is unchanged."""
    from scipy import stats  # noqa: PLC0415
    rp = o.get("repaired", {})
    out = {}
    for arm, yard in (("gbm|frozen", "yardstick_agreement.json"), ("gbm|h5ad", "yardstick_agreement_h5ad.json")):
        reg = J(yard).get("per_method", {})
        if not reg:
            continue
        acs = {m: v["acs"] for m, v in reg.items()}
        acc = {m: v["rho_purity"] for m, v in reg.items()}
        subs = []
        for m, a in rp.get("acs", {}).items():                      # genuine-package anatomy scores
            if m in acs and a.get("acs") is not None and (a.get("n_pairs") or 0) == 57:
                subs.append(f"{m} ACS {acs[m]:.3f} -> {a['acs']:.3f}")
                acs[m] = a["acs"]
        for m, r in (rp.get("tcga", {}).get(arm) or {}).items():     # repaired TCGA accuracy
            t = (r.get("Tumor|ABSOLUTE") or {}).get("rho")
            if t is None:
                continue
            if m in acc:
                subs.append(f"{m} accuracy {acc[m]:.3f} -> {t:.3f} ({r.get('implementation')})")
                acc[m] = t
            elif m in o["acs"] and o["acs"][m]["comparable"] and not o["acs"][m]["is_control"]:
                subs.append(f"{m} added: accuracy {t:.3f} ({r.get('implementation')}), ACS {o['acs'][m]['acs']:.3f}")
                acs[m], acc[m] = o["acs"][m]["acs"], t
        g = (rp.get("genuine_bayesprism_tcga") or {}).get(arm)
        if g and "bayesprism" in acc:
            subs.append(f"bayesprism accuracy {acc['bayesprism']:.3f} -> {g['Tumor|ABSOLUTE']['rho']:.3f} (genuine)")
            acc["bayesprism"] = g["Tumor|ABSOLUTE"]["rho"]
        ms = sorted(set(acs) & set(acc))
        r = stats.spearmanr([acs[m] for m in ms], [acc[m] for m in ms])
        reg_r = J(yard)["arm_2_absolute_purity"]["all_methods"]["spearman"]
        out[arm] = {"n_methods": len(ms), "registered_spearman": reg_r, "repaired_spearman": round(float(r.statistic), 4),
                    "p": round(float(r.pvalue), 4), "meets_registered_bar": bool(r.statistic >= 0.60),
                    "substitutions": subs}
    return out


def fmt(v) -> str:
    return "--" if v is None else f"{v:+.2f}"


def write_doc(o: dict) -> None:
    L = ["# Evaluation matrix: every method against every truth", "",
         f"Generated by `scripts/evaluation_matrix.py` from saved artefacts (`{OUT.relative_to(config.PROJECT_ROOT)}`). "
         "Spearman correlation between each method's estimate and each truth, per cohort and reference build. "
         "No estimate is recomputed; the sample joins are the published ones.", "",
         f"**Join check:** the recomputed tumour-content correlation equals the registered yardstick's in "
         f"{o['n_join_checks'] - o['n_join_disagreements']} of {o['n_join_checks']} method-arms (tolerance {TOL}).", "",
         "Key: `py` this project's implementation; `py-reimpl` this project's Python version of a published "
         "package (never the package itself; the SCDC and BayesPrism versions are stand-ins for different algorithms, "
         "and the legacy DWLS inverts the published dampening: OPEN_DEFECTS D32, D34); `R:x` the published package; "
         "DEGEN the method reduced to another "
         "(reported, never under its own name as a result). ACS exists for glioblastoma only (Ivy GAP), on the "
         "donor-level reference. T>B: yes if the method's mean T exceeds its mean B (direct measurement says it "
         "should, in both tumour types); the count is samples with no lymphoid signal at all.", ""]
    names = {"gbm": "Glioblastoma", "lgg": "Lower-grade glioma"}
    for cohort in ("gbm", "lgg"):
        for ref in ("frozen", "h5ad"):
            cells = o["arms"][f"{cohort}|{ref}"]
            refname = "frozen signature" if ref == "frozen" else "donor-level reference (raw counts)"
            L += [f"## {names[cohort]}, {refname}", ""]
            hdr = "| method | implementation |"
            if cohort == "gbm" and ref == "h5ad":
                hdr += " ACS (Ivy GAP) |"
            hdr += " tumour vs ABSOLUTE | leukocytes vs LF | lymphoid vs EpiDISH | T | B | NK | T>B |"
            if cohort == "gbm":
                hdr += " tumour vs CPTAC WGS |"
            L += [hdr, "|" + "---|" * (hdr.count("|") - 1)]
            for m in sorted(cells, key=lambda k: -(cells[k].get("Tumor|ABSOLUTE", {}).get("rho") or -9)):
                c = cells[m]
                impl = {"python": "py", "python-reimplementation": "py-reimpl"}.get(c["implementation"], c["implementation"] or "--")
                impl += " DEGEN" if c["degenerate"] else ""
                row = f"| {m} | {impl} |"
                if cohort == "gbm" and ref == "h5ad":
                    a = o["acs"].get(m)
                    row += (f" {a['acs']:.3f}{'' if a['comparable'] else ' (not comparable)'} |" if a else " -- |")
                if c["status"] != "ok":
                    row += f" **{c['status'][:70]}** |" + " |" * (6 + (1 if cohort == "gbm" else 0))
                else:
                    row += "".join(f" {fmt(c[f'{k}|{t}']['rho'])} |" for k, t in COMPS)
                    tb = c.get("T>B")
                    row += (f" {'yes' if tb['T_exceeds_B'] else 'no'} ({tb['samples_with_no_lymphoid_signal']} empty) |"
                            if tb else " -- |")
                    if cohort == "gbm":
                        w = o["cptac_wgs"][ref].get(m)
                        row += f" {fmt(w['rho']) if w else '--'} |"
                L.append(row)
            L.append("")
    L += ["## Second methylation instrument: GIMiCC (glioma-specific)", "",
          f"GIMiCC's primary runs ({', '.join(f'{k}: {v}' for k, v in o['gimicc_runs'].items())}), keyed as every "
          "methylation comparison is. Its tumour fraction is a DNA purity (InfiniumPurify, Layer 0); in LGG it places "
          "B above T, against flow cytometry (`docs/WHY_B_OVER_T.md` §7n).", ""]
    for cohort in ("gbm", "lgg"):
        for ref in ("frozen", "h5ad"):
            cells = o["arms"][f"{cohort}|{ref}"]
            L += [f"**{names[cohort]}, {'frozen signature' if ref == 'frozen' else 'donor-level reference'}**", "",
                  "| method | tumour | leukocytes | lymphoid | T | B | NK |", "|---|---|---|---|---|---|---|"]
            for m in sorted(cells, key=lambda k: -(cells[k].get("Tumor|GIMiCC", {}).get("rho") or -9)):
                c = cells[m]
                if c["status"] != "ok":
                    continue
                L.append(f"| {m} |" + "".join(f" {fmt(c[f'{k}|GIMiCC']['rho'])} |" for k in GIM_COMPS))
            L.append("")
    L += ["## Degenerate per-sample solves (OPEN_DEFECTS D7)", "",
          "Samples for which a method returned everything in one cell type, an exactly uniform split, or a row "
          "identical to another sample's. Listed where any occurs.", ""]
    flagged = [(arm, m, c["per_sample_solves"]) for arm, cells in o["arms"].items() for m, c in cells.items()
               if c.get("per_sample_solves") and (c["per_sample_solves"]["all_in_one_type"] or
                                                  c["per_sample_solves"]["exactly_uniform"] or
                                                  c["per_sample_solves"]["duplicate_rows"])]
    L += [f"- `{m}`, {arm}: {d['all_in_one_type']} all-in-one-type, {d['exactly_uniform']} uniform, "
          f"{d['duplicate_rows']} duplicate rows, of {d['n']}" for arm, m, d in flagged] or ["- none"]
    L.append("")
    ia = o.get("instrument_agreement", {})
    L += ["## Do two independent DNA instruments rank the methods alike? (post hoc, descriptive)", "",
          "Spearman across methods between each method's accuracy against one instrument and against the other. "
          "High values mean the methods' accuracy ordering does not depend on which DNA truth is used.", "",
          "| arm | tumour: ABSOLUTE vs GIMiCC | leukocytes: LF vs GIMiCC |", "|---|---|---|"]
    for arm, rec in ia.items():
        L.append(f"| {arm} | " + " | ".join(
            (f"{rec[k]['spearman']:+.2f} (p {rec[k]['p']:.3f}, {rec[k]['n_methods']} methods)" if k in rec else "--")
            for k in ("tumour", "leukocytes")) + " |")
    L.append("")
    ra = o.get("repaired_agreement", {})
    if ra:
        L += ["## The registered ranking test (R2) after the repairs (post hoc sensitivity)", "",
              "ACS against tumour accuracy on TCGA-GBM, the registered test, recomputed with each broken "
              "implementation's value replaced by its repaired or genuine-package value. The registered reading "
              "stands; this asks whether it depended on the broken implementations.", ""]
        for arm, v in ra.items():
            L.append(f"- **{arm}**: registered rho {v['registered_spearman']:+.3f} -> repaired **{v['repaired_spearman']:+.3f}** "
                     f"(p {v['p']:.3f}, {v['n_methods']} methods; bar 0.60 {'met' if v['meets_registered_bar'] else 'not met'}). "
                     f"Substitutions: {'; '.join(v['substitutions']) or 'none yet'}.")
        L.append("")
    rp = o.get("repaired", {})
    L += ["## Repaired and genuine-package results, beside the registered ones", "",
          "Repairs are documented in `docs/METHOD_REPAIRS.md`. Registered values are unchanged; each repaired value "
          "comes from `results/repaired/` or from a genuine-package run that already exists.", "",
          "| arm | method | implementation (registered -> repaired) | tumour vs ABSOLUTE | leukocytes vs LF | "
          "lymphoid vs EpiDISH | T>B |", "|---|---|---|---|---|---|---|"]
    for arm, methods in rp.get("tcga", {}).items():
        for m, r in methods.items():
            reg = o["arms"].get(arm, {}).get(m, {})
            if "failed" in r:
                L.append(f"| {arm} | {m} | {reg.get('implementation') or '--'} -> **failed: {r['failed'][:60]}** | | | | |")
                continue
            def pair(k):
                a = (reg.get(k) or {}).get("rho") if reg.get("status") == "ok" else None
                b = (r.get(k) or {}).get("rho")
                return f"{fmt(a)} -> **{fmt(b)}**"
            tb = r.get("T>B", {})
            L.append(f"| {arm} | {m} | {reg.get('implementation') or reg.get('status', '--')[:40]} -> {r['implementation']} | "
                     f"{pair('Tumor|ABSOLUTE')} | {pair('Leukocytes|LF')} | {pair('Lymphoid|EpiDISH')} | "
                     f"{'yes' if tb.get('T_exceeds_B') else 'no' if 'T_exceeds_B' in tb else '--'} |")
    for arm, g in rp.get("genuine_bayesprism_tcga", {}).items():
        reg = o["arms"].get(arm, {}).get("bayesprism", {})
        def pair(k):
            return f"{fmt((reg.get(k) or {}).get('rho'))} -> **{fmt(g[k]['rho'])}**"
        L.append(f"| {arm} | bayesprism | {reg.get('implementation')} -> {g['implementation']} | {pair('Tumor|ABSOLUTE')} | "
                 f"{pair('Leukocytes|LF')} | {pair('Lymphoid|EpiDISH')} | -- |")
    for ref, methods in rp.get("cptac", {}).items():
        for m, r in methods.items():
            w = r.get("Tumor|WGS") or {}
            if w.get("rho") is None:
                continue
            reg = o["cptac_wgs"].get(ref, {}).get(m, {})
            L.append(f"| cptac|{ref} | {m} | {reg.get('implementation') or '--'} -> {r['implementation']} | "
                     f"vs whole-genome: {fmt(reg.get('rho'))} -> **{fmt(w['rho'])}** | -- | -- | -- |")
    if rp.get("acs"):
        L += ["", "**Anatomy (ACS, Ivy GAP), genuine packages:**", ""]
        for m, a in rp["acs"].items():
            reg = o["acs"].get(m, {})
            L.append(f"- `{m}`: registered {reg.get('acs')} ({reg.get('implementation')}) -> genuine {a['acs']} on "
                     f"{a['n_pairs']} of 57 pairs, {a.get('n_samples_unestimated') or 0} samples unestimated "
                     f"(`{a['source']}`)")
    L.append("")
    L += ["## Every empty cell, and why", ""]
    L += [f"- `{g['method']}`, {g['arm']}, {g['what']}: {g['why']}" for g in o["gaps"]] or ["- none"]
    L += ["", "## Implementation differences between arms (the same name, different code)", ""]
    L += [f"- `{x['method']}`: anatomy arm {x['acs_arm']}; TCGA-GBM donor-level {x['tcga_gbm_donor_level']}; "
          f"CPTAC donor-level {x['cptac_donor_level']}" for x in o["implementation_mismatches"]] or ["- none"]
    bad = [j for j in o["join_checks"] if not j["agree"]]
    if bad:
        L += ["", "## Join disagreements (recomputed vs registered tumour rho)", ""]
        L += [f"- `{j['method']}` {j['arm']}: registered {j['registered']} (n {j['n_registered']}), recomputed "
              f"{j['recomputed']} (n {j['n_recomputed']})" for j in bad]
    DOC.write_text("\n".join(L) + "\n")


def write_supplementary_csv(o: dict) -> Path:
    """Supplementary table: one row per method x arm, every truth's rho and n, registered beside repaired."""
    rows = []
    rp = o.get("repaired", {})
    for arm, cells in o["arms"].items():
        for m, c in cells.items():
            r = {"arm": arm, "method": m, "implementation": c.get("implementation"), "degenerate": c.get("degenerate"),
                 "status": c.get("status")}
            if arm == "gbm|h5ad" and m in o["acs"]:
                r["acs_ivygap"] = o["acs"][m]["acs"]
            for k, v in c.items():
                if "|" in k and isinstance(v, dict) and "rho" in v:
                    r[f"rho {k}"], r[f"n {k}"] = v["rho"], v.get("n")
            tb = c.get("T>B")
            if tb:
                r["T_exceeds_B"], r["samples_no_lymphoid"] = tb.get("T_exceeds_B"), tb.get("samples_with_no_lymphoid_signal")
            if arm.startswith("gbm"):
                w = o["cptac_wgs"].get(arm.split("|")[1], {}).get(m)
                if w:
                    r["rho Tumor|CPTAC_WGS"] = w["rho"]
            rep = (rp.get("tcga", {}).get(arm) or {}).get(m)
            if rep and "failed" not in rep:
                r["repaired_implementation"] = rep.get("implementation")
                for k in ("Tumor|ABSOLUTE", "Leukocytes|LF", "Lymphoid|EpiDISH"):
                    r[f"repaired rho {k}"] = (rep.get(k) or {}).get("rho")
            g = (rp.get("genuine_bayesprism_tcga") or {}).get(arm) if m == "bayesprism" else None
            if g:
                r["repaired_implementation"] = g["implementation"]
                for k in ("Tumor|ABSOLUTE", "Leukocytes|LF", "Lymphoid|EpiDISH"):
                    r[f"repaired rho {k}"] = g[k]["rho"]
            rows.append(r)
    out = config.PROJECT_ROOT / "docs" / "supplementary" / "evaluation_matrix.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return out


def main() -> int:
    o = build()
    OUT.write_text(json.dumps(o, indent=2, default=float))
    write_doc(o)
    csv = write_supplementary_csv(o)
    print(f"wrote {csv.relative_to(config.PROJECT_ROOT)}")
    print(f"join checks: {o['n_join_checks'] - o['n_join_disagreements']} of {o['n_join_checks']} agree")
    print(f"gaps: {len(o['gaps'])}; implementation mismatches: {len(o['implementation_mismatches'])}")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)} and {DOC.relative_to(config.PROJECT_ROOT)}")
    return 0 if o["n_join_disagreements"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
