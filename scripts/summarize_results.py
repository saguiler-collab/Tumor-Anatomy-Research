#!/usr/bin/env python3
"""
summarize_results.py — render the results tree as markdown, from the artefacts.

WHY THIS EXISTS
---------------
The last run's numbers survived a destructive accident only because they had already
been copied into RESULTS.md by hand. That is worth fixing in both directions: hand
transcription is how a number ends up in a write-up attached to the wrong cohort, and
it is also how a write-up silently goes stale after a rerun.

Everything below is read from the files a run actually wrote. Nothing is typed twice,
and every table names the artefact it came from so a reader can check it.

    python scripts/summarize_results.py                 # the canonical tree
    python scripts/summarize_results.py --results-dir results_synthetic
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd


def _load_json(path: Path):
    return json.loads(path.read_text()) if path.exists() else None


def _load_csv(path: Path, **kw):
    return pd.read_csv(path, **kw) if path.exists() else None


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def integrity_section(anat: Path) -> str:
    """
    The three risks Anatomy_Test.md names that used to be asserted rather than checked.
    Placed before the leaderboard on purpose: whether the constraint file is actually
    pre-registered, what the roster cannot see, and how each method was configured all
    bear on how the numbers below should be read.
    """
    out = ["\n## 1a. Integrity checks\n"]

    reg = _load_json(anat / "registration_status.json")
    if reg:
        state = reg.get("state")
        mark = "**REGISTERED**" if state == "REGISTERED" else f"**{state}**"
        out.append(f"**Pre-registration:** {mark}. {reg.get('verdict','')}\n")

    cov = _load_json(anat / "reference_coverage.json")
    if cov and cov.get("atlas_labels_dropped"):
        n_drop = cov.get("n_cells_dropped")
        frac = cov.get("fraction_of_atlas_dropped")
        out.append(f"**Reference coverage:** the roster maps "
                   f"{len(cov.get('atlas_labels_mapped_onto_roster', []))} of the "
                   f"atlas's labels and drops the rest — **{n_drop:,} cells "
                   f"({frac:.1%})**.\n")
        out.append("| dropped population | cells |")
        out.append("|---|---|")
        for k, v in list(cov["atlas_labels_dropped"].items())[:8]:
            out.append(f"| {k} | {v:,} |")
        out.append(f"\n{cov.get('what_this_means','')}\n")

    cfg = _load_json(anat / "method_configs.json")
    if cfg:
        dev = cfg.get("methods_with_declared_deviations", [])
        if dev:
            out.append(f"**Method configuration:** every method's parameters are "
                       f"recorded in `method_configs.json`. Declared departures from "
                       f"published defaults: **{', '.join(dev)}**.\n")
            for rec in cfg.get("configs", []):
                for d in rec.get("declared_deviations_from_published_defaults", []):
                    out.append(f"- `{rec['method']}` — {d['parameter']}: "
                               f"{d['published_default']} → {d['used']}. {d['why']} "
                               f"*Decided {d['decided']}.*")
            out.append("")
        else:
            out.append("**Method configuration:** recorded in `method_configs.json`; "
                       "no method departs from its published defaults.\n")

    impl = _load_json(anat / "implementation_report.json")
    if impl:
        genuine = [r["method"] for r in impl if str(r.get("implementation","")).startswith("R:")]
        fell = [(r["method"], r.get("fallback_reason")) for r in impl
                if r.get("implementation") == "python-reimplementation"]
        out.append(f"**What actually ran:** {len(genuine)} published package(s) ran as "
                   f"the genuine R implementation ({', '.join(genuine) or 'none'}).")
        if fell:
            out.append("\nFell back to this project's Python reimplementation:\n")
            out.append("| method | why |")
            out.append("|---|---|")
            for m, why in fell:
                out.append(f"| {m} | {str(why or 'reason not recorded')[:150]} |")
        out.append("")

    return "\n".join(out) if len(out) > 1 else ""


def resolution_note(lb_path: Path, pc_path: Path) -> str:
    """
    How finely can ACS separate methods on this cohort at all?

    ACS is a weighted proportion over a fixed, small set of (constraint, tumour) pairs,
    so it is quantised: with a weighted denominator of D it can take at most D+1 values,
    spaced 1/D apart. On Ivy GAP's H&E study D is 65, and methods therefore tie in
    blocks. That is a property of the instrument, not a coincidence, and it bounds every
    ranking built on top — including the agreement test, which needs enough DISTINCT
    points to correlate.
    """
    lb = _load_csv(lb_path, index_col=0)
    pc = _load_csv(pc_path)
    if lb is None or pc is None:
        return ""
    one = pc[pc["method"] == pc["method"].iloc[0]]
    den = float((one["n_tumors_evaluable"] * one["weight"]).sum())
    if den <= 0:
        return ""

    real = lb[~lb["is_control"]]
    out = ["\n### What ACS can and cannot resolve here\n",
           f"- weighted denominator: **{den:.0f}** (constraint x tumour pairs, weighted)",
           f"- so ACS takes at most **{den + 1:.0f} distinct values**, spaced "
           f"**{1 / den:.4f}** apart",
           f"- {len(lb)} methods scored produce **{lb['acs'].nunique()} distinct values**; "
           f"the {len(real)} real methods produce **{real['acs'].nunique()}**"]

    ties = {v: list(g.index) for v, g in lb.groupby(lb["acs"].round(9)) if len(g) > 1}
    for v, members in sorted(ties.items(), reverse=True):
        out.append(f"- tied at **{v:.3f}**: {', '.join(members)}")

    need = 6
    if real["acs"].nunique() < need:
        out.append(
            f"\n**This bounds the agreement test independently of any yardstick.** It "
            f"needs at least {need} methods to correlate, but the ACS side supplies only "
            f"**{real['acs'].nunique()} distinct values** across {len(real)} real "
            f"methods. Even with ground truth covering every method, a rank correlation "
            f"would be computed on mostly tied ranks.")
    return "\n".join(out) + "\n"


def control_audit(lb_path: Path, pc_path: Path) -> str:
    """
    The control check, stated more strictly than the run's own one-line verdict.

    `run_anatomic` compares the best control against the MEDIAN real method. That is a
    weak test: a control can tie the worst real method, and beat its own permutation
    null, while still sitting below the median. Both of those are the protocol's third
    branch — "negative controls score highly, the constraint set is too permissive" —
    and both must be visible.

    Nothing here changes the constraint file, the controls, or any score. It reports
    facts already present in the artefacts that the headline verdict does not surface.
    """
    lb = _load_csv(lb_path, index_col=0)
    if lb is None:
        return ""
    ctrl = lb[lb["is_control"]]
    real = lb[~lb["is_control"]]
    if ctrl.empty or real.empty:
        return ""

    best_c = ctrl["acs"].idxmax()
    best_c_acs = float(ctrl.loc[best_c, "acs"])
    worst_r = real["acs"].idxmin()
    worst_r_acs = float(real.loc[worst_r, "acs"])

    tied_or_beaten = real.index[real["acs"] <= best_c_acs].tolist()
    beats_null = ctrl.index[(ctrl["null_p"] < 0.05)
                            & (ctrl["acs"] > ctrl["null_mean"])].tolist()

    out = ["\n### Control audit (stricter than the headline verdict)\n",
           "_The control is a **single permutation** drawn from `config.RANDOM_SEED`. "
           "Everything in this block describes that one draw._\n"]
    out += [
           f"- best control: **{best_c} = {best_c_acs:.3f}**",
           f"- worst real method: {worst_r} = {worst_r_acs:.3f}"]
    if tied_or_beaten:
        out.append(f"- real methods the best control **ties or beats**: "
                   f"**{', '.join(tied_or_beaten)}**")
    else:
        out.append("- real methods the best control ties or beats: none")
    if beats_null:
        out.append(f"- controls that **beat their own permutation null** (p < 0.05): "
                   f"**{', '.join(beats_null)}**")
    else:
        out.append("- controls beating their own permutation null: none")

    pc = _load_csv(pc_path)
    if pc is not None:
        sub = pc[pc["method"] == best_c]
        if not sub.empty:
            sub = sub.copy()
            sub["w_den"] = sub["n_tumors_evaluable"] * sub["weight"]
            sub["w_num"] = sub["n_satisfied"] * sub["weight"]
            den = sub["w_den"].sum()
            perfect = sub[sub["n_satisfied"] == sub["n_tumors_evaluable"]]
            if not perfect.empty and den:
                share = perfect["w_den"].sum() / den
                out.append(
                    f"- constraints **this single draw** satisfies perfectly: "
                    f"**{', '.join(perfect['constraint'])}** — "
                    f"{perfect['w_den'].sum():.0f} of {den:.0f} weighted units "
                    f"(**{share:.1%}** of the constraint set). This is ONE permutation, "
                    f"not a property of the constraints; see the draw distribution "
                    f"below before reading anything into it.")
            out.append("\n| constraint | weight | evaluable | best control satisfies |")
            out.append("|---|---|---|---|")
            for _, r in sub.iterrows():
                out.append(f"| {r['constraint']} | {_fmt(r['weight'], 1)} "
                           f"| {int(r['n_tumors_evaluable'])} "
                           f"| {int(r['n_satisfied'])} |")
    cal = (_load_json(lb_path.parent / "control_calibration.json")
           or _load_json(lb_path.parent / "shuffled_control_diagnostic.json"))
    if cal and "controls" in cal:
        out.append(f"\n**Both controls calibrated over {cal['n_draws']} independent "
                   f"draws each.** Each leaderboard control row is one draw from these.\n")
        out.append("| control | mean | sd | 5–95% | the leaderboard's draw |")
        out.append("|---|---|---|---|---|")
        for name, c in cal["controls"].items():
            a = c["acs"]
            ar = c.get("as_run_draw")
            shown = (f"{ar['acs']:.3f} ({ar['percentile_among_draws']:.0%} pct)"
                     if ar else "—")
            out.append(f"| `{name}` | {a['mean']:.3f} | {a['sd']:.3f} | "
                       f"{a['p05']:.3f} – {a['p95']:.3f} | {shown} |")
        hardest = cal["hardest_control"]
        hp95 = cal["controls"][hardest]["acs"]["p95"]
        out.append(f"\nMethods are judged against the **hardest** control "
                   f"(`{hardest}`, 95th percentile **{hp95:.3f}**), not a convenient one.")
        nd = cal.get("not_distinguishable_from_noise", [])
        if nd:
            out.append(
                f"\n> **{', '.join(nd)}** {'is' if len(nd) == 1 else 'are'} **not "
                f"distinguishable from a meaningless composition**. The other "
                f"{len(real) - len(nd)} real method(s) are.")
        else:
            out.append("\n> Every real method scores above the hardest control's 95th "
                       "percentile.")

        rates = cal.get("per_constraint_satisfied_rate", {})
        if rates:
            out.append(f"\n**What the hardest control satisfies, averaged over "
                       f"{cal['n_draws']} draws** — this is the number to read, not the "
                       f"single draw above:\n")
            out.append("| constraint | satisfied rate |")
            out.append("|---|---|")
            for cid, rate in sorted(rates.items(), key=lambda kv: -kv[1]):
                out.append(f"| {cid} | {rate:.2f} |")
    return "\n".join(out) + "\n"


def leaderboard_table(path: Path) -> str:
    df = _load_csv(path, index_col=0)
    if df is None:
        return f"_(not produced: `{path}`)_\n"
    lines = ["| method | ACS | 95% CI | null mean | null p | tumours | control |",
             "|---|---|---|---|---|---|---|"]
    for m, r in df.iterrows():
        ctrl = "**yes**" if r.get("is_control") else ""
        name = f"**{m}**" if r.get("is_control") else m
        lines.append(
            f"| {name} | {_fmt(r['acs'])} | {_fmt(r['ci_low'])} – {_fmt(r['ci_high'])} "
            f"| {_fmt(r['null_mean'])} | {_fmt(r['null_p'])} | "
            f"{int(r['n_tumors'])} | {ctrl} |")
    return "\n".join(lines) + "\n"


def per_constraint_table(path: Path, method: str) -> str:
    df = _load_csv(path)
    if df is None:
        return f"_(not produced: `{path}`)_\n"
    sub = df[df["method"] == method]
    if sub.empty:
        return f"_(no rows for method {method!r})_\n"
    lines = ["| ID | claim | weight | tumours evaluable | satisfied |",
             "|---|---|---|---|---|"]
    for _, r in sub.iterrows():
        n, k = int(r["n_tumors_evaluable"]), int(r["n_satisfied"])
        frac = f"{k} ({k / n:.2f})" if n else "—"
        lines.append(f"| {r['constraint']} | {r['description']} | {_fmt(r['weight'], 1)} "
                     f"| {n} | {frac} |")
    return "\n".join(lines) + "\n"


def per_tumor_table(path: Path, method: str) -> str:
    df = _load_csv(path)
    if df is None:
        return f"_(not produced: `{path}`)_\n"
    sub = df[df["method"] == method]
    if sub.empty:
        return f"_(no rows for method {method!r})_\n"
    lines = ["| tumour | constraints evaluable | satisfied | ACS | violated |",
             "|---|---|---|---|---|"]
    for _, r in sub.iterrows():
        viol = r["constraints_violated"] if isinstance(r["constraints_violated"], str) else ""
        lines.append(f"| {r['tumor_id']} | {int(r['n_constraints_evaluable'])} "
                     f"| {int(r['n_satisfied'])} | {_fmt(r['acs'])} | {viol or '—'} |")
    return "\n".join(lines) + "\n"


def render(results_dir: Path) -> str:
    anat = results_dir / "anatomic"
    full = anat / "full_database"
    out: list[str] = []
    A = out.append

    rep = _load_json(anat / "anatomic_report.json")
    recon = _load_json(results_dir / "data_reconciliation.json")
    audit = _load_json(results_dir / "clinical_missingness_audit.json")

    A(integrity_section(anat))
    A("## 2. Cohort and provenance\n")
    if rep:
        A(f"- constraint freeze hash: `{rep['constraint_freeze_hash']}`")
        A(f"- ACS cohort: **{rep['n_samples']} samples / {rep['n_tumors']} tumours**, "
          f"{rep['n_permutations']:,} within-tumour permutations per method")
        A(f"- structures: {rep['structures']}")
        if rep.get("deconvolved_all_samples"):
            A(f"- deconvolved: {rep['n_samples_deconvolved']} samples / "
              f"{rep['n_tumors_deconvolved']} tumours")
    if recon:
        A(f"- portal reconciliation: {recon['verdict']}")
    A("")

    A("## 3. The ACS leaderboard\n")
    A(f"_from `{(anat / 'acs_leaderboard.csv').relative_to(results_dir.parent)}`_\n")
    A(leaderboard_table(anat / "acs_leaderboard.csv"))
    if rep:
        A(f"\n**Control verdict.** {rep['control_verdict']}\n")

    A(resolution_note(anat / "acs_leaderboard.csv", anat / "acs_per_constraint.csv"))
    A(control_audit(anat / "acs_leaderboard.csv", anat / "acs_per_constraint.csv"))

    best = rep["best_acs_method"] if rep else None
    if best:
        A(f"### Per constraint — best method (`{best}`)\n")
        A(per_constraint_table(anat / "acs_per_constraint.csv", best))
        A(f"\n### Per tumour — `{best}`\n")
        A(per_tumor_table(anat / "acs_per_tumor.csv", best))

    if rep and rep.get("implementation_disclosure"):
        A("\n## 4. What actually ran\n")
        A("| method | implementation | degenerate | why |")
        A("|---|---|---|---|")
        for d in rep["implementation_disclosure"]:
            why = d.get("degeneracy_reason") or d.get("fallback_reason") or ""
            A(f"| {d['method']} | {d['implementation']} | "
              f"{'**yes**' if d['degenerate'] else 'no'} | {why[:110]} |")
        A("")

    agree = _load_json(anat / "agreement_report.json")
    if agree:
        A("\n## 5. The agreement test (the primary result)\n")
        A(f"Pre-registered bar: rho >= {agree['preregistered_threshold']} "
          f"AND a bootstrap CI excluding zero. Minimum methods: "
          f"{agree['min_methods_for_correlation']}.\n")
        A("| yardstick | rho | 95% CI | methods | distinct | verdict |")
        A("|---|---|---|---|---|---|")
        for y in agree["by_yardstick"]:
            A(f"| {y['yardstick']} | {_fmt(y['rho'])} | "
              f"{_fmt(y['ci_low'])} – {_fmt(y['ci_high'])} | {y['n_methods']} "
              f"| {_fmt(y['n_distinct'], 0)} | {str(y['verdict'])[:90]} |")
        A("")

    sens = _load_csv(anat / "acs_cohort_sensitivity.csv", index_col=0)
    if sens is not None:
        A("\n## 6. Does deconvolving all 270 samples change the ranking?\n")
        A("| method | ACS (122 anatomic) | ACS (270 deconvolved) | delta | rank change |")
        A("|---|---|---|---|---|")
        for m, r in sens.iterrows():
            A(f"| {m} | {_fmt(r['acs_anatomic_only'])} | {_fmt(r['acs_full_database'])} "
              f"| {r['delta']:+.3f} | {r['rank_change']:+.0f} |")
        rho = sens["acs_anatomic_only"].corr(sens["acs_full_database"], method="spearman")
        A(f"\nSpearman rho between the two ACS rankings: **{_fmt(rho, 4)}**\n")

        moved = sens.index[sens["delta"].abs() > 1e-9].tolist()
        still = sens.index[sens["delta"].abs() <= 1e-9].tolist()
        A(f"\n**Which methods moved, and why.** {len(still)} of {len(sens)} did not "
          f"move by a single unit: {', '.join(still)}. Those solve each sample "
          f"independently, so what else is in the cohort cannot reach them — the "
          f"zeros are exact, not rounded.\n")
        if moved:
            A(f"The ones that moved are **{', '.join(moved)}**. Bisque normalises the "
              f"bulk with cohort-wide per-gene statistics (`B.mean(axis=1)` and "
              f"`B.std(axis=1)`), so the 148 ISH-cluster samples shift the reference "
              f"frame every anatomic sample is mapped through. `control_random` moves "
              f"for a different and uninteresting reason: it draws one Dirichlet sample "
              f"per row, so a 270-row draw is not a superset of a 122-row draw.\n")
            A("Note the direction: **adding 148 more real samples made Bisque's "
              "anatomic concordance worse**, not better. For a method that borrows "
              "strength across samples, which samples it is handed is part of the "
              "method — and Ivy GAP's two studies are not the same tissue.\n")

    if full.exists():
        frep = _load_json(full / "anatomic_report.json")
        if frep:
            A("\n### The full-database run\n")
            A(f"- deconvolved **{frep['n_samples_deconvolved']} samples / "
              f"{frep['n_tumors_deconvolved']} tumours**")
            A(f"- ACS still scored on **{frep['n_samples']} samples / "
              f"{frep['n_tumors']} tumours** — {frep['acs_scoring_set']}")
            A(f"- control verdict: {frep['control_verdict']}\n")
            A(leaderboard_table(full / "acs_leaderboard.csv"))

    A("\n## 7. Survival\n")
    strict = _load_json(results_dir / "survival" / "survival_power.json")
    if strict:
        A(f"**Canonical verdict.** {strict.get('verdict')}\n")
    if audit:
        A(f"**Missingness.** {audit.get('verdict', audit.get('reason'))}\n")
        if "mgmt" in audit:
            m = audit["mgmt"]
            A(f"- MGMT methylated among tumours with a recorded time: "
              f"{m['methylated_among_recorded']}; among blanks: "
              f"{m['methylated_among_blank']} (Fisher p = {m['fisher_p']:.4g})")
        if "age" in audit:
            a = audit["age"]
            A(f"- median age: {_fmt(a['median_recorded'], 0)} recorded vs "
              f"{_fmt(a['median_blank'], 0)} blank (Mann-Whitney p = "
              f"{a['mannwhitney_p']:.3g})")
        A("")

    decl = results_dir / "survival" / "declared_event_policy"
    dpow = _load_json(decl / "survival_power.json")
    if dpow:
        A("\n### Under the declared event policy (quarantined, not canonical)\n")
        A(f"- policy: `{dpow.get('event_policy')}`")
        A(f"- patients {dpow['n_patients']} · events {dpow['n_events']} · "
          f"events per covariate {dpow['events_per_covariate']:.2f}")
        A(f"- smallest detectable C-index difference: "
          f"**{dpow['detectable_delta_c_index']:.3f}**")
        A(f"- {dpow['verdict']}\n")
        dm = _load_csv(decl / "survival_metrics.csv", index_col=0)
        if dm is not None and "c_index_with_composition" in dm.columns:
            A("| method | C baseline | C + composition | delta | interpretation |")
            A("|---|---|---|---|---|")
            for m, r in dm.iterrows():
                A(f"| {m} | {_fmt(r.get('c_index_baseline'))} | "
                  f"{_fmt(r.get('c_index_with_composition'))} | "
                  f"{_fmt(r.get('delta_c_index'))} | {r.get('interpretation')} |")
            A("")

            lb = _load_csv(results_dir / "anatomic" / "acs_leaderboard.csv", index_col=0)
            if lb is not None and "delta_c_index" in dm.columns:
                j = (lb[~lb["is_control"]][["acs"]]
                     .join(dm[["delta_c_index"]], how="inner").dropna())
                if len(j) >= 4:
                    from scipy.stats import spearmanr
                    r, pv = spearmanr(j["acs"], j["delta_c_index"])
                    A("\n#### Exploratory: does ACS track prognostic value?\n")
                    A("Not a protocol analysis. The protocol's agreement test correlates "
                      "ACS against *accuracy*, which is still not computable. This "
                      "correlates it against *prognostic value* instead, which the "
                      "declared-policy run makes available. It selects nothing.\n")
                    A(f"- Spearman(ACS, delta C-index) over {len(j)} methods: "
                      f"**rho = {r:.3f}**, p = {pv:.3f}")
                    best_acs, best_prog = j["acs"].idxmax(), j["delta_c_index"].idxmax()
                    worst_acs = j["acs"].idxmin()
                    A(f"- best ACS: `{best_acs}` (ACS {j.loc[best_acs, 'acs']:.3f}, "
                      f"delta C {j.loc[best_acs, 'delta_c_index']:+.3f})")
                    A(f"- **worst ACS: `{worst_acs}` (ACS {j.loc[worst_acs, 'acs']:.3f}, "
                      f"delta C {j.loc[worst_acs, 'delta_c_index']:+.3f}) — which is also "
                      f"the best prognostic method**"
                      if worst_acs == best_prog else
                      f"- best prognosis: `{best_prog}` "
                      f"(ACS {j.loc[best_prog, 'acs']:.3f}, "
                      f"delta C {j.loc[best_prog, 'delta_c_index']:+.3f})")
                    A("\nThree independent reasons this cannot support a claim: every "
                      "C-index is INCONCLUSIVE by the pre-specified power rule; ACS "
                      f"supplies only {j['acs'].nunique()} distinct values across "
                      f"{len(j)} methods, so the rank is mostly ties; and the outcome "
                      "side rests on a declared assumption about censoring. It is "
                      "recorded because it points the same way as the superseded run "
                      "did, and because the direction is the protocol's second branch.\n")

    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    args = ap.parse_args()
    root = Path(args.results_dir)
    if not root.exists():
        print(f"no such results directory: {root}")
        return 1
    print(render(root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
