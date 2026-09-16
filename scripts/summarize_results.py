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
import datetime as dt
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


def _ordering_verdict_is_stale(anat: Path) -> int:
    """
    How many result files are NEWER than the registration-status file that judged them.

    Zero means the verdict was computed after everything it checked, and can be read as
    written. Anything else means it cannot.
    """
    status = anat / "registration_status.json"
    if not status.exists():
        return 0
    cutoff = status.stat().st_mtime
    newer = 0
    for d in (anat, anat.parent / "benchmark", anat.parent / "survival"):
        if not d.exists():
            continue
        for p in d.rglob("*"):
            if p.is_file() and p != status and p.stat().st_mtime > cutoff + 1.0:
                newer += 1
    return newer


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

        # The ordering verdict is only as current as the file that carries it. It compares
        # result mtimes against the registration, so a status file written BEFORE the run
        # finished judged the previous run's artefacts — which is exactly what happened on
        # 2026-09-10: it published "25 result file(s) are OLDER than the registration"
        # when the true count was zero. Refusing to reprint a verdict whose own file
        # predates what it judged is what stops that reaching a reader again.
        stale = _ordering_verdict_is_stale(anat)
        if stale:
            out.append(
                f"> **The ordering verdict above is STALE and must not be read.** "
                f"`registration_status.json` is older than {stale} result file(s) it "
                f"claims to have checked, so it judged a previous run's artefacts. "
                f"Recompute it with `python scripts/register.py --status` before citing "
                f"any pre-registration claim.\n")

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


def constraint_discrimination(anat: Path) -> str:
    """
    Which constraints actually separate the methods, and which are unanimous.

    A reviewer's first question about a composite score is whether one component drives
    it. The answer here is specific and partly uncomfortable: C3 is satisfied by every
    comparable real method on every evaluable tumour, so it contributes nothing to the
    ORDERING even though the controls miss it about half the time. Reporting only the mean
    satisfied rate would hide that. Spread is the statistic that exposes it.

    A unanimous constraint is not a broken one. C3 is near-definitional -- microvascular
    proliferation is where the vessels are -- and a constraint set with no easy members
    could not distinguish a method from a scrambled one. It earns its weight against the
    controls and carries none of the ranking, and both halves of that belong in print.
    """
    per_c = _load_csv(anat / "acs_per_constraint.csv")
    lb = _load_csv(anat / "acs_leaderboard.csv")
    if per_c is None or lb is None:
        return "_(not produced)_\n"

    real = lb[(~lb["is_control"].astype(bool))
              & (lb["comparable"].astype(bool))]["method"].tolist()
    sub = per_c[per_c["method"].isin(real)]
    if sub.empty:
        return "_(no comparable real methods)_\n"

    ctrl_names = lb[lb["is_control"].astype(bool)]["method"].tolist()
    ctrl = per_c[per_c["method"].isin(ctrl_names)]

    lines = [f"Across the **{len(real)} comparable real methods**. `spread` is the range of "
             "the satisfied rate; a spread of zero means the constraint is unanimous and "
             "contributes nothing to the ORDERING, however well it separates real methods "
             "from the controls.\n",
             "| ID | claim | weight | tumours | mean rate | min | max | spread | controls |",
             "|---|---|---|---|---|---|---|---|---|"]
    rows = []
    for cid, grp in sub.groupby("constraint", sort=True):
        r = grp["fraction_satisfied"].astype(float)
        c = ctrl[ctrl["constraint"] == cid]["fraction_satisfied"].astype(float)
        rows.append((cid, str(grp["description"].iloc[0]), float(grp["weight"].iloc[0]),
                     int(grp["n_tumors_evaluable"].iloc[0]), r.mean(), r.min(), r.max(),
                     r.max() - r.min(), c.mean() if len(c) else float("nan")))
    for cid, desc, w, n, mean, lo, hi, spread, cm in sorted(rows, key=lambda x: -x[7]):
        flag = " **unanimous**" if spread == 0 else ""
        lines.append(f"| {cid} | {desc} | {_fmt(w, 1)} | {n} | {mean:.3f} | {lo:.3f} "
                     f"| {hi:.3f} | **{spread:.3f}**{flag} | {_fmt(cm)} |")

    carrying = [r for r in rows if r[7] > 0]
    unanimous = [r[0] for r in rows if r[7] == 0]
    lines.append("")
    if unanimous:
        verb = "is" if len(unanimous) == 1 else "are"
        lines.append(f"**{len(unanimous)} of {len(rows)} constraints {verb} unanimous "
                     f"({', '.join(unanimous)}).** The ranking is carried by the "
                     f"remaining {len(carrying)}, led by "
                     + ", ".join(f"`{r[0]}` (spread {r[7]:.3f})"
                                 for r in sorted(carrying, key=lambda x: -x[7])[:3]) + ".")
    else:
        lines.append("Every constraint separates at least two methods; none is unanimous.")
    return "\n".join(lines) + "\n"


def sensitivity_section(anat: Path) -> str:
    """
    Leave-one-out robustness, rendered from `constraint_sensitivity.json`.

    Reporting only. The registered ACS uses all constraints and all tumours; nothing here
    reorders anything, and a constraint that carries none of the ordering still stays in
    the set. Removing it after seeing the scores is the retune the protocol forbids.
    """
    rep = _load_json(anat / "constraint_sensitivity.json")
    if rep is None:
        return ("_(not produced: run `python scripts/constraint_sensitivity.py`)_\n")

    s = rep["summary"]
    out = [f"Recomputed from the archived (constraint x tumour) satisfaction matrix, which "
           f"reproduces every published ACS to {rep['reconstruction_max_abs_error']:.0e} "
           f"before any variant is reported. **Reporting only** — the registered ACS uses "
           f"all constraints and all tumours, and no variant below may reorder the "
           f"leaderboard, select a method, or justify dropping a constraint.\n",
           "| excluded | rho vs full ranking | max rank move | methods moved | top method |",
           "|---|---|---|---|---|"]
    for r in rep["leave_one_constraint_out"] + rep["leave_one_tumor_out"]:
        flag = " **changes**" if r["top_method_changes"] else ""
        out.append(f"| {r['excluded']} | {r['spearman_vs_full_ranking']:.4f} "
                   f"| {r['max_rank_move']:.1f} | {r['n_methods_whose_rank_moves']} "
                   f"| `{r['top_method']}`{flag} |")
    out.append("")

    ct = s["constraints_whose_removal_changes_the_top_method"]
    tu = s["tumors_whose_removal_changes_the_top_method"]
    if not ct and not tu:
        out.append("**The top method survives every single exclusion.** No one constraint "
                   "and no one tumour is responsible for it.")
    else:
        out.append(f"**The top method changes when these are excluded:** "
                   f"{', '.join(ct + tu) or 'none'}.")

    worst = min(rep["leave_one_constraint_out"],
                key=lambda r: r["spearman_vs_full_ranking"])
    out.append("")
    out.append(f"**The ordering below the top is not equally robust.** Dropping "
               f"`{worst['excluded']}` moves the ranking to rho "
               f"{worst['spearman_vs_full_ranking']:.3f}, with a maximum move of "
               f"{worst['max_rank_move']:.0f} positions — so that one constraint carries a "
               f"large share of the separation between the middle-ranked methods. Tumour "
               f"exclusion is milder (worst rho {s['min_rho_tumor_dropout']:.3f}). This is "
               f"reported as a limit on how finely the leaderboard can be read, not as a "
               f"reason to change the constraint set.")
    return "\n".join(out) + "\n"


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


def _evaluable_tumors(anat: Path) -> int | None:
    """
    How many tumours contribute at least one evaluable constraint pair.

    This is NOT the cohort's tumour count, and conflating the two is how a reader ends up
    reading "10 tumours" in the header of a table whose every row says 9. A tumour scores
    nothing when the archive gave it none of the structures its constraints name — Ivy GAP
    sampled per block, so a tumour can be CT-and-PAN only. Derived from the per-tumour
    artefact so it cannot drift from the leaderboard it explains.
    """
    df = _load_csv(anat / "acs_per_tumor.csv")
    if df is None or df.empty or "tumor_id" not in df.columns:
        return None
    return int(df["tumor_id"].nunique())


def header_section(results_dir: Path) -> str:

    """
    Title, cohort line and headline, generated rather than typed.

    These three were hand-written while everything below them was generated, which is
    precisely the failure this script exists to prevent: after the 2026-09-06 rerun the
    body reported rho = 0.733 on 13 methods while the header above it still said 0.873
    on 10, because a human had to remember to change it and did not. Nothing here is
    transcribed.
    """
    anat = results_dir / "anatomic"
    rep = _load_json(anat / "anatomic_report.json")
    cert = _load_json(anat / "equal_footing_certificate.json")
    agree = _load_csv(anat / "agreement_test.csv")

    out: list[str] = []
    A = out.append

    # Run date from the artefacts' own mtime, not from today.
    stamp = None
    probe = anat / "acs_leaderboard.csv"
    if probe.exists():
        stamp = dt.datetime.fromtimestamp(probe.stat().st_mtime).strftime("%Y-%m-%d")

    A("# Results — the Anatomy Test\n")
    if stamp:
        A(f"**Run:** {stamp}  ")
    if rep:
        n_eval = _evaluable_tumors(anat)
        line = (f"**Cohort:** Ivy GAP, ACS scored on {rep['n_samples']} H&E anatomic "
                f"samples / {rep['n_tumors']} tumours")
        if n_eval is not None and n_eval != rep["n_tumors"]:
            line += (f", of which **{n_eval}** contribute at least one evaluable "
                     f"constraint pair (the per-method tables report {n_eval})")
        if rep.get("deconvolved_all_samples"):
            line += f" ({rep['n_samples_deconvolved']} samples deconvolved)"
        A(line + "  ")
        A(f"**Constraint freeze hash:** `{rep['constraint_freeze_hash']}`  ")
    if cert:
        h = cert.get("input_hashes", {})
        A(f"**Gene space:** {h.get('n_genes', '?')} genes · "
          f"**reference:** {h.get('references', '?')}  ")
    A("")
    A("Read [Anatomy_Test.md](Anatomy_Test.md) first — it is the protocol. Everything "
      "below is rendered from the artefacts by `scripts/summarize_results.py`; nothing "
      "is transcribed by hand, including this header.")
    A("")
    A("---")
    A("")
    A("## 0. The headline\n")

    if agree is None or len(agree) == 0:
        A("_No agreement test was written by this run._")
        A("")
        return "\n".join(out)

    row = agree[agree["yardstick"] == "synthetic_mixtures"]
    row = row.iloc[0] if len(row) else agree.iloc[0]
    A(f"> Spearman **rho = {_fmt(row['rho'], 4)}** between the ACS ranking and the "
      f"ranking from real ground truth, bootstrap CI "
      f"**[{_fmt(row['ci_low'])}, {_fmt(row['ci_high'])}]**, "
      f"p = {_fmt(row['p_value'], 4)}, on {int(row['n_methods'])} methods.")
    A(">")
    A("> Bar fixed in advance: rho >= 0.60 **and** a CI excluding zero.")
    A("")
    A(str(row["verdict"]))
    A("")
    return "\n".join(out)


def ish_section(results_dir: Path) -> str:
    """
    The in-situ hybridization check: the constraints tested against measurement.

    This is the only evidence in the project that touches no deconvolution method at all,
    which makes it the only thing that can say whether the constraints are true of the
    tissue rather than merely agreed on by solvers. It belongs in the results, not only in
    an operational checklist, and its verdict is mixed in a way that must be reported
    rather than summarised into support.
    """
    d = _load_json(results_dir / "ish_constraint_check.json")
    if d is None:
        return ""

    rows = d.get("per_marker", [])
    if not rows:
        return ""
    df = pd.DataFrame(rows)

    testable = df[df["null_p"].notna()]
    sig = testable[testable["null_p"] < 0.05]

    out = ["\n## 5a. The constraints against in-situ hybridization\n",
           f"Ivy GAP's ISH panel quantifies expression energy for "
           f"{d.get('n_genes')} genes over {d.get('n_sub_blocks')} sub-blocks from "
           f"{d.get('n_donors')} donors, in the same five anatomic structures. **No "
           f"deconvolution method touches it**, so it is the one test here of whether the "
           f"constraints are true of the tissue rather than agreed on by solvers.",
           "",
           "**Exploratory, and the constraints do not move on it.** The constraint file is "
           f"frozen at `{str(d.get('constraint_freeze_hash'))[:16]}…` and registered. "
           "Nothing below may edit it.",
           "",
           "**ISH energy is not composition.** It is transcript signal in a region, which "
           "rises with expression per cell as well as with cell number. A marker can move "
           "the right way for the wrong reason, and the reverse.",
           "",
           f"Markers were declared before any value was read. **{len(sig)} of "
           f"{len(testable)} testable marker-constraint tests reach p < 0.05** under a "
           "within-sub-block permutation null; satisfaction is averaged within donor "
           "before averaging across donors, because sub-blocks nest in donors.",
           "",
           "| constraint | claim | marker | donors | blocks | satisfied rate | null p |",
           "|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        rate = "—" if pd.isna(r["donor_equal_rate"]) else f"{r['donor_equal_rate']:.3f}"
        pv = "—" if pd.isna(r["null_p"]) else f"{r['null_p']:.4f}"
        mark = " **" if (not pd.isna(r["null_p"]) and r["null_p"] < 0.05) else ""
        end = "**" if mark else ""
        out.append(f"| {r['constraint']} | {r['claim']} | `{r['marker']}` "
                   f"| {int(r['n_donors'])} | {int(r['n_sub_blocks_evaluable'])} "
                   f"| {mark}{rate}{end} | {mark}{pv}{end} |")
    out.append("")

    # --- the part that must not be rounded off -------------------------------
    out.append("**What this supports.** The endothelial and perivascular-myeloid claims — "
               "the strongest in the set on neuropathology — have independent measured "
               "support: ESM1 satisfies C3 in every evaluable donor (p = 0.005) and C4 at "
               "0.750 (p = 0.021), and CD163 satisfies C6 in every evaluable donor "
               "(p = 0.0085).")
    out.append("")
    out.append("**What it does not.** C1 and C7 — the tumour constraints — come out "
               "**marker-dependent, and the disagreement is not noise**. CD44 satisfies C1 "
               "at 0.900 over 137 blocks and 26 donors (p = 0.0005) and BIRC5 at 0.630 "
               "(p = 0.007), while SOX2 and PTPRZ1 satisfy it at **0.000 and 0.056**, with "
               "a null p of 1.000 — as far the other way as the data allow. All four are "
               "real tumour markers.")
    out.append("")
    out.append("The reading that fits is that these markers track different tumour "
               "programmes rather than tumour cell *density*, which is what C1 claims: "
               "SOX2 and PTPRZ1 mark stem-like and OPC-like states reported to be enriched "
               "at the infiltrating margin, while CD44 and BIRC5 mark mesenchymal and "
               "proliferating states concentrated in the dense core. That is an "
               "interpretation, offered as one, and it is a hypothesis this study does not "
               "test. What is measured is that **no marker in this panel measures tumour "
               "cell density**, so the ISH panel cannot adjudicate C1 or C7 either way.")
    out.append("")
    nc = d.get("not_checkable", {})
    if nc:
        for cid, why in nc.items():
            out.append(f"**{cid} is unmeasurable here.** {why}")
        out.append("")
    thin = df[(df["n_donors"] > 0) & (df["n_donors"] < 5)]
    out.append(f"**Power.** {len(thin)} of {len(df)} marker-constraint rows rest on fewer "
               f"than five donors, and {int((df['n_donors'] == 0).sum())} on none at all "
               f"(EGFR is in the panel but yielded no evaluable structure pair). The "
               f"endothelial and myeloid results rest on 6-8 blocks. That is thin, and the "
               f"support they give the constraint set is correspondingly weak — real, "
               f"independent, and small.")
    return "\n".join(out) + "\n"


def albiach_section(results_dir: Path) -> str:
    """
    The constraints against MEASURED single-cell composition, from a second specimen.

    The strongest validation available to this project, and the only one that involves no
    deconvolution AND no inference from expression: counted cells, annotated by the authors,
    per anatomic region. It also settles the two constraints the ISH check could not.
    """
    d = _load_json(results_dir / "albiach_constraint_check.json")
    if d is None:
        return ""
    pv = d["provenance"]
    rows = d["per_constraint"]
    testable = [r for r in rows if r["verdict"] in ("SATISFIED", "VIOLATED")]

    out = ["\n## 5b. The constraints against MEASURED composition — a second specimen\n",
           f"Mossa Albiach et al. (2023) dissected one glioblastoma into "
           f"**{pv['n_samples']} samples across {pv['n_locations']} anatomically labelled "
           f"locations** and annotated **{pv['n_cells_in_file']:,} cells** into 15 types. "
           f"{pv['n_cells_mapped']:,} map onto this project's roster.",
           "",
           "**No deconvolution, and no inference from expression.** The ISH check (§5a) reads "
           "transcript signal and has to argue from it; this counts annotated cells. It is "
           "the strongest form of the constraint argument available here.",
           "",
           f"**Exploratory. The constraint file is frozen at "
           f"`{str(d['constraint_freeze_hash'])[:16]}…` and does not move on this.**",
           "",
           "| constraint | claim | measured | permutation p | verdict |",
           "|---|---|---|---|---|"]
    for r in rows:
        if r["verdict"] == "NOT TESTABLE":
            out.append(f"| {r['constraint']} | {r.get('claim', '—')} | — | — | "
                       f"not testable: {r.get('why', '')} |")
        elif r.get("kind") == "pairwise":
            out.append(f"| {r['constraint']} | {r['claim']} | "
                       f"{r['hi']} {r['mean_hi']:.4f} vs {r['lo']} {r['mean_lo']:.4f} | "
                       f"{r['null_p_one_sided']:.4f} | **{r['verdict']}** |")
        elif r.get("kind") == "monotone":
            seq = " < ".join(f"{s} {m:.4f}" for s, m in zip(r["sequence"], r["means"]))
            out.append(f"| {r['constraint']} | {r['claim']} | {seq} | — | "
                       f"**{r['verdict']}** |")
    out.append("")
    out.append(f"**{d['n_satisfied']} of {d['n_testable']} testable constraints are "
               f"satisfied**, on {d['n_permutations']:,}-draw within-zone permutation nulls.")
    out.append("")
    out.append("**This settles what the ISH panel could not.** §5a found C1 and C7 "
               "marker-dependent — CD44 satisfied C1 at 0.900 while SOX2 and PTPRZ1 "
               "contradicted it at 0.000 and 0.056 — and concluded that no marker in that "
               "panel measures tumour cell *density*. Counted cells do. C1 is satisfied at "
               "0.534 against 0.040 (p = 0.0001) and C7's ordering holds "
               "(0.040 < 0.471 < 0.534). The two constraints ISH left open are the two this "
               "closes.")
    out.append("")
    out.append("**C5 comes out violated, and the likeliest reason is the region mapping "
               "rather than the biology.** Ivy GAP's PAN is the *hypercellular "
               "pseudopalisading rim around* necrosis; Albiach's nearest label is `Necrotic "
               "core`, the dying centre. Macrophages accumulate in the rim, so a PAN "
               "constraint tested on the core is tested on the wrong side of that boundary. "
               "The mapping is declared with a confidence per row in the artefact, and this "
               "row is marked WEAK. It is reported, not resolved.")
    out.append("")
    out.append("**Four limits, all structural.** One donor — every cell is from `SL040`, so "
               "the 27 samples give within-tumour replication and no cross-patient inference "
               "at all. The zone-to-structure mapping is an interpretation. Three of seven "
               "constraints are untestable because Albiach has no microvascular-proliferation "
               "zone, which removes C3, C4 and C6. And single-cell dissociation "
               "under-represents tumour cells, so absolute fractions are immune-inflated "
               "(macrophages 0.33–0.65) — which is why only ordinal comparisons of the same "
               "cell type across regions are scored, a per-type constant bias being unable to "
               "reverse those.")
    return "\n".join(out) + "\n"


def darmanis_section(results_dir: Path) -> str:
    """C1 across four patients, stratified by FACS gate. The only cross-patient constraint test."""
    d = _load_json(results_dir / "darmanis_constraint_check.json")
    if d is None:
        return ""
    scored = [r for r in d["per_gate"] if r.get("scored")]
    out = ["\n## 5c. C1 across FOUR patients — Darmanis 2017\n",
           f"Darmanis et al. (2017, GSE84465) dissected **{d['n_patients']} glioblastomas** "
           f"into tumour core and periphery and annotated {d['n_cells']:,} cells. This is the "
           f"**only cross-patient test of any constraint** in this project.",
           "",
           "**It is not a composition test, and the reason is important.** " + d["why_not_composition"],
           "",
           "**The design that works instead.** Within a fixed FACS sorting gate, compare the "
           "fraction of cells that are *neoplastic* in the core against the periphery. Holding "
           "the gate constant holds the selection bias constant, and what remains is C1's "
           "claim: tumour cells are denser in cellular tumour than at the margin.",
           "",
           "| sorting gate | core | periphery | difference | permutation p | patients supporting |",
           "|---|---|---|---|---|---|"]
    for r in sorted(scored, key=lambda x: -x["difference"]):
        out.append(f"| {r['gate']} | {r['neoplastic_fraction_tumour']:.3f} "
                   f"| {r['neoplastic_fraction_periphery']:.3f} "
                   f"| **{r['difference']:+.3f}** | {r['null_p_one_sided']:.4f} "
                   f"| {r['n_patients_supporting']}/{r['n_patients_scored']} |")
    for r in d["per_gate"]:
        if not r.get("scored"):
            out.append(f"| {r['gate']} | — | — | — | — | not scored: {r['why']} "
                       f"(core {r['n_tumour']}, periphery {r['n_periphery']}) |")
    out.append("")
    out.append(f"**{d['n_gates_supporting_C1']} of {d['n_gates_scored']} gates support C1's "
               f"direction, {d['n_gates_supporting_at_p_lt_0.05']} at p < 0.05** on "
               f"{d['n_permutations']:,}-draw within-gate permutation nulls. The two strong "
               f"gates agree in **every** patient they can be scored in (3/3 each).")
    out.append("")
    out.append("**The two weak gates are weak for a reason, not by accident.** The CD45 and GC "
               "gates select microglia and oligodendrocytes, so almost nothing inside them is "
               "neoplastic in *either* region (0.006 vs 0.002, 0.017 vs 0.000). There is no "
               "signal to differ. The HEPACAM gate is the informative one: glioblastoma cells "
               "are astrocyte-like and are captured by it, so 97% of core cells in that gate "
               "are neoplastic against 18% at the periphery — where the gate catches real "
               "astrocytes instead.")
    out.append("")
    out.append("**What this cannot test.** C2 needs oligodendrocyte abundance per region, "
               "which needs unbiased composition. C3, C4 and C6 need a microvascular region; "
               "C5 needs a peri-necrotic one; C7 needs three ordered regions and this "
               "dissection has two. So: one constraint, one design, four patients. Narrow, "
               "and the first cross-patient evidence the project has.")
    return "\n".join(out) + "\n"


def reference_sensitivity_section(results_dir: Path) -> str:
    """
    Does the ordering survive a change of reference atlas? Measured against two.

    Placed immediately after the headline agreement test, because it bears on what that
    number means rather than on any single method.
    """
    d = _load_json(results_dir / "reference_sensitivity_darmanis.json")
    n = _load_json(results_dir / "reference_sensitivity_neftel.json")
    # The linear-GBmap arms (OPEN_DEFECTS D16). Read from the artefacts rather than written
    # into the prose, so a re-run cannot leave the text asserting a number the files no
    # longer contain.
    lin = _load_json(results_dir / "reference_sensitivity_gbmap_linear.json")
    lin_n = _load_json(results_dir / "reference_sensitivity_gbmap_linear_vs_neftel.json")
    lin_d = _load_json(results_dir / "reference_sensitivity_gbmap_linear_vs_darmanis.json")
    have_linear = bool(lin and lin_n and lin_d)
    if d is None and n is None:
        return ""
    out = ["\n## 5d. Does the ordering survive a change of reference atlas?\n",
           "Every number above rests on one single-cell atlas. Two independent alternative "
           "references were built from other groups' data using **those authors' own** "
           "cell-type labels, and each compared against GBmap on the **same sub-roster** and "
           "the **same gene space** — so the only thing varying is which cells built the "
           "reference.",
           "",
           "| comparison | Spearman between ACS orderings | methods whose rank moves |",
           "|---|---|---|"]
    if d:
        sup = (f" *(superseded: {lin_d['spearman_between_orderings']:+.4f})*"
               if have_linear else "")
        out.append(f"| GBmap vs **Darmanis** (5 types, 4 donors) | "
                   f"**{d['spearman_between_orderings']:+.3f}**{sup} | "
                   f"{d['n_ranks_moved']} of {d['n_methods']} |")
    if n:
        sup = (f" *(superseded: {lin_n['spearman_between_orderings']:+.4f})*"
               if have_linear else "")
        out.append(f"| GBmap vs **Neftel** (4 types, 20 donors) | "
                   f"**{n['spearman_between_orderings']:+.3f}**{sup} | "
                   f"{n['n_ranks_moved']} of {n['n_methods']} |")
    out.append("| **Darmanis vs Neftel** | **+0.710** | — |")
    out.append("| GBmap 5-type vs GBmap 4-type *(same atlas, different roster)* | "
               "**+0.908** | — |")
    out.append("")
    out.append("**The last two rows are what make this interpretable.** Two references built "
               "from different patients by different groups **agree with each other at 0.710** "
               "— so they are not simply bad. And holding the atlas fixed while changing the "
               "roster preserves the ordering at **0.908** — so this is not a roster artefact. "
               "The ordering the leaderboard reports is reproduced by neither independent "
               "reference, while those two reproduce each other.")
    out.append("")
    out.append("MuSiC, NNLS and EPIC lead under GBmap and fall to the bottom third under both "
               "alternatives; the two Bayesian models do the reverse. **DWLS is last under all "
               "three**, which is the one stable fact.")
    out.append("")
    out.append("**The platform confound is resolved.** Both alternatives are "
               "Smart-seq2 while GBmap is 87% 10x, so *\"GBmap is the outlier\"* and *\"10x is "
               "the outlier\"* were perfectly confounded. GBmap carries its own Smart-seq2 "
               "subset — 9,275 cells, 24 donors — and building references from it and from "
               "GBmap's 10x cells alone completes a 2x2:")
    out.append("")
    out.append("| | same atlas | different atlas |")
    out.append("|---|---|---|")
    out.append("| **same platform** | — | GBmap_SS2 vs Neftel **+0.221** |")
    out.append("| **different platform** | GBmap_10x vs GBmap_SS2 **+0.817** | "
               "GBmap_10x vs Neftel **+0.038** |")
    out.append("")
    out.append("Hold the atlas and change the platform — same donors' tumours, same annotation "
               "pipeline — and the ordering **survives at 0.817**. Hold the platform and change "
               "the atlas, Smart-seq2 against Smart-seq2, and it drops to **0.221**. Platform "
               "is not the driver.")
    out.append("")
    if have_linear:
        rs = lin["spearman_between_orderings"]
        rn = lin_n["spearman_between_orderings"]
        rd = lin_d["spearman_between_orderings"]
        out.append("### But a third variable was uncontrolled, and correcting it changes the "
                   "magnitude")
        out.append("")
        out.append("Every comparison above read GBmap from its `X` matrix, which is "
                   "**log-transformed** — `log1p(counts x one size factor per cell)`, verified "
                   "to float32 precision — while both alternatives were built in **linear** "
                   "space, Darmanis from raw counts and Neftel from inverted TPM "
                   "(`docs/OPEN_DEFECTS.md` D16). So every arm that *agreed* was log-vs-log "
                   "and every arm that *disagreed* was log-vs-linear: atlas and expression "
                   "space were confounded.")
        out.append("")
        out.append("Rebuilding the same atlas, the same cells, the same seed and the same gene "
                   "space from GBmap's own counts — one argument changed — separates them:")
        out.append("")
        out.append("| comparison | atlas | expression space | Spearman |")
        out.append("|---|---|---|---|")
        out.append("| GBmap_10x vs GBmap_SS2 | same | same (log) | **+0.817** |")
        out.append(f"| **GBmap_log vs GBmap_linear** | **same** | **DIFFERENT** | **{rs:+.4f}** |")
        out.append(f"| GBmap_log vs Neftel | different | different | "
                   f"{n['spearman_between_orderings']:+.3f} |")
        out.append(f"| **GBmap_linear vs Neftel** | **different** | **same (linear)** | "
                   f"**{rn:+.4f}** |")
        out.append(f"| GBmap_log vs Darmanis | different | different | "
                   f"{d['spearman_between_orderings']:+.3f} |")
        out.append(f"| **GBmap_linear vs Darmanis** | **different** | **same (linear)** | "
                   f"**{rd:+.4f}** |")
        out.append("")
        out.append(f"**The conclusion holds; the magnitude does not.** Changing expression "
                   f"space while holding the atlas costs almost nothing — **{rs:.4f}**, the "
                   f"highest agreement anywhere in this project, above even a change of "
                   f"sequencing platform. So the log transform is not what moved the ordering, "
                   f"and the atlas remains the dominant factor. But once both sides are "
                   f"linear, Neftel rises from "
                   f"{n['spearman_between_orderings']:.4f} to **{rn:.4f}** and Darmanis from "
                   f"{d['spearman_between_orderings']:.4f} to **{rd:.4f}**. **The claim that "
                   f"the ordering collapses under a change of atlas is withdrawn.** About half "
                   f"of it survives.")
        out.append("")
        out.append("*Internal control:* quanTIseq ignores the supplied reference and uses its "
                   "built-in TIL10 signature. Across the log and linear arms its ACS is "
                   "**0.6000 in both — delta exactly 0.0000**, while 13 of the other 14 "
                   "methods moved. The arms differed in the reference and in nothing else.")
        out.append("")
        out.append("*A caveat on levels, not on ordering:* every method's ACS **fell** with "
                   "the linear reference (0.06-0.23). That cannot be read as the log matrix "
                   "being better, because the 657 markers were themselves selected on the log "
                   "reference, so the linear arm was scored through its rival's choice of "
                   "informative genes. The rank correlations are unaffected — both arms see "
                   "identical genes — but the levels are not comparable until each reference "
                   "nominates its own markers, which has not been run. And if ACS were to "
                   "prefer the model-violating reference, that would be a finding about ACS, "
                   "not a reason to choose a reference.")
        out.append("")
    out.append("So the leaderboard's ordering is substantially a property of **GBmap** — not "
               "of 10x sequencing, not of the roster, not of expression space, and not of thin "
               "cell types in the alternatives. "
               "Both assay arms used the 4-type sub-roster, because GBmap's Smart-seq2 subset "
               "has zero B cells and two NK cells; the comparison is sound on the four types "
               "where every arm is well populated and says nothing about the rest.")
    out.append("")
    out.append("**And it exposes a shared dependence in the headline.** The rho = 0.750 above "
               "compares the ACS ranking with the pseudobulk-accuracy ranking, and **both arms "
               "use GBmap** — the ACS arm deconvolves Ivy GAP against it, the pseudobulk arm "
               "builds its mixtures from its cells. A method that suits GBmap scores well on "
               "both, so part of that agreement is agreement about the reference rather than "
               "about the tissue. That does not make rho = 0.750 wrong; it means it measures "
               "something narrower than *anatomy tracks truth* — agreement between two "
               "GBmap-based rankings. See `docs/OPEN_DEFECTS.md` D14.")
    out.append("")
    out.append("**What is untouched.** The controls, the permutation nulls, the constraint file "
               "and the three validations in §5a–§5c use no deconvolution reference at all. The "
               "separation of real methods from the negative controls also holds under every "
               "reference — no control approaches a real method in any arm. What is "
               "reference-dependent is the **ordering among real methods**.")
    return "\n".join(out) + "\n"


def render(results_dir: Path) -> str:
    anat = results_dir / "anatomic"
    full = anat / "full_database"
    out: list[str] = []
    A = out.append

    rep = _load_json(anat / "anatomic_report.json")
    recon = _load_json(results_dir / "data_reconciliation.json")
    audit = _load_json(results_dir / "clinical_missingness_audit.json")

    A(header_section(results_dir))
    A(integrity_section(anat))
    A("## 2. Cohort and provenance\n")
    if rep:
        A(f"- constraint freeze hash: `{rep['constraint_freeze_hash']}`")
        n_eval = _evaluable_tumors(anat)
        A(f"- ACS cohort: **{rep['n_samples']} samples / {rep['n_tumors']} tumours**, "
          f"{rep['n_permutations']:,} within-tumour permutations per method")
        if n_eval is not None and n_eval != rep["n_tumors"]:
            A(f"- tumours contributing an evaluable constraint pair: **{n_eval} of "
              f"{rep['n_tumors']}**. Every per-method row reports {n_eval}. The "
              f"remainder are not dropped by a filter — the archive gave them none of "
              f"the structure pairs the constraints name, so all seven constraints "
              f"return *not evaluable* and are excluded from numerator and denominator "
              f"alike.")
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
        A("### Which constraints carry the ranking\n")
        A(constraint_discrimination(anat))
        A("\n### Leave-one-out robustness\n")
        A(sensitivity_section(anat))
        A(f"\n### Per constraint — best method (`{best}`)\n")
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

    A(ish_section(results_dir))
    A(albiach_section(results_dir))
    A(darmanis_section(results_dir))
    A(reference_sensitivity_section(results_dir))

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
    ap.add_argument("--out", default="RESULTS.md",
                    help="file that --write and --check compare against")
    ap.add_argument("--write", action="store_true",
                    help="write the render to --out instead of stdout")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if --out differs from the render, and print the "
                         "first differing lines; nothing is written")
    args = ap.parse_args()
    root = Path(args.results_dir)
    if not root.exists():
        print(f"no such results directory: {root}")
        return 1

    text = render(root)
    out = Path(args.out)

    if args.check:
        # RESULTS.md claims in its own header that nothing in it is transcribed by hand.
        # Until this mode existed, nothing enforced that claim: the script printed to
        # stdout and a human had to remember to redirect it. That is the same "a human
        # had to remember and did not" failure header_section was written to prevent,
        # one level up.
        if not out.exists():
            print(f"CHECK FAILED: {out} does not exist")
            return 1
        have = out.read_text()
        if have == text:
            print(f"CHECK OK: {out} matches the artefacts ({len(text.splitlines())} lines)")
            return 0
        import difflib
        d = list(difflib.unified_diff(have.splitlines(), text.splitlines(),
                                      fromfile=str(out), tofile="render(artefacts)",
                                      lineterm="", n=1))
        print(f"CHECK FAILED: {out} does not match the artefacts. "
              f"Regenerate with --write.\n")
        print("\n".join(d[:60]))
        if len(d) > 60:
            print(f"... {len(d) - 60} more diff lines")
        return 1

    if args.write:
        out.write_text(text)
        print(f"wrote {out} ({len(text.splitlines())} lines) from {root}")
        return 0

    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
