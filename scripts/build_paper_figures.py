"""Render the paper's figures from artefacts.

Every value plotted is read from a file under `results/`. A figure showing a number nobody
computed is impossible rather than merely unlikely.

Writes PNG (300 dpi, for drafts) and PDF (vector, for submission) to `results/figures/`.

    python scripts/build_paper_figures.py
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402

from ivygap import config                                          # noqa: E402
from ivygap.paper_style import apply as apply_style                # noqa: E402
from ivygap import paper_style as PS                               # noqa: E402

apply_style()

# TRACKED, unlike results/ which is gitignored. Figures are deliverables: they are
# reviewed, marked up and compared between versions, so they have to be openable by
# someone who did not run the pipeline.
FIG = config.PROJECT_ROOT / "docs" / "figures"
# Journal figures are read in greyscale as often as not, so the palette is distinguishable
# by lightness as well as hue.
C_T, C_NK, C_B = PS.BLUE, PS.LIGHT, PS.RED
C_TRUTH, C_GBM, C_LGG = PS.INK, PS.BLUE, PS.ORANGE


def J(n):
    p = config.RESULTS_DIR / n
    return json.loads(p.read_text()) if p.exists() else {}


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    # derive the path from FIG rather than hardcoding it: this line said
    # "results/figures/" for a while after the output moved to docs/figures/,
    # which is how someone ends up editing a stale copy.
    print(f"  wrote {FIG.relative_to(config.PROJECT_ROOT)}/{name}.png and .pdf")


def fig_lymphoid():
    """THE headline figure: every method against per-type methylation truth, both cohorts."""
    panels = [(J("lymphoid_ordering.json"), "A   Glioblastoma (GBM)"),
              (J("lymphoid_ordering_lgg.json"), "B   Lower-grade glioma (LGG)")]
    panels = [(d, t) for d, t in panels if d]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(15.5, 6.4))
    fig.subplots_adjust(wspace=0.42, top=0.84, bottom=0.16)
    axes = np.atleast_1d(axes)
    for ax, (d, title) in zip(axes, panels):
        ms = d["methods"]
        order = sorted(ms, key=lambda m: -ms[m]["T"])
        labels = order + ["", "DNA methylation"]
        T = [ms[m]["T"] for m in order] + [0, d["truth_mean"]["T_cell"]]
        NK = [ms[m]["NK"] for m in order] + [0, d["truth_mean"]["NK_cell"]]
        B = [ms[m]["B"] for m in order] + [0, d["truth_mean"]["B_cell"]]
        y = np.arange(len(labels))
        ax.barh(y, T, color=C_T, label="T cell")
        ax.barh(y, NK, left=T, color=C_NK, label="NK cell")
        ax.barh(y, B, left=np.array(T) + np.array(NK), color=C_B, label="B cell")
        # "found no lymphocytes at all" is the SECOND failure mode and has to be visible
        # without colliding with the next panel's labels, so it sits inside the axes.
        for i, m in enumerate(order):
            z, n = ms[m].get("n_no_lymphoid_signal", 0), ms[m].get("n", 0)
            if z and n:
                ax.text(0.985, i, f"{z}/{n} zero", va="center", ha="right", fontsize=7.4,
                        color="white" if z / n > .35 else "#333",
                        fontweight="bold" if z / n > .5 else "normal")
        ax.set_yticks(y)
        ax.set_yticklabels([PS.display(l) for l in labels], fontsize=9.5)
        ax.get_yticklabels()[-1].set_fontweight("bold")
        ax.axhline(len(order) + 0.5, color="#999", lw=0.9, ls="--")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Relative composition within {T, NK, B}", fontsize=9.5)
        ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold", pad=10)
        ax.invert_yaxis()
    # one legend for both panels, below the figure so it covers no data
    h, lg = axes[0].get_legend_handles_labels()
    fig.legend(h, lg, loc="lower center", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Relative composition within {T, NK, B}, by method and by methylation",
                 fontsize=12.5, fontweight="regular", y=0.96)
    fig.text(0.5, 0.905, "Grey labels, number of samples returning exactly zero T, NK and B.",
             ha="center", fontsize=9, color="#555")
    save(fig, "Figure_lymphoid_failure")


def fig_model_fit():
    """The bound: what the additive model explains, bracketed by floor and ceiling."""
    mf = J("model_fit_residual.json")
    if not mf:
        return
    co = mf["cohorts"]
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    fig.subplots_adjust(top=0.80, bottom=0.22)
    labs = list(co)
    x = np.arange(len(labs))
    real = [co[k]["controls"]["real_reference"] for k in labs]
    shuf = [co[k]["controls"]["shuffled_profiles"] for k in labs]
    rand = [co[k]["controls"]["random_basis"] for k in labs]
    ceil = [co[k]["controls"]["ceiling_top_k_svd"] for k in labs]
    ax.bar(x, ceil, width=.5, color="#e8e8e8", label="ceiling: best possible 8-dim basis")
    ax.bar(x, real, width=.5, color=C_GBM, label="the reference actually used")
    ax.plot(x, shuf, "v", color=C_B, ms=10, label="floor: gene→cell-type labels shuffled",
            clip_on=False, zorder=5)
    ax.plot(x, rand, "o", color="#666", ms=7, label="floor: random non-negative basis",
            clip_on=False, zorder=5)
    for i, k in enumerate(labs):
        f = co[k]["controls"]["fraction_of_achievable"]
        # value label INSIDE the navy bar, so it cannot collide with the legend above
        ax.text(i, real[i] / 2, f"{real[i]:.3f}\n({f:.0%} of ceiling)", ha="center",
                va="center", fontsize=10, fontweight="bold", color="white")
        ax.text(i + .30, ceil[i], f"{ceil[i]:.3f}", ha="left", va="center", fontsize=9,
                color="#555")
        ax.text(i + .30, shuf[i] - .015, f"{shuf[i]:+.4f}", ha="left", va="top", fontsize=8,
                color=C_B)
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=11)
    ax.set_ylabel("R² of the best non-negative fit to real bulk", fontsize=10)
    ax.set_ylim(-.10, 1.05)
    ax.set_xlim(-.6, len(labs) - .15)
    ax.axhline(0, color="#bbb", lw=.9)
    ax.set_title("Per-sample fit of the non-negative mixing model, against floor and "
                 "ceiling controls",
                 fontsize=12, fontweight="regular", loc="left", pad=12)
    h, lg = ax.get_legend_handles_labels()
    fig.legend(h, lg, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.54, 0.0))
    save(fig, "Figure_model_fit_bound")


def fig_recovery():
    """Tumour-content recovery per method, both cohorts, degenerate rows marked."""
    efg, efl = J("equal_footing_ranking.json"), J("equal_footing_ranking_lgg.json")
    if not efg:
        return
    g, l = efg.get("recovery_frozen", {}), efl.get("recovery_frozen", {})
    nc = set(efg.get("non_comparable_frozen", {}))
    ms = sorted(g, key=lambda m: -g[m])
    y = np.arange(len(ms))
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.barh(y - .2, [g[m] * 100 for m in ms], height=.38, color=C_GBM, label="GBM")
    ax.barh(y + .2, [l.get(m, np.nan) * 100 if m in l else np.nan for m in ms],
            height=.38, color=C_LGG, label="LGG")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{PS.display(m)}  \u2020" if m in nc else PS.display(m)
                        for m in ms], fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="#999", lw=.9)
    ax.set_xlabel("Recovery of true tumour-content variation (%)\n"
                  "0, no information; 100, perfect.", fontsize=9)
    ax.set_title("Recovery of true tumour-content variation, by method and cohort",
                 fontsize=11.5, fontweight="regular", loc="left")
    ax.legend(frameon=False, fontsize=9)
    ax.text(0.99, 0.02, "\u2020 not evaluable \u2014 ran without an input its published "
                        "algorithm requires",
            transform=ax.transAxes, ha="right", fontsize=7.6, color="#666")
    save(fig, "Figure_tumour_recovery")


def fig_detects_not_ranks():
    """Result 1: ACS's ranking predicts accuracy ONLY against the yardstick that shares its atlas."""
    fz, h5 = J("yardstick_agreement.json"), J("yardstick_agreement_h5ad.json")
    if not fz:
        return
    bar = fz.get("registered_bar", 0.60)
    a1 = fz["arm_1_pseudobulk"]
    rows = [("Synthetic pseudobulk\n(built from GBmap, the atlas the ACS arm deconvolves against)",
             a1["rho"], None, None, True)]
    for d, lab in ((fz, "frozen signature"), (h5, "sigma-carrying reference")):
        if not d:
            continue
        a2 = d["arm_2_absolute_purity"]
        for key, tag in (("all_methods", "all methods"),
                         ("excluding_degenerate", "comparable methods only")):
            v = a2.get(key)
            if not v:
                continue
            rows.append((f"DNA tumour purity\n({lab}, {tag}, n={v['n_methods']})",
                         v["spearman"], v["ci"][0], v["ci"][1], False))
    fig, ax = plt.subplots(figsize=(10.2, 4.6))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.86, bottom=0.20)
    y = np.arange(len(rows))[::-1]
    for i, (lab, rho, lo, hi, shares) in zip(y, rows):
        col = C_B if shares else C_GBM
        if lo is not None:
            ax.plot([lo, hi], [i, i], color=col, lw=2.4, alpha=.40, solid_capstyle="round")
        ax.plot(rho, i, "o", color=col, ms=11, zorder=5)
        # Nudge the label off the pre-registered bar line when it would sit on top of it;
        # at rho = 0.637 against a bar of 0.60 the two collided.
        dx, ha = (0.0, "center")
        if abs(rho - bar) < 0.12:
            dx, ha = (-0.055, "right")
        ax.text(rho + dx, i + .22, f"{rho:+.3f}", ha=ha, fontsize=9.5, color=col)
    ax.axvline(bar, color="#111", ls="--", lw=1.4)
    # annotation sits at the BOTTOM of the axes, clear of both the title and the data
    ax.text(bar - .03, -0.62, f"pre-registered bar, \u03c1 \u2265 {bar}", fontsize=9,
            ha="right", va="center")
    ax.axvline(0, color="#ccc", lw=.9)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_ylim(-1.0, len(rows) - 0.4)
    ax.set_xlim(-1.02, 1.02)
    ax.set_xlabel("Spearman correlation between the two rankings\n"
                  "Bars, 95% bootstrap CI over methods.", fontsize=9.5)
    fig.text(0.015, 0.965, "Rank correlation between the ACS ordering and each yardstick's "
                           "ordering", fontsize=12, ha="left", va="top")
    save(fig, "Figure_detects_not_ranks")


def fig_bisque_anchoring():
    """Bisque without overlapping subjects returns its reference's composition (OPEN_DEFECTS D23).

    A: planted truth -- genuine BisqueRNA on synthetic data whose true cohort mean is far from
    the reference's. B: every method's TCGA cohort mean, distance from the reference's donor-mean
    composition, raw/X reference. Both panels read from results/; nothing is typed in.
    """
    syn, real = J("bisque_anchoring_synthetic.json"), J("bisque_anchoring.json")
    if not syn or not real:
        print("  SKIP Figure_bisque_anchoring: artefacts missing"); return
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1, 1.25]})
    far = syn["far"]; types = list(far["planted_mean"])
    prior = syn["reference_prior_donor_mean"]
    x = np.arange(len(types)); w = 0.27
    for k, (lab, vals, col) in enumerate([
            ("Planted (true) cohort mean", far["planted_mean"], C_TRUTH),
            ("Bisque estimate, cohort mean", far["estimated_mean"], PS.RED),
            ("Reference donor-mean composition", prior, PS.LIGHT)]):
        a.bar(x + (k - 1) * w, [vals[t] for t in types], w, label=lab, color=col,
              edgecolor=PS.INK, linewidth=0.6)
    a.set_xticks(x); a.set_xticklabels([f"Type {t}" for t in types], fontsize=10)
    a.set_ylabel("Cohort-mean proportion", fontsize=10)
    a.set_title("A  Planted truth, genuine BisqueRNA (no overlapping subjects)", loc="left",
                fontsize=11, fontweight="bold", pad=10)
    a.text(0.02, 0.97, f"L1 to truth {far['L1_estimate_to_plant']:.3f}\n"
                       f"L1 to reference {far['L1_estimate_to_prior']:.3f}\n"
                       f"per-sample rho {min(far['per_sample_spearman'].values()):.2f}"
                       f"-{max(far['per_sample_spearman'].values()):.2f}",
           transform=a.transAxes, va="top", fontsize=9, color=PS.INK)
    a.set_ylim(0, 0.84)
    a.legend(frameon=False, fontsize=8.6, loc="upper right")
    g, l = real["cohorts"]["gbm"]["methods"], real["cohorts"]["lgg"]["methods"]
    ms = sorted(set(g) & set(l), key=lambda m: g[m]["L1_to_reference_prior_all_types"])
    y = np.arange(len(ms))
    b.scatter([g[m]["L1_to_reference_prior_all_types"] for m in ms], y, color=C_GBM, s=34,
              label="Glioblastoma", zorder=3)
    b.scatter([l[m]["L1_to_reference_prior_all_types"] for m in ms], y, color=C_LGG, s=34,
              marker="s", label="Lower-grade glioma", zorder=3)
    for i, m in enumerate(ms):
        b.plot([g[m]["L1_to_reference_prior_all_types"], l[m]["L1_to_reference_prior_all_types"]],
               [i, i], color=PS.LIGHT, linewidth=1, zorder=2)
    b.set_yticks(y); b.set_yticklabels([PS.display(m) for m in ms], fontsize=9)
    b.invert_yaxis(); b.set_xlim(left=0)
    b.set_xlabel("L1 distance of the TCGA cohort-mean composition from the reference's\n"
                 "donor-mean composition, all eight cell types (raw/X reference)", fontsize=9.5)
    b.set_title("B  TCGA cohorts, every method", loc="left", fontsize=11, fontweight="bold", pad=10)
    b.legend(frameon=False, fontsize=9, loc="lower left", bbox_to_anchor=(0.18, 0.02))
    fig.suptitle("Cohort-mean composition against the single-cell reference's own composition",
                 fontsize=12.5, fontweight="regular", y=1.0)
    fig.tight_layout()
    save(fig, "Figure_bisque_anchoring")


def fig_lymphoid_mechanisms():
    """Reference-side tests of the lymphoid inversion, SVR (CIBERSORT core), frozen reference.

    A: share of samples with B above T, per reference variant. B: where the B column's estimate
    goes when the column is removed. C: correlation of each lymphoid profile with the mean bulk
    tumour profile. Every value from results/; nothing typed in.
    """
    nk, rp = J("nk_reannotation.json"), J("ribosomal_test.json")
    bd, tl = J("b_column_diagnostics.json"), J("b_profile_tissue_likeness.json")
    if not (nk and rp and bd and tl):
        print("  SKIP Figure_lymphoid_mechanisms: artefacts missing"); return
    cohorts = [c for c in ("gbm", "lgg") if c in nk.get("cohorts", {}) and c in rp.get("cohorts", {})
               and "b_removed_svr" in bd.get("cohorts", {}).get(c, {})]
    if len(cohorts) < 2:
        print("  SKIP Figure_lymphoid_mechanisms: LGG not complete yet"); return
    lab = {"gbm": "Glioblastoma", "lgg": "Lower-grade glioma"}
    colr = {"gbm": C_GBM, "lgg": C_LGG}
    ia, iai = J("independent_atlas_test.json"), J("independent_atlas_test_ig_removed.json")
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 9.6))
    (a, b), (c, dd) = axs
    variants = [("as registered", lambda co: nk["cohorts"][co]["variants"]["as_registered"]["svr"]["frac_B_over_T"]),
                ("NK column\nremoved", lambda co: nk["cohorts"][co]["variants"]["nk_removed"]["svr"]["frac_B_over_T"]),
                ("T and NK\nmerged", lambda co: nk["cohorts"][co]["variants"]["t_nk_merged"]["svr"]["frac_B_over_T"]),
                ("ribosomal genes\nremoved", lambda co: rp["cohorts"][co]["svr_rp_removed"]["frac_B_over_T"])]
    x = np.arange(len(variants)); w = 0.36
    for k, co in enumerate(cohorts):
        vals = [100 * f(co) for _, f in variants]
        a.bar(x + (k - 0.5) * w, vals, w, color=colr[co], edgecolor=PS.INK, linewidth=0.5, label=lab[co])
        meth = 100 * (1 - nk["cohorts"][co]["methylation"]["frac_T_over_B"])
        a.axhline(meth, color=colr[co], linestyle="--", linewidth=1)
    a.set_xticks(x); a.set_xticklabels([v for v, _ in variants], fontsize=9)
    a.set_ylim(0, 122); a.set_ylabel("Samples with B above T (%)", fontsize=10)
    a.set_title("A  B above T, by reference variant (SVR)", loc="left", fontsize=11, fontweight="bold", pad=10)
    a.plot([], [], ls="--", color=PS.INK, lw=1, label="DNA methylation (truth)")
    a.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=3)
    order = ["Tumor", "T_cell", "NK_cell"]
    for k, co in enumerate(cohorts):
        sh = bd["cohorts"][co]["b_removed_svr"]["share_of_B_mass_to"]
        left = 0.0
        parts = [(t, max(sh.get(t, 0.0), 0.0)) for t in order]
        other = max(0.0, sum(max(v, 0.0) for t, v in sh.items() if t not in order))
        for t, v in parts + [("other", other)]:
            b.barh(k, 100 * v, left=100 * left, color={"Tumor": PS.INK, "T_cell": C_T, "NK_cell": C_NK, "other": PS.LIGHT}[t],
                   edgecolor="white", linewidth=0.6, label=PS.display(t) if (k == 0 and t != "other") else ("other types" if k == 0 else None))
            if v >= 0.08:
                b.text(100 * (left + v / 2), k, f"{100 * v:.0f}%", ha="center", va="center", fontsize=8.5,
                       color="white" if t in ("Tumor", "T_cell") else PS.INK)
            left += v
    b.set_yticks(range(len(cohorts))); b.set_yticklabels([lab[c] for c in cohorts], fontsize=9.5)
    b.set_xlabel("Share of the removed B estimate re-absorbed (%)", fontsize=9.5)
    b.set_title("B  Where B's estimate goes when B is removed (SVR;\n    under NNLS it goes to T instead)", loc="left", fontsize=11, fontweight="bold", pad=10)
    b.legend(frameon=False, fontsize=8.5, loc="lower right", bbox_to_anchor=(1.0, -0.42), ncol=4)
    types = ["T_cell", "NK_cell", "B_cell"]; xx = np.arange(len(types))
    for k, co in enumerate(cohorts):
        r = tl["cohorts"][co]["r_with_mean_bulk"]
        c.bar(xx + (k - 0.5) * w, [r[t]["all_genes"] for t in types], w, color=colr[co], edgecolor=PS.INK,
              linewidth=0.5, label=f"GBmap, {lab[co]}")
        if ia:
            pp = ia["cohorts"][co]["profile_properties"]
            for j, t in enumerate(types):
                if t in pp:
                    c.scatter(xx[j] + (k - 0.5) * w, pp[t]["r_with_mean_bulk"], marker="D", s=46, color="white",
                              edgecolor=PS.INK, linewidth=1.2, zorder=4,
                              label="independent atlas" if (k == 0 and j == 0) else None)
    c.axhline(0, color=PS.INK, linewidth=0.6)
    c.set_xticks(xx); c.set_xticklabels([PS.display(t) for t in types], fontsize=9.5)
    c.set_ylabel("Correlation with the mean bulk profile\n(log expression)", fontsize=9.5)
    c.set_title("C  Similarity of each lymphoid profile to bulk tumour", loc="left", fontsize=11, fontweight="bold", pad=10)
    c.legend(frameon=False, fontsize=8.3, loc="upper left")
    # panel D: the independent-atlas test
    if ia:
        conds = [("GBmap\n(frozen)", lambda co, m: ia["cohorts"][co]["methods"][m]["gbmap_frozen_for_comparison"]["frac_samples_B_over_T"]),
                 ("independent\natlas", lambda co, m: ia["cohorts"][co]["methods"][m]["frac_samples_B_over_T"])]
        if iai:
            conds.append(("independent atlas,\nIg genes removed", lambda co, m: iai["cohorts"][co]["methods"][m]["frac_samples_B_over_T"]))
        groups = [(co, m) for co in cohorts for m in ("nnls", "svr")]
        gx = np.arange(len(groups)); ww = 0.8 / len(conds)
        shades = [PS.INK, PS.RED, PS.LIGHT]
        for k, (nm, f) in enumerate(conds):
            vals = [100 * (f(co, m) or 0) for co, m in groups]
            dd.bar(gx + (k - (len(conds) - 1) / 2) * ww, vals, ww, color=shades[k], edgecolor=PS.INK, linewidth=0.5,
                   label=nm.replace("\n", " "))
        for j, (co, m) in enumerate(groups):
            meth = 100 * (1 - ia["cohorts"][co]["methylation_frac_T_over_B"])
            dd.plot([j - 0.42, j + 0.42], [meth, meth], color=colr[co], linestyle="--", linewidth=1)
        dd.set_xticks(gx); dd.set_xticklabels([f"{lab[co].split()[0] if co == 'gbm' else 'LGG'}\n{m.upper()}" for co, m in groups], fontsize=9)
        dd.set_ylim(0, 122); dd.set_ylabel("Samples with B above T (%)", fontsize=10)
        dd.set_title("D  B above T: GBmap versus an independent atlas", loc="left", fontsize=11, fontweight="bold", pad=10)
        dd.plot([], [], ls="--", color=PS.INK, lw=1, label="DNA methylation (truth)")
        dd.legend(frameon=False, fontsize=8.0, loc="upper center", ncol=2)
    else:
        dd.axis("off")
    fig.suptitle("Reference-side tests of the lymphoid ordering, TCGA cohorts",
                 fontsize=12.5, fontweight="regular", y=1.01)
    fig.tight_layout()
    save(fig, "Figure_lymphoid_mechanisms")


def fig_acs_vs_auc():
    """ACS and a per-method AUC: one construct on two scales, and neither tracks accuracy."""
    d = J("anatomic_auc.json")
    if not d or not d.get("controls", {}).get("acs_reproduced_for_every_method"):
        return
    ms, cp = d["methods"], d.get("c_purity", {})
    real = {m: r for m, r in ms.items() if not r["is_control"] and r["comparable"]}
    ctrl = {m: r for m, r in ms.items() if r["is_control"]}
    lb = pd.read_csv(config.ANATOMIC_DIR / "acs_leaderboard.csv")
    lb = lb[(~lb["is_control"].astype(bool)) & (lb["comparable"].astype(bool))]
    acs_chance = float(lb["null_mean"].median())
    va = d["agreement_with_auc_in_place_of_acs"]

    fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 5.0))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.84, bottom=0.14, wspace=0.28)
    # A -- the two anatomy scores against each other
    a.plot([0.3, 1.0], [0.3, 1.0], color="#bbb", lw=1, ls=":")
    a.axvline(acs_chance, color=PS.RED, lw=1, ls="--", alpha=.6)
    a.axhline(0.5, color=PS.BLUE, lw=1, ls="--", alpha=.6)
    a.text(acs_chance + .005, 0.31, f"ACS chance\n{acs_chance:.2f}", color=PS.RED, fontsize=8)
    a.text(0.99, 0.505, "AUC chance 0.50", color=PS.BLUE, fontsize=8, va="bottom", ha="right")
    for m, r in real.items():
        a.plot(r["acs"], r["auc_anat"], "o", color=PS.INK, ms=6)
    for m, r in ctrl.items():
        a.plot(r["acs"], r["auc_anat"], "X", color=PS.RED, ms=9)
        above = r["auc_anat"] >= 0.5          # keep each label off the chance line
        a.text(r["acs"] + .012, r["auc_anat"] + (.015 if above else -.015),
               m.replace("control_", "control: "), fontsize=7.5, color=PS.RED,
               va="bottom" if above else "top")
    a.set_xlim(0.3, 1.0)
    a.set_ylim(0.3, 1.0)
    a.set_xlabel("ACS (registered): share of constraint-tumour pairs satisfied")
    a.set_ylabel("Anatomic AUC (exploratory): within-tumour\nMann-Whitney probability")
    a.set_title(f"A  Same construct, different scales\nrank agreement ρ = "
                f"{d['rank_agreement']['acs_vs_auc_anat_spearman_comparable_methods']:.2f} "
                f"({d['rank_agreement']['n_comparable_methods']} methods)", loc="left", fontsize=10)
    # B -- anatomy score against DNA-measured accuracy
    xs = [m for m in real if m in cp]
    b.axhline(0.5, color="#999", lw=1, ls="--")
    # Degenerate methods give byte-identical points (SCDC ENSEMBLE with one reference is
    # SCDC); label the point once with both names rather than letting one hide the other.
    groups: dict[tuple, list[str]] = {}
    for m in xs:
        groups.setdefault((round(real[m]["auc_anat"], 9), round(cp[m]["c"], 9)), []).append(m)
    tcga = (J("absolute_purity_yardstick.json").get("methods") or {})
    same = {"music": "nnls", "scdc_ensemble": "scdc"}
    for (x, y), names in groups.items():
        b.plot(x, y, "o", color=PS.INK, ms=6)
        # A method degenerate on the TCGA frozen signature is never shown under its own name
        # alone (CLAUDE.md): its concordance there is another method's.
        lab = " = ".join(sorted(names))
        if any((tcga.get(n) or {}).get("degenerate") for n in names):
            lab += "\u2020"
        b.text(x + .004, y + .004, lab, fontsize=7, color="#444")
    b.text(0.99, 0.01, "\u2020 degenerate on the frozen TCGA signature: MuSiC = NNLS there; Bisque in its "
           "no-overlap mode", transform=b.transAxes, ha="right", va="bottom", fontsize=7, color="#666")
    b.text(b.get_xlim()[0], 0.502, " chance", fontsize=8, color="#777", va="bottom")
    b.set_xlabel("Anatomic AUC (agreement with anatomy)")
    b.set_ylabel("Concordance with ABSOLUTE DNA purity\n(accuracy, TCGA-GBM)")
    b.set_title(f"B  Neither anatomy score tracks accuracy\nρ = "
                f"{va['vs_rho_purity']['rho']:+.2f} (AUC) vs "
                f"{va['registered_for_comparison']['rho']:+.3f} (ACS, registered); "
                f"n = {va['vs_rho_purity']['n_methods']}", loc="left", fontsize=10)
    save(fig, "Figure_acs_vs_auc")


def fig_identifiability():
    """Extension E2: truth-free stability against agreement with DNA truths, per compartment."""
    d = J("identifiability_diagnostics.json")
    if not d or not all(c.get("pass") for c in d.get("controls", {}).values()):
        return
    U = d["units"]
    names = {"Tumor": "Tumour", "Leukocytes": "Leukocytes", "Lymphoid": "Lymphoid",
             "T_cell": "T", "B_cell": "B", "NK_cell": "NK"}
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.84, bottom=0.17, wspace=0.26)
    panels = (("d1_stability", "truth_rho_unmix", "H1_primary", "D1",
               "Loss-scale stability (DESeq2 unmix, 7 settings):\nmean pairwise Spearman across settings",
               "Agreement with DNA truth, unmix\n(Spearman across samples)",
               "A  Loss-scale stability against agreement with DNA truth"),
              ("d2_agreement", "truth_rho_median_methods", "H1_primary", "D2",
               "Method agreement (registered panel):\nmean pairwise Spearman across methods",
               "Agreement with DNA truth,\nmedian over the registered methods",
               "B  Method agreement against agreement with DNA truth"))
    # Label placement (dx, dy, ha) where two units nearly coincide; everything else sits up-right.
    # Placement only -- every coordinate plotted is read from the artefact.
    place = {("d1_stability", "Tumor|gbm"): (-.012, .035, "right"),
             ("d1_stability", "Leukocytes|gbm"): (-.012, -.06, "right"),
             ("d1_stability", "Leukocytes|lgg"): (-.012, .03, "right"),
             ("d1_stability", "Tumor|lgg"): (0, -.075, "center"),
             ("d2_agreement", "T_cell|lgg"): (-.006, -.055, "right"),
             ("d2_agreement", "NK_cell|lgg"): (.006, -.06, "left"),
             ("d2_agreement", "B_cell|gbm"): (-.01, .02, "right")}
    for ax, (xk, yk, hk, dk, xl, yl, title) in zip(axes, panels):
        ax.axhline(0, color="#999", lw=1, ls="--")
        for k, u in U.items():
            if u["truth_tier"] == "none":
                continue
            col = C_GBM if u["cohort"] == "gbm" else C_LGG
            mk = "o" if u["cohort"] == "gbm" else "s"
            primary = u["truth_tier"] == "primary"
            ax.plot(u[xk], u[yk], mk, ms=8 if primary else 6.5, mec=col,
                    mfc=col if primary else "white", mew=1.6)
            dx, dy, ha = place.get((xk, k), (.008, .012, "left"))
            ax.text(u[xk] + dx, u[yk] + dy, names[u["compartment"]], fontsize=7.5, ha=ha,
                    color=col, fontweight="bold" if primary else "normal")
        h = d[hk][dk]
        hs = d["H1s_secondary"][dk]
        ax.set_title(f"{title}\nfilled, 6 primary units: ρ = {h['rho']:.2f}, exact p = {h['p_one_sided']:.3f};"
                     f" all 12: ρ = {hs['rho']:.2f}, p = {hs['p_one_sided']:.3f}", loc="left", fontsize=9.5)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
    lo = min(u["d1_stability"] for u in U.values()) - .05
    axes[0].set_xlim(lo, 1.0)
    axes[1].set_xlim(0.0, max(u["d2_agreement"] for u in U.values()) + .1)
    for ax in axes:
        ax.set_ylim(-0.4, 0.95)
    notruth = [u for u in U.values() if u["truth_tier"] == "none"]
    fig.text(0.07, 0.015,
             "Circles: TCGA-GBM; squares: TCGA-LGG. Filled: primary units (truths: ABSOLUTE purity, methylation "
             "leukocyte fraction, EpiDISH lymphoid total); open: T, B and NK (secondary). Myeloid, endothelial and "
             f"oligodendrocyte have no DNA truth (stability {min(u['d1_stability'] for u in notruth):.2f}-"
             f"{max(u['d1_stability'] for u in notruth):.2f}). Post-registration, exploratory.",
             fontsize=7, color="#555")
    save(fig, "Figure_identifiability")


def _mean_comp(d: dict) -> tuple[list[float], int]:
    """Mean within-lymphoid composition over the methods that report one (each sums to 1)."""
    ms = [v for v in d.get("methods", {}).values() if all(v.get(k) is not None for k in ("T", "NK", "B"))]
    if not ms:
        return [np.nan] * 3, 0
    return [float(np.mean([m[k] for m in ms])) for k in ("T", "NK", "B")], len(ms)


def _atlas_shares(rows: list[dict], prefix: tuple[str, ...]) -> tuple[list[float], int]:
    """Patient-mean within-lymphoid shares (strict T, NK, B) for the atlas patients of one diagnosis."""
    sh = []
    for r in rows:
        if r["patient"].startswith(prefix):
            v = np.array([r["T_strict"], r["NK"], r["B"]], float)
            if v.sum() > 0:
                sh.append(v / v.sum())
    return (list(np.mean(sh, axis=0)) if sh else [np.nan] * 3), len(sh)


def fig_truth_instruments():
    """The lymphoid truth by every instrument that measures it, beside what the methods report
    (MANUSCRIPT §4.4; WHY_B_OVER_T §7n). Direct measurement first, then DNA methylation, then bulk-RNA
    deconvolution. Every bar is a within-lymphoid composition {T, NK, B}."""
    lo, ll = J("lymphoid_ordering.json"), J("lymphoid_ordering_lgg.json")
    lh, llh = J("lymphoid_ordering_h5ad.json"), J("lymphoid_ordering_lgg_h5ad.json")
    gim, atl, gbmap, kl = (J("gimicc_truth_confirmation.json"), J("atlas_t_vs_b.json"),
                           J("gbmap_t_vs_b.json"), J("klemm_t_vs_b.json"))
    if not all((lo, ll, gim, atl, gbmap, kl)):
        return
    gm = lambda c: [gim["cohorts"][c]["Q1_all_methylation_samples"]["mean_within_lymphoid_share"][k]
                    for k in ("T", "NK", "B")]                                          # noqa: E731
    ep = lambda d: [d["truth_mean"][k] for k in ("T_cell", "NK_cell", "B_cell")]          # noqa: E731
    kw = lambda key: [kl[key]["within_lymphoid_share"][k] for k in ("T", "NK", "B")]      # noqa: E731
    a_gbm, n_ag = _atlas_shares(atl["per_patient"], ("ndGBM", "rGBM"))
    a_lgg, n_al = _atlas_shares(atl["per_patient"], ("LGG",))
    gbs = gbmap["pooled_within_lymphoid_share"]
    m_lo, n_lo = _mean_comp(lo); m_lh, n_lh = _mean_comp(lh)
    m_ll, n_ll = _mean_comp(ll); m_llh, n_llh = _mean_comp(llh)
    panels = [
        ("A   Glioblastoma  /  IDH-wildtype glioma", [
            ("direct", f"Flow cytometry, IDH-wildtype (n = 40) [53]", kw("IDH_wildtype_glioma_n40")),
            ("direct", f"Single-cell atlas, GBM patients (n = {n_ag}) [42]", a_gbm),
            ("direct", f"GBmap atlas, all glioblastoma ({gbmap['n_donors_with_any_lymphoid']} donors)",
             [gbs["T"], gbs["NK"], gbs["B"]]),
            ("meth", "DNA methylation, EpiDISH (blood reference)", ep(lo)),
            ("meth", "DNA methylation, GIMiCC (glioma-specific)", gm("gbm")),
            ("rna", f"Deconvolution, frozen signature (mean of {n_lo})", m_lo),
            ("rna", f"Deconvolution, donor-level reference (mean of {n_lh})", m_lh)]),
        ("B   Lower-grade glioma  /  IDH-mutant glioma", [
            ("direct", f"Flow cytometry, IDH-mutant (n = 17) [53]", kw("IDH_mutant_glioma_n17")),
            ("direct", f"Single-cell atlas, LGG patients (n = {n_al}) [42]", a_lgg),
            ("meth", "DNA methylation, EpiDISH (blood reference)", ep(ll)),
            ("meth", "DNA methylation, GIMiCC (glioma-specific)", gm("lgg")),
            ("rna", f"Deconvolution, frozen signature (mean of {n_ll})", m_ll),
            ("rna", f"Deconvolution, donor-level reference (mean of {n_llh})", m_llh)]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6))
    fig.subplots_adjust(wspace=0.95, top=0.83, bottom=0.15)
    group_name = {"direct": "direct measurement", "meth": "DNA methylation", "rna": "bulk RNA"}
    for ax, (title, bars) in zip(axes, panels):
        y = np.arange(len(bars))
        T = np.array([b[2][0] for b in bars]); NK = np.array([b[2][1] for b in bars]); B = np.array([b[2][2] for b in bars])
        ax.barh(y, T, color=C_T, label="T cell")
        ax.barh(y, NK, left=T, color=C_NK, label="NK cell")
        ax.barh(y, B, left=T + NK, color=C_B, label="B cell")
        for i, v in enumerate(B):
            ax.text(1.01, i, f"B {v:.2f}", va="center", ha="left", fontsize=8, color=PS.RED, clip_on=False)
        ax.set_yticks(y)
        ax.set_yticklabels([b[1] for b in bars], fontsize=8.8)
        for i in range(1, len(bars)):                          # separate the three instrument groups
            if bars[i][0] != bars[i - 1][0]:
                ax.axhline(i - 0.5, color="#888", lw=0.8, ls="--")
        starts = {}
        for i, b in enumerate(bars):
            starts.setdefault(b[0], i)
        for g, i0 in starts.items():
            ax.text(-0.02, i0 - 0.45, group_name[g], transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                    fontsize=7.6, color="#555", style="italic")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Relative composition within {T, NK, B}", fontsize=9.5)
        ax.set_title(title, loc="left", fontsize=11.5, fontweight="bold", pad=10)
        ax.invert_yaxis()
    h, lg = axes[0].get_legend_handles_labels()
    fig.legend(h, lg, loc="lower center", ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Every direct measurement puts T cells far above B cells; bulk-RNA deconvolution does not",
                 fontsize=12.5, y=0.97)
    fig.text(0.5, 0.905, "Direct measurement: flow cytometry and single-cell counts. Methylation instruments disagree in "
             "LGG; direct measurement settles it.", ha="center", fontsize=9, color="#555")
    save(fig, "Figure_truth_instruments")


def fig_accuracy_factors():
    """What deconvolution can be trusted for, and which checks predict it, against real DNA truths.
    A: accuracy by compartment; B-C: the reference build changes accuracy, method by method;
    D: how well each ground-truth-free check predicts real accuracy."""
    idd, gsec, den = J("identifiability_diagnostics.json"), J("gimicc_secondary.json"), J("identifiability_denominators.json")
    yf, yh = J("absolute_purity_yardstick.json"), J("absolute_purity_yardstick_h5ad.json")
    lf, lh = J("absolute_purity_yardstick_lgg.json"), J("absolute_purity_yardstick_lgg_h5ad.json")
    ya, ag = J("yardstick_agreement.json"), J("agreement_selection_test.json")
    if not all((idd, yf, yh, lf, lh, ya, ag)):
        return
    U = idd["units"]
    gq3 = (gsec.get("Q3") or {}).get("units") or {}
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    fig.subplots_adjust(wspace=0.55, hspace=0.42, top=0.9)
    # A: accuracy by compartment
    ax = axs[0, 0]
    comps = [("Tumor", "Tumour content\n(vs ABSOLUTE)"), ("Leukocytes", "All leukocytes\n(vs methylation LF)"),
             ("Lymphoid", "Lymphoid total"), ("T_cell", "T cells"), ("B_cell", "B cells"), ("NK_cell", "NK cells")]
    y = np.arange(len(comps))
    D = den.get("units", {})
    for dy, c, col, mk in ((-0.12, "gbm", C_GBM, "o"), (0.12, "lgg", C_LGG, "s")):
        # lymphoid units: EpiDISH's shares of the immune compartment against each estimate's share of
        # its own leukocyte total (check 7, matched denominators); tumour and leukocytes are tissue
        # fractions on both sides already
        v = [(D.get(f"{k}|{c}") or U[f"{k}|{c}"])["truth_rho_median_methods"] for k, _ in comps]
        ax.scatter(v, y + dy, color=col, marker=mk, s=46, zorder=3,
                   label=f"{c.upper()}, DNA truth")
        g = [gq3.get(f"{k}|{c}", {}).get("gimicc_level3_truth_rho_median_methods", np.nan) if k not in ("Tumor",) else np.nan
             for k, _ in comps]
        ax.scatter(g, y + dy, facecolors="none", edgecolors=col, marker=mk, s=46, zorder=3,
                   label=f"{c.upper()}, GIMiCC truth")
    ax.axvline(0, color="#999", lw=0.8)
    ax.axvline(0.4, color="#bbb", lw=0.8, ls=":")
    ax.text(0.42, len(comps) - 0.62, "tracks (0.40)", fontsize=7.6, color="#777", va="center")
    ax.set_yticks(y); ax.set_yticklabels([lab for _, lab in comps], fontsize=9)
    ax.invert_yaxis(); ax.set_xlim(-0.6, 1.0)
    ax.set_xlabel("Spearman with DNA truth (median across methods)\nDNA truth: ABSOLUTE purity; methylation "
                  "leukocyte fraction; EpiDISH on matched denominators", fontsize=9)
    ax.set_title("A   What bulk RNA determines", loc="left", fontsize=11.5, fontweight="bold")
    ax.set_ylim(len(comps) - 0.35, -0.6)
    ax.legend(fontsize=7.4, frameon=False, loc="lower left", borderaxespad=0.3)
    # B, C: the reference build, method by method
    for ax, (fz, h5, name) in ((axs[0, 1], (yf, yh, "GBM")), (axs[1, 0], (lf, lh, "LGG"))):
        mf, mh = fz["methods"], h5["methods"]
        ms = [m for m in mf if m in mh and mf[m].get("spearman_vs_purity") is not None
              and mh[m].get("spearman_vs_purity") is not None]
        ms.sort(key=lambda m: mh[m]["spearman_vs_purity"])
        yy = np.arange(len(ms))
        a = np.array([mf[m]["spearman_vs_purity"] for m in ms]); b = np.array([mh[m]["spearman_vs_purity"] for m in ms])
        ax.hlines(yy, np.minimum(a, b), np.maximum(a, b), color="#bbb", lw=2)
        ax.scatter(a, yy, facecolors="white", edgecolors=PS.GREY, s=40, zorder=3, label="frozen signature (registered)")
        ax.scatter(b, yy, color=C_GBM if name == "GBM" else C_LGG, s=40, zorder=3, label="donor-level reference (raw counts)")
        ax.set_yticks(yy); ax.set_yticklabels([PS.display(m) for m in ms], fontsize=8.6)
        ax.set_xlim(-0.1, 1.0)
        ax.set_xlabel("Tumour content: Spearman with DNA purity", fontsize=9.5)
        ax.set_title(f"{'B' if name == 'GBM' else 'C'}   The reference build changes accuracy ({name})",
                     loc="left", fontsize=11.5, fontweight="bold")
        ax.legend(fontsize=7.6, frameon=False, loc="upper left")
    # D: does a ground-truth-free check predict real accuracy?
    ax = axs[1, 1]
    h1 = idd["H1_primary"]
    a2 = ya["arm_2_absolute_purity"]["all_methods"]
    p2 = ag["primary_panel"]["P2"]
    rows = [(f"Anatomic concordance (ACS)\nranks {a2['n_methods']} methods", a2["spearman"], a2["p"]),
            ("Agreement between methods\nranks methods (mean of 6 units)", p2["mean_rho"], p2["p_one_sided"]),
            (f"Agreement between methods\nranks {h1['D2']['n_units']} compartments", h1["D2"]["rho"], h1["D2"]["p_one_sided"]),
            (f"Loss-scale stability\nranks {h1['D1']['n_units']} compartments", h1["D1"]["rho"], h1["D1"]["p_one_sided"])]
    yy = np.arange(len(rows))
    ax.barh(yy, [r[1] for r in rows], color=[PS.GREY, PS.GREY, C_GBM, C_GBM], height=0.55)
    for i, r in enumerate(rows):
        ax.text(max(r[1], 0) + 0.02, i, f"rho {r[1]:.2f}, p {r[2]:.3g}", va="center", fontsize=8.6)
    ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
    ax.invert_yaxis(); ax.set_xlim(-0.1, 1.25); ax.axvline(0, color="#999", lw=0.8)
    ax.set_xlabel("Spearman between the check and real accuracy", fontsize=9.5)
    ax.set_title("D   Which ground-truth-free check predicts accuracy?", loc="left", fontsize=11.5, fontweight="bold")
    fig.suptitle("Factors that decide whether deconvolution can be trusted, measured against real DNA truths",
                 fontsize=13, y=0.975)
    save(fig, "Figure_accuracy_factors")


def fig_cptac_per_sample():
    """Per-sample accuracy against single-nucleus counts from the SAME tissue (CPTAC glioblastoma;
    prespecified/cptac_per_sample_truth.md). A: compartments recovered (C1); B: anatomy vs per-sample
    tumour accuracy (C2); C: the reference build (C3); D: lymphoid composition, nuclei vs methods (C4).

    NOT RENDERED (2026-10-06). It was written before any CPTAC result. Every reading it draws is
    INCONCLUSIVE: the primary DNA control failed (rho 0.01), and so did the whole-genome check of the
    nuclei truth (S2.1, rho 0.13). It is kept as written; fig_cptac_wgs() is the CPTAC figure."""
    d = J("cptac_per_sample_truth.json")
    lb = config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv"
    if not d or not lb.exists():
        return
    acs = dict(pd.read_csv(lb)[["method", "acs"]].itertuples(index=False))
    R, S1 = d["registered_truth"], d["s1_truth_nk_without_ncam1"]
    comps = [("Tumor", "Tumour"), ("Macrophage_Microglia", "Macrophage / microglia"), ("Oligodendrocyte", "Oligodendrocyte"),
             ("Endothelial", "Endothelial"), ("Leukocytes", "All leukocytes"), ("Lymphoid", "Lymphoid total"),
             ("T_cell", "T cells"), ("NK_cell", "NK cells"), ("B_cell", "B cells")]
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    fig.subplots_adjust(wspace=0.5, hspace=0.42, top=0.9)
    # A
    ax = axs[0, 0]
    y = np.arange(len(comps))
    for dy, ref, col, mk in ((-0.15, "frozen", PS.GREY, "o"), (0.15, "h5ad", C_GBM, "s")):
        if ref not in R:
            continue
        for i, (k, _) in enumerate(comps):
            vals = [v[k]["rho"] for v in R[ref]["methods"].values() if not v["_is_control"] and v[k].get("rho") is not None]
            ax.scatter(vals, [i + dy] * len(vals), color=col, s=10, alpha=0.45, zorder=2)
        med = [R[ref]["compartment_median_rho"].get(k) for k, _ in comps]
        ax.scatter([m if m is not None else np.nan for m in med], y + dy, color=col, marker=mk, s=70, zorder=3,
                   edgecolors="black", linewidths=0.6,
                   label=("frozen signature" if ref == "frozen" else "donor-level reference") + " (median; dots = methods)")
    ax.axvline(0, color="#999", lw=0.8); ax.axvline(0.4, color="#bbb", lw=0.8, ls=":")
    ax.set_yticks(y); ax.set_yticklabels([lab for _, lab in comps], fontsize=9); ax.invert_yaxis()
    ax.set_xlim(-1, 1)
    ax.set_xlabel(f"Spearman across tumours: bulk estimate vs single-nucleus share (n = {R['n_cases_scored']})", fontsize=9)
    ax.set_title("A   Which cell types the methods recover, per sample", loc="left", fontsize=11.5, fontweight="bold")
    ax.legend(fontsize=7.4, frameon=False, loc="lower left")
    # B
    ax = axs[0, 1]
    for ref, col, mk in (("frozen", PS.GREY, "o"), ("h5ad", C_GBM, "s")):
        if ref not in R or "C2_acs_vs_tumour_accuracy" not in R[ref]:
            continue
        ms = [m for m, v in R[ref]["methods"].items() if not v["_is_control"] and m in acs and v["Tumor"].get("rho") is not None]
        xs, ys = [acs[m] for m in ms], [R[ref]["methods"][m]["Tumor"]["rho"] for m in ms]
        c2 = R[ref]["C2_acs_vs_tumour_accuracy"]
        ax.scatter(xs, ys, color=col, marker=mk, s=40, zorder=3,
                   label=f"{'frozen' if ref == 'frozen' else 'donor-level'}: rho = {c2['spearman']:.2f} ({c2['n_methods']} methods)")
        for m, xv, yv in zip(ms, xs, ys):
            ax.annotate(PS.display(m), (xv, yv), fontsize=6.5, color="#555", xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("Anatomic concordance (ACS, Ivy GAP)", fontsize=9.5)
    ax.set_ylabel("Tumour share accuracy per sample (Spearman)", fontsize=9.5)
    ax.set_title("B   Does the anatomy-chosen method get each sample right?", loc="left", fontsize=11.5, fontweight="bold")
    ax.legend(fontsize=7.6, frameon=False, loc="lower left")
    # C
    ax = axs[1, 0]
    c3 = d["registered_truth"].get("C3_reference_effect_tumour_rho", {})
    ms = [m for m, v in c3.items() if v["frozen"] is not None and v["h5ad"] is not None and not m.startswith("control_")]
    ms.sort(key=lambda m: c3[m]["h5ad"])
    yy = np.arange(len(ms))
    a = np.array([c3[m]["frozen"] for m in ms]); b = np.array([c3[m]["h5ad"] for m in ms])
    ax.hlines(yy, np.minimum(a, b), np.maximum(a, b), color="#bbb", lw=2)
    ax.scatter(a, yy, facecolors="white", edgecolors=PS.GREY, s=40, zorder=3, label="frozen signature")
    ax.scatter(b, yy, color=C_GBM, s=40, zorder=3, label="donor-level reference")
    ax.set_yticks(yy); ax.set_yticklabels([PS.display(m) for m in ms], fontsize=8.6)
    ax.axvline(0, color="#999", lw=0.8); ax.set_xlim(-1, 1)
    ax.set_xlabel("Tumour share accuracy per sample (Spearman)", fontsize=9.5)
    ax.set_title("C   The reference build, per sample", loc="left", fontsize=11.5, fontweight="bold")
    ax.legend(fontsize=7.6, frameon=False, loc="upper left")
    # D
    ax = axs[1, 1]
    bars = []
    for lab, o in (("Single nuclei (registered rule)", R), ("Single nuclei (NK panel without NCAM1)", S1)):
        sh = o["C4_snrna_lymphoid"].get("pooled_within_lymphoid_share")
        if sh:
            bars.append((lab, [sh["T_cell"], sh["NK_cell"], sh["B_cell"]]))
    for ref, lab in (("frozen", "Deconvolution, frozen signature (mean)"), ("h5ad", "Deconvolution, donor-level (mean)")):
        if ref in R:
            ws = [v["_within_lymphoid_mean"] for v in R[ref]["methods"].values()
                  if not v["_is_control"] and all(k in v["_within_lymphoid_mean"] for k in ("T_cell", "NK_cell", "B_cell"))]
            if ws:
                bars.append((lab, [float(np.mean([w[k] for w in ws])) for k in ("T_cell", "NK_cell", "B_cell")]))
    yy = np.arange(len(bars))
    T = np.array([b_[1][0] for b_ in bars]); NK = np.array([b_[1][1] for b_ in bars]); B = np.array([b_[1][2] for b_ in bars])
    ax.barh(yy, T, color=C_T, label="T cell"); ax.barh(yy, NK, left=T, color=C_NK, label="NK cell")
    ax.barh(yy, B, left=T + NK, color=C_B, label="B cell")
    for i, v in enumerate(B):
        ax.text(1.01, i, f"B {v:.2f}", va="center", fontsize=8, color=PS.RED, clip_on=False)
    ax.set_yticks(yy); ax.set_yticklabels([b_[0] for b_ in bars], fontsize=8.8); ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel("Relative composition within {T, NK, B}", fontsize=9.5)
    ax.set_title("D   Lymphoid composition: nuclei vs methods", loc="left", fontsize=11.5, fontweight="bold")
    ax.legend(fontsize=7.6, frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.32))
    fig.suptitle("Per-sample truth from the same tissue: CPTAC glioblastoma, bulk RNA vs single-nucleus counts",
                 fontsize=13, y=0.975)
    save(fig, "Figure_cptac_per_sample")


def fig_cptac_wgs():
    """CPTAC glioblastoma against whole-genome DNA purity (prespecified/cptac_wgs_purity_secondary.md, S2;
    SECONDARY, registered after the primary DNA control failed). Top row: which per-sample truth agrees with
    DNA (S2.1, S2.2, S2.6). Bottom row: the bulk methods against DNA (S2.3), the TCGA ranking (S2.4), anatomy (S2.5).

    Replaces fig_cptac_per_sample() above, which was written before any result and is not rendered: every
    reading it draws is INCONCLUSIVE (the primary control failed, and so did S2.1)."""
    d = J("cptac_wgs_purity.json")
    tab_f = config.RESULTS_DIR / "cptac" / "wgs_purity_per_sample.csv"
    lb = config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv"
    if not d or not tab_f.exists() or not lb.exists():
        return
    t = pd.read_csv(tab_f, index_col=0)
    acs = dict(pd.read_csv(lb)[["method", "acs"]].itertuples(index=False))
    a = d["all_cases"]
    REF = (("frozen", PS.GREY, "o", "frozen signature"), ("h5ad", C_GBM, "s", "donor-level reference"))
    fig, axs = plt.subplots(2, 3, figsize=(17, 11.2))
    fig.subplots_adjust(wspace=0.36, hspace=0.55, top=0.89)

    def stat(ax, txt, loc="upper left"):
        y, va = (0.97, "top") if loc.startswith("upper") else (0.03, "bottom")
        x, ha = (0.03, "left") if loc.endswith("left") else (0.97, "right")
        ax.text(x, y, txt, transform=ax.transAxes, va=va, ha=ha, fontsize=8.4,
                bbox=dict(fc="white", ec="#ddd", alpha=0.9, boxstyle="round,pad=0.3"))

    def unit(ax, ylab):
        ax.plot([0, 1], [0, 1], ls="--", color="#aaa", lw=0.9, zorder=1)
        ax.set_xlim(0.3, 1.0); ax.set_ylim(0, 1.04)
        ax.set_xlabel("Whole-genome DNA purity (AscatNGS)", fontsize=9.5); ax.set_ylabel(ylab, fontsize=9.5)

    def p_(v):
        return "p < 0.001" if v["perm_p"] < 0.001 else f"p = {v['perm_p']:.2f}"

    # A  nuclei
    ax = axs[0, 0]
    s21 = a["S2_1_nuclei_truth_vs_wgs"]["registered"]["all_scored"]
    sc = t.dropna(subset=["snrna_tumour_registered"])
    sc = sc.loc[a["S2_1_nuclei_truth_vs_wgs"]["registered"]["scored_cases"]]
    for pooled, face, lab in ((False, PS.INK, "one tumour piece"), (True, "white", "pooled pieces (nuclei from one)")):
        q = sc[sc["pooled_pieces"] == pooled]
        ax.scatter(q["wgs_purity"], q["snrna_tumour_registered"], facecolors=face, edgecolors=PS.INK, s=42, zorder=3, label=lab)
    unit(ax, "Tumour share of nuclei (single-nucleus RNA)")
    stat(ax, f"ρ = {s21['rho']:+.2f}, {p_(s21)}, n = {s21['n']}\nfails the 0.40 bar", loc="lower left")
    ax.legend(fontsize=7.6, frameon=False, loc="lower right")
    ax.set_title("A   Single-nucleus counts do not track DNA", loc="left", fontsize=11.5, fontweight="bold")
    # B  methylation
    ax = axs[0, 1]
    s22 = a["S2_2_methylation_instrument"]
    low = t["methylation_library_coverage"] < 0.90
    ax.scatter(t.loc[~low, "wgs_purity"], t.loc[~low, "gimicc_769cpg"], facecolors="white", edgecolors=PS.GREY, s=36,
               zorder=2, label=f"all 18 samples, 769 CpGs: ρ = {s22['a_complete_case_769_cpgs']['rho']:+.2f}")
    ax.scatter(t.loc[low, "wgs_purity"], t.loc[low, "gimicc_769cpg"], marker="x", color=PS.RED, s=40, zorder=3,
               label="  of which, the 5 samples missing a third or more of the CpGs")
    ok = t["gimicc_cov90"].notna()
    ax.scatter(t.loc[ok, "wgs_purity"], t.loc[ok, "gimicc_cov90"], color=C_GBM, s=42, zorder=4,
               label=f"13 complete samples, 3,775 CpGs: ρ = {s22['b_coverage_ge_90pct']['rho']:+.2f}")
    unit(ax, "Tumour fraction from methylation (GIMiCC)")
    ax.legend(fontsize=7.4, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    ax.set_title("B   Methylation tracks DNA once inputs are complete", loc="left", fontsize=11.5, fontweight="bold")
    # C  pathology
    ax = axs[0, 2]
    s26 = d["S2_6_pathology_descriptive"]
    ax.scatter(t["wgs_purity"], t["pathology_tumor_nuclei"], color=PS.ORANGE, s=42, zorder=3)
    unit(ax, "Tumour nuclei on the slide (pathologist, mean)")
    v = s26["vs_wgs_purity"]
    stat(ax, f"ρ = {v['rho']:+.2f}, {p_(v)}, n = {v['n']}\nrange {s26['range'][0]:.0%}-{s26['range'][1]:.0%} (CPTAC selected\n"
             "samples on tumour content); descriptive", loc="lower right")
    ax.set_title("C   Pathology, for reference", loc="left", fontsize=11.5, fontweight="bold")
    # D  methods vs DNA
    ax = axs[1, 0]
    ms = sorted(a["frozen"]["methods"], key=lambda m: a["frozen"]["methods"][m]["rho"])
    ms = [m for m in a["h5ad"]["methods"] if m not in ms] + ms          # donor-level-only methods at the bottom
    yy = {m: i for i, m in enumerate(ms)}
    for ref, col, mk, lab in REF:
        M = a[ref]["methods"]
        s23 = a[ref]["S2_3"]
        dy = -0.14 if ref == "frozen" else 0.14
        ax.scatter([M[m]["rho"] for m in M], [yy[m] + dy for m in M], color=col, marker=mk, s=34, zorder=3,
                   label=f"{lab}: median {s23['median_rho']:.2f} (TCGA-GBM, same methods, {s23['tcga_gbm_median_same_methods']:.2f})")
        ax.axvline(s23["median_rho"], color=col, lw=1.2, alpha=0.7)
    ax.axvline(0.40, color="#bbb", lw=0.9, ls=":"); ax.axvline(0, color="#999", lw=0.8)
    ax.set_yticks(range(len(ms))); ax.set_yticklabels([PS.display(m) for m in ms], fontsize=8.2)
    ax.set_xlim(-0.6, 1.0)
    ax.set_xlabel("Bulk tumour estimate vs DNA purity (Spearman, n = 18)", fontsize=9.5)
    ax.legend(fontsize=7.2, frameon=False, loc="lower center", bbox_to_anchor=(0.45, -0.27))
    ax.set_title("D   Bulk methods recover tumour content", loc="left", fontsize=11.5, fontweight="bold")
    # E  ranking transfer
    ax = axs[1, 1]
    impl_t = {ref: {m: v.get("implementation") for m, v in J(y)["methods"].items()}
              for ref, y in (("frozen", "absolute_purity_yardstick.json"), ("h5ad", "absolute_purity_yardstick_h5ad.json"))}
    for ref, col, mk, lab in REF:
        M = a[ref]["methods"]
        ms_ = [m for m in M if M[m]["tcga_gbm_rho_vs_absolute"] is not None]
        tr = a[ref]["S2_4_ranking_transfer"]
        ax.scatter([M[m]["tcga_gbm_rho_vs_absolute"] for m in ms_], [M[m]["rho"] for m in ms_], color=col, marker=mk,
                   s=38, zorder=3, label=f"{lab}: ρ = {tr['spearman']:+.2f} ({tr['n_methods']} methods), {tr['reading']}")
        groups: dict = {}
        for m in ms_:                                  # degenerate methods share a point: one label, "A = B"
            dag = "\u2020" if impl_t[ref].get(m) != M[m]["implementation"] else ""
            groups.setdefault((M[m]["tcga_gbm_rho_vs_absolute"], M[m]["rho"]), []).append(PS.display(m) + dag)
        labels: list = []                              # near-coincident points (within 0.03) share one label
        for (xv, yv), names in sorted(groups.items()):
            near = [L for L in labels if abs(L[0] - xv) < 0.03 and abs(L[1] - yv) < 0.03]
            if near:
                near[0][2].append(" = ".join(names))
            else:
                labels.append([xv, yv, [" = ".join(names)]])
        for xv, yv, names in labels:
            if ref == "frozen":
                ax.annotate(", ".join(names), (xv, yv), fontsize=6.4, color="#555", xytext=(3, 2), textcoords="offset points")
            elif abs(yv - xv) > 0.3:
                ax.annotate(", ".join(names) + " (donor-level)", (xv, yv), fontsize=6.4, color=C_GBM,
                            xytext=(4, -3), textcoords="offset points")
    ax.plot([-0.4, 1], [-0.4, 1], ls="--", color="#aaa", lw=0.9, zorder=1)
    ax.set_xlim(-0.4, 1); ax.set_ylim(-0.4, 1)
    ax.set_xlabel("Accuracy in TCGA-GBM (vs ABSOLUTE, n = 154)", fontsize=9.5)
    ax.set_ylabel("Accuracy in CPTAC (vs AscatNGS, n = 18)", fontsize=9.5)
    ax.legend(fontsize=7.2, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.27))
    ax.set_title("E   The best methods stay best (frozen signature only)", loc="left", fontsize=11.5, fontweight="bold")
    if any(impl_t[r].get(m) != a[r]["methods"][m]["implementation"] for r in ("frozen", "h5ad") for m in a[r]["methods"]
           if a[r]["methods"][m]["tcga_gbm_rho_vs_absolute"] is not None):
        ax.text(0.5, -0.43, "\u2020 a different implementation ran in the two cohorts", transform=ax.transAxes,
                ha="center", fontsize=7.2, color="#555")
    # F  anatomy vs accuracy
    ax = axs[1, 2]
    for ref, col, mk, lab in REF:
        s25 = a[ref]["S2_5_acs_vs_wgs_accuracy"]
        ms_ = s25["methods"]
        ax.scatter([acs[m] for m in ms_], [a[ref]["methods"][m]["rho"] for m in ms_], color=col, marker=mk, s=38,
                   zorder=3, label=f"{lab}: ρ = {s25['spearman']:+.2f} ({s25['n_methods']} methods)")
    ax.set_xlabel("Anatomic concordance (ACS, Ivy GAP)", fontsize=9.5)
    ax.set_ylabel("Accuracy in CPTAC (vs AscatNGS, n = 18)", fontsize=9.5)
    ax.legend(fontsize=7.2, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.27))
    ax.set_title("F   Anatomy does not pick the accurate method", loc="left", fontsize=11.5, fontweight="bold")
    fig.suptitle("CPTAC glioblastoma: the per-sample truths and the bulk methods against whole-genome DNA purity "
                 "(secondary analysis S2)", fontsize=13, y=0.975)
    save(fig, "Figure_cptac_wgs")


def fig_acs_constraints():
    """What the anatomic score separates: each method's share of tumours satisfying each registered constraint
    (results/anatomic/acs_per_constraint.csv), ordered by ACS, beside its tumour-content accuracy against DNA
    (TCGA-GBM donor-level, the reference the ACS arm uses; CPTAC whole-genome). Descriptive; selects nothing."""
    pc_f, lb_f = config.RESULTS_DIR / "anatomic" / "acs_per_constraint.csv", config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv"
    yh, cw = J("absolute_purity_yardstick_h5ad.json").get("methods", {}), J("cptac_wgs_purity.json").get("all_cases", {})
    if not (pc_f.exists() and lb_f.exists()):
        return
    pc, lb = pd.read_csv(pc_f), pd.read_csv(lb_f).set_index("method")
    order = list(lb.sort_values(["is_control", "acs"], ascending=[True, False]).index)
    order = [m for m in order if lb.loc[m, "comparable"]] + [m for m in order if not lb.loc[m, "comparable"]]
    cons = sorted(pc["constraint"].unique(), key=lambda c: int(c[1:]))
    desc = pc.drop_duplicates("constraint").set_index("constraint")["description"]
    frac = pc.pivot(index="method", columns="constraint", values="fraction_satisfied").reindex(index=order, columns=cons)
    nsat = pc.pivot(index="method", columns="constraint", values="n_satisfied").reindex(index=order, columns=cons)
    nev = pc.pivot(index="method", columns="constraint", values="n_tumors_evaluable").reindex(index=order, columns=cons)
    fig = plt.figure(figsize=(15.5, 8.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[7, 1.1, 2.6], wspace=0.08)
    ax = fig.add_subplot(gs[0])
    ax.imshow(frac.to_numpy(float), cmap="Blues", vmin=0, vmax=1, aspect="auto")      # darker = larger share
    for i, m in enumerate(order):
        for j, c in enumerate(cons):
            n = nev.loc[m, c]
            txt = "--" if (pd.isna(n) or n == 0) else f"{int(nsat.loc[m, c])}/{int(n)}"
            v = frac.loc[m, c]
            ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                    color="white" if (pd.notna(v) and v > 0.6) else PS.INK)
    ax.set_xticks(range(len(cons)))
    ax.set_xticklabels([f"{c}\n" + desc[c].replace("Macrophage_Microglia", "Myeloid").replace(": ", ":\n")
                        for c in cons], fontsize=7.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([PS.display(m) + ("  (control)" if lb.loc[m, "is_control"]
                                         and not PS.display(m).lower().startswith("control") else "")
                        + ("  (not comparable)" if not lb.loc[m, "comparable"] else "") for m in order], fontsize=8.6)
    n_real = int(((~lb["is_control"]) & lb["comparable"]).sum())
    ax.axhline(n_real - 0.5, color=PS.RED, lw=1.2)
    ax.set_title("A   Share of tumours satisfying each registered constraint (Ivy GAP)", loc="left",
                 fontsize=11.5, fontweight="bold")
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    ax2.barh(range(len(order)), [lb.loc[m, "acs"] for m in order], color=[PS.GREY if lb.loc[m, "is_control"] else C_GBM
                                                                         for m in order], height=0.6)
    for i, m in enumerate(order):
        ax2.text(lb.loc[m, "acs"] + 0.02, i, f"{lb.loc[m, 'acs']:.2f}", va="center", fontsize=7.6)
    ax2.set_xlim(0, 1.25); ax2.set_xticks([0, 0.5, 1]); ax2.tick_params(labelleft=False)
    ax2.set_title("B   ACS", loc="left", fontsize=11.5, fontweight="bold")
    ax3 = fig.add_subplot(gs[2], sharey=ax)
    for i, m in enumerate(order):
        t = (yh.get(m) or {}).get("spearman_vs_purity")
        c = ((cw.get("h5ad") or {}).get("methods", {}).get(m) or {}).get("rho")
        if t is not None:
            ax3.scatter(t, i, color=C_GBM, marker="s", s=30, zorder=3)
        if c is not None:
            ax3.scatter(c, i, facecolors="white", edgecolors=PS.ORANGE, s=34, zorder=3)
    ax3.scatter([], [], color=C_GBM, marker="s", label="TCGA-GBM vs ABSOLUTE (n = 154)")
    ax3.scatter([], [], facecolors="white", edgecolors=PS.ORANGE, label="CPTAC vs whole-genome (n = 18)")
    ax3.axvline(0, color="#999", lw=0.8); ax3.set_xlim(-0.6, 1.0); ax3.tick_params(labelleft=False)
    ax3.set_xlabel("Tumour accuracy (Spearman),\ndonor-level reference", fontsize=8.6)
    ax3.legend(fontsize=7, frameon=False, loc="lower left")
    ax3.set_title("C   Accuracy against DNA", loc="left", fontsize=11.5, fontweight="bold")
    fig.suptitle("What the anatomic score separates: the controls fail across constraints; working methods differ "
                 "mostly on C6 and C7", fontsize=12.5, y=0.99)
    save(fig, "Figure_acs_constraints")


def main() -> int:
    print("rendering figures from artefacts...")
    fig_detects_not_ranks()
    fig_lymphoid()
    fig_model_fit()
    fig_recovery()
    fig_bisque_anchoring()
    fig_lymphoid_mechanisms()
    fig_acs_vs_auc()
    fig_identifiability()
    fig_truth_instruments()
    fig_accuracy_factors()
    fig_cptac_wgs()                 # fig_cptac_per_sample() is not rendered: its readings are INCONCLUSIVE
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
