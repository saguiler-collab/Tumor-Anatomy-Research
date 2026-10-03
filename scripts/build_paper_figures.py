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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
