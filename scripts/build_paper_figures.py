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


def main() -> int:
    print("rendering figures from artefacts...")
    fig_detects_not_ranks()
    fig_lymphoid()
    fig_model_fit()
    fig_recovery()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
