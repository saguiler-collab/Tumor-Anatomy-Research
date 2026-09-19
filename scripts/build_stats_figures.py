"""Statistical figures: the scatters the Spearman coefficients were computed for, and
paired comparisons with error bars and annotated p-values.

WHY THESE EXIST. The project reports Spearman correlations between each method's estimate and
DNA-measured purity across 154 (GBM) and 510 (LGG) samples, and never plotted them. A
correlation coefficient without its scatter asks the reader to take the relationship on trust,
and it hides the two things a scatter shows immediately: whether the relationship is driven by
a handful of points, and whether the estimates occupy the same RANGE as the truth.

The paired panels follow the convention of a results figure in an experimental paper: bars are
mean +/- SD, the test is named, the p-value is printed, and n is stated.

    python scripts/build_stats_figures.py
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402
from scipy import stats                                            # noqa: E402

from ivygap import config                                          # noqa: E402
from ivygap.paper_style import apply as apply_style                # noqa: E402
from ivygap import paper_style as PS                               # noqa: E402

apply_style()

# TRACKED, unlike results/ which is gitignored. Figures are deliverables: they are
# reviewed, marked up and compared between versions, so they have to be openable by
# someone who did not run the pipeline.
FIG = config.PROJECT_ROOT / "docs" / "figures"
C_PT, C_FIT, C_ID = PS.BLUE, PS.RED, PS.GREY
C_T, C_B = PS.BLUE, PS.RED


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    # derive the path from FIG rather than hardcoding it: this line said
    # "results/figures/" for a while after the output moved to docs/figures/,
    # which is how someone ends up editing a stale copy.
    print(f"  wrote {FIG.relative_to(config.PROJECT_ROOT)}/{name}.png and .pdf")


def p_text(p: float) -> str:
    """Journals print small p-values as a bound, not as a spuriously precise number."""
    if p < 1e-4:
        return "p < 0.0001"
    if p < 0.001:
        return f"p = {p:.1e}".replace("e-0", "e-")
    return f"p = {p:.4f}".rstrip("0").rstrip(".")


def fig_purity_scatter(tag: str, label: str, n_panels: int = 6):
    """Estimated tumour fraction vs DNA purity, per method. THE scatter behind the Spearman."""
    p = config.RESULTS_DIR / f"absolute_purity_per_sample{tag}.csv"
    if not p.exists():
        return
    d = pd.read_csv(p, index_col=0)
    tcol = next(c for c in ("absolute_purity", "purity") if c in d.columns)
    truth = d[tcol].to_numpy(float)
    methods = [c for c in d.columns if c != tcol]
    # rank by |rho| so the panel shows the best and the worst rather than an arbitrary slice
    rho = {}
    for m in methods:
        e = d[m].to_numpy(float)
        ok = np.isfinite(e) & np.isfinite(truth)
        if ok.sum() >= 50 and np.std(e[ok]) > 0:
            rho[m] = stats.spearmanr(e[ok], truth[ok]).statistic
    if not rho:
        return
    ordered = sorted(rho, key=lambda m: -rho[m])
    show = ordered[:n_panels - 2] + ordered[-2:] if len(ordered) > n_panels else ordered
    ncol = 3
    nrow = int(np.ceil(len(show) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.7 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, m in zip(axes, show):
        e = d[m].to_numpy(float)
        ok = np.isfinite(e) & np.isfinite(truth)
        x, y = truth[ok], e[ok]
        r = stats.spearmanr(x, y)
        ax.plot([0, 1], [0, 1], "--", color=C_ID, lw=1.1, zorder=1)
        ax.scatter(x, y, s=13, alpha=.45, color=C_PT, edgecolors="none", zorder=2)
        if np.std(x) > 0:
            b, a = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, a + b * xs, color=C_FIT, lw=2, zorder=3)
        ax.set_title(m, fontsize=10.5, fontweight="bold", loc="left")
        ax.text(.03, .97, f"ρ = {r.statistic:+.3f}\n{p_text(r.pvalue)}\nn = {ok.sum()}",
                transform=ax.transAxes, va="top", fontsize=8.6,
                bbox=dict(fc="white", ec="#ddd", alpha=.85, boxstyle="round,pad=0.3"))
        ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
        ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
    for ax in axes[len(show):]:
        ax.axis("off")
    fig.supxlabel("DNA-measured tumour purity (ABSOLUTE)", fontsize=11)
    fig.supylabel("estimated tumour fraction", fontsize=11)
    # The accurate claim is about SLOPE, not position. Several methods cross the identity
    # line rather than sitting below it: they over-call at low purity and under-call at high,
    # which is range compression and is exactly what the recovery statistic measures. An
    # earlier draft said "sit far below the identity line", which is true of `bayesian` and
    # false of `svr`, `cibersortx` and `dwls`.
    fig.suptitle(f"{label}: every fitted slope is shallower than identity — "
                 f"estimates compress the true range",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.text(.5, -0.015, "dashed grey = identity (a perfect estimate); red = least-squares "
                         "fit; ρ = Spearman. A slope below 1 means the method under-calls "
                         "high-purity samples and over-calls low-purity ones.",
             ha="center", fontsize=8.6, color="#666")
    fig.tight_layout()
    save(fig, f"Figure_purity_scatter{tag or '_gbm'}")


def fig_paired_lymphoid():
    """T vs B per method with mean +/- SD and a paired test, against methylation truth."""
    def load(tag):
        est = config.RESULTS_DIR / f"estimates_full{tag}.csv"
        meth = config.RESULTS_DIR / f"methylation_celltypes{tag}.csv"
        if not (est.exists() and meth.exists()):
            return None
        k4s = lambda s: "-".join(str(s).split("-")[:3] + [str(s).split("-")[3][:2]]) \
            if len(str(s).split("-")) > 3 else str(s)
        M = pd.read_csv(meth, index_col=0)
        ros = pd.DataFrame({"T_cell": M[["CD4T", "CD8T"]].sum(axis=1),
                            "B_cell": M["B"], "NK_cell": M["NK"]})
        ros = ros.div(ros.sum(axis=1).replace(0, np.nan), axis=0)
        ros.index = [k4s(i) for i in ros.index]
        ros = ros[~ros.index.duplicated()].dropna()
        E = pd.read_csv(est)
        sc = next(c for c in ("sample", "sample_id") if c in E.columns)
        E["k"] = E[sc].map(k4s)
        return ros, E
    panels = []
    for tag, lab in (("", "GBM"), ("_lgg", "LGG")):
        got = load(tag)
        if got:
            panels.append((got[0], got[1], lab))
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(7.6 * len(panels), 5.6))
    axes = np.atleast_1d(axes)
    C = ["T_cell", "B_cell", "NK_cell"]
    for panel_letter, (ax, (ros, E, lab)) in zip("ABCD", zip(axes, panels)):
        names, Tm, Ts, Bm, Bs, ps = [], [], [], [], [], []
        n_matched = 0
        for m, g in E.groupby("method"):
            g = g.drop_duplicates("k").set_index("k")
            idx = [k for k in g.index if k in ros.index]
            if len(idx) < 50:
                continue
            sub = g.loc[idx, C]
            sub = sub.div(sub.sum(axis=1).replace(0, np.nan), axis=0).dropna()
            if len(sub) < 20:
                continue
            names.append(m)
            n_matched = max(n_matched, len(sub))
            Tm.append(sub.T_cell.mean()); Ts.append(sub.T_cell.std())
            Bm.append(sub.B_cell.mean()); Bs.append(sub.B_cell.std())
            # paired within sample: the honest test for "does this method put B above T"
            ps.append(stats.wilcoxon(sub.T_cell, sub.B_cell).pvalue
                      if (sub.T_cell != sub.B_cell).any() else np.nan)
        # truth bar last
        names.append("DNA\nmethylation")
        Tm.append(ros.T_cell.mean()); Ts.append(ros.T_cell.std())
        Bm.append(ros.B_cell.mean()); Bs.append(ros.B_cell.std())
        ps.append(stats.wilcoxon(ros.T_cell, ros.B_cell).pvalue)
        x = np.arange(len(names))
        ax.bar(x - .19, Tm, .36, yerr=Ts, color=C_T, capsize=2.5, label="T cell",
               error_kw=dict(lw=.9, ecolor="#555"))
        ax.bar(x + .19, Bm, .36, yerr=Bs, color=C_B, capsize=2.5, label="B cell",
               error_kw=dict(lw=.9, ecolor="#555"))
        for i, pv in enumerate(ps):
            if np.isnan(pv):
                continue
            # stagger alternate labels: adjacent bars of similar height put their p-values at
            # the same y and overlap, which made two of them unreadable in the first draft
            top = max(Tm[i] + Ts[i], Bm[i] + Bs[i]) + .04 + (.10 if i % 2 else 0)
            ax.plot([i - .19, i - .19, i + .19, i + .19],
                    [top, top + .025, top + .025, top], lw=.9, color="#333")
            ax.text(i, top + .035, p_text(pv), ha="center", fontsize=7.0)
        ax.set_xticks(x)
        ax.set_xticklabels([n.replace("_", " ") for n in names], rotation=45, ha="right",
                           fontsize=9)
        ax.get_xticklabels()[-1].set_fontweight("bold")
        ax.set_ylim(0, 1.52)
        ax.set_ylabel("relative composition within {T, NK, B}", fontsize=10)
        # TWO DIFFERENT n's, and conflating them misleads. `len(ros)` is how many samples
        # have methylation truth; `n_matched` is how many of those also have an estimate from
        # this cohort's deconvolution run, which is what every method bar is computed on.
        # In GBM those are 154 and 56 -- most TCGA-GBM methylation is HM27, not HM450.
        ax.set_title(f"{panel_letter}   {lab}   "
                     f"(truth n = {len(ros)}; per-method n = {n_matched})",
                     fontsize=11.5, fontweight="bold", loc="left")
        ax.axvline(len(names) - 1.5, color="#999", ls="--", lw=.9)
    axes[0].legend(frameon=False, fontsize=10, loc="upper left")
    fig.suptitle("Methods place B cells above T cells; DNA methylation places T above B",
                 fontsize=13.5, fontweight="bold", y=1.0)
    fig.text(.5, -0.06, "Bars are mean ± SD. p-values are two-sided Wilcoxon signed-rank, "
                        "paired within sample.", ha="center", fontsize=9, color="#555")
    fig.tight_layout()
    save(fig, "Figure_lymphoid_paired")


def main() -> int:
    print("rendering statistical figures...")
    fig_purity_scatter("", "Glioblastoma")
    fig_purity_scatter("_lgg", "Lower-grade glioma")
    fig_paired_lymphoid()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
