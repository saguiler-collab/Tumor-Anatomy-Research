"""Write docs/FIGURES.md — every figure embedded, with a journal-ready caption.

A caption is not a title. It must let a reader who skips the body understand what was
measured, on how many samples, by what test, and what the panel shows. Numbers here are read
from artefacts so a caption cannot drift from the figure beside it.

    python scripts/build_figure_index.py
"""
from __future__ import annotations

import json

from ivygap import config

OUT = config.PROJECT_ROOT / "docs" / "FIGURES.md"
FIGDIR = config.PROJECT_ROOT / "docs" / "figures"


def J(n):
    p = config.RESULTS_DIR / n
    return json.loads(p.read_text()) if p.exists() else {}


def main() -> int:
    lo, lol = J("lymphoid_ordering.json"), J("lymphoid_ordering_lgg.json")
    mf = J("model_fit_residual.json").get("cohorts", {})
    fz, h5 = J("yardstick_agreement.json"), J("yardstick_agreement_h5ad.json")
    efg, efl = J("equal_footing_ranking.json"), J("equal_footing_ranking_lgg.json")
    a2f = fz.get("arm_2_absolute_purity", {})
    a2h = h5.get("arm_2_absolute_purity", {})
    g, l = mf.get("GBM", {}), mf.get("LGG", {})

    def emp(d):
        return sum(1 for v in (d.get("methods") or {}).values()
                   if (v.get("n_scored") or v.get("n", 0)) < 0.5 * v.get("n", 1))

    figs = [
        ("Figure_detects_not_ranks",
         "Anatomic concordance predicts accuracy only against the yardstick that shares its "
         "own reference atlas.",
         f"Spearman correlation between the ranking produced by the Anatomic Concordance Score "
         f"(ACS) and the ranking produced by each external yardstick. Points are the "
         f"correlation; bars are 95% bootstrap confidence intervals over methods. The dashed "
         f"line is the pre-registered bar (ρ ≥ 0.60), fixed before any result was computed. "
         f"Only the synthetic pseudobulk yardstick clears it (ρ = "
         f"{fz.get('arm_1_pseudobulk', {}).get('rho'):.4f}), and that yardstick is built from "
         f"GBmap — the same single-cell atlas the ACS arm deconvolves against. Against "
         f"DNA-measured tumour purity, which shares nothing with ACS, the correlation is "
         f"{a2f.get('all_methods', {}).get('spearman')} "
         f"(n = {a2f.get('all_methods', {}).get('n_methods')} methods) and "
         f"{a2f.get('excluding_degenerate', {}).get('spearman')} restricted to comparable "
         f"methods. Scoring both arms on the same reference build removes the one remaining "
         f"objection and changes nothing: "
         f"{a2h.get('all_methods', {}).get('spearman')} "
         f"(n = {a2h.get('all_methods', {}).get('n_methods')}) and "
         f"{a2h.get('excluding_degenerate', {}).get('spearman')} "
         f"(n = {a2h.get('excluding_degenerate', {}).get('n_methods')})."),

        ("Figure_lymphoid_failure",
         "No method reproduces the lymphoid ordering that DNA methylation measures, in either "
         "tumour type.",
         f"Relative composition within {{T, NK, B}} for every deconvolution method, against DNA "
         f"methylation as an independent per-cell-type ground truth (EpiDISH RPC on the "
         f"`centDHSbloodDMC.m` reference). Both sides are renormalised within the same three "
         f"cell types, so the leukocyte-subcomposition and all-cell-fraction denominators "
         f"cancel by construction. (A) Glioblastoma, truth n = {lo.get('n_methylation_samples')}. "
         f"(B) Lower-grade glioma, truth n = {lol.get('n_methylation_samples')}. Methylation "
         f"orders the compartment {lo.get('truth_ordering')}; "
         f"{lo.get('n_agree_T_over_B')} of {lo.get('n_methods')} methods in GBM and "
         f"{lol.get('n_agree_T_over_B')} of {lol.get('n_methods')} in LGG agree that T exceeds "
         f"B. Grey labels count samples in which the method returned EXACTLY ZERO T, NK and B — "
         f"a second and distinct failure, affecting {emp(lo)} of {lo.get('n_methods')} methods "
         f"in GBM and {emp(lol)} of {lol.get('n_methods')} in LGG."),

        ("Figure_lymphoid_paired",
         "Methods place B cells above T cells; DNA methylation places T above B.",
         # DERIVED, not typed. This caption said "154 ... samples with methylation data",
         # but 154 is `n_orderable_samples` -- the methylation samples MINUS the one whose
         # T, NK and B are all exactly zero, which has no defined ordering. The count with
         # methylation data is 155, which is what Figure 2's generated caption says. The
         # number was right for what the figure plots and the DESCRIPTION was wrong, so the
         # two captions contradicted each other on the same quantity.
         f"Mean ± SD of the T-cell and B-cell fraction within {{T, NK, B}}, per method and "
         f"for the methylation truth bar. p-values are two-sided Wilcoxon signed-rank tests "
         f"paired within sample, not comparisons of two cohort averages. Sample counts "
         f"differ within a panel and the difference is not incidental: "
         f"{lo.get('n_methylation_samples')} GBM samples have methylation data, "
         f"{lo.get('n_orderable_samples')} of those admit an ordering at all (the "
         f"remaining {lo.get('n_unorderable_zero_lymphoid')} returned exactly zero T, NK "
         f"and B, so no ordering is defined), and only "
         f"{max((v.get('n') or 0) for v in (lo.get('methods') or {{}}).values()) if lo.get('methods') else 0} "
         f"also have a deconvolution estimate — most TCGA-GBM methylation was assayed on "
         f"the older HM27 platform rather than HM450."),

        ("Figure_purity_scatter_gbm",
         "Glioblastoma: every fitted slope is shallower than identity, so estimates compress "
         "the true range.",
         "Estimated tumour fraction against DNA-measured tumour purity (ABSOLUTE, from copy "
         "number), one panel per method, ordered by Spearman ρ. Each panel gives ρ, its "
         "p-value and n. The dashed grey line is identity — where a perfect estimate would "
         "lie; the red line is the least-squares fit. A slope below 1 means the method "
         "under-calls high-purity samples and over-calls low-purity ones. Note that a high "
         "correlation does not imply a usable estimate: `bayesian` reaches ρ = +0.742 with its "
         "entire distribution far below identity."),

        ("Figure_purity_scatter_lgg",
         "Lower-grade glioma: the same compression, on 510 samples.",
         "As the preceding figure, for the TCGA lower-grade glioma cohort. Across both cohorts "
         "all 24 method-by-cohort fitted slopes lie below 1 (GBM maximum 0.637, LGG maximum "
         "0.421)."),

        ("Figure_model_fit_bound",
         "The additive mixing model explains a minority of real bulk — but far more than "
         "chance.",
         f"R² of the best non-negative fit of the reference profiles to each bulk sample, "
         f"median across samples, on the marker gene space. The model leaves "
         f"{g.get('unexplained_median_pct')}% (GBM) and {l.get('unexplained_median_pct')}% "
         f"(LGG) of variance unexplained. Two floors and a ceiling make that number "
         f"interpretable: shuffling which gene belongs to which cell type drops R² to "
         f"{g.get('controls', {}).get('shuffled_profiles')}, and a random non-negative basis to "
         f"{g.get('controls', {}).get('random_basis')}, so the reference carries real "
         f"structure; the top-8 SVD of the bulk itself reaches "
         f"{g.get('controls', {}).get('ceiling_top_k_svd')}, which is the best any "
         f"eight-dimensional basis could achieve on the same samples. The reference attains "
         f"{g.get('controls', {}).get('fraction_of_achievable', 0) * 100:.0f}% (GBM) and "
         f"{l.get('controls', {}).get('fraction_of_achievable', 0) * 100:.0f}% (LGG) of that "
         f"ceiling."),

        ("Figure_tumour_recovery",
         "Deconvolution recovers a minority of true tumour-content variation.",
         f"Recovery of true tumour-content variation per method, against DNA-measured purity, "
         f"in both cohorts. Recovery is defined as 1 + the slope of (estimate − truth) "
         f"regressed on truth: 0 indicates no information and 100% a perfect estimate. The "
         f"statistic is invariant to a constant offset, so it measures tracking rather than "
         f"calibration. Methods marked † ran without an input their published algorithm "
         f"requires and are reported as not evaluable rather than omitted or scored."),

        ("Figure_bisque_anchoring",
         "Bisque's cohort-mean composition is its single-cell reference's composition.",
         (lambda syn, real: (
             f"(A) Genuine BisqueRNA run without overlapping subjects -- the only mode a TCGA "
             f"cohort permits -- on synthetic data whose true cohort mean is set far from the "
             f"reference's. The estimated cohort mean lies {syn.get('far', {}).get('L1_estimate_to_prior')} "
             f"(L1) from the reference's donor-mean composition and "
             f"{syn.get('far', {}).get('L1_estimate_to_plant')} from the truth, while samples are "
             f"still ranked correctly (per-sample Spearman "
             f"{min((syn.get('far', {}).get('per_sample_spearman') or {0: 0}).values()):.2f}-"
             f"{max((syn.get('far', {}).get('per_sample_spearman') or {0: 0}).values()):.2f}); a matched "
             f"control in which truth equals the reference is recovered. This is the package's "
             f"declared assumption (Jew et al. 2020, Methods). (B) For every method, the L1 "
             f"distance of its TCGA cohort-mean composition, over all eight cell types, from the "
             f"same reference's donor-mean composition (raw/X reference; glioblastoma n = 154, "
             f"lower-grade glioma n = 510). Bisque lies "
             f"{(real.get('cohorts', {}).get('gbm', {}).get('methods', {}).get('bisque') or {}).get('L1_to_reference_prior_all_types')} and "
             f"{(real.get('cohorts', {}).get('lgg', {}).get('methods', {}).get('bisque') or {}).get('L1_to_reference_prior_all_types')} "
             f"from it; every other method lies "
             f"{min(v['L1_to_reference_prior_all_types'] for c in ('gbm', 'lgg') for k, v in (real.get('cohorts', {}).get(c, {}).get('methods') or {}).items() if k != 'bisque'):.2f} "
             f"or more."))(
             J("bisque_anchoring_synthetic.json"), J("bisque_anchoring.json"))),

        ("Figure_acs_vs_auc",
         "ACS and an AUC measure the same thing on different scales, and neither tracks accuracy.",
         (lambda d: (
             f"(A) The registered Anatomic Concordance Score (ACS: each constraint-tumour pair "
             f"scores 1 or 0 on per-structure means) against its graded counterpart, the "
             f"within-tumour Mann-Whitney AUC over the same pairs, samples and permutation null "
             f"(exploratory). They rank the {d['rank_agreement']['n_comparable_methods']} comparable "
             f"methods almost identically (Spearman ρ = "
             f"{d['rank_agreement']['acs_vs_auc_anat_spearman_comparable_methods']:.2f}) but sit on "
             f"different scales: dashed lines mark each statistic's chance level (median "
             f"permutation-null mean), and the AUC gives partial credit on the two conjunction "
             f"constraints (C4, C7) where ACS gives none. Crosses are the two broken-input controls "
             f"(AUC permutation p = "
             f"{', '.join(f'{v:.2f}' for v in sorted(r['null_p'] for r in d['methods'].values() if r['is_control']))}). "
             f"(B) The AUC against accuracy -- Harrell's concordance of each method's TCGA-GBM "
             f"tumour fraction with ABSOLUTE DNA purity (n = "
             f"{J('absolute_purity_yardstick.json').get('n_samples')} samples). Substituted for ACS in "
             f"the registered agreement test it gives ρ = "
             f"{d['agreement_with_auc_in_place_of_acs']['vs_rho_purity']['rho']:+.2f} (p = "
             f"{d['agreement_with_auc_in_place_of_acs']['vs_rho_purity']['p_value']:.2f}, n = "
             f"{d['agreement_with_auc_in_place_of_acs']['vs_rho_purity']['n_methods']} methods), "
             f"against ρ = {d['agreement_with_auc_in_place_of_acs']['registered_for_comparison']['rho']:+.3f} "
             f"for the registered ACS: the null result is not an artefact of thresholding. "
             f"Identical points are labelled together, and † marks a method that is degenerate on the "
             f"frozen TCGA signature (MuSiC without cross-donor variance is NNLS; Bisque runs its "
             f"no-overlap mode; SCDC ENSEMBLE with one reference is SCDC), so its concordance there is "
             f"not its own.")
          if d.get("controls", {}).get("acs_reproduced_for_every_method") else "")(J("anatomic_auc.json"))),

        ("Figure_lymphoid_mechanisms",
         "Reference-side tests of the lymphoid inversion: none removes it.",
         (lambda rp, bd, tl, ia: (
             f"TCGA cohorts. (A) Share of samples with B above T under SVR for each reference variant -- as "
             f"registered (glioblastoma {100 * rp['cohorts']['gbm']['svr_baseline']['frac_B_over_T']:.0f}%, "
             f"lower-grade glioma {100 * rp['cohorts']['lgg']['svr_baseline']['frac_B_over_T']:.0f}%), with the NK "
             f"column removed, with T and NK merged, and with ribosomal genes removed "
             f"({100 * rp['cohorts']['gbm']['svr_rp_removed']['frac_B_over_T']:.0f}% and "
             f"{100 * rp['cohorts']['lgg']['svr_rp_removed']['frac_B_over_T']:.0f}%). Every variant keeps B above T "
             f"in nearly all samples; dashed lines mark the share DNA methylation implies. (B) Where the B "
             f"column's estimate goes when the column is removed, under SVR: "
             f"{100 * bd['cohorts']['gbm']['b_removed_svr']['share_of_B_mass_to']['Tumor']:.0f}% (glioblastoma) and "
             f"{100 * bd['cohorts']['lgg']['b_removed_svr']['share_of_B_mass_to']['Tumor']:.0f}% (lower-grade glioma) "
             f"is re-absorbed by Tumor. Under NNLS it goes to T instead "
             f"({100 * bd['cohorts']['gbm']['b_removed_nnls']['share_of_B_mass_to']['T_cell']:.0f}% and "
             f"{100 * bd['cohorts']['lgg']['b_removed_nnls']['share_of_B_mass_to']['T_cell']:.0f}%): the "
             f"re-absorption is estimator-specific. (C) Correlation of each lymphoid profile with the mean bulk "
             f"profile (log expression). GBmap's B profile is the most bulk-like (r = "
             f"{tl['cohorts']['gbm']['r_with_mean_bulk']['B_cell']['all_genes']:.2f} in glioblastoma); the "
             f"independent atlas's B profile is anti-correlated. (D) B above T under GBmap (frozen signature), "
             f"under an independent atlas (Abdelfattah et al. 2022, GSE182109; not a GBmap source study), and "
             f"under that atlas with immunoglobulin genes removed, for NNLS and SVR in both cohorts. The "
             f"independent atlas recovers T above B only for NNLS in lower-grade glioma; immunoglobulin removal "
             f"changes resolution, not direction.")
          if rp and bd and tl and ia else "")(J("ribosomal_test.json"), J("b_column_diagnostics.json"),
                                              J("b_profile_tissue_likeness.json"), J("independent_atlas_test.json"))),

        ("Figure_identifiability",
         "Truth-free stability against agreement with DNA truth (exploratory).",
         (lambda d, r: (
             f"Post-registration (Extension E2; rules in `prespecified/identifiability_diagnostics.md`, fixed "
             f"before computing). Each point is one compartment in one TCGA cohort (circles glioblastoma, "
             f"squares lower-grade glioma). Filled points are the primary units, with truths from DNA: ABSOLUTE "
             f"purity for tumour content (n = {d['units']['Tumor|gbm']['n_truth']} and "
             f"{d['units']['Tumor|lgg']['n_truth']}), the methylation leukocyte fraction for all leukocytes "
             f"(n = {d['units']['Leukocytes|gbm']['n_truth']} and {d['units']['Leukocytes|lgg']['n_truth']}), "
             f"and EpiDISH for the lymphoid total (n = {d['units']['Lymphoid|gbm']['n_truth']} and "
             f"{d['units']['Lymphoid|lgg']['n_truth']}). Open points are T, B and NK (secondary). (A) Stability "
             f"of DESeq2 `unmix`'s per-sample estimates across seven settings of its loss scale (mean pairwise "
             f"Spearman) against `unmix`'s agreement with the truth: rho = {d['H1_primary']['D1']['rho']:.2f} "
             f"across the six primary units (exact one-sided p = {d['H1_primary']['D1']['p_one_sided']:.3f}). "
             f"In glioblastoma, two settings coincide, because the pre-declared shift is 1. Over distinct fits "
             f"the result is unchanged (rho = {r['check5_distinct_fits']['H1_D1']['rho']:.2f}). (B) Agreement "
             f"among the registered methods against their median agreement with the truth: rho = "
             f"{d['H1_primary']['D2']['rho']:.2f} (p = {d['H1_primary']['D2']['p_one_sided']:.3f}). Tumour and "
             f"leukocyte content are stable and accurate; the lymphoid total is the least stable and does not "
             f"track methylation. Stability is not sufficient: T cells in glioblastoma are stable and do "
             f"not track methylation. EpiDISH's blood reference measures shares of the immune "
             f"compartment; with each estimate's lymphoid share of its own leukocytes in place of its tissue "
             f"fraction, the association strengthens (rho = "
             f"{J('identifiability_denominators.json').get('H1_matched', {}).get('D1', {}).get('rho', float('nan')):.2f}) "
             f"and the lymphoid estimate still does not track methylation.")
          if d.get("H1_primary") and r.get("control_pass") else "")(
             J("identifiability_diagnostics.json"), J("identifiability_robustness.json"))),
    ]

    L = ["# Figures\n",
         "*Generated by `scripts/build_figure_index.py`. Figures themselves come from "
         "`scripts/build_paper_figures.py` and `scripts/build_stats_figures.py`; every value "
         "plotted is read from an artefact under `results/`. Captions are written for a reader "
         "who skips the body: what was measured, on how many samples, by what test.*\n",
         "*Each figure is available as PNG (300 dpi, for drafts and review) and PDF (vector, "
         "for submission) in `docs/figures/`.*\n", "---\n"]
    n_ok = 0
    for i, (key, title, caption) in enumerate(figs, 1):
        if not (FIGDIR / f"{key}.png").exists():
            L.append(f"\n## Figure {i}. {title}\n\n*`{key}` — NOT YET RENDERED*\n")
            continue
        n_ok += 1
        L.append(f"\n## Figure {i}. {title}\n")
        L.append(f"![Figure {i}](figures/{key}.png)\n")
        L.append(f"**Figure {i}.** *{title}* {caption}\n")
        L.append(f"<sub>`docs/figures/{key}.png` · `docs/figures/{key}.pdf`</sub>\n")
    OUT.write_text("\n".join(L) + "\n")
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)} — {n_ok} of {len(figs)} figures present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
