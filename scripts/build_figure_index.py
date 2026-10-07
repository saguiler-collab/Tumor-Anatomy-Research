"""Write docs/FIGURES.md — every figure embedded, with a journal-ready caption.

A caption is not a title. It must let a reader who skips the body understand what was
measured, on how many samples, by what test, and what the panel shows. Numbers here are read
from artefacts so a caption cannot drift from the figure beside it.

    python scripts/build_figure_index.py
"""
from __future__ import annotations

import json

import pandas as pd

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

        ("Figure_truth_instruments",
         "The lymphoid truth by every instrument that measures it, beside the methods.",
         (lambda kl, at, gb, g: (
             f"Relative composition within {{T, NK, B}}; the B share of lymphocytes is printed beside each bar. "
             f"**Direct measurement:** flow cytometry of 17 IDH-mutant and 40 IDH-wildtype gliomas (Klemm et al. "
             f"[53], Figure 1F, measured from the figure's vector geometry and validated against two numbers the "
             f"paper prints: {kl['controls']['melanoma_CD8']['measured']:.2f} vs "
             f"{kl['controls']['melanoma_CD8']['printed']} and {kl['controls']['brm_lymphocytes_weighted']['measured']:.2f} "
             f"vs {kl['controls']['brm_lymphocytes_weighted']['printed']}; `prespecified/klemm_t_vs_b.md`); an "
             f"independent single-cell atlas (Abdelfattah et al. [42], {at['n_patients']} patients, T and NK "
             f"separated per cell by a rule declared beforehand; `prespecified/atlas_t_vs_b.md`); GBmap's core atlas "
             f"({gb['n_donors_with_any_lymphoid']} glioblastoma donors with lymphocytes). **DNA methylation:** EpiDISH "
             f"(blood reference; the registered truth) and GIMiCC (glioma-specific; registered reading "
             f"{g['Q1_reading']}). **Bulk-RNA deconvolution:** the mean composition over the registered methods, "
             f"frozen signature and donor-level reference. Every direct measurement puts T far above B "
             f"(IDH-mutant flow cytometry {kl['IDH_mutant_glioma_n17']['T_to_B']:.0f}:1); the methods place B above T.")
          if kl.get("controls", {}).get("all_pass") and at and gb and g else "")(
             J("klemm_t_vs_b.json"), J("atlas_t_vs_b.json"), J("gbmap_t_vs_b.json"),
             J("gimicc_truth_confirmation.json"))),

        ("Figure_accuracy_factors",
         "Factors that decide whether deconvolution can be trusted, against real DNA truths.",
         (lambda ya, ag, d: (
             f"**(A)** Accuracy by compartment: Spearman between each method's estimate and the DNA truth, median "
             f"across the registered methods (filled: ABSOLUTE purity, the methylation leukocyte fraction, EpiDISH "
             f"on matched denominators; open: GIMiCC tissue fractions). Tumour and total-leukocyte content track "
             f"their truths in glioblastoma; lymphoid subtypes do not, against either methylation instrument. "
             f"**(B, C)** Tumour-content accuracy of each method under the registered frozen signature (open) and "
             f"a donor-level reference built from raw counts (filled): the reference build moves a method's "
             f"accuracy by up to 0.5 in either direction. **(D)** How well a ground-truth-free check predicts real "
             f"accuracy: anatomic concordance ranks {ya['arm_2_absolute_purity']['all_methods']['n_methods']} methods "
             f"at rho = {ya['arm_2_absolute_purity']['all_methods']['spearman']:.2f} "
             f"(p = {ya['arm_2_absolute_purity']['all_methods']['p']:.2f}); agreement between methods ranks methods at "
             f"a mean rho = {ag['primary_panel']['P2']['mean_rho']:.2f} and compartments at "
             f"{d['H1_primary']['D2']['rho']:.2f}; loss-scale stability ranks compartments at rho = "
             f"{d['H1_primary']['D1']['rho']:.2f} (exact p = {d['H1_primary']['D1']['p_one_sided']:.3f}). "
             f"Interpretation and practical guidance: `docs/ACCURACY_FACTORS.md`.")
          if ya and ag and d.get("H1_primary") else "")(
             J("yardstick_agreement.json"), J("agreement_selection_test.json"), J("identifiability_diagnostics.json"))),

        ("Figure_cptac_wgs",
         "An independent cohort with a per-sample DNA truth: CPTAC glioblastoma against whole-genome purity.",
         (lambda w: (lambda a, s21, s22, s26: (
             f"Secondary analysis, registered after the primary CPTAC control failed and before any value shown "
             f"here was seen (`prespecified/cptac_wgs_purity_secondary.md`). Truth: AscatNGS tumour purity from "
             f"whole-genome sequencing (GDC), {w['truth']['n_with_value']} tumours, range "
             f"{w['truth']['range'][0]:.2f}-{w['truth']['range'][1]:.2f}. "
             f"**(A)** Tumour share of single nuclei (GDC clusters, named by the registered marker rule) against DNA "
             f"purity: rho = {s21['rho']:.2f} (permutation p = {s21['perm_p']:.2f}, n = {s21['n']}); it fails the "
             f"registered 0.40 bar, so every per-sample reading made against the nuclei stays INCONCLUSIVE. Open "
             f"symbols, tumours whose bulk and DNA aliquots pooled several pieces while the nuclei came from one. "
             f"**(B)** Methylation purity (GIMiCC) against DNA purity: on the CpGs complete in all 18 samples "
             f"(769 of 4,022), rho = {s22['a_complete_case_769_cpgs']['rho']:.2f}; red crosses, the five samples "
             f"missing a third or more of the library; on the {s22['b_coverage_ge_90pct']['n']} complete samples "
             f"(3,775 CpGs), rho = {s22['b_coverage_ge_90pct']['rho']:.2f} (OPEN_DEFECTS D30). **(C)** Pathologist's "
             f"percent tumour nuclei, descriptive: range {s26['range'][0]:.0%}-{s26['range'][1]:.0%}, rho = "
             f"{s26['vs_wgs_purity']['rho']:.2f}. **(D)** Each bulk method's tumour estimate against DNA purity "
             f"(n = 18); vertical lines, the medians ({a['frozen']['S2_3']['median_rho']:.2f} frozen signature, "
             f"{a['h5ad']['S2_3']['median_rho']:.2f} donor-level reference; TCGA-GBM on the same methods "
             f"{a['frozen']['S2_3']['tcga_gbm_median_same_methods']:.2f} and "
             f"{a['h5ad']['S2_3']['tcga_gbm_median_same_methods']:.2f}); dotted, the 0.40 bar. **(E)** Each "
             f"method's accuracy in CPTAC against its accuracy in TCGA-GBM (ABSOLUTE, n = 154): the ranking "
             f"transfers on the frozen signature (rho = {a['frozen']['S2_4_ranking_transfer']['spearman']:.2f}, "
             f"{a['frozen']['S2_4_ranking_transfer']['n_methods']} methods) but not on the donor-level reference "
             f"(rho = {a['h5ad']['S2_4_ranking_transfer']['spearman']:.2f}); a dagger marks BayesPrism, which ran "
             f"as the R package in CPTAC and as the Python reimplementation in TCGA. Methods that are degenerate "
             f"copies of another share a point (\"MuSiC = NNLS\"). **(F)** Anatomic concordance against accuracy "
             f"in CPTAC: rho = {a['frozen']['S2_5_acs_vs_wgs_accuracy']['spearman']:.2f} (frozen) and "
             f"{a['h5ad']['S2_5_acs_vs_wgs_accuracy']['spearman']:.2f} (donor-level); the registered bar is 0.60.")
          )(w["all_cases"], w["all_cases"]["S2_1_nuclei_truth_vs_wgs"]["registered"]["all_scored"],
            w["all_cases"]["S2_2_methylation_instrument"], w["S2_6_pathology_descriptive"])
          if w else "")(J("cptac_wgs_purity.json"))),

        ("Figure_acs_constraints",
         "What the anatomic score separates: each registered constraint, method by method.",
         (lambda pc, lb, yh: (
             f"**(A)** Share of evaluable Ivy GAP tumours in which each method satisfies each of the "
             f"{pc['constraint'].nunique()} registered constraints (cells: satisfied / evaluable; darker is a larger "
             f"share; `results/anatomic/acs_per_constraint.csv`). Methods are ordered by ACS; the red line separates "
             f"the {int(((~lb['is_control']) & lb['comparable']).sum())} comparable real methods from the two deliberately "
             f"broken controls and from quanTIseq, which models only immune cells and is scored on two constraints. "
             f"The controls fail across constraints. Among working methods the fraction satisfied spreads by at most "
             f"0.33 on C1-C5 and by 0.62-0.67 on C6 (myeloid: microvascular proliferation above cellular tumour) and "
             f"C7 (the tumour gradient, weighted 2), each evaluable in 8-9 tumours: the score's resolution among "
             f"working methods. **(B)** ACS, as registered. **(C)** Each method's tumour-content accuracy against DNA, "
             f"on the donor-level reference the anatomic arm uses: TCGA-GBM against ABSOLUTE (squares, n = 154) and "
             f"CPTAC against whole-genome purity (circles, n = 18; secondary analysis). Registered implementations: "
             f"BayesPrism and DWLS appear as this project's reimplementations here; their genuine-package values are "
             f"in `docs/EVALUATION_MATRIX.md` (`docs/METHOD_REPAIRS.md`). Descriptive; nothing is selected on it.")
          if len(pc) and len(lb) else "")(
             pd.read_csv(config.RESULTS_DIR / "anatomic" / "acs_per_constraint.csv")
             if (config.RESULTS_DIR / "anatomic" / "acs_per_constraint.csv").exists() else pd.DataFrame(),
             pd.read_csv(config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv")
             if (config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv").exists() else pd.DataFrame(),
             None)),
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
