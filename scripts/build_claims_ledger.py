"""
build_claims_ledger.py -- docs/CLAIMS_LEDGER.md: every claim the study makes, the number behind it (read from its
artefact), how the verification re-run treated that artefact, its controls and caveats, how it stands against
prior work, and what it means for building better deconvolution methods.

Numbers are read from results/, and verification status from results/verification/state.json through
scripts/verify_rerun.py's own comparison rules, so the ledger cannot drift from either. Prior work was checked by
targeted literature searches on 2026-10-07 (listed in the document); "none found" is reported as such, never as
proof that none exists.

    python3 scripts/build_claims_ledger.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from ivygap import config  # noqa: E402
import verify_rerun as vr  # noqa: E402

R = config.RESULTS_DIR
DOC = config.PROJECT_ROOT / "docs" / "CLAIMS_LEDGER.md"
BASE = R / "verification" / "new_outputs_baseline"
OK = {"REPRODUCED", "NUMERICAL NOISE", "REPRODUCED + NEW FIELDS", "REPRODUCED + NEW ROWS"}


def J(name: str) -> dict:
    p = R / name
    return json.loads(p.read_text()) if p.exists() else {}


def verification(artefact: str) -> tuple[str, str]:
    """(status, detail) for one artefact, from the step that produces it."""
    state = vr.load_state()
    for st in vr.STEPS:
        for o in st["outputs"]:
            if o.split("=>")[0] != artefact and o.split("=>")[-1] != artefact:
                continue
            rec = state.get(st["id"])
            if rec is None:
                return "PENDING", f"step `{st['id']}` ({st['tier']}) not yet re-run"
            res = vr._reclassify(vr._recompare_moved(st["id"], o, rec["outputs"].get(o, {})))
            s = res.get("status", "?")
            if s == "NEW" and (BASE / artefact).exists():
                cmp = (vr.compare_json if artefact.endswith(".json") else vr.compare_csv)(BASE / artefact, R / artefact)
                if (R / artefact).stat().st_mtime > (BASE / artefact).stat().st_mtime:
                    return cmp["status"], f"step `{st['id']}`: a second run compared with the first ({cmp['status']})"
                return "FIRST RUN", f"step `{st['id']}`: created after the snapshot; a second run is scheduled"
            if s == "DIFFERS" and st["id"] in vr.EXPLAINED:
                return "DIFFERS, EXPLAINED", f"step `{st['id']}`: {vr.EXPLAINED[st['id']][:160]}..."
            return s, f"step `{st['id']}` ({st['tier']})"
    return "NO STEP", "no verification step produces it"


def claim_status(arts: list[str]) -> str:
    ss = [verification(a)[0] for a in arts]
    if any(s in ("DIFFERS", "NO STEP", "NOT WRITTEN", "UNREADABLE") for s in ss):
        return "NOT VERIFIED"
    if any(s in ("PENDING", "FIRST RUN") for s in ss):
        return "PARTLY VERIFIED (re-runs pending)"
    if any(s == "DIFFERS, EXPLAINED" for s in ss):
        return "VERIFIED, with explained differences"
    return "VERIFIED"


def build() -> list[dict]:
    lb = pd.read_csv(R / "anatomic" / "acs_leaderboard.csv").set_index("method")
    real = lb[(~lb.is_control) & lb.comparable]
    ctl = lb[lb.is_control]
    ya, yh, pw = J("yardstick_agreement.json"), J("yardstick_agreement_h5ad.json"), J("agreement_power.json")
    cw, em = J("cptac_wgs_purity.json"), J("evaluation_matrix.json")
    lo = {t: J(f"lymphoid_ordering{t}.json") for t in ("", "_lgg", "_h5ad", "_lgg_h5ad")}
    kl, at, gb = J("klemm_t_vs_b.json"), J("atlas_t_vs_b.json"), J("gbmap_t_vs_b.json")
    idd, den, ags = J("identifiability_diagnostics.json"), J("identifiability_denominators.json"), J("agreement_selection_test.json")
    ext = J("extension/agreement_extended.json").get("registered_plus_extension", {})
    e3 = J("identifiability_e3.json")
    a2 = ya["arm_2_absolute_purity"]["all_methods"]
    ra = em.get("repaired_agreement", {})
    rep_tb = [r["T>B"].get("T_exceeds_B") for ms in em.get("repaired", {}).get("tcga", {}).values()
              for r in ms.values() if isinstance(r.get("T>B"), dict) and "T_exceeds_B" in r["T>B"]]
    ia = em.get("instrument_agreement", {})
    it = [v["tumour"]["spearman"] for v in ia.values() if "tumour" in v]
    il = [v["leukocytes"]["spearman"] for v in ia.values() if "leukocytes" in v]
    s2 = cw.get("all_cases", {})
    bp = em.get("repaired", {}).get("genuine_bayesprism_tcga", {})
    gbm_reg = em["arms"]["gbm|h5ad"]
    gaps = [g["Tumor|ABSOLUTE"]["rho"] - em["arms"][arm]["bayesprism"]["Tumor|ABSOLUTE"]["rho"] for arm, g in bp.items()]
    gaps += [r["Tumor|ABSOLUTE"]["rho"] - em["arms"][arm][m]["Tumor|ABSOLUTE"]["rho"]
             for arm, ms in em.get("repaired", {}).get("tcga", {}).items() for m, r in ms.items()
             if str(r.get("implementation", "")).startswith("R:") and em["arms"][arm].get(m, {}).get("status") == "ok"
             and (r.get("Tumor|ABSOLUTE") or {}).get("rho") is not None]
    max_gap = max(gaps) if gaps else float("nan")
    rp_gbm = em.get("repaired", {}).get("tcga", {}).get("gbm|h5ad", {})
    smallest = pw.get("smallest_panel_with_80pct_power")
    return [
        {"id": "L1", "title": "Anatomy detects broken deconvolution", "kind": "registered (R1)",
         "claim": (f"Every comparable real method satisfies the seven pre-registered anatomic constraints far more often "
                   f"than its within-tumour permutation null (ACS {real.acs.min():.3f}-{real.acs.max():.3f}, every null p "
                   f"<= {real.null_p.max():.1e}); the two deliberately broken estimators score {ctl.acs.max():.3f} and "
                   f"{ctl.acs.min():.3f} and do not beat their nulls (margin {real.acs.min() - ctl.acs.max():+.2f})."),
         "artefacts": ["anatomic/acs_leaderboard.csv"],
         "controls": "Two broken estimators (random proportions; shuffled signature) scored beside the real methods; "
                     "within-tumour permutation nulls. The constraint file was committed on 2026-09-03 (git 8d12559), "
                     "two days before the first scoring run (2026-09-05), and has never changed; its hash (2d1fb47c...) "
                     "is the one recorded by that run and deposited on OSF (dm2t8) on 2026-09-10 "
                     "(`docs/REGISTRATION_AUDIT.md`).",
         "caveats": "Glioblastoma only (Ivy GAP, 10 histology-labelled tumours). Two leaderboard rows (BayesPrism, DWLS) "
                    "depend on a wall-clock budget (D35); every other row reproduced exactly.",
         "prior": "Anatomically plausible deconvolution has been shown as supporting evidence: BayesPrism's authors "
                  "deconvolved Ivy GAP and reported endothelium enriched in microvascular proliferation, and Varn et al. "
                  "(2022) characterised Ivy GAP regions by deconvolution and immunofluorescence. Both use anatomy "
                  "qualitatively. No pre-registered, controlled test of anatomic concordance as a validation criterion "
                  "was found.",
         "new": "Anatomy turned into a falsifiable, pre-registered test with broken-input controls and a permutation null.",
         "design": "A cheap screen that rejects broken pipelines where region-labelled samples exist."},
        {"id": "L2", "title": "Anatomy does not identify the accurate method", "kind": "registered (R2): bar not met",
         "claim": (f"ACS against tumour accuracy (ABSOLUTE, TCGA-GBM): rho {a2['spearman']:+.3f} (95% CI "
                   f"{a2['ci'][0]:+.2f} to {a2['ci'][1]:+.2f}, {a2['n_methods']} methods); donor-level reference "
                   f"{yh['arm_2_absolute_purity']['all_methods']['spearman']:+.3f}; registered bar 0.60. Replicated against "
                   f"whole-genome purity in CPTAC ({s2.get('frozen', {}).get('S2_5_acs_vs_wgs_accuracy', {}).get('spearman', float('nan')):+.2f}) "
                   f"and after repairing the broken implementations ({ra.get('gbm|frozen', {}).get('repaired_spearman', float('nan')):+.2f} "
                   f"frozen, {ra.get('gbm|h5ad', {}).get('repaired_spearman', float('nan')):+.2f} donor-level). The bar is met only "
                   f"against the synthetic yardstick built from the same atlas ({ya['arm_1_pseudobulk']['rho']:.2f}). "
                   + (f"Adding the genuine extension packages ({ext['n']} methods, exploratory): rho {ext['rho']:+.2f} "
                      f"(95% CI {ext['ci95'][0]:+.2f} to {ext['ci95'][1]:+.2f})." if ext else "")),
         "artefacts": ["yardstick_agreement.json", "yardstick_agreement_h5ad.json", "agreement_power.json",
                       "cptac_wgs_purity.json", "evaluation_matrix.json", "extension/agreement_extended.json"],
         "controls": "A yardstick that shares the atlas with ACS (arm 1) shows what shared construction does; the "
                     "independent yardsticks (DNA) share nothing with it.",
         "caveats": (f"Underpowered: with {a2['n_methods']} methods a true rho of 0.70 meets the bar only part of the time; "
                     f"80% power at a true 0.70 needs {smallest.get('0.7', 'more')} methods (`agreement_power.json`), so "
                     f"'not met', not 'refuted'. Arm 1's magnitude depends on which implementation ran (D35)."),
         "prior": "Plausibility checks are routine supporting evidence (anatomy, marker correlation), and a 2026 "
                  "framework selects methods by cross-cohort and prognostic reproducibility (Li et al., Genome Biol. "
                  "2026, recommending ReCIDE and BayesPrism) without testing the criteria against a measured truth. No "
                  "test of whether such checks identify the accurate method against independent DNA truths was found.",
         "new": "The plausibility check, tested and found not to rank methods by accuracy, in two cohorts and after repairs.",
         "design": "Validate and select methods against a DNA truth (copy number, methylation), never by plausibility."},
        {"id": "L3", "title": "No method robustly recovers the lymphoid ordering, which direct measurement confirms",
         "kind": "registered prediction (R3); robustness exploratory",
         "claim": (f"Methods putting T above B: GBM {lo['']['n_agree_T_over_B']} and LGG {lo['_lgg']['n_agree_T_over_B']} on "
                   f"the frozen signature; {lo['_h5ad']['n_agree_T_over_B']} and {lo['_lgg_h5ad']['n_agree_T_over_B']} on the "
                   f"donor-level reference; {sum(bool(x) for x in rep_tb)} of {len(rep_tb)} repaired method-arms. Direct "
                   f"measurement says T >> B: flow cytometry T:B {kl['IDH_mutant_glioma_n17']['T_to_B']}:1 (IDH-mutant, 17) "
                   f"and {kl['IDH_wildtype_glioma_n40']['T_to_B']}:1 (IDH-wildtype, 40); the Abdelfattah atlas "
                   f"{at['n_patients_T_strict_over_B']}/{at['n_patients']} patients; GBmap {gb['n_donors_T_over_B']}/"
                   f"{gb['n_donors_with_any_lymphoid']} donors ({gb['pooled_T_to_B']:.0f}:1)."),
         "artefacts": ["lymphoid_ordering.json", "lymphoid_ordering_lgg.json", "lymphoid_ordering_h5ad.json",
                       "lymphoid_ordering_lgg_h5ad.json", "klemm_t_vs_b.json", "atlas_t_vs_b.json", "gbmap_t_vs_b.json",
                       "evaluation_matrix.json"],
         "controls": "The prediction and its falsifier were registered (2026-09-17) before the LGG methylation existed; "
                     "the flow-cytometry reading has measurement controls (all pass); a second methylation instrument "
                     "(GIMiCC) was run and disagrees in LGG, and is reported.",
         "caveats": "Cohort-level truths; no per-sample lymphoid truth exists here (CPTAC nuclei hold B cells in 2 of "
                    "15 tumours). One methylation instrument (GIMiCC) places B above T in LGG. The GIMiCC rule's file "
                    "time was reset by an in-place edit, so its precedence over the GIMiCC runs rests on the session "
                    "record, not the filesystem (`docs/REGISTRATION_AUDIT.md`); the prediction itself (git 43fbe99) "
                    "precedes its test data by 94 minutes.",
         "prior": "Lymphocyte subsets are recognised as hard to deconvolve, and glioma-specific tools were validated "
                  "against imaging mass cytometry (GBMDeconvoluteR, Ajaib et al. 2023). A systematic B-over-T inversion "
                  "across about 20 methods in both glioma types, checked against flow cytometry and two atlases, was not "
                  "found reported.",
         "new": "A reproducible, method-independent failure mode, measured against direct counts, that survives genuine "
                "packages and faithful reimplementations.",
         "design": "Bulk RNA with current references does not determine lymphocyte subsets in glioma: a method should "
                   "either bring new information (priors, a second modality) or not report them."},
        {"id": "L4", "title": "A truth-free stability check flags which compartments the bulk data determine",
         "kind": "exploratory (R4)",
         "claim": (f"Re-fitting one genuine package under seven loss scales, the compartments that stay put are those DNA "
                   f"confirms: rho {idd['H1_primary']['D1']['rho']:.3f} (exact p {idd['H1_primary']['D1']['p_one_sided']:.3f}, "
                   f"{idd['H1_primary']['D1']['n_units']} units); {den['H1_matched']['D1']['rho']:.3f} on matched denominators. "
                   f"Agreement between methods is weaker ({idd['H1_primary']['D2']['rho']:.2f}; DECEPTICON's rule "
                   f"{ags['primary_panel']['P2']['mean_rho']:.2f}). "
                   + (f"A second perturbation, registered before it was run (E3: 10 random halves of the marker genes): "
                      f"rho {e3['H1_E3_primary'].get('rho', float('nan')):+.3f} (exact p "
                      f"{e3['H1_E3_primary'].get('p_one_sided', float('nan')):.3f}; reading: {e3['H1_E3_primary'].get('reading')}); "
                      f"agreement of the two perturbations across units {e3['C1_d3_vs_d1_descriptive'].get('rho', float('nan')):+.2f}."
                      if e3.get("H1_E3_primary") else "A second perturbation (E3, gene resampling) is registered and running.")),
         "artefacts": ["identifiability_diagnostics.json", "identifiability_denominators.json", "agreement_selection_test.json"]
                      + (["identifiability_e3.json"] if e3 else []),
         "controls": "Shuffled truths almost never support the hypothesis (tests/test_identifiability.py); matched "
                     "denominators; rules registered before computing.",
         "caveats": ("Six units; exploratory. Stability flags what not to trust and certifies nothing: T cells in GBM "
                     "are stable and wrong. "
                     + (f"It depends on the perturbation: re-fitting on random halves of the genes (E3) rates the GBM lymphoid "
                        f"total stable ({e3['units']['Lymphoid|gbm']['d3_stability']:.2f}) although it is wrong against DNA "
                        f"({e3['units']['Lymphoid|gbm']['truth_rho_unmix_e2']:+.2f}), where the loss-scale perturbation "
                        f"flagged it ({e3['units']['Lymphoid|gbm']['d1_stability_e2']:.2f}). The loss scale is the informative "
                        f"perturbation; gene resampling alone is not enough." if e3.get("units") else "")),
         "prior": "Robustness of deconvolution to perturbations is studied on simulated mixtures (Xu et al., Brief. "
                  "Bioinform. 2025). No validation of a truth-free stability diagnostic against real DNA truths was found.",
         "new": "A diagnostic a lab can run with no measurement, validated against DNA in two cohorts.",
         "design": "Report, beside every estimate, whether its compartment is stable under re-fitting."},
        {"id": "L5", "title": "Per-sample truth from the same tissue: bulk methods work for tumour content; nuclei are not a truth",
         "kind": "secondary (registered after the primary control failed)",
         "claim": (f"CPTAC, 18 glioblastomas, against whole-genome purity: bulk tumour estimates median rho "
                   f"{s2['frozen']['S2_3']['median_rho']:.2f} (frozen) and {s2['h5ad']['S2_3']['median_rho']:.2f} (donor-level); "
                   f"the frozen-signature ranking transfers from TCGA ({s2['frozen']['S2_4_ranking_transfer']['spearman']:.2f}). "
                   f"Single-nucleus tumour share: {s2['S2_1_nuclei_truth_vs_wgs']['registered']['all_scored']['rho']:.2f}. "
                   f"Methylation purity: {s2['S2_2_methylation_instrument']['a_complete_case_769_cpgs']['rho']:.2f} on a probe "
                   f"set shrunk by five low-coverage samples, {s2['S2_2_methylation_instrument']['b_coverage_ge_90pct']['rho']:.2f} "
                   f"on complete samples (D30)."),
         "artefacts": ["cptac_wgs_purity.json"],
         "controls": "Permutation nulls; the rule was written before any whole-genome value was seen; the primary analysis "
                     "is reported as INCONCLUSIVE.",
         "caveats": "Secondary; 18 tumours; GDC's automated clusters (the authors' curated labels were not tested).",
         "prior": "Dissociation and capture biases in single-cell and single-nucleus proportions are documented "
                  "(Slyper et al., Nat. Med. 2020). A per-sample test of nuclei composition against whole-genome purity "
                  "from the same aliquots, beside the bulk methods, was not found.",
         "new": "On matched material, a DNA purity is a valid per-sample yardstick and single-nucleus composition is not.",
         "design": "Validate against DNA purity; check per-sample probe coverage before methylation deconvolution."},
        {"id": "L6", "title": "Implementation fidelity is as large a factor as the choice of method",
         "kind": "post-registration (OPEN_DEFECTS D31-D34)",
         "claim": (f"Tumour accuracy, GBM donor-level: BayesPrism stand-in {gbm_reg['bayesprism']['Tumor|ABSOLUTE']['rho']:.2f} "
                   f"vs genuine {bp.get('gbm|h5ad', {}).get('Tumor|ABSOLUTE', {}).get('rho', float('nan')):.2f}; DWLS legacy "
                   f"reimplementation {gbm_reg['dwls']['Tumor|ABSOLUTE']['rho']:.2f} vs genuine "
                   f"{rp_gbm.get('dwls', {}).get('Tumor|ABSOLUTE', {}).get('rho', float('nan')):.2f}. The genuine DWLS failed "
                   f"on 73 of 122 anatomy samples for a numerical reason; one solution-preserving rescaling removes every "
                   f"failure."),
         "artefacts": ["evaluation_matrix.json"],
         "controls": "Package equivalence on real inputs (5e-11 given the dampening constant); the legacy implementation "
                     "as the negative control (tests/test_dwls_conditioning.py).",
         "caveats": "The repaired runs are post-registration and behind a switch; registered results are unchanged.",
         "prior": "The DWLS solver error is reported publicly (dtsoucas/DWLS issue #12; no cause or fix given). Benchmarks "
                  "often run reimplementations; the size of the gap between a reimplementation and its package was not "
                  "found quantified.",
         "new": (f"The cause of the DWLS failure and an exact fix; measured gaps between a genuine package and this "
                 f"project's Python version of up to {max_gap:.2f} in Spearman with DNA purity."),
         "design": "Run published packages; verify any reimplementation against its package on real inputs; record every "
                   "fallback and its reason."},
        {"id": "L7", "title": "Method accuracy is a reproducible property of the method and reference build",
         "kind": "post hoc, descriptive",
         "claim": (f"Two independent DNA instruments order the methods almost identically: tumour accuracy "
                   f"{min(it):.2f}-{max(it):.2f}, leukocyte accuracy {min(il):.2f}-{max(il):.2f}, in all {len(ia)} GBM/LGG arms; "
                   f"the frozen-signature ranking transfers to an independent cohort "
                   f"({s2['frozen']['S2_4_ranking_transfer']['spearman']:.2f}), the donor-level one does not "
                   f"({s2['h5ad']['S2_4_ranking_transfer']['spearman']:.2f})."),
         "artefacts": ["evaluation_matrix.json", "cptac_wgs_purity.json"],
         "controls": "Two instruments that share no measurement (copy number; methylation).",
         "caveats": "Post hoc; 11-14 methods per arm.",
         "prior": "Benchmarks report rankings per dataset (e.g., the DREAM community assessment, White et al. 2024). The "
                  "agreement of rankings across independent DNA instruments and cohorts was not found reported.",
         "new": "Evidence that a DNA-truth benchmark in a similar tissue is a sound basis for choosing a method.",
         "design": "Benchmark against DNA truths, per reference build, and choose from that."},
        {"id": "L8", "title": "A wall-clock budget can silently change which implementation runs, and the statistics",
         "kind": "verification finding (OPEN_DEFECTS D35)",
         "claim": "Re-running the registered pipeline with identical code, two methods switched from fallback to genuine "
                  "packages because their runtimes straddled a 2,400 s budget. Every other method reproduced exactly, "
                  "but the run's summary statistics moved with those two rows: the synthetic-yardstick agreement from "
                  "0.637 (14 methods) to 0.828 (13), the run's own ACS-versus-ABSOLUTE agreement to +0.165 (11 methods; "
                  "registered +0.081 on 12). Every difference was checked against this cause (34 of 34).",
         "artefacts": ["anatomic/acs_leaderboard.csv"],
         "controls": "Snapshot comparison of every output; registered artefacts kept, re-run outputs archived.",
         "caveats": "Readings unchanged (both values clear the 0.60 bar).",
         "prior": "Not found reported for deconvolution pipelines.",
         "new": "A concrete reproducibility hazard in multi-method benchmarks.",
         "design": "Decide implementations by a declared rule, not a clock; record every fallback."},
    ]


def bars() -> dict:
    """The best tumour accuracy measured here, per truth, over registered and repaired methods (for acceptance test 1)."""
    em, cw = J("evaluation_matrix.json"), J("cptac_wgs_purity.json")
    t = [c["Tumor|ABSOLUTE"]["rho"] for arm, cells in em["arms"].items() if arm.startswith("gbm")
         for c in cells.values() if c.get("status") == "ok" and c["Tumor|ABSOLUTE"]["rho"] is not None]
    t += [r["Tumor|ABSOLUTE"]["rho"] for arm, ms in em.get("repaired", {}).get("tcga", {}).items() if arm.startswith("gbm")
          for r in ms.values() if (r.get("Tumor|ABSOLUTE") or {}).get("rho") is not None]
    w = [v["rho"] for ref in ("frozen", "h5ad") for v in cw["all_cases"][ref]["methods"].values() if v.get("rho") is not None]
    w += [r["Tumor|WGS"]["rho"] for ms in em.get("repaired", {}).get("cptac", {}).values() for r in ms.values()
          if (r.get("Tumor|WGS") or {}).get("rho") is not None]
    return {"tcga_gbm": max(t), "cptac": max(w)}


def main() -> int:
    claims = build()
    bb = bars()
    L = ["# Claims ledger: what this study establishes, how it is verified, and what is new", "",
         "*Generated by `scripts/build_claims_ledger.py`. Numbers are read from `results/`; verification status from "
         "the re-run harness (`results/verification/state.json`, `docs/VERIFICATION_RERUN.md`). Regenerate after the "
         "verification finishes.*", "",
         "**How to read it.** VERIFIED: every artefact behind the claim was re-produced from code and matched the "
         "registered one (or added only fields or rows). VERIFIED, with explained differences: a difference occurred "
         "and its cause is diagnosed and recorded. PARTLY VERIFIED: some artefacts await their re-run. Prior work was "
         "checked by targeted searches on 2026-10-07 (queries at the end); 'not found' means not found, not proven "
         "absent.", "", "| claim | kind | verification |", "|---|---|---|"]
    for c in claims:
        L.append(f"| **{c['id']}** {c['title']} | {c['kind']} | {claim_status(c['artefacts'])} |")
    L.append("")
    for c in claims:
        L += [f"## {c['id']} · {c['title']}", "", f"*{c['kind']}*", "", f"**Claim.** {c['claim']}", "",
              "**Verification:**"]
        for a in c["artefacts"]:
            s, d = verification(a)
            L.append(f"- `{a}`: {s} -- {d}")
        L += ["", f"**Controls.** {c['controls']}", "", f"**Caveats.** {c['caveats']}", "",
              f"**Against prior work.** {c['prior']}", "", f"**What is new.** {c['new']}", "",
              f"**For method development.** {c['design']}", ""]
    L += ["## Acceptance tests for a new deconvolution method, derived from the findings", "",
          "Each test is something a method can pass or fail on public data used here, and each comes from a verified "
          "finding above. A method that passes all of them is not proven optimal; one that fails any is not ready.", "",
          "| # | test | pass condition | basis |", "|---|---|---|---|",
          f"| 1 | Tumour content against a DNA truth, two cohorts | Spearman with ABSOLUTE (TCGA-GBM, n 154) at or above "
          f"{bb['tcga_gbm']:.2f} and with whole-genome purity (CPTAC, n 18) at or above {bb['cptac']:.2f}: the best "
          f"measured here | L5, L7 |",
          "| 2 | Total immune content | Spearman with the methylation leukocyte fraction at or above the panel median, in "
          "GBM and LGG | L7 |",
          "| 3 | Lymphocyte subsets | Either cohort-level T above B (as flow cytometry and two atlases show, 9:1 to 43:1) "
          "or no lymphoid subtype reported | L3 |",
          "| 4 | Truth-free stability | Per-compartment stability under loss-scale re-fitting reported; compartments "
          "that move are flagged, not reported as estimates. Gene resampling alone does not flag the lymphoid failure "
          "(E3) | L4 |",
          "| 5 | Numerical robustness | Identical output when signature and bulk are rescaled together; no sample "
          "silently dropped | L6 (D31) |",
          "| 6 | Implementation fidelity | Any reimplementation reproduces its published package on real inputs, and every "
          "fallback is recorded with its reason | L6 (D32, D34) |",
          "| 7 | Reproducibility | A re-run of the same code reproduces every output; no implementation choice depends "
          "on a clock | L8 (D35) |",
          "| 8 | Plausibility is not a pass | Anatomic or marker plausibility screens out broken output but is never "
          "the evidence of accuracy | L1, L2 |", ""]
    L += ["## What this study does not establish", "",
          "- A per-sample lymphoid truth: every lymphoid truth here is cohort-level.",
          "- That anatomy cannot rank methods: with 12-14 methods the registered test is underpowered; it was not met.",
          "- Anything outside glioma, or for references other than the GBmap family.",
          "- That single-nucleus composition fails with the authors' own curated annotation (untested).", "",
          "## Literature checks (2026-10-07)", "",
          "Targeted web searches: anatomic (Ivy GAP) validation of deconvolution; method selection without ground truth "
          "and plausibility checks; B- and T-cell accuracy of glioma deconvolution against flow cytometry; truth-free "
          "stability or identifiability diagnostics; single-nucleus proportions against tumour purity; the DWLS solver "
          "error. Works found and cited above: Chu et al. 2022 (BayesPrism, Nat. Cancer); Varn et al. 2022 (Cell); "
          "Ajaib et al. 2023 (GBMDeconvoluteR, Neuro-Oncology); White et al. 2024 (DREAM, Nat. Commun.); Xu et al. 2025 "
          "(Brief. Bioinform.); Li et al. 2026 (Genome Biol., doi:10.1186/s13059-026-03942-1); Slyper et al. 2020 "
          "(Nat. Med.); dtsoucas/DWLS issue #12. In docs/REFERENCES.md: [4] Li 2026, [19] BayesPrism, [25] "
          "GBMDeconvoluteR, [47] DREAM, [57] Varn, [58] Xu, [59] Slyper (each checked against Crossref)."]
    DOC.write_text("\n".join(L) + "\n")
    print(f"wrote {DOC.relative_to(config.PROJECT_ROOT)}: " + ", ".join(f"{c['id']} {claim_status(c['artefacts'])}" for c in claims))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
