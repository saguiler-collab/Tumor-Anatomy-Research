# Confirming the lymphoid truth with a glioma-specific methylation method (GIMiCC)

**Written 2026-10-03.** The time is this file's birth time (`stat`). At writing, GIMiCC had run only on
its authors' three example samples (the reproduction control below) and on no TCGA sample.

## Why this is the most important check left

The study's headline is that no deconvolution method robustly reproduces methylation's T > B
lymphoid ordering. That "truth" comes from EpiDISH run on brain tumour DNA with a **blood**
reference (`centDHSbloodDMC.m`):
- its seven columns are shares of blood cell types;
- brain and tumour DNA have nowhere to go but into them;
- its per-sample split failed a reference-free positive control (D23).

If the truth is wrong, the headline is wrong. GIMiCC (Pike et al., Acta Neuropathol Commun 2024; Salas
lab) is built for glioma:
- InfiniumPurify for tumour against non-tumour;
- then neuronal, glial, angiogenic and immune;
- then myeloid against lymphoid;
- then NK, B and T;
- then CD4 against CD8.

Each layer is weighted by its parent, so the outputs are fractions of the tissue.

## Software (verified before use)

- GIMiCC 0.99.1. The user's copy (`pipeline_packages / repos/GIMiCC-main`) is byte-identical to
  github.com/SalasLab/GIMiCC @ 26cb8a15. It was installed with `R CMD INSTALL`. The library comes from
  ExperimentHub EH9483.
- Dependencies: ExperimentHub 3.2.2, FlowSorted.Blood.EPIC 2.16.0, InfiniumPurify 1.3.1, minfi
  1.58.0, plus 32 supporting packages. **36 packages added; no existing package changed version**
  (versions recorded before and after).
- **Reproduction control: PASSED.** The authors' example (EH9482, three GBM samples, `tumor.type =
  "GBM"`, `h = 5`) reproduces the vignette's printed table: 54 of 54 values within print precision.

## Data

- TCGA-GBM and TCGA-LGG Illumina 450K betas: Xena, D07/D08, local and hash-checked in
  `docs/DATA_INVENTORY.md`.
- Only GIMiCC's library CpGs are extracted, the union over every layer and tumour-type library.
- **Missing values are never imputed.** A CpG missing in any sample of a cohort is dropped for that
  cohort -- the complete-case rule EpiDISH used. GIMiCC's own messages report the probes missing per
  layer, and they are recorded.
- **Samples:** every sample in each 450K matrix (primary); the RNA-matched subsets used in the
  registered comparisons, 56 GBM and 510 LGG (secondary).

## Tumour type

GIMiCC's `tumor.type` changes only Layer 0, the tumour/non-tumour split.
- **Primary:** each cohort is run with one default -- "GBM" for TCGA-GBM, "AST" for TCGA-LGG.
  1p/19q status is not available locally, so oligodendrogliomas cannot be assigned "OLG".
- **Sensitivity:** both cohorts are run under all four types.
- **Structural check:** within-lymphoid shares must be identical under every type (max |difference| <
  1e-9), as the code implies. Tissue fractions are reported under each type.

## Controls (read before the primary result)

1. **Reproduction:** passed (above).
2. **Positive, tumour layer:** GIMiCC Tumor against ABSOLUTE purity, Spearman >= 0.40 in both cohorts.
3. **Positive, immune total:** GIMiCC Immune + Microglia (all CD45+ cells; GIMiCC places microglia in
   its glial layer) against the Thorsson methylation leukocyte fraction, Spearman >= 0.40 in both
   cohorts. Immune alone is reported beside it.
4. **Negative:** CpG labels shuffled within the extracted matrix must break control 2 (|rho| < 0.20).
5. **Structural:** the invariance above.

**If control 2 or 3 fails,** GIMiCC's lymphoid output is still reported, with that failure stated next
to it; the confirmation reading below is then INCONCLUSIVE, never CONFIRMED.

## The primary question (Q1): does GIMiCC confirm methylation's T > B ordering?

- **Definition**, the registered criterion's own computation:
  - per sample, renormalise within {T, NK, B}, with T = CD4T + CD8T (GIMiCC level 4);
  - average over the samples with any lymphoid estimate;
  - compare the cohort-level T and B.
  - Samples with no lymphoid estimate (GIMiCC zeroes projections under 1e-5) are counted, not
    zero-filled.
- **CONFIRMED:** T > B in both cohorts, with controls 2 and 3 passing.
- **NOT CONFIRMED:** B >= T in either cohort. The study's headline would then rest on a truth a
  glioma-specific method contradicts, and the manuscript says so before anything else.
- **INCONCLUSIVE:** otherwise, for example a control failing.
- **Reported alongside:**
  - the full ordering of T, NK and B, against EpiDISH's T > NK > B;
  - the share of samples with T > B;
  - CD4 against CD8;
  - per-sample Spearman agreement with EpiDISH's within-lymphoid T share and B share.

## Secondary

- **Q2.** Is the headline re-scored against GIMiCC still "no method robustly reproduces T > B"? That
  is, do the registered and extension methods' cohort orderings match GIMiCC's?
- **Q3.** E2's lymphoid units against GIMiCC tissue fractions, whose denominators match by
  construction: does "the lymphoid estimate does not track methylation" hold against a second
  methylation truth?

Exploratory, post-registration. Nothing here is selected on an outcome, and no registered statistic
changes; the registered criterion was scored against EpiDISH and stays as recorded.

## Addendum 1 -- probe coverage and its robustness check (written after the first TCGA run's diagnostics, before any GIMiCC result was read)

The first run (TCGA-GBM, "GBM", h = 4) printed only GIMiCC's diagnostics; its estimates were not
opened.
- After the complete-case rule, 406 of the 4,022 library CpGs are missing. They are missing in
  essentially every sample: about 127 of 155 per probe, consistent with probes masked in TCGA's
  processing.
- GIMiCC reports the per-layer loss: L1 15/145, L2A 27/215, L2B 39/259, L2C 14/143, L2D 28/190,
  **L3A (the T / NK / B split) 46/266**, L3B 30/187, L4 (CD4 / CD8) 50/282, L5A 59/334, L5B 51/316,
  L5C 59/304.
- That is 10-21% per layer: a deviation from the published library, reported as such.
- GIMiCC also flagged noisy samples ("NaN ... removed"). These are counted, never zero-filled.

**Robustness check, fixed now.** 20 repeats, each dropping a further random 20% of the available
library CpGs across all layers (seed `config.RANDOM_SEED` + repeat). Each repeat recomputes the Q1
quantity.
- Reported: the share of repeats in which cohort-level T > B holds.
- **Reading:** robust if the Q1 direction holds in >= 18 of 20 repeats in each cohort; fragile
  otherwise. A fragile result cannot be read as CONFIRMED.

## Addendum 2 -- three fixes from reading GIMiCC's source (written 2026-10-03 13:03 EDT -- the file's modification time after appending, 13:03:51 -- before the analysis stage ran; `results/gimicc_truth_confirmation.json` did not exist, and no GIMiCC estimate had been opened -- only column names and run logs)

**1. Level-4 T can be lost or silently zeroed while Layer 3A's T is defined.**
The source (`R/GIMiCC_Deconvo.R`, lines 486-505) handles a layer whose projections are all under
1e-5 as follows:
- The renormalisation divides 0 by 0 and gives NaN. GIMiCC sets NaN to 0, then re-sums the sample.
- If the sum rounds below 1 (lost mass >= 0.005 of the tissue), the sample's row is restored with
  its NaNs (the "noisy ... removed" message).
- If the lost mass is smaller, the zeros stay.
- So a sample whose CD4/CD8 projection (Layer 4) fails has CD4 and CD8 at level 4 that are either
  NaN or exactly 0, even though Layer 3A gave it a T fraction.

The registered T (CD4 + CD8 at level 4) equals Layer 3A's T only when Layer 4 is defined. A silent
zero biases T downward, which is toward B >= T.

**Added check:** each cohort is also run at level 3, the default tumour type. GIMiCC's own output
there gives T, NK and B straight from Layer 3A, without Layer 4. Q1 is recomputed on it.
- Reported: the number of samples where level-4 CD4 + CD8 is 0 or NaN while level-3 T is positive.
- The registered reading stays the level-4 one.
- If level 3 gives the opposite direction in a cohort, the reading is stated together with that
  disagreement and its cause (samples with an undefined CD4/CD8 projection). It is not replaced.

**2. Control 3 must not zero-fill.**
The analysis summed the immune columns with pandas' default, which skips NaN, and skipping NaN in a
sum is zero-filling. Now:
- Samples with a NaN in any immune column are excluded from control 3, and counted.
- The level-3 version (Tcell + Bcell + NK + Mono + Neu + Microglia, which does not depend on Layer 4)
  is reported beside it.

**3. The reading follows this document's text.**
The analysis code applied the robustness requirement in both directions: a fragile B >= T became
INCONCLUSIVE. The text above restricts only CONFIRMED ("a fragile result cannot be read as
CONFIRMED"), and NOT CONFIRMED is "B >= T in either cohort".
- The reading now follows the text: a fragile B >= T is NOT CONFIRMED, stated with its fragility
  (k of 20).
- The symmetric reading is reported beside it.
- A failing control still makes the reading INCONCLUSIVE, as above.

## Addendum 3 -- exact definitions of Q2 and Q3 (clock read at 13:09:07 EDT, before writing; neither `results/gimicc_truth_confirmation.json` nor `results/gimicc_secondary.json` existed, and no GIMiCC estimate had been opened)

Q2 and Q3 above were stated in words only. Their computations are fixed here, in
`scripts/gimicc_secondary.py`, before GIMiCC's lymphoid estimates are seen.

**GIMiCC's lymphoid quantities for both questions:**
- Cohort ordering: Q1's level-4 registered computation, on the RNA-matched samples.
- Per-sample tissue fractions: level 3 (T, NK and B straight from Layer 3A, by Addendum 2), with
  level 4 reported beside them.
- NaN rows are excluded and counted, never zero-filled.

**Q2 -- the headline re-scored against GIMiCC.**
- Inputs: the methods' cohort orderings exactly as recorded in `results/lymphoid_ordering{,_lgg}{,_h5ad}.json`
  (the registered panel on both reference builds, both cohorts). They are not recomputed.
- Each method's recorded `T_exceeds_B` is compared with GIMiCC's cohort-level `T_exceeds_B`.
- Reported per cohort and reference: the number of methods matching GIMiCC's direction, beside the
  number matching EpiDISH's.
- Per-sample discordance is also reported, as in the registered per-sample table: the share of
  samples where the method puts B above T while GIMiCC puts T above B, and the reverse.
- **Reading:** the headline ("no method robustly reproduces the methylation T > B ordering")
  - is **unchanged** if GIMiCC's direction is T > B in both cohorts;
  - is **reversed for a cohort** if GIMiCC puts B >= T there. The methods that put B above T then
    agree with the glioma-specific truth, and the manuscript says so.

**Q3 -- E2's lymphoid units against a second methylation truth with matched denominators.**
- D1 and D2 are truth-free. They are taken unchanged from `results/identifiability_diagnostics.json`.
- Only the truth column changes:
  - Lymphoid: GIMiCC T + NK + B, tissue fraction.
  - T_cell, B_cell, NK_cell: each GIMiCC tissue fraction.
  - Leukocytes: GIMiCC Immune + Microglia, beside the Thorsson LF.
  - Tumor: ABSOLUTE, unchanged.
- Estimates: the saved per-sample estimates. These are the `unmix_{cohort}_predeclared.csv` E2
  wrote, and the registered panel's `estimates_full*.csv` with the same degenerate duplicates
  dropped.
- **Reproduction control first:** the saved predeclared unmix estimates must reproduce E2's
  recorded Tumor rho against ABSOLUTE (|difference| < 1e-4). If they do not, Q3 is not reported.
- **Reading:**
  - "The lymphoid estimate does not track methylation" **holds** if the Lymphoid unit's truth
    rho against GIMiCC is below 0.40 for both unmix and the methods' median, in both cohorts
    (0.40 is E2's own "tracks" bar, H2).
  - It **fails** if either reaches 0.40 in either cohort.
- H1 is recomputed with GIMiCC's Lymphoid truth, using the same exact one-sided test, and reported.
  It does not replace the registered H1.

## Addendum 4 -- D23's positive control, re-run with GIMiCC's per-sample split (written 2026-10-03 13:17:11 EDT -- the file's modification time after appending, copied from stat, before the analysis stage ran; neither GIMiCC output JSON existed and no GIMiCC estimate had been opened)

D23 found that EpiDISH's per-sample T/(T+B) **fails a reference-free positive control**:
- the control is a T-minus-B marker index from the same bulk RNA;
- GBM rho +0.138, p 0.27, n 66; LGG rho -0.123, p 0.004, n 530;
- so methylation was usable only as a cohort-level truth.

GIMiCC removes brain lineages before splitting the lymphocytes, so its per-sample split may carry
signal EpiDISH's does not. The same control is run on it, unchanged:
- `scripts/lymphoid_tracking.py`'s `marker_index` and `perm_spearman`;
- panels fixed 2026-09-30: T = CD3D CD3E CD3G CD5 CD6 TRAT1; B = CD19 MS4A1 CD79A CD79B CD22 BLK
  PAX5 FCRLA;
- 10,000 permutations, seed `config.RANDOM_SEED`.

**The GIMiCC quantity:**
- per-sample T/(T+B), level 3 (Layer 3A's T, by Addendum 2), with level 4 beside it;
- samples with a NaN, or with T + B = 0, are excluded and counted;
- the PTPRC-tertile stratification is reported as in the original (exploratory).

**Reading (Q4):**
- **PASSES** if rho > 0 and permutation p < 0.05 in both cohorts.
- **FAILS** otherwise.
- If it passes, the registered and extension methods' per-sample tracking of GIMiCC's T/(T+B) is
  reported, with the same permutation test, as an informative per-sample comparison.
- If it fails, that tracking is reported as INCONCLUSIVE, as D23 did for EpiDISH.

## Addendum 5 -- an execution defect, the registered reading, and diagnostics of LGG's B > T (written 2026-10-03 13:36:28 EDT -- the file's modification time after appending, copied from stat; the shuffled run's lymphoid output, Q2-Q4 and every diagnostic below had not been computed)

**Execution defect, fixed; no rule changed.** The first analysis run (13:34:00) matched 0 samples in
control 2 (tumour layer against ABSOLUTE) and in control 4 (shuffled CpGs).
- Cause: ABSOLUTE's sample IDs carry the vial letter (`TCGA-xx-xxxx-01A-...`); the methylation IDs
  do not (`TCGA-xx-xxxx-01`). The key function kept the letter.
- An n = 0 control was scored as a failure, which made that run's reading INCONCLUSIVE for no
  scientific reason.
- Fix: both sides keyed without the vial letter, as every other methylation comparison here is. A
  guard now blocks the analysis if any control matches fewer than 100 samples.
- The defective run's output is kept (`results/gimicc_truth_confirmation.superseded_keybug_133400.json`).
- **Its Q1 values were visible before the fix.** The fix touches only sample matching. Q1 does not
  use purity, and its values are identical in both runs.

**The registered reading, after the fix: INCONCLUSIVE.**
- Passing controls: tumour layer (GBM 0.935, LGG 0.627); immune total (0.959, 0.921); structural
  invariance; probe-drop robustness (GBM T > B in 18/20, LGG B > T in 19/20).
- **Failing control:** LGG's negative control. Shuffled CpG labels give a tumour layer that still
  tracks ABSOLUTE at rho 0.646 (GBM -0.128).
- By the rule, a failing control makes the reading INCONCLUSIVE, never CONFIRMED, and the control is
  not retuned.
- Reported as measured: GBM T > B > NK (0.462 / 0.277 / 0.261); **LGG B > T > NK (0.493 / 0.405 /
  0.102)**. EpiDISH gives T > NK > B in both.

**Diagnostics, fixed now (exploratory; they explain, they never change the reading):**
- **D-a, is LGG's B > T carried by lymphoid-specific CpGs?** The lymphoid split of the shuffled-CpG
  runs (already computed, not yet opened), by Q1's computation.
  - If the shuffled run also gives B > T in LGG, the ordering does not depend on which CpGs are
    which: it is a property of the sample's overall methylation distribution, not of lymphoid signal.
  - If shuffling removes or reverses it, the B > T is carried by the library's lymphoid CpGs.
- **D-b, global methylation.** Per LGG sample, Spearman of GIMiCC's level-3 B share within T + B
  against the sample's mean beta over the complete-case library CpGs (a CIMP-like summary). Same for
  GBM.
- **D-c, why the negative control fails in LGG.** Spearman of each sample's mean beta over the
  library CpGs against ABSOLUTE purity, per cohort. If the global summary tracks purity in LGG and not
  in GBM, label shuffling cannot break the tumour layer in LGG, and the control is uninformative
  there by construction.
- **D-d, the RNA side of the same samples.** Q4's marker index is already registered (Addendum 4).
  Reported beside D-a to D-c: whether LGG samples GIMiCC calls B-dominant are B-dominant by the
  reference-free RNA index too.

## Addendum 6 -- two more diagnostics (written 2026-10-03 13:40:27 EDT -- the file's modification time after appending, copied from stat, after Addendum 5's results, before these were computed)

Addendum 5's results:
- Shuffling CpG labels **reverses** the lymphoid ordering in both cohorts: GBM T > B > NK becomes
  B > NK > T, and LGG B > T > NK becomes T > NK > B. Both orderings are carried by which CpGs are
  which.
- The LGG negative control fails because GIMiCC's IDH-mutant purity libraries are 98-99.7%
  tumour-hypermethylated CpGs (GBM: 465 of 1,000), so a shuffled LGG tumour layer reads global
  methylation (rho 0.922 with mean beta), and global methylation tracks purity.
- GIMiCC's per-sample dominance agrees weakly with the RNA marker index in both cohorts.

Two questions follow.

**D-e, is GIMiCC's B an overflow for tumour signal in LGG?** Layer 3A has no tumour column. In the RNA
arm, B already behaves as an overflow for unmodelled tumour signal (WHY_B_OVER_T §7h). If the same
happens in methylation, GIMiCC's within-lymphoid B share should **rise with tumour content**.
- Computed per cohort: Spearman of GIMiCC's level-3 B/(T+B) against ABSOLUTE purity, and against
  GIMiCC's own Tumor.
- **Reading:** consistent with overflow in a cohort if rho > 0 and p < 0.05 against ABSOLUTE there.
  Not consistent otherwise.
- Either way it is a diagnostic. It does not decide the truth.

**D-f, what does Q3's tracking reflect?** Q3 found `unmix`'s lymphoid total tracking GIMiCC's lymphoid
tissue fraction (rho 0.647 GBM, 0.526 LGG). A lymphoid tissue fraction rises with the whole leukocyte
compartment.
- Computed: the partial Spearman of `unmix` lymphoid against GIMiCC lymphoid, controlling for the
  Thorsson leukocyte fraction (rank residuals); the same for the methods' median.
- **Reading:** if the partial rho falls below 0.20 in a cohort, the tracking there is the leukocyte
  level, not lymphoid-specific.

## Note (not a rule change) -- a second sample-key defect, D27 (written 2026-10-03 13:59:52 EDT, the file's modification time after appending, copied from stat)

`tests/test_gimicc_independent.py` recomputes the controls through ABSOLUTE's own `array` column and
matched more samples than the analysis did.
- Cause: 795 ABSOLUTE rows carry a Broad-internal ID in `sample` (for example
  `GBM-TCGA-06-5416-Tumor-SM-1QETM`), which no key function can parse.
- The analysis now keys purity by `array`, the TCGA barcode in exactly the methylation IDs' form.
- Controls 2 and 4 move from n 142 / 510 to **146 / 519**:
  - tumour layer: GBM 0.935, LGG 0.627 -> 0.630;
  - shuffled: GBM -0.128 -> -0.116, LGG 0.646 -> 0.648.
- **No reading changes.** The values quoted in Addendum 5 were correct for the join used then; the
  superseded output is `results/gimicc_truth_confirmation.superseded_samplekey_135000.json`.
- The same defect drops 1 GBM and 9 LGG RNA samples from the registered ABSOLUTE yardstick
  (OPEN_DEFECTS D27). Its statistics are reproduced exactly, and the bound on the effect leaves the
  conclusion unchanged.
