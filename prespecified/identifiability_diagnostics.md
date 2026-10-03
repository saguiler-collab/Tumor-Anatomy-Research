# Extension E2 -- truth-free identifiability diagnostics: binding rules

**Written 2026-10-02, before any computation below.** Rationale, literature and the full list of
what has already been seen: `docs/EXTENSION_IDENTIFIABILITY.md` §4. Exploratory, post-registration;
selects nothing.

## Inputs (fixed)

- **D1 estimates:** genuine DESeq2 `unmix` on the TCGA frozen arm (`scripts/extension_tcga.py`'s
  inputs, unchanged), GBM and LGG. Seven settings:
  - the pre-declared shift (per-cohort rule, `prespecified/deseq2_unmix_config.md`);
  - forced shift 1, 10, 100, 1000 and 10000 at power 1;
  - the pre-declared shift at power 2.

  Per-sample cell fractions saved per setting.
- **D2 estimates:** `results/estimates_full.csv` (GBM) and `results/estimates_full_lgg.csv` (LGG),
  frozen arm. `music` (= NNLS) and `scdc_ensemble` (= SCDC) are dropped as degenerate duplicates.
- **Truths:**
  - ABSOLUTE purity (`data/raw/tcga/...mastercalls...`, 'Cancer DNA fraction', as the registered
    yardstick);
  - Thorsson leukocyte fraction (`Immune_Fraction/TCGA_all_leuk_estimate.masked.20170107.tsv`);
  - EpiDISH per-type fractions (`results/methylation_celltypes{,_lgg}.csv`).

  Samples are joined on 15-character barcodes as in the registered scripts.

## Compartments

Each is a sum of cell-fraction columns:
- Tumor;
- Leukocytes = Macrophage_Microglia + T_cell + NK_cell + B_cell;
- Myeloid = Macrophage_Microglia;
- Lymphoid = T_cell + NK_cell + B_cell;
- T_cell, B_cell, NK_cell;
- Endothelial; Oligodendrocyte.

Astrocyte is excluded (unidentifiable sidecar).

**Truth for a compartment:**
- Tumor: purity;
- Leukocytes: the leukocyte fraction;
- Lymphoid: EpiDISH (CD4T + CD8T + NK + B), with the EpiDISH column names mapped exactly as
  `scripts/methylation_celltypes.py` maps them;
- T, B, NK separately: the EpiDISH type (secondary);
- the others: none.

## Statistics (per compartment x cohort)

- **D1 stability** = mean of the 21 pairwise Spearman correlations of per-sample estimates between
  the 7 settings.
- **D2 agreement** = mean pairwise Spearman correlation between the remaining methods.
- **Truth agreement** = Spearman rho of estimate against truth across matched samples:
  - D1: `unmix` at the pre-declared setting;
  - D2: the median rho over the remaining methods.
- **Not estimable:** a Spearman correlation needs variation. A setting or method returning a
  constant (for example all zero) for a compartment is excluded from that compartment's pairs, and
  the number excluded is reported. A compartment with fewer than 3 usable settings (D1) or methods
  (D2) is "not estimable".
- Zero estimates are reported per compartment (the share of samples with exactly 0), because ties
  shape a rank correlation.

## Hypotheses and reading rules (fixed now)

- **H1 (primary).** Across the 6 primary units -- Tumor, Leukocytes and Lymphoid, each in GBM and
  LGG -- truth-free stability ranks the units as truth agreement does. Tested separately for D1 and
  D2: Spearman across units, with an exact one-sided permutation p over 6! = 720 orderings.
  - **SUPPORTED** if, for D1 or D2, rho > 0 and p < 0.05 (with 6 units this needs rho >= 0.83).
  - **NOT SUPPORTED** if rho <= 0 for both.
  - **INCONCLUSIVE** otherwise.
- **H1s (secondary).** The same with the 12 units that add T, B and NK in both cohorts (exact
  permutation, Monte Carlo 100,000).
- **H2 (a specific prediction).** All-leukocyte content is identifiable even though its lymphoid split
  is not: D1 stability >= 0.8 **and** `unmix` truth agreement with the leukocyte fraction >= 0.4, in
  both cohorts. Each cohort is reported pass or fail.
- **S1 (data scale versus loss scale).** NNLS on log1p(bulk) against log1p(signature), with
  per-sample proportions normalised to sum to 1 (this project's code, labelled), on the frozen arm.
  Prediction (from Avila Cobos 2020): its purity rho is lower than the registered linear NNLS's
  in both cohorts. Each cohort is reported pass or fail.
- **Controls:**
  - `unmix`'s Tumor rho at the pre-declared setting must reproduce the reported rows (GBM 0.7801,
    LGG 0.5916);
  - the registered NNLS purity rho computed here from `estimates_full*.csv` must reproduce the
    registered yardstick's value.

  Either failing aborts.

**Power.** Six units allow only a strong association to reach significance. An absent signal is
INCONCLUSIVE, never "no relation".

---

## Addendum 1 -- robustness of H1 (written 2026-10-02 after H1 was seen -- E2's output is dated 08:47:00 -- and before any of this was computed: the first robustness output is dated 09:00:36. File times; clock times stated in an earlier draft of this addendum were estimates and are removed)

H1 came out SUPPORTED on D1 (rho 0.886, exact p 0.0167). The two lymphoid units are the two least
stable, while the other four units' D1 values lie within 0.033 of each other (0.927 to 0.960). The
registered reading stands as recorded. The checks below **qualify** it and change nothing; they are
post hoc.

1. **Contrast versus fine ordering (exact, no data).** Under the registered null, the probability
   that the two units lowest on truth are also the two lowest on D1 is 1/15 = 0.067. This is reported
   next to p = 0.0167, so a reader can see how much the order inside the top four contributes.
2. **Near-tie sensitivity (exact).** The GBM Tumor and GBM Leukocytes truth values differ by 0.0003.
   H1 is recomputed with the two swapped, and the exact p is reported.
3. **Sampling robustness (bootstrap).** Resample each cohort's samples with replacement, B = 2000,
   seed `config.RANDOM_SEED`. Recompute D1, D2 and both truth agreements for the 6 primary units,
   then the two H1 Spearman correlations. Report:
   - the 2.5th, 50th and 97.5th percentiles of each H1 rho;
   - the share of replicates in which the two lymphoid units are the two lowest on D1 (and on D2).

   **Reading:**
   - **contrast robust** if that share is >= 0.95;
   - **fine ordering not robust** if the 95% interval of the H1 rho reaches 0.
4. **H2 wording check.** Each registered method's all-leukocyte truth agreement is listed, so the
   text can say whether `unmix` exceeds all of them or only their median.
5. **Duplicate setting (a design flaw, found from file sizes and confirmed by sha256 before the first robustness run, whose output is
   dated 09:00:36, i.e. before any D1 value was recomputed).**
   - In GBM the pre-declared shift is 1, so the settings `predeclared` and `s1` are the same fit:
     `results/identifiability/unmix_gbm_predeclared.csv` and `unmix_gbm_s1.csv` are byte-identical
     (sha256 837911747e5e79ae...).
   - GBM's registered D1 therefore averages 21 pairs over 6 distinct fits. One pair is exactly 1, and
     five pairs are counted twice.
   - LGG's pre-declared shift is 0.5, so its 7 fits are distinct. The registered D1 favours the GBM
     units in a comparison across cohorts.

   **Check:** recompute D1 over distinct fits only (GBM 6, LGG 7) and re-test H1 and H1s. Report it
   next to the registered values.

   **Reading:** if H1's D1 reading changes, the paper reports the corrected value as the one to
   believe and the registered value as superseded by a design error. It does not retune or replace
   the rule.
6. **Degraded methods in D2 (a second design inconsistency, found while cross-checking check 4 against
   `results/immune_arm.json`, between the first robustness output (09:00:36) and the run that computed
   this check (started 09:03:34) -- before any D2 value was recomputed without them).**
   - E2's D2 dropped MuSiC and SCDC ENSEMBLE but kept Bisque (no subjects shared between bulk and
     single cells) and EPIC (no per-gene reference variance).
   - Both run in a degraded mode here. The registered immune arm excluded them, with MuSiC and SCDC
     ENSEMBLE, from its comparison (its `comparable_methods`, 8 methods).

   **Check:** recompute D2 and the median-of-methods truth agreement on that comparable panel, read
   from the artefact rather than typed. Re-test H1 and H1s on D2 with the same exact and Monte Carlo
   permutation tests, and report the result next to the registered values.

   **Reading:**
   - H1's registered reading is decided by D1 and does not change.
   - If D2's own result changes, both versions are reported; the comparable panel is labelled post
     hoc and named as the cleaner panel.
   - Bisque and EPIC are never listed under their own names without their degraded-mode label.
7. **Denominator mismatch in the EpiDISH units (a third design flaw, found while reading DREAM's
   detection limits -- its text was fetched at 09:21:26 -- and before the script that computes it was
   created (09:24:44), so before any matched value was computed).**
   - `scripts/methylation_celltypes.py` runs EpiDISH with the blood reference
     (`centDHSbloodDMC.m`, RPC), so its seven columns sum to 1. CD4T + CD8T + NK + B is therefore
     lymphocytes' share of the immune (blood-type) compartment.
   - E2 compared it with deconvolution's lymphoid share of all cells: T + NK + B over every cell type.
     Across samples the second quantity also carries the overall leukocyte level. Even a perfect
     estimate would agree only partly, so the Lymphoid, T, B and NK truth agreements are attenuated by
     design.

   **Check:** recompute each EpiDISH unit with matched denominators, within the immune compartment:
   - deconvolution: (T + NK + B) / (Macrophage_Microglia + T + NK + B), and T, B, NK each over the same
     leukocyte total;
   - EpiDISH: the corresponding share, unchanged;
   - D1 stability: for the same ratio across the 7 settings (an all-zero leukocyte total leaves the
     ratio undefined; such samples are dropped and counted).

   Then re-test H1 (Tumor, Leukocytes, matched Lymphoid) and H1s, exact and Monte Carlo as registered.

   **Reading:**
   - If the matched Lymphoid truth agreement of `unmix` is >= 0.30 in either cohort, the statement
     "the lymphoid estimate does not track methylation" is withdrawn as a denominator artefact, and
     the matched H1 is reported as the one to believe.
   - If it stays below 0.30 in both, that statement stands on matched denominators.
   - Either way the registered values are reported beside the matched ones, with the reason.

---

## Addendum 2 -- S2, the anatomy arm (written 2026-10-02 12:55 -- this file's last-modified time -- before the S2 code change (12:56:04) and run (12:56-13:03:59), so before any of it was computed)

S2 was described in `docs/EXTENSION_IDENTIFIABILITY.md` §3 but never entered these rules, and was not
run. It is declared here, post hoc and exploratory, after E2's TCGA results were seen.

- **Inputs.** `scripts/remeasure_method.py`'s exact anatomic inputs for `deseq2_unmix`, built once:
  - Ivy GAP anatomic samples, the leaderboard's gene space, training donors only, raw/X, all eight
    equivalence conditions;
  - DESeq2 `unmix` at E2's seven settings, through the same environment overrides.
- **Outputs** go to `results/identifiability/` only. The reported `unmix` anatomic artefacts are never
  written.
- **Control.** The pre-declared setting must reproduce the reported row: ACS 0.9077, and the estimates
  in `results/extension/ivygap_deseq2_unmix.csv` to max |diff| < 1e-9. Failing aborts.
- **Duplicates.** The reported run's flatness rule chose shift 1, so `predeclared` and `s1` are expected
  to be the same fit. Settings whose estimates are identical are collapsed for the stability statistic
  (as in Addendum 1, check 5); ACS is reported for all seven.
- **Statistics.**
  - Per setting: ACS, its permutation null p, and each constraint's fraction satisfied.
  - Across settings: the ACS range.
  - Per cell type: the stability of per-sample estimates (mean pairwise Spearman over distinct fits)
    for the four constrained cell types and for T, B and NK.
- **Prediction: anatomy does not see the loss-scale non-identifiability.**
  - **SUPPORTED** if the ACS range across the seven settings is <= 0.05 and every setting beats its
    null (p < 0.05);
  - **NOT SUPPORTED** if the range exceeds 0.10 or any setting fails its null;
  - **INCONCLUSIVE** otherwise.
- **Descriptive, untested (seven cell types):** whether the four constrained cell types are more stable
  on Ivy GAP than T, B and NK, as E2 found on TCGA.
