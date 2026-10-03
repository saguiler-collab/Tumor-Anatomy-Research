# Does agreement between methods pick the accurate estimate? DECEPTICON's criterion against DNA truth

**Written 2026-10-02 09:32 (file time; an earlier draft said 11:05, an estimate, wrong), before any computation below.** Post-registration and exploratory. It
selects nothing for this study and changes no registered statistic. It is part of Extension E2
(`docs/EXTENSION_IDENTIFIABILITY.md`).

## Why

DECEPTICON (Deng et al., *Brief. Bioinform.* 2025;26:bbaf234, doi:10.1093/bib/bbaf234; full text read
via Europe PMC, PMC12107245) picks deconvolution estimates **without ground truth** by agreement
between methods.

- Its premise, quoted: "the most effective analysis strategies for deconvolution demonstrate
  significant concordance".
- Its rule, quoted: "Pearson correlation coefficients are computed between these prediction
  outcomes, and the two pairs of methods with the highest correlation for each cell type are
  selected as the optimal analytical strategy. ... the selected methods are initially assigned
  equal weights (0.25). If a method is selected multiple times, its weight is proportionally
  increased. The final estimate for each cell type is derived from the weighted sum of these
  methods."
- It was validated on simulated mixtures, on flow-cytometry truth in blood and non-brain tumours,
  and in TCGA through prognosis, **not against a per-sample orthogonal molecular truth in tumours.**
  This study has such truths in two glioma cohorts.

E2 already found that agreement *across compartments* points the right way but is weaker than
loss-scale stability (D2). This test asks the method-level question DECEPTICON's rule depends on.

## Inputs (fixed)

- **Panel (primary):** the registered immune arm's `comparable_methods`, read from
  `results/immune_arm*.json` (8 methods; the degraded-mode methods left out, as in E2's Addendum 1,
  check 6).
- **Panel (secondary):** the 10 non-duplicate registered methods (adds Bisque and EPIC, labelled
  degraded).
- **Estimates:** frozen arm, `results/estimates_full{,_lgg}.csv`, as E2.
- **Truths and units:** exactly E2's, with check 7's matched denominators.
  - Primary units:
    - Tumor (tissue fraction vs ABSOLUTE purity);
    - Leukocytes (Macrophage_Microglia + T + NK + B vs the Thorsson leukocyte fraction);
    - Lymphoid share of leukocytes vs EpiDISH (CD4T + CD8T + NK + B).
  - Secondary units: T, B and NK shares of leukocytes vs EpiDISH.
  - Each unit in GBM and LGG.

## The rule as implemented (this project's code; NOT the DECEPTICON package)

- **Per cell type of the frozen signature** (Tumor, Macrophage_Microglia, T_cell, NK_cell, B_cell,
  and the rest):
  - compute Pearson correlations between every pair of methods' per-sample estimates;
  - select the two pairs with the highest correlation;
  - weight each selected method 0.25 per selection, so a method selected twice gets 0.5;
  - the cell type's ensemble estimate is the weighted sum.
- **Zero handling:** a method whose estimate for a cell type is constant (for example all zero) is
  excluded for that cell type. This follows DECEPTICON's rule 5 against "artificial correlations
  arising from the simultaneous presence of zero values".
- **No step-4 normalisation**: every method here already returns proportions.
- **Compartments** are formed from the ensemble's cell types exactly as from any method's, and
  scored against the truths as in E2. A sample whose ensemble leukocyte total is 0 has no lymphoid
  share, and is dropped and counted.

## Tests and reading rules (fixed now)

### P1 · Does the agreement-selected ensemble beat chance selection?

- **Null:** for each cell type independently, select two distinct pairs uniformly at random from all
  method pairs. Build the ensemble the same way, 2000 draws, seed `config.RANDOM_SEED`.
- **Per unit:** the DECEPTICON ensemble's truth Spearman, its percentile in the null, and whether it
  beats the null median.

**Reading over the 6 primary units:**
- **SUPPORTED** if it beats the null median in all 6 (one-sided sign test p = 0.016);
- **NOT SUPPORTED** if it beats it in 3 or fewer;
- **INCONCLUSIVE** at 4 or 5.

### P2 · Do the methods that agree most with the others agree most with the truth?

- **Per unit:** each method's centrality is its mean Pearson correlation with the other methods, for
  that unit's compartment quantity. Compute Spearman across methods between centrality and truth
  agreement.
- **Statistic:** the mean of that Spearman over the 6 primary units.
- **Null:** shuffle truth agreements among methods within each unit, 10,000 permutations, seed
  `config.RANDOM_SEED`.

**Reading:**
- **SUPPORTED** if the mean is > 0 and p < 0.05;
- **NOT SUPPORTED** if the mean is <= 0;
- **INCONCLUSIVE** otherwise.

### Reported alongside, descriptive only

- The median single-method truth agreement, and the best single method's, with its implementation
  label.
- The secondary units and the 10-method panel.

**Power.** Six units and eight methods: only a consistent pattern can reach significance. An absent
signal is INCONCLUSIVE, never "agreement does not matter".
