# Results — the Anatomy Test, first valid full run

**Run:** 2026-09-04 · **Cohort:** Ivy GAP, 270 archive samples deconvolved, ACS scored on
the 122 H&E anatomic samples / 10 tumours (9 evaluable).
**Reference:** GBmap Core, 15,311 cells across
110 donors, 16,758 genes.
**Constraint freeze hash:** `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`

Read [Anatomy_Test.md](Anatomy_Test.md) first — it is the protocol. Every table below is
rendered from the artefacts by `scripts/summarize_results.py`; nothing is transcribed by
hand.

---

## 0. The headline

**The primary result meets its pre-registered bar.**

> Spearman **rho = 0.873** between the ACS ranking and the ranking from real
> ground truth, bootstrap CI **[0.396, 0.990]**, p = 0.0010,
> on 10 methods (9 distinct).
> Bar fixed in advance: rho >= 0.60 **and** a CI excluding zero. Both hold.

This is the protocol's first branch: *anatomic concordance is usable to choose a
deconvolution method in tissue with an anatomical atlas, using no ground truth and no
patient outcomes.*

Four supporting facts, each of which had to hold for that sentence to mean anything:

1. **The controls behave.** The hardest negative control sits at
   0.386 (sd 0.066) over 200 draws, and every real method —
   0.769 at worst —
   clears its 95th percentile of 0.492. Neither control beats its own
   permutation null.
2. **The biology is right.** Tumour fraction rises 0.33 → 0.65 → 0.83 across leading
   edge → infiltrating → cellular tumour; endothelium reaches 0.64 in microvascular
   proliferation; oligodendrocyte peaks at the leading edge. Nobody told the methods
   where the structures were.
3. **All four published packages ran as themselves** — `R:MuSiC`, `R:DWLS`,
   `R:BisqueRNA`, `R:SCDC` — not as this project's reimplementations.
4. **Nothing selected on the thing being tested.** The method was frozen by an
   outcome-blind, ACS-blind pseudobulk benchmark (88 training / 22 held-out donors,
   200 mixtures) before the anatomy was scored. That the benchmark and the anatomy
   independently agree on MuSiC is the result, not the method.

### What this is not

- **Not pre-registered.** The constraint file is hashed and committed, but no public
  registration receipt exists, so `registration_status.json` reports UNREGISTERED. The
  hash proves the constraints have not changed; it does not prove they predate the
  results. Until that is fixed this is a strong internal result, not a pre-registered one.
- **Not a prognostic finding.** See §7: 29 events against the ~80 needed.
- **Not transferable yet.** One tissue is a case study; the protocol's step 7 is untouched.

---


## 1a. Integrity checks

**Pre-registration:** **UNREGISTERED**. UNREGISTERED: the constraint file is hashed and committed, but no public registration receipt exists. The hash proves the constraints have not changed; it does not prove they were written before the results, which is what the pre-registration claim rests on. Register 2d1fb47c98832adf... publicly (OSF takes minutes) and record the receipt in REGISTRATION.json. Until then this project must not describe itself as pre-registered.

**Reference coverage:** the roster maps 12 of the atlas's labels and drops the rest — **23,864 cells (7.0%)**.

| dropped population | cells |
|---|---|
| Mono | 14,215 |
| DC | 3,961 |
| RG | 2,807 |
| Mural cell | 1,418 |
| Plasma B | 572 |
| OPC | 496 |
| Mast | 373 |
| Neuron | 22 |

Populations the atlas labels but the roster drops do not disappear from the tissue, only from the model. A solver with eight columns must still explain the expression those cells contribute, so it is absorbed by whichever retained type is closest. Largest dropped populations: Mono (14,215 cells), DC (3,961 cells), RG (2,807 cells) — 7.0% of the atlas in total. The dominant one is myeloid (DC, Mast, Mono), and the nearest retained column is Macrophage_Microglia — which carries C5 and C6. So the absorbed signal is concentrated on two of the seven constraints rather than spread evenly. Note the frozen constraint file lists neuronal content as untestable for want of a neuron column. That holds, but not for the reason it gives: this atlas contains only 22 neurons, so a neuron column could not have been estimated from it either. The limitation is the reference, not just the roster.

**Method configuration:** every method's parameters are recorded in `method_configs.json`. Declared departures from published defaults: **dwls**.

- `dwls` — buildSignatureMatrixMAST(diff.cutoff, pval.cutoff): diff.cutoff = 0.5, pval.cutoff = 0.01 → diff.cutoff = 0, pval.cutoff = 1 (filter fully open). Equal footing requires every method to receive the same pre-filtered gene space, so the informative genes have already been chosen identically for all methods. DWLS's own differential-expression pass is then a second filter, not a second opinion: at the published cutoffs it left 31 genes for 8 cell types and the driver refused to deconvolve on them. *Decided from the gene count, before any DWLS score existed.*
- `dwls` — cells used for the signature build: all cells in the reference → <= 250 cells per cell type. buildSignatureMatrixMAST costs ~30 minutes per gene set on this hardware regardless of the DE method. A signature is a per-cell-type summary, so type proportions are irrelevant to it and capping per type spends the budget where it buys accuracy. *Decided from measured runtime, before any DWLS score existed.*
- `dwls` — per-sample failure handling: an error aborts the call → a sample that fails to solve becomes NA and is counted. solveDampenedWLS fails on particular mixtures with 'NA/NaN argument'. One such sample was sending the whole cohort to the Python fallback. NA is already what this pipeline means by a failed sample and it is reported in n_failed_samples. *Decided from the failure mode, before any DWLS score existed.*

**What actually ran:** 5 published package(s) ran as the genuine R implementation (music, dwls, bisque, scdc, scdc_ensemble).

## 2. Cohort and provenance

- constraint freeze hash: `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`
- ACS cohort: **122 samples / 10 tumours**, 10,000 within-tumour permutations per method
- structures: {'CT': 30, 'MVP': 25, 'PAN': 24, 'IT': 24, 'LE': 19}
- portal reconciliation: RECONCILED: 122 anatomic samples across 5 structures agree exactly with the portal's own table ({'CT': 30, 'IT': 24, 'LE': 19, 'MVP': 25, 'PAN': 24}), and the study assignment agrees on all 270 shared samples. 9 sample(s) exist on the live portal but not in the 2014-11-25 archive — a release skew, all in the 'Cancer Stem Cells RNA Seq' study, not a parsing error. Step 3's gate is closed.

## 3. The ACS leaderboard

_from `results/anatomic/acs_leaderboard.csv`_

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| music | 1.000 | 1.000 – 1.000 | 0.370 | 0.000 | 9 |  |
| nnls | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| elastic_net | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| svr | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| scdc_ensemble | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| scdc | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| bisque | 0.923 | 0.866 – 0.983 | 0.380 | 0.000 | 9 |  |
| dwls | 0.877 | 0.742 – 0.985 | 0.372 | 0.000 | 9 |  |
| bayesian | 0.769 | 0.627 – 0.921 | 0.380 | 0.000 | 9 |  |
| bayesian_hierarchical | 0.769 | 0.627 – 0.921 | 0.381 | 0.000 | 9 |  |
| **control_random** | 0.400 | 0.288 – 0.530 | 0.382 | 0.431 | 9 | **yes** |
| **control_shuffled_signature** | 0.138 | 0.059 – 0.228 | 0.349 | 0.999 | 9 | **yes** |


**Control verdict.** CONTROLS BEHAVE: best control 0.400 sits below the median real method 0.954, ties no real method, and does not beat its own permutation null. The constraint set discriminates.


### What ACS can and cannot resolve here

- weighted denominator: **65** (constraint x tumour pairs, weighted)
- so ACS takes at most **66 distinct values**, spaced **0.0154** apart
- 12 methods scored produce **8 distinct values**; the 10 real methods produce **6**
- tied at **0.985**: nnls, elastic_net, svr
- tied at **0.954**: scdc_ensemble, scdc
- tied at **0.769**: bayesian, bayesian_hierarchical


### Control audit (stricter than the headline verdict)

_The control is a **single permutation** drawn from `config.RANDOM_SEED`. Everything in this block describes that one draw._

- best control: **control_random = 0.400**
- worst real method: bayesian = 0.769
- real methods the best control ties or beats: none
- controls beating their own permutation null: none

| constraint | weight | evaluable | best control satisfies |
|---|---|---|---|
| C1 | 1.0 | 8 | 5 |
| C2 | 1.0 | 8 | 3 |
| C3 | 1.0 | 9 | 7 |
| C4 | 1.0 | 9 | 4 |
| C5 | 1.0 | 6 | 2 |
| C6 | 1.0 | 9 | 5 |
| C7 | 2.0 | 8 | 0 |

**Both controls calibrated over 200 independent draws each.** Each leaderboard control row is one draw from these.

| control | mean | sd | 5–95% | the leaderboard's draw |
|---|---|---|---|---|
| `control_shuffled_signature` | 0.386 | 0.066 | 0.292 – 0.492 | 0.138 (0% pct) |
| `control_random` | 0.385 | 0.072 | 0.262 – 0.492 | 0.400 (63% pct) |

Methods are judged against the **hardest** control (`control_shuffled_signature`, 95th percentile **0.492**), not a convenient one.

> Every real method scores above the hardest control's 95th percentile.

**What the hardest control satisfies, averaged over 200 draws** — this is the number to read, not the single draw above:

| constraint | satisfied rate |
|---|---|
| C3 | 0.52 |
| C1 | 0.52 |
| C5 | 0.52 |
| C6 | 0.50 |
| C2 | 0.46 |
| C4 | 0.23 |
| C7 | 0.18 |

### Per constraint — best method (`music`)

| ID | claim | weight | tumours evaluable | satisfied |
|---|---|---|---|---|
| C1 | Tumor: CT > LE | 1.0 | 8 | 8 (1.00) |
| C2 | Oligodendrocyte: LE > CT | 1.0 | 8 | 8 (1.00) |
| C3 | Endothelial: MVP > CT | 1.0 | 9 | 9 (1.00) |
| C4 | Endothelial: MVP is the maximum | 1.0 | 9 | 9 (1.00) |
| C5 | Macrophage_Microglia: PAN > LE | 1.0 | 6 | 6 (1.00) |
| C6 | Macrophage_Microglia: MVP > CT | 1.0 | 9 | 9 (1.00) |
| C7 | Tumor: LE < IT < CT | 2.0 | 8 | 8 (1.00) |


### Per tumour — `music`

| tumour | constraints evaluable | satisfied | ACS | violated |
|---|---|---|---|---|
| 292163427 | 6 | 6 | 1.000 | — |
| 703393 | 7 | 7 | 1.000 | — |
| 703493 | 7 | 7 | 1.000 | — |
| 705757 | 7 | 7 | 1.000 | — |
| 705758 | 7 | 7 | 1.000 | — |
| 705803 | 7 | 7 | 1.000 | — |
| 705859 | 7 | 7 | 1.000 | — |
| 711547 | 6 | 6 | 1.000 | — |
| 711560 | 3 | 3 | 1.000 | — |


## 4. What actually ran

| method | implementation | degenerate | why |
|---|---|---|---|
| nnls | python | no |  |
| svr | python | no |  |
| elastic_net | python | no |  |
| bayesian | python | no |  |
| bayesian_hierarchical | python | no |  |
| music | R:MuSiC | no |  |
| dwls | R:DWLS | no |  |
| bisque | R:BisqueRNA | **yes** | no subjects are assayed both as bulk and as single cells, so this runs Bisque's documented no-overlap mode (us |
| scdc | R:SCDC | no |  |
| scdc_ensemble | R:SCDC | **yes** | one reference supplied ('gbmap'), so there is nothing to weight across: SCDC ENSEMBLE reduces exactly to SCDC. |
| control_random | python | no |  |
| control_shuffled_signature | python | no |  |


## 5. The agreement test (the primary result)

Pre-registered bar: rho >= 0.6 AND a bootstrap CI excluding zero. Minimum methods: 6.

| yardstick | rho | 95% CI | methods | distinct | verdict |
|---|---|---|---|---|---|
| synthetic_mixtures | 0.873 | 0.396 – 0.990 | 10 | 9 | HEADLINE: ACS ranking tracks true accuracy at the pre-registered bar. Anatomic concordance |
| absolute_purity | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| sc_pseudobulk | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| simulated_donor_mismatch | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |


## 6. Does deconvolving all 270 samples change the ranking?

| method | ACS (122 anatomic) | ACS (270 deconvolved) | delta | rank change |
|---|---|---|---|---|
| music | 1.000 | 1.000 | +0.000 | +0 |
| nnls | 0.985 | 0.985 | +0.000 | +0 |
| elastic_net | 0.985 | 0.985 | +0.000 | +0 |
| svr | 0.985 | 0.985 | +0.000 | +0 |
| scdc | 0.954 | 0.954 | +0.000 | +1 |
| scdc_ensemble | 0.954 | 0.954 | +0.000 | +1 |
| bisque | 0.923 | 0.969 | +0.046 | -2 |
| dwls | 0.877 | 0.738 | -0.138 | +2 |
| bayesian | 0.769 | 0.769 | +0.000 | -0 |
| bayesian_hierarchical | 0.769 | 0.785 | +0.015 | -2 |
| control_random | 0.400 | 0.338 | -0.062 | +0 |
| control_shuffled_signature | 0.138 | 0.138 | +0.000 | +0 |

Spearman rho between the two ACS rankings: **0.9554**


**Which methods moved, and why.** 8 of 12 did not move by a single unit: music, nnls, elastic_net, svr, scdc, scdc_ensemble, bayesian, control_shuffled_signature. Those solve each sample independently, so what else is in the cohort cannot reach them — the zeros are exact, not rounded.

The ones that moved are **bisque, dwls, bayesian_hierarchical, control_random**. Bisque normalises the bulk with cohort-wide per-gene statistics (`B.mean(axis=1)` and `B.std(axis=1)`), so the 148 ISH-cluster samples shift the reference frame every anatomic sample is mapped through. `control_random` moves for a different and uninteresting reason: it draws one Dirichlet sample per row, so a 270-row draw is not a superset of a 122-row draw.

Note the direction: **adding 148 more real samples made Bisque's anatomic concordance worse**, not better. For a method that borrows strength across samples, which samples it is handed is part of the method — and Ivy GAP's two studies are not the same tissue.


### The full-database run

- deconvolved **270 samples / 37 tumours**
- ACS still scored on **122 samples / 10 tumours** — H&E-selected anatomic study only; ISH-cluster samples are excluded from scoring by construction because their structure labels were assigned using expression
- control verdict: CONTROLS BEHAVE: best control 0.338 sits below the median real method 0.962, ties no real method, and does not beat its own permutation null. The constraint set discriminates.

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| music | 1.000 | 1.000 – 1.000 | 0.370 | 0.000 | 9 |  |
| nnls | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| elastic_net | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| svr | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| bisque | 0.969 | 0.930 – 1.000 | 0.379 | 0.000 | 9 |  |
| scdc_ensemble | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| scdc | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| bayesian_hierarchical | 0.785 | 0.648 – 0.925 | 0.380 | 0.000 | 9 |  |
| bayesian | 0.769 | 0.627 – 0.921 | 0.381 | 0.000 | 9 |  |
| dwls | 0.738 | 0.600 – 0.867 | 0.374 | 0.000 | 9 |  |
| **control_random** | 0.338 | 0.217 – 0.486 | 0.388 | 0.779 | 9 | **yes** |
| **control_shuffled_signature** | 0.138 | 0.059 – 0.228 | 0.349 | 0.999 | 9 | **yes** |


## 7. Survival

**Canonical verdict.** BLOCKED: tumor_details.csv is present and joins correctly, but it publishes survival_days and NO vital-status column. Who was censored is not knowable from the release, so a C-index cannot be computed. Nothing is imputed and no prognostic claim is made.

**Missingness.** INFORMATIVE MISSINGNESS: a blank survival time is strongly associated with MGMT methylation (Fisher p = 0.002134), the strongest favourable prognostic factor in GBM. The blanks are concentrated among the patients most likely to have been alive at the data freeze, so dropping them biases the cohort toward short survival. Any complete-case survival estimate on this table is biased downward and must be reported as such.

- MGMT methylated among tumours with a recorded time: 9/31; among blanks: 8/9 (Fisher p = 0.002134)
- median age: 61 recorded vs 52 blank (Mann-Whitney p = 0.0153)


### Under the declared event policy (quarantined, not canonical)

- policy: `observed-time-is-death`
- patients 29 · events 29 · events per covariate 3.62
- smallest detectable C-index difference: **0.515**
- UNDERPOWERED: 29 events supports detecting a C-index difference of about 0.515 at best. Differences smaller than that are reported as INCONCLUSIVE and are not used to rank methods.

| method | C baseline | C + composition | delta | interpretation |
|---|---|---|---|---|
| svr | 0.484 | 0.564 | 0.080 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| nnls | 0.484 | 0.561 | 0.077 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bisque | 0.484 | 0.521 | 0.037 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bayesian_hierarchical | 0.484 | 0.514 | 0.029 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| dwls | 0.484 | 0.510 | 0.026 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| elastic_net | 0.484 | 0.505 | 0.021 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bayesian | 0.484 | 0.503 | 0.019 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| music | 0.484 | 0.486 | 0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc | 0.484 | 0.482 | -0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc_ensemble | 0.484 | 0.482 | -0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |


#### Exploratory: does ACS track prognostic value?

Not a protocol analysis. The protocol's agreement test correlates ACS against *accuracy*, which is still not computable. This correlates it against *prognostic value* instead, which the declared-policy run makes available. It selects nothing.

- Spearman(ACS, delta C-index) over 10 methods: **rho = 0.031**, p = 0.932
- best ACS: `music` (ACS 1.000, delta C +0.002)
- best prognosis: `svr` (ACS 0.985, delta C +0.080)

Three independent reasons this cannot support a claim: every C-index is INCONCLUSIVE by the pre-specified power rule; ACS supplies only 6 distinct values across 10 methods, so the rank is mostly ties; and the outcome side rests on a declared assumption about censoring. It is recorded because it points the same way as the superseded run did, and because the direction is the protocol's second branch.



---

## Provenance of this document

Generated from `results/` by `scripts/summarize_results.py` and
`scripts/figure_data.py`. The figures in
[`results/figures/`](results/figures/) read the same artefacts.

Superseded runs are preserved under `results_superseded/` with the reason each was
retired. Two are worth knowing about, because both produced plausible-looking numbers:

- A run whose reference was **scrambled** — cell expression attached to the wrong
  cell's metadata. It put every real method below a random control and SCDC at exactly
  0.000. Column order, shapes and 161 tests all passed; it surfaced only as biology
  (23 of 24 marker genes in the wrong column).
- A run whose **stage-3 artefacts were overwritten by the test suite**, which called
  `run_benchmark.run()` against the real results tree. Stage 4 reads
  `benchmark_summary.csv` as the ground-truth yardstick, so the primary result would
  have been computed against a 10-mixture fixture.

Both are fixed, both have regression tests, and both are described in
`docs/METHODS.md`.
