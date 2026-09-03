# Results — full-database run, 2026-09-03

**Cohort:** Ivy GAP. 270 archive samples deconvolved, ACS scored on the 122 H&E anatomic
samples / 10 tumours (9 evaluable).
**Constraint freeze hash:** `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`
**Gene space:** 1,614 of 16,007 shared genes.

Read [Anatomy_Test.md](Anatomy_Test.md) first — it is the protocol. This file reports
what the pipeline produced and, more importantly, what it does not support.

Every table below is rendered from the artefacts by
`python scripts/summarize_results.py`. Nothing is transcribed by hand. That is
deliberate: the previous run's numbers survived a destructive accident *only* because
someone had copied them into a markdown file, and a number copied by hand is a number
that can be copied wrong.

---

## 0. The headline, in five lines

1. **ACS is at its ceiling on this instrument, and the ceiling is low.** The weighted
   denominator is 65 and the run reached all 65. Ten real methods produce **five
   distinct scores**; five tie at exactly 0.800. Nothing about how a method is run
   changes this.
2. **Every real method clears the negative controls** — but only once the controls are
   read as distributions. Against the hardest control's 95th percentile
   (0.492), all ten are above. Against the *single pre-registered draw* (0.585),
   one of them ties. Both facts are reported; §8 explains why the second is an artefact
   of the seed.
3. **Deconvolving all 270 samples does not change the ACS ranking** (Spearman 0.982).
   Only methods that borrow strength across samples move, and Bisque moves *down*.
4. **The agreement test is not computable, now for two independent reasons.** No
   yardstick covers more than two methods, *and* the ACS side supplies only five
   distinct values. Obtaining a yardstick is necessary but no longer sufficient.
5. **Prognosis is BLOCKED, and would be INCONCLUSIVE even if it were not.** Ivy GAP
   publishes no vital-status column. Under the most generous reading available the
   cohort reaches 29 events against the ~80 needed, detecting at best a C-index
   difference of 0.515 — wider than the entire interpretable range.

---

## 1. What ran

| Protocol step | Status |
|---|---|
| 1. Restore the frozen inputs | **done** — `reference_frozen/`, hashes recorded |
| 2. Write + timestamp the constraint file | **done in repo**, hash above. *Still not publicly registered (OSF).* |
| 3. Ingest Ivy GAP | **done, and the gate is closed externally** — counts reconciled sample-by-sample against the Allen Institute's own live metadata export, not against the archive describing itself |
| 4. Build + prove the scorer | **done** — known-answer fixtures pass; permutation null centred at chance |
| 5. Score every method + controls | **done** — 12 methods, 10,000 permutations each, two cohorts, both controls calibrated over 200 draws |
| 6. The agreement test | **NOT COMPUTABLE** — for two independent reasons, not one |
| 7. Prove it transfers | not started |

## 2. Cohort and provenance

- constraint freeze hash: `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`
- ACS cohort: **122 samples / 10 tumours**, 10,000 within-tumour permutations per method
- structures: {'CT': 30, 'MVP': 25, 'PAN': 24, 'IT': 24, 'LE': 19}
- portal reconciliation: RECONCILED: 122 anatomic samples across 5 structures agree exactly with the portal's own table ({'CT': 30, 'IT': 24, 'LE': 19, 'MVP': 25, 'PAN': 24}), and the study assignment agrees on all 270 shared samples. 9 sample(s) exist on the live portal but not in the 2014-11-25 archive — a release skew, all in the 'Cancer Stem Cells RNA Seq' study, not a parsing error. Step 3's gate is closed.

## 3. The ACS leaderboard

_from `results/anatomic/acs_leaderboard.csv`_

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| scdc | 0.877 | 0.761 – 1.000 | 0.373 | 0.000 | 9 |  |
| scdc_ensemble | 0.877 | 0.761 – 1.000 | 0.373 | 0.000 | 9 |  |
| nnls | 0.800 | 0.721 – 0.870 | 0.380 | 0.000 | 9 |  |
| svr | 0.800 | 0.719 – 0.883 | 0.367 | 0.000 | 9 |  |
| music | 0.800 | 0.721 – 0.870 | 0.380 | 0.000 | 9 |  |
| elastic_net | 0.800 | 0.721 – 0.870 | 0.381 | 0.000 | 9 |  |
| bisque | 0.800 | 0.726 – 0.875 | 0.383 | 0.000 | 9 |  |
| dwls | 0.738 | 0.625 – 0.868 | 0.371 | 0.000 | 9 |  |
| bayesian | 0.615 | 0.525 – 0.721 | 0.382 | 0.001 | 9 |  |
| bayesian_hierarchical | 0.585 | 0.522 – 0.677 | 0.382 | 0.003 | 9 |  |
| **control_shuffled_signature** | 0.585 | 0.482 – 0.662 | 0.372 | 0.002 | 9 | **yes** |
| **control_random** | 0.400 | 0.288 – 0.530 | 0.384 | 0.437 | 9 | **yes** |


**Control verdict.** CONTROLS PARTLY BEHAVE: the best control (0.585) is below the median real method (0.800), but it ties or beats bayesian_hierarchical and control_shuffled_signature beats its own permutation null. The constraint set separates most methods from noise and does NOT separate all of them. Reported as-is; the constraint file is NOT retuned. Note the control is a SINGLE permutation — see the control calibration artefact before reading a tie as a property of the constraints.


### What ACS can and cannot resolve here

- weighted denominator: **65** (constraint x tumour pairs, weighted)
- so ACS takes at most **66 distinct values**, spaced **0.0154** apart
- 12 methods scored produce **6 distinct values**; the 10 real methods produce **5**
- tied at **0.877**: scdc, scdc_ensemble
- tied at **0.800**: nnls, svr, music, elastic_net, bisque
- tied at **0.585**: bayesian_hierarchical, control_shuffled_signature

**This bounds the agreement test independently of any yardstick.** It needs at least 6 methods to correlate, but the ACS side supplies only **5 distinct values** across 10 real methods. Even with ground truth covering every method, a rank correlation would be computed on mostly tied ranks.


### Control audit (stricter than the headline verdict)

_The control is a **single permutation** drawn from `config.RANDOM_SEED`. Everything in this block describes that one draw._

- best control: **control_shuffled_signature = 0.585**
- worst real method: bayesian_hierarchical = 0.585
- real methods the best control **ties or beats**: **bayesian_hierarchical**
- controls that **beat their own permutation null** (p < 0.05): **control_shuffled_signature**
- constraints **this single draw** satisfies perfectly: **C1, C7** — 24 of 65 weighted units (**36.9%** of the constraint set). This is ONE permutation, not a property of the constraints; see the draw distribution below before reading anything into it.

| constraint | weight | evaluable | best control satisfies |
|---|---|---|---|
| C1 | 1.0 | 8 | 8 |
| C2 | 1.0 | 8 | 3 |
| C3 | 1.0 | 9 | 6 |
| C4 | 1.0 | 9 | 0 |
| C5 | 1.0 | 6 | 4 |
| C6 | 1.0 | 9 | 1 |
| C7 | 2.0 | 8 | 8 |

**Both controls calibrated over 200 independent draws each.** Each leaderboard control row is one draw from these.

| control | mean | sd | 5–95% | the leaderboard's draw |
|---|---|---|---|---|
| `control_shuffled_signature` | 0.373 | 0.073 | 0.277 – 0.492 | 0.585 (100% pct) |
| `control_random` | 0.389 | 0.068 | 0.277 – 0.492 | 0.400 (57% pct) |

Methods are judged against the **hardest** control (`control_shuffled_signature`, 95th percentile **0.492**), not a convenient one.

> Every real method scores above the hardest control's 95th percentile.

**What the hardest control satisfies, averaged over 200 draws** — this is the number to read, not the single draw above:

| constraint | satisfied rate |
|---|---|
| C5 | 0.51 |
| C1 | 0.50 |
| C3 | 0.48 |
| C2 | 0.48 |
| C6 | 0.47 |
| C4 | 0.23 |
| C7 | 0.17 |

### Per constraint — best method (`scdc`)

| ID | claim | weight | tumours evaluable | satisfied |
|---|---|---|---|---|
| C1 | Tumor: CT > LE | 1.0 | 8 | 8 (1.00) |
| C2 | Oligodendrocyte: LE > CT | 1.0 | 8 | 8 (1.00) |
| C3 | Endothelial: MVP > CT | 1.0 | 9 | 9 (1.00) |
| C4 | Endothelial: MVP is the maximum | 1.0 | 9 | 9 (1.00) |
| C5 | Macrophage_Microglia: PAN > LE | 1.0 | 6 | 6 (1.00) |
| C6 | Macrophage_Microglia: MVP > CT | 1.0 | 9 | 7 (0.78) |
| C7 | Tumor: LE < IT < CT | 2.0 | 8 | 5 (0.62) |


### Per tumour — `scdc`

| tumour | constraints evaluable | satisfied | ACS | violated |
|---|---|---|---|---|
| 292163427 | 6 | 6 | 1.000 | — |
| 703393 | 7 | 5 | 0.625 | C6,C7 |
| 703493 | 7 | 5 | 0.625 | C6,C7 |
| 705757 | 7 | 6 | 0.750 | C7 |
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
| music | python-reimplementation | **yes** | reference 'gbmap_frozen' carries no cross-donor variance, so MuSiC's gene weighting is constant and the result |
| dwls | python-reimplementation | no | dwls unavailable via R: R package DWLS is not installed |
| bisque | python-reimplementation | **yes** | no subjects are assayed both as bulk and as single cells, so this runs Bisque's documented no-overlap mode (us |
| scdc | python-reimplementation | no | scdc unavailable via R: R package SCDC is not installed |
| scdc_ensemble | python-reimplementation | **yes** | one reference supplied ('gbmap_frozen'), so there is nothing to weight across: SCDC ENSEMBLE reduces exactly t |
| control_random | python | no |  |
| control_shuffled_signature | python | no |  |


## 5. The agreement test (the primary result)

Pre-registered bar: rho >= 0.6 AND a bootstrap CI excluding zero. Minimum methods: 6.

| yardstick | rho | 95% CI | methods | distinct | verdict |
|---|---|---|---|---|---|
| synthetic_mixtures | — | — – — | 2 | 2 | NOT COMPUTABLE: only 2 methods have both scores. |
| absolute_purity | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| sc_pseudobulk | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| simulated_donor_mismatch | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |


## 6. Does deconvolving all 270 samples change the ranking?

| method | ACS (122 anatomic) | ACS (270 deconvolved) | delta | rank change |
|---|---|---|---|---|
| scdc | 0.877 | 0.877 | +0.000 | +0 |
| scdc_ensemble | 0.877 | 0.877 | +0.000 | +0 |
| nnls | 0.800 | 0.800 | +0.000 | -0 |
| svr | 0.800 | 0.800 | +0.000 | -0 |
| music | 0.800 | 0.800 | +0.000 | -0 |
| elastic_net | 0.800 | 0.800 | +0.000 | -0 |
| bisque | 0.800 | 0.738 | -0.062 | +2 |
| dwls | 0.738 | 0.708 | -0.031 | +0 |
| bayesian | 0.615 | 0.615 | +0.000 | +0 |
| bayesian_hierarchical | 0.585 | 0.585 | +0.000 | +0 |
| control_shuffled_signature | 0.585 | 0.585 | +0.000 | +0 |
| control_random | 0.400 | 0.338 | -0.062 | +0 |

Spearman rho between the two ACS rankings: **0.9816**


**Which methods moved, and why.** 9 of 12 did not move by a single unit: scdc, scdc_ensemble, nnls, svr, music, elastic_net, bayesian, bayesian_hierarchical, control_shuffled_signature. Those solve each sample independently, so what else is in the cohort cannot reach them — the zeros are exact, not rounded.

The ones that moved are **bisque, dwls, control_random**. Bisque normalises the bulk with cohort-wide per-gene statistics (`B.mean(axis=1)` and `B.std(axis=1)`), so the 148 ISH-cluster samples shift the reference frame every anatomic sample is mapped through. `control_random` moves for a different and uninteresting reason: it draws one Dirichlet sample per row, so a 270-row draw is not a superset of a 122-row draw.

Note the direction: **adding 148 more real samples made Bisque's anatomic concordance worse**, not better. For a method that borrows strength across samples, which samples it is handed is part of the method — and Ivy GAP's two studies are not the same tissue.


### The full-database run

- deconvolved **270 samples / 37 tumours**
- ACS still scored on **122 samples / 10 tumours** — H&E-selected anatomic study only; ISH-cluster samples are excluded from scoring by construction because their structure labels were assigned using expression
- control verdict: CONTROLS PARTLY BEHAVE: the best control (0.585) is below the median real method (0.800), but it ties or beats bayesian_hierarchical and control_shuffled_signature beats its own permutation null. The constraint set separates most methods from noise and does NOT separate all of them. Reported as-is; the constraint file is NOT retuned. Note the control is a SINGLE permutation — see the control calibration artefact before reading a tie as a property of the constraints.

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| scdc | 0.877 | 0.761 – 1.000 | 0.373 | 0.000 | 9 |  |
| scdc_ensemble | 0.877 | 0.761 – 1.000 | 0.373 | 0.000 | 9 |  |
| nnls | 0.800 | 0.721 – 0.870 | 0.380 | 0.000 | 9 |  |
| svr | 0.800 | 0.719 – 0.883 | 0.367 | 0.000 | 9 |  |
| music | 0.800 | 0.721 – 0.870 | 0.380 | 0.000 | 9 |  |
| elastic_net | 0.800 | 0.721 – 0.870 | 0.381 | 0.000 | 9 |  |
| bisque | 0.738 | 0.657 – 0.817 | 0.366 | 0.000 | 9 |  |
| dwls | 0.708 | 0.607 – 0.820 | 0.371 | 0.000 | 9 |  |
| bayesian | 0.615 | 0.525 – 0.721 | 0.381 | 0.001 | 9 |  |
| bayesian_hierarchical | 0.585 | 0.522 – 0.677 | 0.381 | 0.003 | 9 |  |
| **control_shuffled_signature** | 0.585 | 0.482 – 0.662 | 0.372 | 0.002 | 9 | **yes** |
| **control_random** | 0.338 | 0.217 – 0.486 | 0.389 | 0.788 | 9 | **yes** |


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
| bayesian_hierarchical | 0.484 | 0.684 | 0.200 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| nnls | 0.484 | 0.652 | 0.168 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| music | 0.484 | 0.652 | 0.168 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc | 0.484 | 0.644 | 0.160 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc_ensemble | 0.484 | 0.644 | 0.160 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bisque | 0.484 | 0.616 | 0.131 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bayesian | 0.484 | 0.589 | 0.105 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| svr | 0.484 | 0.585 | 0.101 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| elastic_net | 0.484 | 0.585 | 0.101 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| dwls | 0.484 | 0.437 | -0.047 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |


#### Exploratory: does ACS track prognostic value?

Not a protocol analysis. The protocol's agreement test correlates ACS against *accuracy*, which is still not computable. This correlates it against *prognostic value* instead, which the declared-policy run makes available. It selects nothing.

- Spearman(ACS, delta C-index) over 10 methods: **rho = 0.052**, p = 0.886
- best ACS: `scdc` (ACS 0.877, delta C +0.160)
- **worst ACS: `bayesian_hierarchical` (ACS 0.585, delta C +0.200) — which is also the best prognostic method**

Three independent reasons this cannot support a claim: every C-index is INCONCLUSIVE by the pre-specified power rule; ACS supplies only 5 distinct values across 10 methods, so the rank is mostly ties; and the outcome side rests on a declared assumption about censoring. It is recorded because it points the same way as the superseded run did, and because the direction is the protocol's second branch.



---

## 8. The controls are single draws, and that turns out to matter

Both negative controls are specified as one draw from `config.RANDOM_SEED`. The
leaderboard's `control_shuffled_signature` row scored **0.585** — tying
`bayesian_hierarchical` exactly and beating its own permutation null at p = 0.002.

Read alone that looks like the protocol's third branch, *"negative controls score
highly, the constraint set is too permissive."* It is not. Scored over
**200 independent draws** on the same cohort, through the pipeline's own
normalisation:

| control | mean | sd | 5–95% | full range | the leaderboard's draw |
|---|---|---|---|---|---|
| `control_shuffled_signature` | 0.373 | 0.073 | 0.277 – 0.492 | 0.185 – 0.600 | 0.585 (100% pct) |
| `control_random` | 0.389 | 0.068 | 0.277 – 0.492 | 0.200 – 0.600 | 0.400 (57% pct) |

**A tie worth naming.** Both controls have a 95th percentile of exactly 0.492. That is
not a duplicated computation — their draw sequences correlate at −0.07 and their means
differ (0.373 vs 0.389) — it is the 1/65 lattice: ACS can only take multiples of 0.0154,
so two genuinely different distributions centred near chance land on the same lattice
point. This run resolved the tie by dict order and reported
`control_shuffled_signature`; the tie-break is now explicit (p95, then mean, then name),
which would name `control_random` instead. **The threshold, and therefore every
conclusion below, is identical either way.**

The pre-registered shuffled-signature draw is the **most extreme of 200**.
The typical draw scores 0.373 — essentially its own permutation null
(0.372). Judged against the hardest control's 95th percentile (0.492), **every real
method is above it**, including `bayesian_hierarchical` at 0.585.

So the tie in the leaderboard is a property of one seed, not of the constraint set.

**What a meaningless signature satisfies, averaged over draws** — the number to read,
rather than the single draw:

| constraint | weight | satisfied rate |
|---|---|---|
| C5 | 1 | 0.51 |
| C1 | 1 | 0.50 |
| C3 | 1 | 0.48 |
| C2 | 1 | 0.48 |
| C6 | 1 | 0.47 |
| C4 | 1 | 0.23 |
| C7 | 2 | 0.17 |

No constraint is systematically satisfiable without gene-level information, and C7 — the
highest-weighted claim — is among the hardest for a meaningless signature, which is what
a three-step monotone conjunction should be.

### Two readings this replaces

**Reading A, from the single draw:** *"C1 and C7 are structurally satisfiable by a
signature carrying no gene-to-cell-type information (36.9% of the constraint weight)."*
Wrong. That draw satisfied both perfectly; across draws neither is.

**Reading B, from a standalone diagnostic script:** *"`bayesian` and
`bayesian_hierarchical` are not distinguishable from a randomly permuted signature."*
Also wrong. That script reported sd 0.156 and a 95th percentile of 0.646; the calibration
module reports 0.073 and 0.492, and was verified
draw-for-draw against an independent from-scratch implementation. The script has been
replaced by one that delegates to the module — one implementation, so there is nothing
left to disagree.

**What survives both corrections** is the methodological point, and it is the reason the
calibration exists: a control with sd 0.073, reported as a single number,
cannot establish that the controls behave in *either* direction. Whichever way the seed
lands, the reading is an artefact of it.

Nothing here changed the constraint file, the controls' definitions, or the leaderboard
rows. Replacing a point estimate with a distribution makes a control **harder** to beat,
not easier — the only direction it is safe to move an instrument after seeing what it
scored.

---

## 9. Defects found and fixed in this session

Each would have produced plausible-looking wrong numbers.

| defect | consequence had it shipped |
|---|---|
| **`results/` was a mixture of two runs** — synthetic-fixture scores (`IVY00`–`IVY23`, 212 samples, 500 permutations) on top of an interrupted real run | every scored artefact was synthetic while being read as real. Quarantined with evidence to `results_superseded/`, not deleted |
| **`config.USE_SIGNATURE_GENE_SUBSET` declared, documented as load-bearing, and never read** on the frozen-signature path — the only path this project runs on | every method received all 16,007 shared genes instead of 1,614. **Every NuSVR fit hit the 200,000-iteration cap**, so no sample got a converged SVR solution, and the two cohorts would have taken 11.2 hours |
| **Bisque set neither `degenerate_` nor `degeneracy_reason_`** while its own docstring said the caveat must travel | `implementation_report.json` recorded `degenerate: false` for a method degraded twice over — no-overlap mode, *and* a donor-profile fallback substituting spread across cell types for spread across donors |
| SCDC ENSEMBLE set `degenerate_` with no reason | every report surfaces degeneracy by printing the reason, so `scdc_ensemble` would have appeared as though an ensemble had happened |
| **`acs_leaderboard.csv` carried no degeneracy at all** | the headline artefact showed ten independent real methods where there are eight |
| **`control_verdict` compared against the MEDIAN real method** | a control could tie the *worst* real method and beat its own null while the run announced "CONTROLS BEHAVE" — which is exactly what it did |
| **both controls reported as single draws** | a seed-dependent tie read as a property of the constraint set |
| clinical URLs all dead; one returned HTTP 200 + `text/csv` with an HTML body | survival BLOCKED for a reason that was fixable by reading the portal's own download page |
| `tumor_details.csv` keyed on `donor_id`, expression on `tumor_id`, **neither file containing the other's key** | clinical data could not be joined to anything at all |
| `age_in_years` ships as `"61 yrs"`; `to_numeric` yields all-NaN | the only baseline covariate in the survival model silently dropped |
| release bundle pointed at `benchmark/` for artefacts that moved to `anatomic/` | four files silently reported as "not produced this run" |

---

## 10. What would move this forward, in order

1. **Publicly register the constraint file** (OSF). The hash is above; the timestamp must
   precede any result anyone cites. Cheapest credibility in the project.
2. **Install the R packages** (`Rscript R/install_deps.R`). All four published tools
   currently run as this project's Python reimplementations. Every artefact says so, but
   they are not MuSiC, DWLS, Bisque or SCDC.
3. **Obtain a real yardstick covering ≥6 methods** — a GBM single-cell atlas at
   `data/reference/gbmap_core.h5ad` enables both `synthetic_mixtures` and
   `sc_pseudobulk`. Necessary but **not sufficient**: see §3's resolution note.
4. **A second tissue (protocol step 7) is the only way to widen ACS's evidence base.**
   Ivy GAP has ten H&E-annotated tumours and cannot supply more.

---

## 11. What was deliberately not done

- **`constraints.py` was not edited.** Not when a control tied a real method, not
  afterwards. Freeze hash unchanged and verified.
- **No method was tuned to raise its ACS.** One shared signature, documented defaults,
  no per-method tuning. "Optimising ACS" in the sense of making methods score better is
  what this study exists to *test*, not to perform.
- **ACS selected nothing.** `method_selection_decision.json` records that no method was
  frozen by this run, and why.
- **No outcome selected anything.** Survival ran last, after the selection record was
  written, and `assert_selection_frozen` verified the recorded criterion mentions no
  outcome.
