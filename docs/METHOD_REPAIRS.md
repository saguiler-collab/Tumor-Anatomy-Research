# Method repairs: what was broken, how it was fixed, and how the fix was checked

Started 2026-10-07, at the user's request: "fix the broken methods and ensure the LGG and glioblastoma data
... is effectively applied so the deconvolution methods are tested against them."

**Rules that do not move:**
- **The two deliberately broken controls stay broken.** The random-proportion and shuffled-signature
  estimators exist to show the anatomic score can fail (0.400 and 0.369, against 0.708-0.969 for real
  methods).
- **The anatomic constraints do not change.** `constraints.py` is hashed and registered. Re-tuning it after
  seeing which methods win would void the main result.
- **Registered artefacts are never overwritten.** Repaired results go to `results/repaired/` and are
  reported beside the registered ones.
- **Every repair is behind `IVYGAP_REPAIRED=1`** (`ivygap/config.py`, `REPAIRED_METHODS`). This lets the
  2026-10-06 verification re-run reproduce the archived results with the code that made them.
  - Checked: that run's genuine-DWLS process started at 23:32:10, before the first repair edit (23:49:45).

**Where every method stands against every truth:** `docs/EVALUATION_MATRIX.md`
(`scripts/evaluation_matrix.py`). Its recomputed tumour correlations equal the registered ones in 51 of 51
method-arms.

---

## The repairs

| method | what was broken | evidence | repair | validation | status |
|---|---|---|---|---|---|
| **DWLS (genuine package)** | quadprog fails in DWLS's unscaled first step (OPEN_DEFECTS D31) | 73 of 122 Ivy GAP samples failed; anatomy covered 26 of 57 pairs | signature and bulk divided by one constant in `R/run_dwls.R` | 73 → 0 failures on the real inputs; identical to 3e-16 where both solve | **GBM donor-level, genuine package, all 154 samples: tumour 0.16 → 0.63, leukocytes 0.31 → 0.72. Ivy GAP anatomy, genuine package, all 122 samples estimated (was 73 failing): ACS 0.723 [0.59, 0.87] on all 57 pairs, the registered Python version's value to three decimals, from different estimates (its CI was [0.61, 0.85]). CPTAC donor-level, genuine package: 0.740 → 0.754 against whole-genome purity.** |
| **DWLS (this project's reimplementation)** | dampening capped relative to the largest weight, not the smallest (D32) | GBM donor-level: 37% endothelium on average, tumour accuracy 0.16; LGG donor-level: 1 of 510 finite | the published algorithm reimplemented exactly | equals the package to 5e-11 given its dampening constant (8 real samples) | **GBM frozen: tumour 0.69 → 0.77, leukocytes 0.67 → 0.76. LGG frozen: tumour 0.52 → 0.59, leukocytes 0.61. LGG donor-level: from 1 of 510 samples finite to all 510; tumour 0.45, leukocytes 0.56. CPTAC frozen: 0.851 → 0.854** |
| **quanTIseq** | not broken: dropped from every TCGA arm because it models no tumour (D33) | "only 0 finite estimates" in all four arms; immune estimates never saved | score every compartment a method models | -- | **Genuine quantiseqr. GBM frozen: leukocytes 0.11, lymphoid 0.22. LGG frozen: leukocytes 0.46, lymphoid 0.16. Within lymphocytes NK ~74%, T ~6% in both** |
| **Bayesian hierarchical** | missing from both frozen arms (D20: `KeyError: 'patient_id'`, before the manifest was supplied) | absent from two rankings | re-run with the manifest | **GBM frozen: tumour 0.75, leukocytes 0.51. LGG frozen: tumour 0.54, leukocytes -0.09** | done |
| **BayesPrism** | the registered arms ran this project's reimplementation (the genuine package hung or timed out; D18) | TCGA tumour accuracy 0.16 (GBM) and -0.03 (LGG) | **already measured:** the genuine package, run unbudgeted on 2026-10-01/02 | genuine: tumour 0.72 (GBM), 0.59 (LGG); anatomy 0.815 on all 57 pairs | report the genuine results |
| **SCDC ENSEMBLE** | one reference, so it reduces to SCDC (D2) | identical to SCDC everywhere | give it a second single-cell reference | -- | **needs a cell-level export** (the stored Abdelfattah, Neftel and Darmanis references are averaged profiles; SCDC ENSEMBLE needs cells). The Abdelfattah raw data (GSE182109) and its cluster-naming rule are on disk; memory-heavy, so after the verification |
| **Bisque** | no tumour has both bulk and single cells, so it runs degraded | declared in every artefact | CPTAC has both for 17 tumours, so paired mode could run there | -- | **deferred, for a reason:** paired mode learns its reference from CPTAC's nuclei labels, which failed against whole-genome purity (S2.1, rho 0.13), so it would inherit that error. Revisit only with the authors' curated labels |
| **MuSiC, SCDC, Bisque, S-mode on the frozen signature** | not applicable: the frozen signature has no cells | MuSiC = NNLS there; S-mode cannot run | report "not applicable", not a number | -- | reporting change |

## What the repairs show (all runs complete, 2026-10-07 08:38)

- **Repair raises tumour and immune accuracy.**
  - DWLS on the frozen signature: GBM tumour 0.69 → 0.77, leukocytes 0.67 → 0.76.
  - BayesPrism as the genuine package: tumour 0.16 → 0.72 (GBM) and -0.03 → 0.59 (LGG); leukocytes
    0.34 → 0.71 and -0.18 → 0.58.
- **Every previously broken or missing method-arm now has estimates.**
  - LGG donor-level DWLS went from 1 of 510 samples to all 510 (tumour 0.45).
  - Genuine DWLS on the anatomy cohort went from 49 of 122 samples to all 122 (ACS 0.723).
- **No repair recovers the lymphoid ordering.**
  - Every repaired method still puts T below B, or below NK in quanTIseq's case.
  - Every one correlates with the methylation lymphoid truth at between -0.28 and +0.22.
  - The lymphoid failure (R3) is therefore not an artefact of a broken implementation. It survives
    the genuine packages and the faithful reimplementation.

## Stand-ins (OPEN_DEFECTS D34)

Two of the Python versions are not reimplementations at all.
- **SCDC:** Huber residual weights instead of SCDC's cross-subject-variance weights.
- **BayesPrism:** a reference-shrinking heuristic instead of the Gibbs sampler.
- **On the frozen signature, none of the single-cell methods can run as published** (MuSiC, SCDC, SCDC
  ENSEMBLE, Bisque, BayesPrism). Their frozen-arm numbers are stand-ins', and MuSiC's is arithmetically
  NNLS.
- `docs/METHODS.md` now says which version is which.

## Does the registered ranking test (R2) depend on the broken implementations? (post hoc)

The test is ACS against tumour accuracy on TCGA-GBM, bar 0.60. It was recomputed with each broken value
replaced by its repaired or genuine one (`docs/EVALUATION_MATRIX.md`, "after the repairs"):

| arm | registered | with the repairs available so far |
|---|---|---|
| frozen signature | 0.081 | -0.17 (13 methods) |
| donor-level reference | 0.314 | 0.11 (14 methods; DWLS and BayesPrism accuracy repaired) |

- Anatomy still does not rank methods by accuracy once the implementations are fixed. The repairs move
  the correlation further from the bar, because the repaired methods are among the most accurate and have
  middling anatomic scores.
- The genuine DWLS's anatomy score is in: 0.723 on all 57 pairs, unchanged from the registered value. So
  the donor-level re-check stays at 0.11. DWLS's accuracy against DNA quadrupled (0.16 → 0.63) while its
  anatomy score did not move: anatomy does not see accuracy.

## What the anatomy score can and cannot separate

Per-constraint results (`results/anatomic/acs_per_constraint.csv`, the registered run):

| constraint | spread across real methods (fraction of tumours satisfied) |
|---|---|
| C1 tumour: cellular tumour > leading edge | 0.25 |
| C2 oligodendrocyte: leading edge > cellular tumour | 0.12 |
| C3, C4 endothelium in microvascular proliferation | 0.11 |
| C5 macrophages: pseudopalisading > leading edge | 0.33 |
| **C6 macrophages: microvascular proliferation > cellular tumour** | **0.67** |
| **C7 tumour gradient: leading edge < infiltrating < cellular tumour** (weight 2) | **0.62** |

- **The score separates working methods from broken ones on every constraint.** The controls fail C1, C2,
  C5 and C7 broadly.
- **Among working methods, the ordering rests on two constraints (C6, C7), each measured on 8-9
  tumours.** That is the resolution limit of the score on Ivy GAP's 10 histology-labelled tumours. More
  region-labelled tumours would sharpen it; a different constraint set chosen now would not be a test.
- **The two constraints pull against DNA accuracy.**
  - The methods most accurate against DNA purity (the Bayesian pair: 0.70-0.75 in TCGA-GBM, 0.63-0.66
    in CPTAC) satisfy the tumour gradient (C7) in 3 of 8 tumours.
  - MuSiC, nu-SVR and SCDC satisfy it in 8 of 8.
  - This is descriptive and reported, not used to select anything.
- **Methods are separable by accuracy, robustly. ACS just does not do that separating.**
  - Two independent DNA instruments rank the methods almost identically (`docs/EVALUATION_MATRIX.md`,
    post hoc):
    - by tumour accuracy, ABSOLUTE against GIMiCC: Spearman 0.80-0.87;
    - by leukocyte accuracy, LF against GIMiCC: 0.82-0.99;
    - in all four arms, GBM and LGG.
  - On the frozen signature, the ranking also transfers to CPTAC (0.68).
- **There is no anatomic score for LGG.** Ivy GAP is glioblastoma only. LGG methods are judged against DNA
  (ABSOLUTE purity, the methylation leukocyte fraction, EpiDISH) and against flow cytometry for the
  lymphoid ordering.

## Runs (as executed)

| arm | DWLS | quanTIseq | Bayesian hierarchical |
|---|---|---|---|
| GBM and LGG, frozen signature | repaired reimplementation (the package needs cells to build its signature) | genuine | re-run (was missing) |
| GBM, donor-level | **genuine package**, rescaled, unbudgeted | genuine (identical to frozen) | re-run (reproduces 0.70) |
| LGG, donor-level | repaired reimplementation: the genuine package on 510 samples would take about 8 h on this machine, and the reimplementation equals it to 5e-11 given its dampening constant | not re-run: its estimates do not depend on the reference build | registered value used |
| CPTAC, frozen / donor-level | repaired reimplementation / **genuine package** | genuine / not re-run | re-run / registered |
| Ivy GAP anatomy | **genuine package**, rescaled (`results/repaired/ivygap_dwls_genuine_rescaled.json`) | registered (partial coverage) | registered |

- quanTIseq gives identical numbers on both reference builds: it uses its own TIL10 signature.
- The order was changed on 2026-10-07 so these ran first. The verification re-run pauses after its
  registered-pipeline step and resumes afterwards (`docs/STS_WORKLOG.md` §2.33).
