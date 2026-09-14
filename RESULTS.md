# Results — the Anatomy Test

**Run:** 2026-09-10  
**Cohort:** Ivy GAP, ACS scored on 122 H&E anatomic samples / 10 tumours, of which **9** contribute at least one evaluable constraint pair (the per-method tables report 9)  
**Constraint freeze hash:** `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`  
**Gene space:** 657 genes · **reference:** gbmap:6fa50a03ee163760  

Read [Anatomy_Test.md](Anatomy_Test.md) first — it is the protocol. Everything below is rendered from the artefacts by `scripts/summarize_results.py`; nothing is transcribed by hand, including this header.

---

## 0. The headline

> Spearman **rho = 0.7501** between the ACS ranking and the ranking from real ground truth, bootstrap CI **[0.328, 0.957]**, p = 0.0020, on 14 methods.
>
> Bar fixed in advance: rho >= 0.60 **and** a CI excluding zero.

HEADLINE: ACS ranking tracks true accuracy at the pre-registered bar. Anatomic concordance is usable to choose a deconvolution method in tissue with an anatomical atlas, using no ground truth and no outcomes. NOTE: 14 methods but only 13 distinct score pairs — scdc=scdc_ensemble are degenerate duplicates. Judge the correlation on 13 points, not 14.


## 1a. Integrity checks

**Pre-registration:** **REGISTERED**. REGISTERED: constraint file 2d1fb47c98832adf... registered at 2026-09-10T03:32:03Z with OSF, and every result file postdates it.

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

**Method configuration:** every method's parameters are recorded in `method_configs.json`. Declared departures from published defaults: **bayesprism, bisque, cibersortx, dwls, music, quantiseq, scdc, scdc_ensemble**.

- `cibersortx` — batch-correction mode (B-mode vs S-mode): no single default; Newman et al. (2019) choose per configuration. Their Supplementary Table 1d records the choice for every deconvolution in the paper. → B-mode. Supplementary Table 1d shows a consistent split: every signature derived from 10x Chromium and applied to a bulk RNA-seq mixture is deconvolved in S-MODE, while B-mode is used for SMART-Seq2- and microarray-derived signatures (LM22) against bulk. The split is mechanistic — droplet 10x data carries a strong 3' bias and UMI counting, so it sits further from bulk RNA-seq than full-length SMART-Seq2 does. This project's reference is GBmap, measured at 87.1% 10x (214,284 cells 3' v2, 49,262 3' v3, 31,316 5' v1 of 338,564) against 2.7% Smart-seq2, and the mixtures are Ivy GAP bulk RNA-seq. By the paper's own practice that is an S-mode configuration, and this implementation runs B-mode. S-mode's per-cell-type expression adjustment is not implemented here, so the row is B-mode and says so rather than claiming to be the mode the authors would have used. *Decided from Supplementary Table 1d and the atlas assay counts, both external to any score this project computed.*
- `music` — the scale of the single-cell matrix this HARNESS exports, and hence whether the package estimates cell size itself: raw per-cell counts, from which the package derives its own per-cell-type mRNA content — `music_basis` computes `M.S`, the mean library size per cell type, and `SCDC_basis` derives one the same way → cells normalised to 1e6 each before export, so every cell type's mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The export is then restricted to the bulk's marker gene space, over which per-type library sizes are NOT identical — measured at 2.17x spread on the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether this package's own cell-size estimate is neutralised is UNRESOLVED. CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central conversion is the IDENTITY on the reference every score was computed against, so it substitutes NOTHING and converts nothing. run_benchmark rebuilds the reference from the already-normalised matrix without passing cell_totals, so build_reference falls back to expression.sum(axis=0) and cell_size comes out 1e6 for all eight types (spread 1.0000). An earlier version of this entry claimed the factors were 'applied once, centrally, to every method'. They are applied nowhere.. Equal footing WAS the intent: one set of factors, computed one way, applied identically to all 16 methods rather than each package deriving its own. D12 shows the intent was never realised — no factors are applied. What the normalisation DOES still do is keep a deeply sequenced cell from dominating its type's mean profile, which is a real effect and the reason the deviation is still declared. MEASURED, not assumed: on a probe where two types differ 3x in mRNA and are mixed 50/50 by cell count (true cell fraction 0.500, true mRNA fraction 0.750), this package returns 0.7501 when handed the probe's FULL gene space and 0.5001 when handed raw per-cell counts — so it does convert cell size when its input lets it estimate one. The probe's full-space export is NOT what production writes: production exports a marker subset, over which per-type library sizes are measured to span 2.17x on the real atlas. So this deviation is DECLARED but its consequence is UNRESOLVED, and no claim is made here that the package's own conversion is neutralised. See scripts/verify_cell_size_semantics.py, results/cell_size_semantics_{normalized,raw}.json, and docs/OPEN_DEFECTS.md D1. *Decided from the measurement above, which is independent of any ACS or benchmark score — the probe is synthetic and its truth is known by construction.*
- `dwls` — buildSignatureMatrixMAST(diff.cutoff, pval.cutoff): diff.cutoff = 0.5, pval.cutoff = 0.01 → diff.cutoff = 0, pval.cutoff = 1 (filter fully open). Equal footing requires every method to receive the same pre-filtered gene space, so the informative genes have already been chosen identically for all methods. DWLS's own differential-expression pass is then a second filter, not a second opinion: at the published cutoffs it left 31 genes for 8 cell types and the driver refused to deconvolve on them. *Decided from the gene count, before any DWLS score existed.*
- `dwls` — cells used for the signature build: all cells in the reference → <= 250 cells per cell type. buildSignatureMatrixMAST costs ~30 minutes per gene set on this hardware regardless of the DE method. A signature is a per-cell-type summary, so type proportions are irrelevant to it and capping per type spends the budget where it buys accuracy. *Decided from measured runtime, before any DWLS score existed.*
- `dwls` — per-sample failure handling: an error aborts the call → a sample that fails to solve becomes NA and is counted. solveDampenedWLS fails on particular mixtures with 'NA/NaN argument'. One such sample was sending the whole cohort to the Python fallback. NA is already what this pipeline means by a failed sample and it is reported in n_failed_samples. *Decided from the failure mode, before any DWLS score existed.*
- `bisque` — the scale of the single-cell matrix this HARNESS exports, and hence whether the package estimates cell size itself: raw per-cell counts → cells normalised to 1e6 each before export, the same treatment every cell-consuming method receives. NO cell-size factors are applied on top: docs/OPEN_DEFECTS.md D12 measures the central conversion as the IDENTITY, because run_benchmark rebuilds the reference from the already-normalised matrix without passing cell_totals, so cell_size comes out 1e6 for every type (spread 1.0000). An earlier version of this entry said the factors were 'applied once and centrally'; they are applied nowhere. Equal footing: one set of factors applied identically to all 16 methods. Declared for this method even though its own convention is UNMEASURED — the cell-size probe cannot read it, because on a roster with six types absent from the mixture this package assigns ~75% of its mass to the absent types and the remaining ratio measures nothing about a cell-size convention. So this entry records what the harness does, and explicitly does NOT claim the package would have corrected cell size itself. An earlier version of docs/OPEN_DEFECTS.md D1 asserted this package was double-corrected; that claim is withdrawn as unmeasured. See scripts/verify_cell_size_semantics.py. *Decided from the harness's own code path, which is a fact about this pipeline and not about any score.*
- `scdc` — the scale of the single-cell matrix this HARNESS exports, and hence whether the package estimates cell size itself: raw per-cell counts, from which the package derives its own per-cell-type mRNA content — `music_basis` computes `M.S`, the mean library size per cell type, and `SCDC_basis` derives one the same way → cells normalised to 1e6 each before export, so every cell type's mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The export is then restricted to the bulk's marker gene space, over which per-type library sizes are NOT identical — measured at 2.17x spread on the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether this package's own cell-size estimate is neutralised is UNRESOLVED. CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central conversion is the IDENTITY on the reference every score was computed against, so it substitutes NOTHING and converts nothing. run_benchmark rebuilds the reference from the already-normalised matrix without passing cell_totals, so build_reference falls back to expression.sum(axis=0) and cell_size comes out 1e6 for all eight types (spread 1.0000). An earlier version of this entry claimed the factors were 'applied once, centrally, to every method'. They are applied nowhere.. Equal footing WAS the intent: one set of factors, computed one way, applied identically to all 16 methods rather than each package deriving its own. D12 shows the intent was never realised — no factors are applied. What the normalisation DOES still do is keep a deeply sequenced cell from dominating its type's mean profile, which is a real effect and the reason the deviation is still declared. MEASURED, not assumed: on a probe where two types differ 3x in mRNA and are mixed 50/50 by cell count (true cell fraction 0.500, true mRNA fraction 0.750), this package returns 0.7501 with the production export and 0.5336 with a raw export, so it does convert when it can and production prevents it ONLY when the export spans the full gene space, which production's does not. UNRESOLVED for production. See scripts/verify_cell_size_semantics.py, results/cell_size_semantics_{normalized,raw}.json, and docs/OPEN_DEFECTS.md D1. *Decided from the measurement above, which is independent of any ACS or benchmark score — the probe is synthetic and its truth is known by construction.*
- `scdc_ensemble` — the scale of the single-cell matrix this HARNESS exports, and hence whether the package estimates cell size itself: raw per-cell counts, from which the package derives its own per-cell-type mRNA content — `music_basis` computes `M.S`, the mean library size per cell type, and `SCDC_basis` derives one the same way → cells normalised to 1e6 each before export, so every cell type's mean library size is identical ACROSS THE FULL TRANSCRIPTOME. The export is then restricted to the bulk's marker gene space, over which per-type library sizes are NOT identical — measured at 2.17x spread on the real atlas (Tumor 29,999 to Oligodendrocyte 65,003) — so whether this package's own cell-size estimate is neutralised is UNRESOLVED. CORRECTED 2026-09-14, docs/OPEN_DEFECTS.md D12: this project's central conversion is the IDENTITY on the reference every score was computed against, so it substitutes NOTHING and converts nothing. run_benchmark rebuilds the reference from the already-normalised matrix without passing cell_totals, so build_reference falls back to expression.sum(axis=0) and cell_size comes out 1e6 for all eight types (spread 1.0000). An earlier version of this entry claimed the factors were 'applied once, centrally, to every method'. They are applied nowhere.. Equal footing WAS the intent: one set of factors, computed one way, applied identically to all 16 methods rather than each package deriving its own. D12 shows the intent was never realised — no factors are applied. What the normalisation DOES still do is keep a deeply sequenced cell from dominating its type's mean profile, which is a real effect and the reason the deviation is still declared. MEASURED, not assumed: on a probe where two types differ 3x in mRNA and are mixed 50/50 by cell count (true cell fraction 0.500, true mRNA fraction 0.750), this package returns 0.7501 when handed the probe's FULL gene space and 0.5336 when handed raw per-cell counts — so it does convert cell size when its input lets it estimate one. The probe's full-space export is NOT what production writes: production exports a marker subset, over which per-type library sizes are measured to span 2.17x on the real atlas. So this deviation is DECLARED but its consequence is UNRESOLVED, and no claim is made here that the package's own conversion is neutralised. See scripts/verify_cell_size_semantics.py, results/cell_size_semantics_{normalized,raw}.json, and docs/OPEN_DEFECTS.md D1. *Decided from the measurement above, which is independent of any ACS or benchmark score — the probe is synthetic and its truth is known by construction.*
- `quantiseq` — the bulk gene space this HARNESS hands the method: the full expression matrix → the full shared space, NOT the shared marker subset. Every other method solves against this project's signature matrix, so the marker subset — top-N genes per cell type ranked against that signature, plus its markers — is the shared gene space and equal footing holds. quanTIseq does not solve against it: it ships TIL10 and ignores the signature it is handed. Ranking genes against a reference it never reads kept 34 of TIL10's 138 signature genes (24.6%), and quanTIseq answered with macrophages at 100% of cells in the median sample and T cells identically zero in all 122. On the full space it finds 136 of 138 (98.6%) and returns a GBM-plausible 15.7% median myeloid fraction. The subset was not equal footing for this method; it was mutilation. EPIC is also signature-only but DOES read the supplied signature, so it keeps the subset and is untouched by this. *Decided from the signature-gene recovery rate, before any quanTIseq ACS existed — the 2026-09-05T2154 run scored it NaN.*
- `quantiseq` — cell-size (mRNA) correction: scale_mRNA = TRUE, quanTIseq's own mRNA scaling → quanTIseq's own scaling only; this project's is not applied on top. The pipeline's cell-size correction and quanTIseq's scale_mRNA are the same correction. Applying both would double-correct, and this project's version renormalises to sum 1, which on a method covering four immune types would assert the tumour is entirely immune. Keeping the published default and skipping ours leaves quanTIseq on its own documented scale. *Decided from the two definitions, before any quanTIseq ACS existed.*
- `bayesprism` — the scale of the single-cell matrix this HARNESS exports, and hence whether the package estimates cell size itself: raw per-cell counts → cells normalised to 1e6 each before export, the same treatment every cell-consuming method receives. NO cell-size factors are applied on top: docs/OPEN_DEFECTS.md D12 measures the central conversion as the IDENTITY, because run_benchmark rebuilds the reference from the already-normalised matrix without passing cell_totals, so cell_size comes out 1e6 for every type (spread 1.0000). An earlier version of this entry said the factors were 'applied once and centrally'; they are applied nowhere. Equal footing, as above. This package's own convention is UNMEASURED under the production export — the probe shows ~31% leakage onto roster types absent from the mixture, which makes the reading unusable. Its PYTHON reimplementation measures 0.7500 on the same probe, i.e. mRNA share, which if it carries over to the R package would make this project's conversion the first and correct one. Not assumed. See scripts/verify_cell_size_semantics.py and docs/OPEN_DEFECTS.md D1. *Decided from the harness's own code path, which is a fact about this pipeline and not about any score.*

**What actually ran:** 6 published package(s) ran as the genuine R implementation (music, bisque, scdc, scdc_ensemble, epic, quantiseq).

Fell back to this project's Python reimplementation:

| method | why |
|---|---|
| dwls | exceeded its 2400s budget for the genuine R package and was stopped; this is a wall-clock limit, not a failure of the method |
| bayesprism | exceeded its 2400s budget for the genuine R package and was stopped; this is a wall-clock limit, not a failure of the method |

## 2. Cohort and provenance

- constraint freeze hash: `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`
- ACS cohort: **122 samples / 10 tumours**, 10,000 within-tumour permutations per method
- tumours contributing an evaluable constraint pair: **9 of 10**. Every per-method row reports 9. The remainder are not dropped by a filter — the archive gave them none of the structure pairs the constraints name, so all seven constraints return *not evaluable* and are excluded from numerator and denominator alike.
- structures: {'CT': 30, 'MVP': 25, 'PAN': 24, 'IT': 24, 'LE': 19}
- portal reconciliation: RECONCILED: 122 anatomic samples across 5 structures agree exactly with the portal's own table ({'CT': 30, 'IT': 24, 'LE': 19, 'MVP': 25, 'PAN': 24}), and the study assignment agrees on all 270 shared samples. 9 sample(s) exist on the live portal but not in the 2014-11-25 archive — a release skew, all in the 'Cancer Stem Cells RNA Seq' study, not a parsing error. Step 3's gate is closed.

## 3. The ACS leaderboard

_from `results/anatomic/acs_leaderboard.csv`_

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| music | 1.000 | 1.000 – 1.000 | 0.370 | 0.000 | 9 |  |
| nnls | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| svr | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| elastic_net | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| epic | 0.985 | 0.953 – 1.000 | 0.377 | 0.000 | 9 |  |
| cibersortx | 0.969 | 0.930 – 1.000 | 0.374 | 0.000 | 9 |  |
| cibersortx_smode | 0.969 | 0.906 – 1.000 | 0.376 | 0.000 | 9 |  |
| scdc | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| scdc_ensemble | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| bisque | 0.923 | 0.866 – 0.983 | 0.380 | 0.000 | 9 |  |
| bayesprism | 0.877 | 0.786 – 0.968 | 0.374 | 0.000 | 9 |  |
| bayesian | 0.769 | 0.627 – 0.921 | 0.380 | 0.000 | 9 |  |
| bayesian_hierarchical | 0.769 | 0.627 – 0.921 | 0.381 | 0.000 | 9 |  |
| dwls | 0.738 | 0.600 – 0.867 | 0.374 | 0.000 | 9 |  |
| **control_random** | 0.400 | 0.288 – 0.530 | 0.382 | 0.431 | 9 | **yes** |
| **control_shuffled_signature** | 0.138 | 0.059 – 0.228 | 0.349 | 0.999 | 9 | **yes** |
| quantiseq | 0.600 | 0.333 – 0.882 | 0.502 | 0.309 | 9 |  |


**Control verdict.** CONTROLS BEHAVE: best control 0.400 sits below the median real method 0.962, ties no real method, and does not beat its own permutation null. The constraint set discriminates.


### What ACS can and cannot resolve here

- weighted denominator: **65** (constraint x tumour pairs, weighted)
- so ACS takes at most **66 distinct values**, spaced **0.0154** apart
- 17 methods scored produce **11 distinct values**; the 15 real methods produce **9**
- tied at **0.985**: nnls, svr, elastic_net, epic
- tied at **0.969**: cibersortx, cibersortx_smode
- tied at **0.954**: scdc, scdc_ensemble
- tied at **0.769**: bayesian, bayesian_hierarchical


### Control audit (stricter than the headline verdict)

_The control is a **single permutation** drawn from `config.RANDOM_SEED`. Everything in this block describes that one draw._

- best control: **control_random = 0.400**
- worst real method: quantiseq = 0.600
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

### Which constraints carry the ranking

Across the **14 comparable real methods**. `spread` is the range of the satisfied rate; a spread of zero means the constraint is unanimous and contributes nothing to the ORDERING, however well it separates real methods from the controls.

| ID | claim | weight | tumours | mean rate | min | max | spread | controls |
|---|---|---|---|---|---|---|---|---|
| C6 | Macrophage_Microglia: MVP > CT | 1.0 | 9 | 0.770 | 0.222 | 1.000 | **0.778** | 0.278 |
| C7 | Tumor: LE < IT < CT | 2.0 | 8 | 0.875 | 0.500 | 1.000 | **0.500** | 0.000 |
| C1 | Tumor: CT > LE | 1.0 | 8 | 0.929 | 0.625 | 1.000 | **0.375** | 0.375 |
| C5 | Macrophage_Microglia: PAN > LE | 1.0 | 6 | 0.917 | 0.667 | 1.000 | **0.333** | 0.250 |
| C2 | Oligodendrocyte: LE > CT | 1.0 | 8 | 0.991 | 0.875 | 1.000 | **0.125** | 0.375 |
| C4 | Endothelial: MVP is the maximum | 1.0 | 9 | 0.992 | 0.889 | 1.000 | **0.111** | 0.278 |
| C3 | Endothelial: MVP > CT | 1.0 | 9 | 1.000 | 1.000 | 1.000 | **0.000** **unanimous** | 0.556 |

**1 of 7 constraints is unanimous (C3).** The ranking is carried by the remaining 6, led by `C6` (spread 0.778), `C7` (spread 0.500), `C1` (spread 0.375).


### Leave-one-out robustness

Recomputed from the archived (constraint x tumour) satisfaction matrix, which reproduces every published ACS to 1e-16 before any variant is reported. **Reporting only** — the registered ACS uses all constraints and all tumours, and no variant below may reorder the leaderboard, select a method, or justify dropping a constraint.

| excluded | rho vs full ranking | max rank move | methods moved | top method |
|---|---|---|---|---|
| C1 | 1.0000 | 0.0 | 0 | `music` |
| C2 | 1.0000 | 0.0 | 0 | `music` |
| C3 | 1.0000 | 0.0 | 0 | `music` |
| C4 | 1.0000 | 0.0 | 0 | `music` |
| C5 | 0.9728 | 2.0 | 10 | `music` |
| C6 | 0.6780 | 6.0 | 13 | `music` |
| C7 | 0.9287 | 5.0 | 9 | `music` |
| tumor:292163427 | 1.0000 | 0.0 | 0 | `music` |
| tumor:703393 | 0.9261 | 5.0 | 11 | `music` |
| tumor:703493 | 0.9682 | 2.0 | 9 | `music` |
| tumor:705757 | 1.0000 | 0.0 | 0 | `music` |
| tumor:705758 | 1.0000 | 0.0 | 0 | `music` |
| tumor:705803 | 0.9705 | 2.0 | 9 | `music` |
| tumor:705859 | 0.9966 | 1.0 | 3 | `music` |
| tumor:711547 | 1.0000 | 0.0 | 0 | `music` |
| tumor:711560 | 1.0000 | 0.0 | 0 | `music` |

**The top method survives every single exclusion.** No one constraint and no one tumour is responsible for it.

**The ordering below the top is not equally robust.** Dropping `C6` moves the ranking to rho 0.678, with a maximum move of 6 positions — so that one constraint carries a large share of the separation between the middle-ranked methods. Tumour exclusion is milder (worst rho 0.926). This is reported as a limit on how finely the leaderboard can be read, not as a reason to change the constraint set.


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
| cibersortx | python | no |  |
| cibersortx_smode | python | no |  |
| elastic_net | python | no |  |
| bayesian | python | no |  |
| bayesian_hierarchical | python | no |  |
| music | R:MuSiC | no |  |
| dwls | python-reimplementation | no | exceeded its 2400s budget for the genuine R package and was stopped; this is a wall-clock limit, not a failure |
| bisque | R:BisqueRNA | **yes** | no subjects are assayed both as bulk and as single cells, so this runs Bisque's documented no-overlap mode (us |
| scdc | R:SCDC | no |  |
| scdc_ensemble | R:SCDC | **yes** | one reference supplied ('gbmap'), so there is nothing to weight across: SCDC ENSEMBLE reduces exactly to SCDC. |
| epic | R:EPIC | no |  |
| quantiseq | R:quantiseqr | no |  |
| bayesprism | python-reimplementation | no | exceeded its 2400s budget for the genuine R package and was stopped; this is a wall-clock limit, not a failure |
| control_random | python | no |  |
| control_shuffled_signature | python | no |  |


## 5. The agreement test (the primary result)

Pre-registered bar: rho >= 0.6 AND a bootstrap CI excluding zero. Minimum methods: 6.

| yardstick | rho | 95% CI | methods | distinct | verdict |
|---|---|---|---|---|---|
| synthetic_mixtures | 0.750 | 0.328 – 0.957 | 14 | 13 | HEADLINE: ACS ranking tracks true accuracy at the pre-registered bar. Anatomic concordance |
| absolute_purity | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| sc_pseudobulk | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |
| simulated_donor_mismatch | — | — – — | 0 | — | UNAVAILABLE: this yardstick produced no scores in this run. |


## 5a. The constraints against in-situ hybridization

Ivy GAP's ISH panel quantifies expression energy for 480 genes over 899 sub-blocks from 42 donors, in the same five anatomic structures. **No deconvolution method touches it**, so it is the one test here of whether the constraints are true of the tissue rather than agreed on by solvers.

**Exploratory, and the constraints do not move on it.** The constraint file is frozen at `2d1fb47c98832adf…` and registered. Nothing below may edit it.

**ISH energy is not composition.** It is transcript signal in a region, which rises with expression per cell as well as with cell number. A marker can move the right way for the wrong reason, and the reverse.

Markers were declared before any value was read. **8 of 17 testable marker-constraint tests reach p < 0.05** under a within-sub-block permutation null; satisfaction is averaged within donor before averaging across donors, because sub-blocks nest in donors.

| constraint | claim | marker | donors | blocks | satisfied rate | null p |
|---|---|---|---|---|---|---|
| C1 | Tumor: CT > LE | `SOX2` | 5 | 5 | 0.000 | 1.0000 |
| C1 | Tumor: CT > LE | `PTPRZ1` | 18 | 22 | 0.056 | 1.0000 |
| C1 | Tumor: CT > LE | `EGFR` | 0 | 0 | — | — |
| C1 | Tumor: CT > LE | `CD44` | 26 | 137 |  **0.900** |  **0.0005** |
| C1 | Tumor: CT > LE | `BIRC5` | 26 | 134 |  **0.630** |  **0.0070** |
| C1 | Tumor: CT > LE | `TOP2A` | 5 | 5 | 0.600 | 0.4963 |
| C3 | Endothelial: MVP > CT | `ESM1` | 6 | 8 |  **1.000** |  **0.0050** |
| C3 | Endothelial: MVP > CT | `CD34` | 2 | 2 | 0.000 | — |
| C3 | Endothelial: MVP > CT | `KDR` | 3 | 3 | 0.667 | 0.4848 |
| C3 | Endothelial: MVP > CT | `CAV1` | 4 | 4 | 0.500 | 0.6837 |
| C4 | Endothelial: MVP is the maximum | `ESM1` | 6 | 8 |  **0.750** |  **0.0205** |
| C4 | Endothelial: MVP is the maximum | `CD34` | 2 | 2 | 0.000 | — |
| C4 | Endothelial: MVP is the maximum | `KDR` | 3 | 3 | 0.667 | 0.2159 |
| C4 | Endothelial: MVP is the maximum | `CAV1` | 4 | 4 | 0.250 | 0.7291 |
| C5 | Macrophage_Microglia: PAN > LE | `CD163` | 1 | 1 | 1.000 | — |
| C5 | Macrophage_Microglia: PAN > LE | `LAPTM5` | 1 | 1 | 1.000 | — |
| C6 | Macrophage_Microglia: MVP > CT | `CD163` | 6 | 7 |  **1.000** |  **0.0085** |
| C6 | Macrophage_Microglia: MVP > CT | `LAPTM5` | 1 | 1 | 1.000 | — |
| C7 | Tumor: LE < IT < CT | `SOX2` | 4 | 4 | 0.000 | 1.0000 |
| C7 | Tumor: LE < IT < CT | `PTPRZ1` | 13 | 15 | 0.000 | 1.0000 |
| C7 | Tumor: LE < IT < CT | `EGFR` | 0 | 0 | — | — |
| C7 | Tumor: LE < IT < CT | `CD44` | 24 | 112 |  **0.529** |  **0.0005** |
| C7 | Tumor: LE < IT < CT | `BIRC5` | 24 | 109 |  **0.260** |  **0.0160** |
| C7 | Tumor: LE < IT < CT | `TOP2A` | 4 | 4 |  **0.750** |  **0.0165** |

**What this supports.** The endothelial and perivascular-myeloid claims — the strongest in the set on neuropathology — have independent measured support: ESM1 satisfies C3 in every evaluable donor (p = 0.005) and C4 at 0.750 (p = 0.021), and CD163 satisfies C6 in every evaluable donor (p = 0.0085).

**What it does not.** C1 and C7 — the tumour constraints — come out **marker-dependent, and the disagreement is not noise**. CD44 satisfies C1 at 0.900 over 137 blocks and 26 donors (p = 0.0005) and BIRC5 at 0.630 (p = 0.007), while SOX2 and PTPRZ1 satisfy it at **0.000 and 0.056**, with a null p of 1.000 — as far the other way as the data allow. All four are real tumour markers.

The reading that fits is that these markers track different tumour programmes rather than tumour cell *density*, which is what C1 claims: SOX2 and PTPRZ1 mark stem-like and OPC-like states reported to be enriched at the infiltrating margin, while CD44 and BIRC5 mark mesenchymal and proliferating states concentrated in the dense core. That is an interpretation, offered as one, and it is a hypothesis this study does not test. What is measured is that **no marker in this panel measures tumour cell density**, so the ISH panel cannot adjudicate C1 or C7 either way.

**C2 is unmeasurable here.** Oligodendrocyte: LE > CT cannot be checked. Ivy GAP's ISH panel contains no myelin or oligodendrocyte-lineage marker — MBP, PLP1, MOG, MAG, CNP, SOX10, MOBP and CLDN11 are all absent. OLIG2 IS present and is deliberately not used: in glioma it is expressed by the tumour cells themselves and by OPCs, so scoring an oligodendrocyte claim on it would measure tumour content and report it as oligodendrocyte content.

**Power.** 11 of 24 marker-constraint rows rest on fewer than five donors, and 2 on none at all (EGFR is in the panel but yielded no evaluable structure pair). The endothelial and myeloid results rest on 6-8 blocks. That is thin, and the support they give the constraint set is correspondingly weak — real, independent, and small.


## 6. Does deconvolving all 270 samples change the ranking?

| method | ACS (122 anatomic) | ACS (270 deconvolved) | delta | rank change |
|---|---|---|---|---|
| music | 1.000 | 1.000 | +0.000 | +0 |
| svr | 0.985 | 0.985 | +0.000 | +0 |
| nnls | 0.985 | 0.985 | +0.000 | +0 |
| elastic_net | 0.985 | 0.985 | +0.000 | +0 |
| epic | 0.985 | 0.985 | +0.000 | +0 |
| cibersortx_smode | 0.969 | 0.969 | +0.000 | +0 |
| cibersortx | 0.969 | 0.969 | +0.000 | +0 |
| scdc_ensemble | 0.954 | 0.954 | +0.000 | +1 |
| scdc | 0.954 | 0.954 | +0.000 | +1 |
| bisque | 0.923 | 0.969 | +0.046 | -3 |
| bayesprism | 0.877 | 0.877 | +0.000 | +0 |
| bayesian | 0.769 | 0.769 | +0.000 | +0 |
| bayesian_hierarchical | 0.769 | 0.785 | +0.015 | -0 |
| dwls | 0.738 | 0.738 | +0.000 | +0 |
| quantiseq | 0.600 | 0.600 | +0.000 | +0 |
| control_random | 0.400 | 0.338 | -0.062 | +0 |
| control_shuffled_signature | 0.138 | 0.138 | +0.000 | +0 |

Spearman rho between the two ACS rankings: **0.9850**


**Which methods moved, and why.** 14 of 17 did not move by a single unit: music, svr, nnls, elastic_net, epic, cibersortx_smode, cibersortx, scdc_ensemble, scdc, bayesprism, bayesian, dwls, quantiseq, control_shuffled_signature. Those solve each sample independently, so what else is in the cohort cannot reach them — the zeros are exact, not rounded.

The ones that moved are **bisque, bayesian_hierarchical, control_random**. Bisque normalises the bulk with cohort-wide per-gene statistics (`B.mean(axis=1)` and `B.std(axis=1)`), so the 148 ISH-cluster samples shift the reference frame every anatomic sample is mapped through. `control_random` moves for a different and uninteresting reason: it draws one Dirichlet sample per row, so a 270-row draw is not a superset of a 122-row draw.

Note the direction: **adding 148 more real samples made Bisque's anatomic concordance worse**, not better. For a method that borrows strength across samples, which samples it is handed is part of the method — and Ivy GAP's two studies are not the same tissue.


### The full-database run

- deconvolved **270 samples / 37 tumours**
- ACS still scored on **122 samples / 10 tumours** — H&E-selected anatomic study only; ISH-cluster samples are excluded from scoring by construction because their structure labels were assigned using expression
- control verdict: CONTROLS BEHAVE: best control 0.338 sits below the median real method 0.969, ties no real method, and does not beat its own permutation null. The constraint set discriminates.

| method | ACS | 95% CI | null mean | null p | tumours | control |
|---|---|---|---|---|---|---|
| music | 1.000 | 1.000 – 1.000 | 0.370 | 0.000 | 9 |  |
| nnls | 0.985 | 0.953 – 1.000 | 0.375 | 0.000 | 9 |  |
| svr | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| elastic_net | 0.985 | 0.953 – 1.000 | 0.374 | 0.000 | 9 |  |
| epic | 0.985 | 0.953 – 1.000 | 0.377 | 0.000 | 9 |  |
| cibersortx | 0.969 | 0.930 – 1.000 | 0.374 | 0.000 | 9 |  |
| cibersortx_smode | 0.969 | 0.906 – 1.000 | 0.377 | 0.000 | 9 |  |
| bisque | 0.969 | 0.930 – 1.000 | 0.379 | 0.000 | 9 |  |
| scdc | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| scdc_ensemble | 0.954 | 0.911 – 0.986 | 0.374 | 0.000 | 9 |  |
| bayesprism | 0.877 | 0.786 – 0.968 | 0.374 | 0.000 | 9 |  |
| bayesian_hierarchical | 0.785 | 0.648 – 0.925 | 0.380 | 0.000 | 9 |  |
| bayesian | 0.769 | 0.627 – 0.921 | 0.381 | 0.000 | 9 |  |
| dwls | 0.738 | 0.600 – 0.867 | 0.374 | 0.000 | 9 |  |
| **control_random** | 0.338 | 0.217 – 0.486 | 0.388 | 0.779 | 9 | **yes** |
| **control_shuffled_signature** | 0.138 | 0.059 – 0.228 | 0.349 | 0.999 | 9 | **yes** |
| quantiseq | 0.600 | 0.333 – 0.882 | 0.502 | 0.309 | 9 |  |


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
| cibersortx_smode | 0.484 | 0.512 | 0.027 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| dwls | 0.484 | 0.510 | 0.026 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bayesprism | 0.484 | 0.509 | 0.025 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| elastic_net | 0.484 | 0.505 | 0.021 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| epic | 0.484 | 0.504 | 0.020 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| bayesian | 0.484 | 0.503 | 0.019 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| cibersortx | 0.484 | 0.495 | 0.011 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| music | 0.484 | 0.486 | 0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc | 0.484 | 0.482 | -0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| scdc_ensemble | 0.484 | 0.482 | -0.002 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |
| quantiseq | 0.484 | 0.466 | -0.019 | INCONCLUSIVE (cohort supports detecting ~0.515; not used for ranking) |


#### Exploratory: does ACS track prognostic value?

Not a protocol analysis. The protocol's agreement test correlates ACS against *accuracy*, which is still not computable. This correlates it against *prognostic value* instead, which the declared-policy run makes available. It selects nothing.

- Spearman(ACS, delta C-index) over 15 methods: **rho = 0.168**, p = 0.549
- best ACS: `music` (ACS 1.000, delta C +0.002)
- best prognosis: `svr` (ACS 0.985, delta C +0.080)

Three independent reasons this cannot support a claim: every C-index is INCONCLUSIVE by the pre-specified power rule; ACS supplies only 9 distinct values across 15 methods, so the rank is mostly ties; and the outcome side rests on a declared assumption about censoring. It is recorded because it points the same way as the superseded run did, and because the direction is the protocol's second branch.
