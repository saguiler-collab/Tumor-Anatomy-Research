# Why do the methods put B cells above T cells?

> **Updated 2026-10-01 — read §7 first.** Sections 1-6 describe the frozen-signature result as it
> stood on 2026-09-28. §7 adds what was measured since: the raw/X reference arm, why Bisque
> appears to fix the inversion and does not, and why the per-sample framing has to go. The
> claims table in §6 is superseded by the one in §7.


The study's strongest result is that **0 of 12 methods reproduce the T > B ordering that DNA
methylation measures**, in two cohorts, against a criterion registered before the data existed.
The obvious next question is *why*, and it is also the question the paper must be most careful
about, because a wrong mechanism is worse than no mechanism.

**Short answer: nobody knows, five candidate explanations have been tested and rejected, and
two structural properties of the reference are strong suspects that have not been proven.**

This file sets out what is established, what is excluded, and what remains — with the epistemic
status of each marked, because that distinction is the whole value of the document.

---

## 1 · What the finding actually is — precisely

It is easy to state this too broadly. The precise claim is:

> Within the lymphoid compartment {T, NK, B}, renormalised on both sides so the denominators
> cancel, DNA methylation places **T above B** in 94.8% of glioblastomas and 93.2% of
> lower-grade gliomas. **No method reproduces that ordering in either cohort.**

Two failure modes hide inside that zero and must not be merged:

| mode | what happens | who |
|---|---|---|
| **ABSENCE** | returns *exactly* 0.0000 for T, NK and B in most samples | 4 of 12 in each cohort — `music`, `nnls`, `bayesprism`, `elastic_net` |
| **MISASSIGNMENT** | reports lymphocytes, then orders them wrongly | 6 of 8 (GBM) and 8 of 8 (LGG) of the rest |

A method that reports no lymphocytes has not *inverted* the compartment; it has declined to
estimate it. The mechanism question really applies to the second group.

**And the inflation is specific to B, not to "anything but T".** Ten methods in each cohort
return a non-zero T and therefore have a computable B/T ratio; **all ten exceed the true
ratio, in both cohorts**. The remaining two — `music` and `nnls` — return T exactly zero, so
their ratio is undefined and they belong to the absence mode rather than this one. Over those
same ten, NK/T exceeds truth in only **5 (GBM) and 7 (LGG)**. Whatever is happening favours B
in particular.

*(A note on counting, since this document exists to insist on the distinction: an earlier
draft of this paragraph wrote "10 of 12", which silently counted the two absence-mode methods
in a denominator that describes misassignment. Ten of ten is the correct statement and it is
the stronger one.)*

---

## 2 · Why it looks contradictory — and why that intuition is right

Three things make this result feel impossible, and each is a fact worth stating in the paper:

1. **The reference itself holds ~5× more T than B.** The atlas these methods solve against
   knows perfectly well that T cells outnumber B cells in glioma. The methods are not being
   misled by a reference that says otherwise.
2. **The signature can separate them.** Planting a known T:B = 2.333 in mixtures built from
   the reference, NNLS recovers it *exactly* — and keeps recovering it under 100% added noise
   and after deleting an entire cell type. Condition number 5.29.
3. **The matrix is well-conditioned.** Over the 634 genes actually solved on, the full 8-type
   signature has condition number **30.3**, and the lymphoid 3-type block only **3.78**. This
   is not a numerically unstable system.

So: a reference that holds the right ratio, a signature demonstrably able to recover it, and a
well-posed linear system — and every method still gets it backwards on real tissue. **The
contradiction is real, and it localises the problem to the gap between the reference and real
bulk, not to the algebra.**

---

## 3 · Five mechanisms tested and RULED OUT

Each was tested with a working control (`docs/supplementary/S5_rejected_mechanisms.csv`).

| # | proposed explanation | test | outcome |
|---|---|---|---|
| 1 | The signature cannot separate T from B | NNLS on mixtures built from the reference, planted T:B = 2.333 | **recovered exactly**; survives 100% noise and deletion of a whole type |
| 2 | B absorbs macrophage/microglia spillover | corr(B, Macrophage) > corr(T, Macrophage) across samples | **5 of 10** methods; sign-test p = 0.62 |
| 3 | High tumour purity explains the absence mode | one-sided Mann–Whitney, purity in zero-lymphoid vs other samples | GBM 2/8, LGG 1/8; **several run the opposite way** |
| 4 | Per-sample model fit is a usable trust signal | corr(R², \|estimate − purity\|) per method | GBM 1/12, LGG 6/12; **2 LGG methods significant in the wrong direction** |
| 5 | Misfit in Tumor/Astrocyte markers propagates to lymphoid coefficients | refit on lymphoid-owned markers only, 596 genes | **T > B in 0.0% of samples** — and the mass moves to **NK** |

**Mechanism 5 is the most informative failure.** Strip the problem down to genes that only
lymphoid types own — removing any possibility of contamination from tumour, astrocyte or
macrophage columns — and the methods *still* do not put T first. The lymphoid failure survives
in isolation, which means it originates **inside the lymphoid block itself**.

---

## 4 · Two structural facts about the reference — HYPOTHESIS, not established

These are properties of the GBmap-derived signature, measured here. Neither has been shown to
*cause* the ordering failure, and the paper must not say they do.

### 4a · The `NK_cell` column carries the defining T-cell markers

CD3 is *the* T-cell marker; NK cells are CD3-negative by definition. In this reference:

| marker | `T_cell` | `NK_cell` | highest in | NK/T |
|---|---|---|---|---|
| CD3D | 1,489.5 | **3,411.3** | NK_cell | 2.29 |
| CD3E | 579.1 | **1,271.7** | NK_cell | 2.20 |
| CD3G | 400.4 | **1,720.4** | NK_cell | 4.30 |
| CD2 | 1,349.4 | **1,918.6** | NK_cell | 1.42 |
| LCK | 656.5 | **3,827.2** | NK_cell | 5.83 |
| IL7R | **1,127.2** | 43.1 | T_cell | 0.04 |
| CD8A | **1,193.6** | 20.4 | T_cell | 0.02 |

**Five of nine canonical pan-T markers are strongest in the `NK_cell` column.** The most likely
reading is that this atlas's "NK" population contains CD3⁺ cells — NKT, γδ T, or cytotoxic T
cells annotated as NK — which is a known and common difficulty in single-cell annotation, not
a defect unique to this atlas.

**Why this is a suspect.** If real T cells in bulk tissue express CD3 strongly, their signal is
better explained by a column that *also* expresses CD3 strongly. T-cell signal would then be
attributed to NK, depressing the T estimate. **This is exactly what mechanism 5 observed** when
the fit was restricted to lymphoid markers: the mass moved to NK.

**Why it is not sufficient.** It predicts T→NK confusion. Among the ten methods that return a
non-zero T, NK/T exceeds truth in only **5 (GBM) and 7 (LGG)**, while **B/T exceeds truth in
all ten, in both cohorts**. The CD3 story explains part of the T under-call; it does not
explain why the freed mass lands on **B**.

### 4b · The `B_cell` profile is closest to `Macrophage_Microglia`

Correlations between profiles, log scale, over the solved gene space:

| | closest neighbour | r |
|---|---|---|
| `T_cell` | `NK_cell` | 0.497 |
| `NK_cell` | `T_cell` | 0.497 |
| **`B_cell`** | **`Macrophage_Microglia`** | **0.497** |

B's nearest neighbour is the population that *dominates* the glioma immune compartment —
methylation puts monocyte/macrophage lineage at ~0.58 of leukocytes in GBM. A profile that
shares structure with a hugely abundant population is exposed to absorbing its residual.

**Note the tension with mechanism 2**, which tested macrophage spillover and rejected it. The
two are not the same test: mechanism 2 correlated the *estimates* across samples; this is a
property of the *profiles*. A reviewer will notice, and the paper should say which was tested.

---

## 5 · What has NOT been tested, and would settle it

Honest list, in descending order of how decisive each would be.

1. **Re-annotate NK.** Rebuild the reference with CD3⁺ cells removed from the NK population,
   or with NK and T merged into a single "T/NK" class, and re-run. If T > B appears, 4a is the
   mechanism. This is the single most informative experiment available and it is not expensive.
2. **Swap the atlas.** The reference-sensitivity arm already holds Neftel, Darmanis and Albiach.
   If the B > T inversion is a property of GBmap's annotation it should weaken or vanish under
   a reference annotated by different people.
3. **Drop `Macrophage_Microglia` from the roster** and re-fit. If B's inflation falls, 4b is
   implicated — a cleaner test than correlating estimates.
4. **Simulate the failure.** Build mixtures with a *known* T:B using T-cell profiles drawn from
   a source other than GBmap's T column. If the methods still invert, the problem is the
   reference; if they recover, the problem is the real-bulk mismatch.

---

## 6 · What the paper should claim

| claim | status |
|---|---|
| No method reproduces the registered T > B ordering, in either cohort | **state plainly** |
| The failure has two distinct modes, absence and misassignment, counted separately | **state plainly** |
| The reference holds ~5× more T than B, and the signature recovers a planted ratio exactly | **state plainly** — it is what makes the result surprising |
| Five candidate mechanisms were tested and rejected, each with a control | **state plainly** |
| The NK column carries pan-T markers, and B's profile is closest to macrophage | **state as measured properties of the reference** |
| Either of those *causes* the inversion | **must not be claimed** — untested |
| The cause is unknown | **state plainly.** It is the honest position and it is not a weakness |

A negative result with five excluded mechanisms and two named suspects is a stronger
contribution than a mechanism asserted on a correlation. The experiments in §5 are what a
follow-up study does — and naming them precisely is how this paper earns that follow-up.

---

## 7 · What changed on 2026-09-30 / 10-01 — and what the paper may now claim

### 7a · The "0 of 12" is one reference's answer

Everything above was measured against the **frozen** signature. The same comparison run against
GBmap rebuilt from genuine counts with donor structure (the h5ad arm, `raw/X`) gives:

| | frozen | h5ad (raw/X) |
|---|---|---|
| GBM: methods ordering T above B (cohort mean) | 0 of 12 | **3 of 14** — bisque, cibersortx_smode, elastic_net |
| LGG: methods ordering T above B (cohort mean) | 0 of 12 | **2 of 13** — bisque, elastic_net |

This was already in `docs/IMMUNE_ARM.md`; it never reached the manuscript (OPEN_DEFECTS D23).

### 7b · Bisque's success is its reference talking

Bisque on the h5ad reference orders T > NK > B in both cohorts, with B above T in only 3.6% (GBM)
and 1.0% (LGG) of samples. It is not measuring the lymphoid compartment:

- **The package says so.** With no subjects shared between single-cell and bulk data — always
  the case for TCGA — Bisque *"assume[s] that the observed mean of Yj is the true mean of our goal
  distribution"* (Jew et al. 2020, *Nat Commun*, Methods). The cohort mean is set by the reference.
- **Planted truth proves it.** Synthetic bulk whose true mean composition is far from the
  reference's comes back as the reference's composition **to four decimals** (L1 to the
  reference 0.000; to the truth 1.091), while samples are still ranked almost perfectly
  (rho 0.88-0.97). `results/bisque_anchoring_synthetic.json`.
- **Real data agrees.** Bisque's TCGA cohort means sit on its reference's donor-mean composition
  across all eight cell types (L1 0.087 GBM, 0.081 LGG; every other method 0.68-1.40). It reports
  glioblastoma as 27% tumour — the reference's number, identical in LGG. `results/bisque_anchoring.json`.

So one of the three, and the only one that agrees in both cohorts robustly, is the reference's
own T > B handed back. The other two are weak (CIBERSORTx S-mode: GBM only, B above T in 41% of
samples) or fragile (elastic net: no lymphocytes in 30-39% of samples, and B at 0.94-0.98 of the
compartment under the frozen signature).

### 7c · Methylation is a cohort-level truth here, not a per-sample one

A reference-free check — the mean expression of six T-cell genes minus eight B-cell genes, from
the same bulk RNA, panels fixed before running — should rise and fall with methylation's
per-sample T/(T+B) if both measure the same thing. It does not: rho +0.14 (GBM, p 0.27) and
**-0.12 (LGG, p 0.004)**, and a pre-declared stratification by total leukocyte RNA does not
rescue it. The positive control failed. Per-sample statements ("T > B in 94.8% of samples") are
true *of methylation* but are not a validated per-sample truth, and no method can be faulted for
failing to track them. `results/lymphoid_tracking.json`.

The cohort-level ordering is unaffected and independently supported: GBmap itself holds far more
T than B (pooled lymphoid T:B about 43:1).

### 7d · The extension panel inverts too

Two top-tier methods from Nguyen et al. 2024 added after registration, on the frozen reference:

| method | GBM: B above T | LGG: B above T | anatomic ACS |
|---|---|---|---|
| FARDEEP | **100%** of samples | **100%** | 0.9846 — above every registered method |
| ARIC | 87.5% | *pending* | *pending* |

The highest-ACS method in the whole study inverts the lymphoid compartment in every sample.

### 7f · Pooling T and NK does not rescue the ordering — suspect 4a is not sufficient

If the NK column's CD3 markers were draining T-cell signal (§4a), T + NK together should still
exceed B. Methylation puts B above T + NK in only 1.3% (GBM) and
1.5% (LGG) of samples. Methods placing B above T + NK in the majority of samples
(`results/lymphoid_pooled_tnk.json`, `scripts/lymphoid_pooled_tnk.py`):

| | frozen | raw/X (h5ad) |
|---|---|---|
| GBM, registered methods | **6 of 8** | 4 of 8 |
| LGG, registered methods | **9 of 12** | 3 of 10 |

FARDEEP places B above T + NK in 100% of GBM samples; ARIC in 77%.

**Reading.** The T-to-NK leak can contribute to the T under-call, but it cannot be the main
cause of the inversion: even crediting T with everything NK received, B wins for most methods.
The excess lands on **B specifically** — consistent with §1, where B/T exceeds truth in all ten
methods with a computable ratio while NK/T does so in only five to seven.

**What this is not.** Summing two fitted columns is not the same as refitting with one merged
T/NK column; collinearity between T and NK could change how much a merged column takes from B.
That refit is §5.1 — **now run, §7g: it does not fix it either.**

### 7g · The decisive test of suspect 4a was run — REJECTED

§5.1's experiment: re-annotate NK and refit. `scripts/nk_reannotation.py` →
`results/nk_reannotation.json`. The frozen reference, the registered NNLS and SVR, the TCGA
arm's own gene space (held fixed across variants) — and three rosters:

| | GBM SVR: B > T | LGG SVR: B > T | LGG NNLS: B > T | T share / B share (LGG SVR) |
|---|---|---|---|---|
| as registered | 90.2% | 96.5% | 98.5% | 0.001 / 0.936 |
| **NK column removed** | **100%** | **100%** | **100%** | 0.001 / 0.999 |
| **T and NK merged** | **100%** | **100%** | **100%** | 0.001 / 0.999 |

Methylation: T above B in 95.5% (GBM) and 93.2% (LGG); T + NK above B in
98.7% and 98.5%.

**The control passed first:** the as-registered refit reproduces the registered estimates exactly
(same samples scored, same B > T fraction, same shares to four decimals), so the variants differ
from it by the roster alone. Merge weights came from the atlas's training cells (T 40,784,
NK 1,746).

**Reading.** If the NK column were draining T-cell signal, deleting it would hand that signal back
to T. It does the opposite: the NK share goes to **B**, and the T profile gets almost nothing in
every variant. Suspect 4a is **rejected as a cause** — the sixth rejected mechanism. GBM NNLS is
uninformative here (lymphoid signal in 1 of 56 samples in every variant); LGG NNLS (66–67 samples)
agrees with SVR.

**What is left.** In every roster the B profile out-competes the T profile for the lymphoid
signal in real bulk, even though the same signature recovers a planted T:B exactly (§2). Suspect
4b — B's profile is closest to `Macrophage_Microglia`, the most abundant immune population — is the
one structural candidate left untested. §5.2 (swap the atlas) and §5.4 (planted T:B using T profiles
from another source) are the experiments that would test it.

### 7h · Three more candidate mechanisms in the B profile -- each tested, none survives

| # | candidate | test | result | artefact |
|---|---|---|---|---|
| 7 | Immunoglobulin (plasma-cell antibody) transcripts inflate B | count Ig genes in GBmap's solved marker space | **excluded by construction**: one Ig gene (JCHAIN) is in the space, carrying 0.057% of the B column's mass | `b_profile_tissue_likeness.json` |
| 8 | Ribosomal genes admitted as B markers -- the B column is 62% ribosomal mass in the solved space (13 RP genes) | refit with the 13 RP genes removed, same samples | **rejected in all four cohort-method cells**: SVR GBM B>T 90.2% -> 100% (51 -> 56 of 56 samples with lymphoid signal); SVR LGG 96.5% -> 99.8% (482 -> 510 of 510); NNLS GBM 100% -> 100% (1 -> 23); NNLS LGG 98.5% -> 98.6% (67 -> 352). Removing them *spreads* the inversion | `ribosomal_test.json` |
| 9 | Ambient myelin RNA captured with B cells | myelin-gene content per cell, raw/X counts | **rejected**: fewer B cells than T cells carry any myelin transcript (23.9% vs 39.2%; mean myelin fraction 0.014% vs 0.023%) | `atlas_cell_fractions.json` |

The ribosomal hypothesis was formed after the fact (B's ribosomal share was noticed first); it was
then tested with matched and any-gene random-removal nulls, and rejected by the test itself. The SVR
nulls were dropped only after the test had failed, because they guard against false positives
alone; the reason is recorded in the script.

**What is measured about GBmap's B profile -- properties, not causes.**

- It is the most bulk-tissue-like lymphoid profile: correlation with the mean bulk profile 0.345
  (T 0.152) in GBM, 0.265 (T 0.090) in LGG (`b_profile_tissue_likeness.json`).
- GBmap's B cells carry a median 42% of their molecules in ribosomal genes (T 26%), and 78% of them
  come from 3 of 31 donors (T: 47% from 3 of 77; NK: 87% from 3 of 39).
- **Where B's estimate goes when the column is removed depends on the estimator**
  (`b_column_diagnostics.json`). Under SVR (GBM, 124 samples with B) 92% of it moves to Tumor.
  Under NNLS it moves to T instead -- 182% of the removed B mass in GBM (8 samples with B; Tumor
  falls) and 109% in LGG (66 samples). Under SVR in LGG (465 samples with B, run 2026-10-02) 99% goes
  to Tumor, replicating GBM. **The pattern holds in both cohorts:** under SVR B is an overflow for
  tumour signal; under NNLS it competes with T.

### 7i · An independent atlas, pre-specified -- the result is MIXED

The reference was rebuilt from Abdelfattah *et al.* 2022 (GEO GSE182109), which is not one of
GBmap's 16 source studies: all 201,986 barcodes carry one of the series' own 44 GSM prefixes, none
from Neftel's GSE131928. Clusters were named by a marker rule fixed before any expression was seen
(`prespecified/abdelfattah_cluster_mapping.md`); 4,062 cells from 18 patients, B cells from 16 of
them. Methylation orders T > NK > B, T above B in 95.5% (GBM) and 93.2% (LGG) of samples.

| cohort | method | GBmap (frozen) | independent atlas | samples with lymphoid signal |
|---|---|---|---|---|
| GBM | NNLS | B>NK>T, B>T 100% | B>T>NK, B>T 69.2% | 13 / 56 |
| GBM | SVR | B>NK>T, B>T 90.2% | B>T>NK, B>T 89.5% | 19 / 56 |
| LGG | NNLS | B>NK>T, B>T 98.5% | **T>B>NK, B>T 5.2%** | 249 / 510 |
| LGG | SVR | B>NK>T, B>T 96.5% | B>T>NK, B>T 83.7% | 49 / 510 |

Read as pre-specified, **mixed**. Under SVR the inversion persists with an independent atlas in
both cohorts. Under NNLS in LGG -- the cell with by far the most lymphoid signal -- the independent
atlas recovers methylation's T above B. The atlas changes two things at once (B's bulk-likeness:
its profile is *anti*-correlated with mean bulk, -0.04 GBM / -0.13 LGG, where GBmap's is the most
bulk-like lymphoid profile; and B's content: plasma cells with immunoglobulin), so this result does
not isolate a mechanism. `results/independent_atlas_test.json`.

**Seed sensitivity (Addendum 2, declared after the main result; can only qualify it).** The build was
repeated with cell-sampling seeds 1, 2 and 3 (`results/independent_atlas_seed_sensitivity.json`):

| cell | registered | seeds 1 / 2 / 3 | reading |
|---|---|---|---|
| GBM NNLS | B>T>NK | T>B>NK / B>T>NK / B>T>NK | **seed-dependent**: loses its reading (7-17 samples with signal) |
| GBM SVR | B>T>NK | B>T>NK ×3 | seed-stable |
| **LGG NNLS** | **T>B>NK** | **T>B>NK ×3** (B > T in 2.9-10.1% of 208-248 samples) | **seed-stable** |
| LGG SVR | B>T>NK | B>T>NK ×3 | seed-stable |

The cell counts per type are identical across builds (the 50-per-patient-and-type cap binds only
where a group is larger), so the check probes the sampling of the abundant types more than B, most of
whose cells enter every build. The reading MIXED stands, and its two informative cells are both
seed-stable: SVR inverts in both cohorts, and NNLS recovers T > B in LGG.

### 7j · Immunoglobulin-removed refit, declared before the main result -- AGREES with it

The independent B column has 62 immunoglobulin genes in its marker space carrying 86.7% of its
mass. Declared at 11:53 on 2026-10-01, before the main result was seen: remove them (the same rule
used for GBmap) and refit the same samples, methods and comparison; if the two agree,
immunoglobulin content is not driving the result.

| cohort | method | main | immunoglobulin removed (1,145 genes) |
|---|---|---|---|
| GBM | NNLS | B>T>NK, B>T 69.2% (13 / 56) | B>T>NK, B>T 40.0% (5 / 56) |
| GBM | SVR | B>T>NK, B>T 89.5% (19 / 56) | B>T>NK, B>T 85.7% (7 / 56) |
| LGG | NNLS | T>B>NK, B>T 5.2% (249 / 510) | T>B>NK, B>T 19.8% (248 / 510) |
| LGG | SVR | B>T>NK, B>T 83.7% (49 / 510) | B>T>NK, B>T 64.3% (28 / 510) |

**The cohort-level ordering agrees in all four cells**, so by the declared reading immunoglobulin
content does not drive the independent-atlas result. What removal does change is resolution: fewer
samples keep any lymphoid signal (GBM SVR 19 -> 7, LGG SVR 49 -> 28), and the per-sample B > T share
moves by 4 to 29 points, in both directions. At 5 to 7 samples the GBM cells are too thin to read
per sample. `results/independent_atlas_test_ig_removed.json`.

### 7k · BayesPrism in its authors' configuration -- the pre-declared prediction is SUPPORTED, and the inversion stays

Genuine BayesPrism, unbudgeted, on TCGA-GBM's raw/X arm in two configurations. The rule was fixed
before either produced output (`prespecified/bayesprism_authors_prediction.md`):
- the **authors'** tutorial configuration: key = "Tumor", their gene cleanup, protein-coding genes;
- this project's **pipeline** driver: key = NULL, no cleanup.

| | authors' configuration | pipeline driver |
|---|---|---|
| lymphoid B share (the pre-declared quantity) | **0.582** | 0.704 |
| B above T, per sample | 80.4% | 87.5% |
| lymphoid ordering | B > NK > T | B > NK > T |
| purity rho (DNA, n = 154) | 0.718 | 0.719 |
| bias, tumour minus purity | -0.518 | -0.555 |
| R time | 12,177 s | 20,429 s |

**Verdict: SUPPORTED.** The authors' configuration lowers the B share as predicted. Per sample, it
is lower in 35 of the 56 samples with signal in both, higher in 21 (median -0.083; Wilcoxon p =
0.045). Modelling each tumour's own malignant expression takes part of B's excess away, which fits
B absorbing tumour signal.

**It does not fix the ordering.** B still exceeds T in 80% of samples, where methylation says
T > NK > B. So the tumour-overflow account explains part of the B excess under BayesPrism, not the
inversion. Its premise was SVR-specific (§7h), and that is how the paper should state it.

Two more findings from this pair:
- **The genuine package ranks tumour content far better than the reimplementation** in the
  registered TCGA rows: 0.72 against 0.16 on the same arm.
- **Both configurations under-call the level of tumour content by about half.** They rank well and
  calibrate badly -- the contrast-versus-level distinction again (§7b).

`results/extension/bayesprism_verdict.json`, `results/extension/bayesprism_tcga_gbm_*.json`.

**LGG, authors' configuration (2026-10-02; `results/extension/bayesprism_tcga_lgg_authors.json`).**
- Purity rho 0.626 on 510 samples, second-best of the study's 22 LGG tumour-content rows.
- Lymphoid B > NK > T, **B above T in 93% of samples** (B share 0.746), against methylation's
  T > NK > B.
- The authors' configuration does not remove the inversion in either cohort.
- **LGG pipeline driver (2026-10-02 21:48, 21,667 s):** purity rho 0.589, B > NK > T, B above T in
  73% (B share 0.467).
- **Pre-declared LGG verdict: NOT SUPPORTED.** The authors' configuration RAISES B's share in LGG (0.746
  against 0.467). Per sample it is higher in 390 of 510 and lower in 120 (median +0.215; Wilcoxon p ~ 0).
  GBM went the other way (0.582 against 0.704; p 0.045).
- **The prediction does not replicate.** Modelling each tumour's own malignant expression is not a
  general account of B's excess. Read GBM's support as cohort-specific. B stays above T under both
  configurations in both cohorts. `results/extension/bayesprism_verdict.json`.

### 7l · The loss decides the lymphoid ordering -- DESeq2 `unmix` and a post hoc ablation (2026-10-02)

**DESeq2 `unmix`** -- one of Nguyen 2024's seven top-10-in-all-scenarios methods, new to this study --
was run genuinely with a configuration fixed before output. It passed a planted-truth gate. On the
registered **frozen** signature it is the first method in the study to order **T above B**:

| | GBM frozen | LGG frozen | GBM raw/X | LGG raw/X |
|---|---|---|---|---|
| lymphoid ordering | **T > B > NK** | **T > B > NK** | B > NK > T | B > NK > T |
| B > T per sample | 23.8% | 48.2% | 73.7% | 83.3% |
| purity rho | **0.780** (best in the study) | 0.592 | 0.598 | 0.452 |
| samples with no lymphoid estimate | 14 / 56 | 261 / 510 | 37 / 56 | 414 / 510 |

So it meets the registered criterion (cohort-level T > B) on the frozen signature in both cohorts --
robustly in GBM, at a coin-toss per-sample split in LGG -- and loses it on the raw/X donor-level
reference in both cohorts. "No method robustly recovers T > B across reference builds" survives, now with one more
method that recovers it on one build only.

**Why it differs: an ablation inside the genuine package** (post hoc, rule fixed first:
`prespecified/unmix_loss_scale_ablation.md`, `results/unmix_loss_ablation_gbm.json`). `unmix`'s loss
is computed on log(x + shift). As shift grows the loss becomes linear, like NNLS and SVR.

| shift (power 1) | 1 | 10 | 100 | 1000 | 10000 | 1, power 2 |
|---|---|---|---|---|---|---|
| ordering | **T>B** | **T>B** | B>T | B>T | B>T | B>T |
| B > T per sample | 24% | 33% | 91% | 100% | 91% | 49% |
| purity rho | 0.780 | 0.784 | 0.805 | 0.780 | 0.748 | 0.775 |

**LGG, the same ablation** (`results/unmix_loss_ablation_lgg.json`). Its control sets no override, so the
per-cohort rule picks shift 0.5 exactly as in the reported row; that control reproduces the row exactly.

| shift (power 1) | 0.5 (pre-declared) | 1 | 10 | 100 | 1000 | 10000 | 1, power 2 |
|---|---|---|---|---|---|---|---|
| ordering | **T>B** | **T>B** | B>T | B>T | B>T | B>T | B>T |
| B > T per sample | 48% | 45% | 64% | 98% | 99.6% | 97% | 88% |
| purity rho | 0.592 | 0.602 | 0.613 | 0.620 | 0.612 | 0.592 | 0.612 |

In LGG the ordering flips sooner (by shift 10), and tumour-content tracking moves even less
(0.59-0.62).

Both controls reproduced the reported rows exactly. **Reading, in both cohorts: loss scale implicated** -- a strongly
stabilised (near-log) loss gives T > B and a near-linear loss gives B > T, with the same package,
samples, genes and reference. L1 matters too: L2 at the same scale loses T > B. Tumour-content
tracking barely moves (rho 0.75-0.81) across all six settings.

**What this means.** The bulk data pin down tumour content; they do not pin down the T-versus-B
split, which follows the analyst's choice of loss. That is an identifiability statement, and it
explains why no method recovers the ordering robustly: there is no single answer in the data for a
method to recover. For method development the lever is concrete -- how a loss weights low-abundance
genes against the dominant tumour genes -- and it can only be judged against an orthogonal truth,
never against anatomy.

### 7m · The lymphoid split is not identified by the bulk data -- Extension E2 (2026-10-02)

§7l showed, inside one package, that the loss scale flips T versus B while tumour content holds.
Extension E2 (`docs/EXTENSION_IDENTIFIABILITY.md`; rules written before computing) asks whether that
generalises across compartments, with DNA truths.
- **It does.** Loss-scale stability ranks compartments the way agreement with DNA does: rho 0.886
  (p 0.017) as registered, and 0.943 (p 0.008) with truth and estimate on matched denominators.
- **Tumour and all-leukocyte content are stable and track their truths.** The lymphoid compartment is
  the least stable in both cohorts. Its T/B/NK split tracks nothing against EpiDISH; against GIMiCC
  (tissue fractions), `unmix`'s T and lymphoid totals do track, the panel's median does not (§7n;
  EXTENSION_IDENTIFIABILITY §7.4).
- **Robust:** to bootstrap resampling, and to three design flaws disclosed and checked (Addendum 1).

**What this changes here.** The B-above-T inversion is best read as **non-identifiability**: the bulk
data do not determine the lymphoid split under this reference, so each method's own modelling choices
decide it (the loss scale for `unmix`, re-absorption differing by estimator in §7h). That is
consistent with nine mechanism tests finding no single cause (§7a-§7k); it is not shown to be the
reason. It does not say why most methods land on B rather than T; that remains open.

### 7n · Is the truth itself right? A second methylation instrument, and direct cell counts (2026-10-03)

Every section above takes methylation's T > B as the truth. That truth came from EpiDISH with a
**blood** reference, into which brain and tumour DNA can only be forced, and its per-sample split
had already failed a reference-free control (§7c). It was tested three ways. Rules were written first
(`prespecified/gimicc_truth_confirmation.md` with Addenda 1-6; `prespecified/atlas_t_vs_b.md`).

**1. GIMiCC [52], a methylation method built for glioma.** It separates tumour, neuronal, glial and
angiogenic DNA before it splits the lymphocytes.
- The genuine package reproduces its authors' example (54/54 values).
- **Registered reading: INCONCLUSIVE.** Every control passes except one: in LGG, shuffled CpG labels
  do not break the tumour layer (rho 0.648).
  - The cause is construction, not GIMiCC's estimate. Its IDH-mutant purity libraries are 98-99.7%
    tumour-hypermethylated CpGs, so a shuffled LGG sample reads its own global methylation, which
    tracks purity.
  - That explains the failure; it does not rescue it.

| | GBM | LGG |
|---|---|---|
| tumour layer vs ABSOLUTE | **0.935** | **0.630** |
| immune + microglia vs methylation leukocyte fraction | **0.959** | **0.921** |
| shuffled CpGs (must break the tumour layer) | -0.116 (breaks) | 0.648 (**does not**) |
| within-lymphoid T : NK : B | **0.462 : 0.261 : 0.277** | **0.405 : 0.102 : 0.493** |
| direction holds when 20% of CpGs are dropped | T > B in 18/20 | B > T in 19/20 |
| shuffled CpGs, lymphoid ordering | B > NK > T (reversed) | T > NK > B (reversed) |

- **GBM: T > B, with every control passing.** Shuffling the CpG labels reverses the ordering, so it
  is carried by the lymphoid CpGs themselves.
- **LGG: B > T.** It is robust to probe dropping, and shuffling reverses it, so it too is
  CpG-specific.
  - It is **not** B absorbing tumour signal: GIMiCC's B share falls as purity rises (rho -0.145,
    p 0.003), where an overflow would rise.
  - Whatever produces it, GIMiCC and EpiDISH disagree on T versus B in LGG.

**2. Direct cell counts: no methylation, no deconvolution.**
- **Abdelfattah et al. [42] (18 patients, not a GBmap source study).** The T/NK cluster was split per
  cell by the rule declared on 2026-10-01.
  - **T exceeds B in 18 of 18 patients**: 17,311 T, 1,302 NK and 1,205 B cells pooled.
  - The patient-mean shares are 0.733 : 0.213 : 0.054, **EpiDISH's T > NK > B**.
  - Its B cluster includes plasma cells, which can only inflate B.
- **GBmap's own core atlas (all glioblastoma):** 54,257 T against 1,250 B cells.
- **Flow cytometry of brain tumours [53]:** "the lymphocyte compartment was mostly composed of T cells
  with fewer NK cells and B cells".

**3. What follows.**
- **GBM, the discovery cohort: the truth is corroborated three ways.** A glioma-specific methylation
  instrument with every control passing, two single-cell atlases, and flow cytometry. The headline
  "no method reproduces T > B" stands on it.
- **LGG, the registered replication: the methylation truth depends on the instrument.**
  - EpiDISH says T > B; GIMiCC says B > T, and 12 of 12 frozen-arm methods (11 of 13 raw/X) agree
    with GIMiCC.
  - The direct counts favour T > B, but rest on only 2 LGG patients in the atlas and on pooled flow
    data.
  - **The LGG replication of the lymphoid finding rests on a contested truth.**
- **Both methylation instruments put B above every direct count; GIMiCC several-fold.** B share of lymphocytes: GIMiCC
  0.28 / 0.49, EpiDISH 0.11 / 0.20, atlas 0.05, GBmap about 0.02.
  - T > B is a fact about glioma tissue on every direct count. The methylation shares are not
    measurements of it at that precision.
- GIMiCC's per-sample split also fails the reference-free RNA control, as EpiDISH's did (GBM n.s.).
  In LGG it agrees weakly where EpiDISH disagreed (+0.10, p 0.035, against -0.12). Methylation stays a
  cohort-level truth.
- **What would settle LGG:** direct counts in IDH-mutant gliomas. Klemm et al. [53] Table S2
  (per-patient flow cytometry with IDH status) is the external item that would do it.

### 7e · What the paper should claim now (supersedes §6)

| claim | status |
|---|---|
| Under the frozen signature, 0 of 12 order T above B in either cohort | **state plainly** |
| Under the donor-level raw/X reference, 3 of 14 (GBM) and 2 of 13 (LGG) do | **state plainly — must be reported** |
| Bisque's agreement is its reference's composition, by its own declared assumption | **state plainly**, with the planted-truth test |
| No method orders T above B robustly across both references and both cohorts | **state plainly** |
| ACS cannot detect a level-replacement failure (it scores contrasts) | **state plainly** — it is a direct instance of the thesis |
| Methylation validates the cohort-level ordering only; per-sample tracking is INCONCLUSIVE | **state plainly** |
| The GBM truth (T > B) is corroborated by a glioma-specific methylation method and by direct cell counts | **state plainly** (§7n): GIMiCC passes every GBM control; atlas 16/16 GBM patients; GBmap 54,257 T vs 1,250 B; flow cytometry [53] |
| The LGG truth is settled | **must not be claimed** -- GIMiCC puts B above T in LGG (registered reading INCONCLUSIVE: its LGG negative control is uninformative by construction); direct counts favour T > B on 2 patients (§7n) |
| Methylation measures the lymphoid shares accurately | **must not be claimed** -- both instruments put B above every direct count, GIMiCC several-fold (§7n) |
| Suspect 4a (NK carries T markers) explains the inversion | **rejected** — removing or merging the NK column sends the mass to B (§7g); state it as the sixth rejected mechanism |
| Suspect 4b (B's profile is closest to macrophage) causes it | **must not be claimed** — untested |
| Immunoglobulin, ribosomal or myelin content of GBmap's B profile causes it | **rejected / excluded** (§7h) -- three more failed explanations, each with its control |
| B absorbs tumour signal | **estimator-specific**: true of SVR (92% of removed B goes to Tumor), not of NNLS (it goes to T) -- §7h |
| The inversion is specific to GBmap | **not established** -- an independent atlas gives a mixed result (§7i): SVR still inverts in both cohorts; NNLS recovers T > B in LGG |
| Immunoglobulin content drives the independent-atlas result | **rejected** by the sensitivity declared before the result (§7j) |
| Modelling tumour-specific expression (BayesPrism, key = Tumor) removes the inversion | **no** -- in GBM it lowers B's share as predicted (0.70 -> 0.58) but B still exceeds T in 80%; **in LGG it raises B's share (0.47 -> 0.75; the pre-declared prediction does not replicate)**, §7k |
| The inversion has one mechanism | **not supported** -- the split is not identified by the bulk data (§7m, E2); each method's modelling choices decide it. Why most land on B stays open |

The honest headline is stronger than "0 of 12": it survives the reviewer who finds the h5ad arm,
and it explains the one method that seems to get it right.

