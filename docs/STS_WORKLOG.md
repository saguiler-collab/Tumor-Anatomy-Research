# STS work log — read this first, every session

Started 2026-09-30 when the project moved from CJSJ to **Regeneron STS**. Maintained so that
nothing is lost when the conversation is compacted. **Update it after every major result.**
Last updated: 2026-10-01 12:15 EDT.

---

## 0 · How to resume (machine state)

Work runs in **two lanes**, one heavy job per lane, on an 8 GB / 2-core Mac.

| lane | what | status file / log |
|---|---|---|
| **A** | genuine BayesPrism on TCGA, **no time budget**: GBM authors → GBM pipeline → LGG authors → LGG pipeline (`scripts/bayesprism_tcga.py`) | `/tmp/bp_tcga.status`, `/tmp/bp_tcga_<cohort>_<variant>.log`, `results/extension/bayesprism_tcga_*_r.log` |
| **B** | strictly sequential: ① independent-atlas test (done) → ①b Ig-removed sensitivity → ② raw/X TCGA extension arm (paused job resumes) → ③ Ivy GAP re-runs of LinDeconSeq then RNA-Sieve → ④ ReCIDE gate then real run → ⑤ resume the two paused LGG SVR diagnostics | `/tmp/lane_b.status`, `/tmp/indep_test*.log`, `/tmp/queue.status`, `/tmp/q_*.log`, `/tmp/recide_*.log` |

Driver scripts live in the session scratchpad
(`/private/tmp/claude-501/.../scratchpad/`: `lane_b.sh`, `bp_tcga_chain.sh`, `lindeconseq_rerun.sh`,
`rnasieve_rerun.sh`, `recide_real_after_gate.sh`). They are detached (`nohup … & disown`); the
harness kills `run_in_background` jobs at 30 min, so long jobs must be detached.

Paused processes (SIGSTOP; resumed by lane B): the raw/X TCGA extension job, and
`b_column_diagnostics.py` / `ribosomal_test.py` (their LGG SVR fits).

**State at 2026-10-01 18:05 (supersedes the 16:10 note below).** Lane A unchanged (GBM authors, final
stage). Lane-B slot (`scratchpad/slot_b_now.sh`): old-analysis re-runs done (all REPRODUCED);
LinDeconSeq Ivy GAP re-run done (ACS 0.9692); **RNA-Sieve Ivy GAP re-run PAUSED** mid-build (pid
67082, memory thrash) -> `scratchpad/sieve_resume.sh` resumes it at "END gbm/authors"; the slot
chain resumes ARIC LGG only after RNA-Sieve finishes. `aric_resume.sh`, `lindeconseq_rerun.sh`,
`rnasieve_rerun.sh` were stopped (superseded by the slot chain; their log markers are the same).

**State at 2026-10-01 16:10.** Lane A: GBM authors' BayesPrism in its FINAL Gibbs stage (snowfall
worker, ETA 20:06); its Python parent still holds the atlas (2.2 GB, started before the free-atlas
fix). Lane B: the raw/X arm finished GBM (FARDEEP +0.524, LinDeconSeq +0.688, ARIC +0.418, RNA-Sieve
+0.719 purity rho; all B over T) and LGG FARDEEP (+0.351) and LinDeconSeq (+0.557); **ARIC LGG is
PAUSED at 265/510** (swap thrash: 100 s/sample vs 3) and resumes automatically when lane A's GBM
authors run ends (`scratchpad/aric_resume.sh`). Queued after lane B: the independent-atlas seed
sensitivity (`scratchpad/lane_b_tail_seeds.sh`, Addendum 2) -> then run
`scripts/independent_atlas_seeds.py`. Verdict watcher: `scratchpad/bp_verdict_watch.sh` runs
`scripts/bayesprism_verdict.py` after each cohort's authors + pipeline pair.

Check everything:
```
cat /tmp/lane_b.status /tmp/bp_tcga.status /tmp/queue.status
ps -axo pid,stat,%cpu,time,command | grep -E "bayesprism|extension_tcga|independent_atlas|ribosomal|b_column|recide|remeasure"
sysctl -n vm.swapusage
```
Known competitors for the machine: Microsoft Edge (user closed it 2026-10-01), **Google Drive**
(was using ~70% CPU; user asked to pause it; still at ~40% CPU on 2026-10-01 15:20), Safari tabs, and the 7.2 GB duplicate
`GEO_GSE182109/Processed_matrix.mtx.gz.download/matrix.mtx` (safe to delete; the build reads the .gz).

---

## 0.5 · What is reportable now -- re-assess after every result (user: "always try to see if we have reached a meaningful conclusion that can be originally reported")

*Assessed 2026-10-01 14:30. Every number is in an artefact named in §2.*

**Registered, complete, both cohorts -- this is the paper, and it is reportable today:**
1. **Anatomy detects broken estimators.** Both broken-input controls fail their within-tumour
   permutation nulls; every real estimator beats its null (p < 1e-4); ACS margin 0.31.
2. **Anatomy did not rank estimators by accuracy.** ACS vs DNA-measured purity: rho = 0.081
   (n = 12, p = 0.80; registered bar 0.60 not met). Say it the manuscript's way: the rank
   correlation is UNDERPOWERED (95% CI -0.64 to 0.70; with 12 methods a true 0.7 would meet the
   criterion only 65% of the time, §2.11), so it is "not met", not a refutation; the second clause rests on
   the registered lymphoid criterion (3). Consistent across instruments: a graded within-tumour
   AUC gives rho = -0.06 (§2.10), and a protein-level test (IMC) saw ACS reverse the ordering.
3. **Getting the anatomy right is not getting the biology right.** No estimator robustly reproduces
   methylation's T above B: 0/12 on the registered signature (GBM and LGG); 3/14 and 2/13 on the
   raw/X rebuild, one of them by returning its reference's composition.

**Original, post-registration, labelled exploratory -- reportable as mechanism:**
4. Bisque's apparent biological success *is* its reference composition (planted truth: L1 to the
   prior 0.000, to the truth 1.091); ACS cannot see this by construction (it scores contrasts).
5. The B-over-T inversion survives 9 tested explanations; reabsorption is estimator-specific (SVR ->
   Tumor, NNLS -> T); an independent atlas gives a MIXED answer (SVR still inverts in both cohorts,
   NNLS recovers T > B in LGG); immunoglobulin removal changes resolution, not direction.

**Re-assessed 2026-10-02, after E2 -- a sixth item, original and post-registration:**
6. **A truth-free check identifies which compartments the bulk data determine** (§2.19). Re-fitting
   one genuine package under different loss scales, the compartments that stay put are the ones DNA
   confirms (rho 0.886, p 0.017, 6 units; robust to sampling and to both design flaws).
   - Tumour and total leukocyte content are identified; the lymphoid split is not.
   - Stability flags what NOT to trust; it certifies nothing: T in GBM is stable and wrong.
   - This is the most original result of the extension. It is exploratory and rests on 6 units.
   - Robust to sampling and to three design flaws (duplicate fit, degraded methods, truth denominators;
     the matched version is stronger, rho 0.943).
   - **Agreement between methods (DECEPTICON's criterion) is a weaker signal:** it tracks accuracy
     on average (0.28), but its top pairs are algorithmic siblings (9/10), so it can select shared
     construction over correctness.

**What can still change, and what cannot.** Still running: BayesPrism authors' configuration
(pre-declared prediction about B), ReCIDE, the LGG SVR diagnostics, the extension agreement at
n = 16. These can add or qualify exploratory detail in 4-5; none can move 1-3, which are
registered and closed. ~~Open ask: the SCP accession~~ -- resolved 2026-10-01 by verifying the processed matrix against
GEO's raw archive (§2.13).

## 1 · Standing user directives

- **Regeneron STS**, not CJSJ. "Make me a top researcher." Never stop unless the project needs a
  different direction.
- **No time budgets.** BayesPrism must actually finish. Use the official tutorial the user added
  (`pipeline_packages / repos/BayesPrism-main/tutorial_deconvolution.pdf`).
- **Always ask for external items** (data, packages, repos, papers, a faster machine) — now a section
  in `CLAUDE.md` and a memory. Lead with the single most valuable item.
- **One by one, efficiently** — two lanes, never more than ~2 heavy jobs.
- The **prose is the user's**: the manuscript generator writes facts and `> WRITE:` beats, not sentences
  to keep.
- **Nothing committed or pushed** this session (CLAUDE.md: only when asked).

---

## 2 · Findings, in order (every number is in an artefact)

### 2.1 D22 — every genuine-package re-measurement ran on GBmap's LOG layer (FIXED)
`scripts/remeasure_method.py` never passed `matrix=`; `build_from_h5ad` defaulted to `X` (log1p).
Bisque 0.9231 and EPIC 0.9846 were bit-for-bit their `acs_log`. Fix: `--matrix` read from
`results/run_provenance.json` + an **8th equivalence condition** that blocks on a mismatch. Old
artefacts kept in `results/superseded_D22/`. Re-runs on raw/X:
- bisque **0.7077**, epic **0.7692** — reproduce their (already-genuine) leaderboard rows exactly.
- dwls **0.7000** on 26 pairs: **73/122 samples fail inside DWLS's quadratic program**
  ("constraints are inconsistent, no solution!"), identical on an independent re-run
  (`results/dwls_remeasured_with_stdout.json`). Not comparable to the reimplementation.
- **bayesprism 0.8154 [0.6999, 0.9286]**, all 122 samples, 4.6 h of R time — vs reimplementation
  0.8000 (**+0.0154**, inside the CI). Finished only because the Python parent was SIGSTOPped so the
  4 h timer could not fire (Popen checks child exit before the clock).
- r_bridge now keeps R's stdout on success (`results/diagnostics/<method>_r_stdout.log`).

### 2.2 D23 — the lymphoid headline is reference-specific (reporting defect)
- "0 of 12" is the **frozen** signature. On the raw/X donor-level (h5ad) reference **3/14 (GBM)**
  and **2/13 (LGG)** order T>B — the manuscript had omitted this.
- **Bisque's agreement is its reference's composition, by its own declared assumption** (Jew 2020,
  Nat Commun 11:1971, Methods p.10, quoted). Truth-free test `scripts/bisque_anchoring.py`: cohort-mean
  L1 to reference prior 0.087/0.081 (others 0.68–1.40). Planted truth
  `scripts/bisque_anchoring_synthetic.R`: L1 to prior **0.000**, to truth 1.091, per-sample rho
  0.88–0.97. **ACS is blind to this by construction** (it scores contrasts).
- **Methylation's per-sample lymphoid split failed a reference-free positive control**
  (`scripts/lymphoid_tracking.py`): marker index vs methylation T/(T+B) LGG rho −0.123, p 0.004;
  pre-declared PTPRC strata do not rescue. Methylation is a **cohort-level** truth only.
- Abstract, findings, Discussion and title notes now qualified ("no method robustly…").

### 2.3 D24 — `build_from_h5ad` silent default (FIXED)
No default layer now (raises). Five legacy callers fixed (`run_benchmark.main`,
`reference_sensitivity`, `build_gbmap_assay_reference`, `remeasure_dwls` → raw/X;
`reconstruct_gene_space` → `X` deliberately). Registered run unaffected (`run_all.py` always passed it).
Negative-control test added.

### 2.4 Why B above T — mechanism tests (WHY_B_OVER_T §7)
| # | explanation | verdict | artefact |
|---|---|---|---|
| 1–5 | (pre-existing) | rejected | S5 |
| 6 | NK column absorbs T signal (CD3) | **REJECTED**: NK removed or T+NK merged → B>T in 100% (SVR GBM & LGG, NNLS LGG) | `nk_reannotation.json` |
| — | pooled T+NK vs B | B > T+NK in most samples for 6/8 (GBM), 9/12 (LGG) registered methods | `lymphoid_pooled_tnk.json` |
| 7 | immunoglobulin genes | **EXCLUDED BY CONSTRUCTION** in GBmap: only JCHAIN, 0.057% of B | `b_profile_tissue_likeness.json` |
| 8 | ribosomal genes admitted as B markers (B column 62% RP in solved space) | **REJECTED**: removal → SVR GBM B>T 90%→100% | `ribosomal_test.json` |
| 9 | ambient myelin RNA | **REJECTED**: B cells carry less myelin than T (24% vs 39% of cells) | `atlas_cell_fractions.json` |

Measured properties (not causes): GBmap's B profile is the **most bulk-tissue-like** lymphoid profile
(r with mean bulk 0.345 vs T 0.152 GBM); removing B under SVR sends **92% of its estimate to Tumor**
(GBM); GBmap B cells: median RP fraction 0.42 vs T 0.26; **78% of B cells from 3 donors** (NK 87%).
Working hypothesis: B is an overflow for patient-specific tumour signal → tested by BayesPrism
`key="Tumor"` (pre-declared prediction: authors' config LOWERS the B share).

### 2.5 Independent atlas — Abdelfattah 2022, GSE182109 (pre-specified)
- User supplied Single Cell Portal files in `GEO_GSE182109/`: `Processed_matrix.mtx.gz` (38,224 genes
  × 201,986 cells, 542.7 M nonzeros), `genes.tsv`, `barcodes.tsv`, `Cluster_GBM.txt` (C1–C12), plus
  `gsm_to_sample.tsv` (from GEO). **Not one of GBmap's 16 source studies.**
- Matrix is **Seurat log1p(count/lib × 1e4)** → inverted exactly to counts (max dev 2.2e-15; totals vs
  library 7.9e-16). Streamed twice (`scripts/build_abdelfattah_reference.py`).
- Clusters named by a rule fixed **before any data** (`prespecified/abdelfattah_cluster_mapping.md`):
  C1/C4/C7 macrophage, C2/C6/C8/C9 tumour, **C3 T (NK inside)**, **C11 B (+plasma)**, C5 oligo,
  C10 endothelial, C12 dropped. 4,062 cells, 18 patients (11 ndGBM, 5 rGBM, 2 LGG); B from 16 patients.
- B column **87% immunoglobulin** in its marker space; 0 RP genes. Ig-removed sensitivity declared in
  an addendum **before** the main result was seen.
- **Main result** (`results/independent_atlas_test.json`): GBM NNLS B>T 69% (13 samples), GBM SVR 89.5%
  (19), **LGG NNLS T>B>NK, B>T in only 5.2% of 249 samples** (GBmap 98.5%), LGG SVR 83.7% (49).
  Independent B profile is *anti*-correlated with bulk (−0.04/−0.13). Pre-specified reading: **MIXED**.
- **Ig-removed sensitivity: DONE — AGREES** (`results/independent_atlas_test_ig_removed.json`): same
  cohort-level ordering in all 4 cells (GBM NNLS 40.0% B>T on 5/56; GBM SVR 85.7% on 7/56; LGG NNLS
  T>B, 19.8% on 248/510; LGG SVR 64.3% on 28/510). Pre-declared reading: Ig content does not drive it;
  removal costs resolution, not direction. Written up: WHY_B_OVER_T §7h–7j; manuscript auto-sentence.
- **Correction of an earlier summary:** B-removal reabsorption is estimator-specific — SVR sends 92% to
  Tumor (GBM), NNLS sends it to T (182% GBM, 109% LGG; Tumor falls). "B = tumour overflow" is SVR-only.
  **LGG SVR (run 2026-10-02 06:30): 99% to Tumor (465 samples)** -- replicates in both cohorts.

### 2.6 Post-registration extension panel (`ivygap/deconv/extension.py`, MANUSCRIPT §4.10)
Genuine packages only, no fallback, never in a registered statistic; each passed a planted-truth gate.
| method | route | anatomic ACS | TCGA frozen: purity rho / B>T (GBM; LGG) | raw/X GBM |
|---|---|---|---|---|
| FARDEEP | CRAN 1.0.1 (permn=0, blocked-parallel; both verified identical) | **0.9846** | +0.700/100%; +0.543/100% | +0.524, B>T 100%, 49/56 zero-lymphoid |
| ARIC | PyPI 1.0.1 | 0.8154 | +0.685/87.5%; +0.575/77.4% | +0.418, 98.1% |
| LinDeconSeq | GitHub @20f1aec (deconSeq w/ our signature) | **0.9692** [0.906, 1.000] (re-run 17:28, after the renaming fix) | +0.732/93.6%; +0.553/97.9% | +0.688, 88.9%, 47/56 zero |
| RNA-Sieve | PyPI 0.1.4 + disclosed SciPy x0 shim, serial | **0.8769** [0.786, 0.971] (re-run 19:08, after the integer-check fix) | n/a (needs cells) | running |
| ReCIDE | vendored @31bbd6b + dplyr summarize_each shim + v3 assays | **gate PASSED** 2026-10-01 21:05 (12 planted mixtures × 4 types: max abs error 0.059, per-type Pearson 0.997–0.999, 2,692 s); real run **first real run KILLED 06:05 2026-10-02 by the driver's 6-h default budget (re-run queued, unbounded)**; it had started 21:05 on 1 core; donor rule chose **k = 20** of 58 eligible training donors (the k = 10 and 15 draws left a roster type at <= 20 cells), reference 45,934 cells, 651 genes | — | — |
| DESeq2 unmix | DESeq2 1.52.0 `unmix` (shift = 1 by the pre-declared flatness rule, |rho| 0.027) | **0.9077** [0.836, 0.982] (2026-10-02 06:55; gate PASSED first) | TCGA runs queued | — |
| BayesPrism (authors) | tutorial config, key="Tumor", cleanup.genes | — | lane A running | — |
Citations: 3 of 4 first written from memory were wrong; all now Crossref-verified with DOIs (test).
ReCIDE's benchmark (Li 2026) is from ReCIDE's own lab — note it.
Exploratory agreement (`scripts/extension_agreement.py`): registered reproduced (rho 0.081, n=12);
+FARDEEP rho 0.182, n=13; **+FARDEEP, ARIC, LinDeconSeq rho 0.237, n=15 (perm p 0.39, CI -0.35..0.72;
2026-10-01 17:30)** -- still below the bar; RNA-Sieve cannot enter (no frozen-arm result: needs cells).
LinDeconSeq anatomic re-run (genuine, raw/X): **ACS 0.9692 [0.906, 1.000]**, 57 pairs, ties MuSiC.

### 2.7 Recovered analyses ("lost things") → MANUSCRIPT §4.11 + Limitations
Data was never lost (only two shelved raw files, intact as .gz). What was lost were **completed analyses
that never reached the paper** (47/71 artefacts unused). Restored:
- Constraints validated with no deconvolution: **Albiach 3/4** (C5 tested on the wrong side of the
  necrosis boundary), **Darmanis C1 4/4 gates**, **ISH** CD44/BIRC5 support C1/C7.
- **Pre-registered IMC protein test**: ACS *reversed* the protein ordering (MCP_GBM 0.5128 vs MCP_GBmap
  0.6667; IMC r 0.37 vs 0.06) — a third instrument says ACS does not rank (uses Ajaib's published r).
- **CDSeq reference-free** ACS 0.846 / 0.683.
- Published benchmarks (Avila Cobos top tier): Mann-Whitney p = 0.125.
- Thorsson immune arm: pre-registered over-call P1/P2 **failed** (2/8; median error negative).
- Reference dependence: rho 0.51 (Neftel-only) / 0.37 (Darmanis-only) — **both are GBmap
  constituents**, so not independent (wording corrected).
- Limitations: survival BLOCKED (no vital status; MGMT-informative missingness, Fisher p 0.0021);
  B-cell donor concentration; methylation cohort-level only; extension post-registration.
- Re-run on raw/X: `benchmark_concordance`, `design_factors` (old copies in `results/superseded_pre_rawX/`).

### 2.8 BayesPrism (user's explicit request)
Official tutorial (TCGA-GBM example) prescribes: `cleanup.genes` (Rb, Mrp, other_Rb, chrM, MALAT1, chrX,
chrY), protein-coding only, `key = "tumor"`, default Gibbs. Our driver did none of these.
`R/run_bayesprism_authors.R` implements it; smoke-tested. GBM authors run: 1,601 → 1,251 genes (78
removed: 1 Rb, 20 other_Rb, 58 chrX, 5 chrY), Gibbs ETA 6h14m (first stage) from 11:23; second stage
follows. LGG ≈ 2–3 days per variant on this machine. The authors' variant now frees the atlas
before the R call (saves ~2 GB; applies from the GBM pipeline run on).

**GBM authors' configuration DONE (2026-10-01 18:09, exit 0, 12,177 s):** purity rho **+0.718**
(the Python reimplementation on the same raw/X arm: +0.164 — the genuine package ranks tumour content
far better), bias **-0.518** (it under-calls tumour content badly in level), lymphoid **B>NK>T, B > T in
80.4%** of 56 matched samples, no zero-lymphoid samples, L1 to the reference prior 0.919 (not
anchored). Methylation: T>NK>B. The pre-declared verdict needs the pipeline variant (running from 18:09).

**GBM pipeline DONE (2026-10-02 04:53, exit 0, 20,429 s):** purity rho +0.719, bias -0.555, B>NK>T,
B > T in 87.5%. **Pre-declared verdict (GBM): SUPPORTED** -- B share authors 0.582 vs pipeline 0.704;
per sample 35 lower / 21 higher (median -0.083, Wilcoxon p 0.045). The inversion persists (B>T 80%).
Written up: WHY_B_OVER_T §7k. LGG authors' run started 04:53 (days; frees the atlas before R).

**LGG authors' configuration DONE (2026-10-02 12:49, exit 0, 20,618 s):**
- **Purity rho +0.626**, second of the 22 LGG tumour-content rows in the study (RNA-Sieve raw/X 0.637
  is first). Bias -0.472.
- **Lymphoid B>NK>T, B > T in 93.3%** of 510 matched samples (lymphoid mean T 0.062, B 0.746, NK
  0.192), against methylation's T>NK>B. No zero-lymphoid samples. L1 to the reference prior 0.810, so
  not anchored.
- **The inversion persists in both cohorts** under the authors' configuration.
- **LGG pipeline DONE 21:48** (21,667 s reported; wall clock 12:49–21:48, which includes ~2 h 20 m
  paused while the lid was shut): purity rho 0.589,
  B>NK>T, B>T 72.6% (B share 0.467).
  - **Pre-declared LGG verdict: NOT SUPPORTED.** B share is 0.746 under the authors' configuration
    against 0.467 under the pipeline; higher in 390/510 samples (median +0.215, p ~ 0).
  - **The prediction holds in GBM and reverses in LGG.** It does not replicate.
  - Written up: WHY_B_OVER_T §7k and its claims table; MANUSCRIPT §4.10 (from the artefact).

**Extended agreement re-run with `unmix` (2026-10-02 12:52; `results/extension/agreement_extended.json`).**
`unmix`'s anatomic ACS (0.908) post-dated the n = 15 artefact.
- Registered control reproduced (0.081).
- **n = 16: rho 0.2375 (prints +0.237), permutation p 0.37, 95% CI -0.32 to 0.69** (n = 15:
  0.2366).
- Exploratory, reported and not tested. RNA-Sieve has an ACS but no frozen-arm purity, so it is excluded.

### 2.9 Data inventory — every public dataset, measured (user request 2026-10-01)
`scripts/build_data_inventory.py` → `docs/DATA_INVENTORY.md` + `docs/supplementary/S7_data_sources.csv`
(+ cache `docs/supplementary/data_inventory_cache.json`). 30 entries: 22 used, 1 derived, 3 on disk
unused, 1 duplicate, 3 considered. Per file: bytes, sha256, recorded-hash match; per dataset: the
scripts that read it (token found in code), Crossref check of the citation, re-check vs the public
server (16 files: 12 identical size incl. GBmap 8,127,057,404 B; 4 identical decompressed sha256).
Large files hashed 2026-10-01 17:00 (`--hash-large`, duplicates skipped): GBmap h5ad sha256
`459444da84ef…` (full value in `docs/supplementary/data_inventory_cache.json`).
**0 open problems.** Found and fixed on the way:
- **D25: GBmap was cited to the wrong paper** (REFERENCES [10] = Ravi 2022's DOI under GBmap's authors)
  **and the wrong CELLxGENE collection** (Siletti's `283d65eb…`). True: collection `999f2a15…`,
  dataset `861acfd8…` (read from the h5ad's own uns/citation); paper Ruiz-Moreno 2025 Neuro-Oncology
  27:2281, doi 10.1093/neuonc/noaf113 (preprint 10.1101/2022.08.27.505439). Fixed + negative control.
- Albiach check's hard-coded path broke when the vendored tree was renamed → `config.VENDORED_REPOS_DIR`.
- Unused on disk: Thorsson CIBERSORT table (DATA_SOURCES claimed a use that was never built), ABSOLUTE
  segtabs, LGG GDC clinical/biospecimen. Duplicates: 6.7 GB `matrix.mtx` (safe to delete), etc.
- **Abdelfattah processed files are NOT on GEO** (GSE182109 suppl = RAW tar only) → Single Cell Portal;
  the SCP study accession was never recorded → **ask the user**.
- `tcga_lgg_bulk_provenance.json` says "TCGA-GBM" in what_this_is (label slip; cohort field says lgg).

### 2.10 ACS vs AUC — differentiated, measured (user request 2026-10-01)
No per-method AUC existed (predecessor's AUCs are 6/12-month survival AUCs per model variant A–D,
0.42–0.60; this project's only AUC was 0.831, MES score vs Verhaak call). `scripts/anatomic_auc.py`
→ `results/anatomic_auc.json`: AUC_anat = within-tumour Mann–Whitney AUC on the same constraint–tumour
pairs, samples and permutation null as ACS (fast path == readable definition, tested). Controls: ACS
reproduced for all 17 rows; registered agreement rho 0.081 reproduced. Results: ACS vs AUC rank rho
**0.94** (14 comparable); chance **0.50 vs ~0.37**; AUC gives partial credit on C4/C7 (conjunctions in
ACS); negative controls 0.533 (p 0.20) / 0.413 (p 0.99). **AUC in ACS's place: rho −0.060 vs ABSOLUTE
(p 0.85)**, −0.004 vs C_purity → the null result is not an artefact of thresholding. C_purity
(Harrell's C vs ABSOLUTE): svr 0.779, cibersortx 0.780 … bisque 0.508 (chance; Bisque in its
no-overlap mode). On TCGA's frozen signature MuSiC = NNLS and SCDC-ENSEMBLE = SCDC (degenerate) —
S8 and Figure 9 label them so. Manuscript §4.11
paragraph auto-built. 12 tests in `tests/test_anatomic_auc.py` (incl. inverted signal, Simpson's paradox).
A test caught my own wrong claim (identity holds for pairwise constraints only, not C4/C7) — fixed.

### 2.11 Power of the registered agreement test — "underpowered" now has a number
`scripts/agreement_power.py` → `results/agreement_power.json` (sensitivity power analysis for effect
sizes fixed in advance — NOT observed power). Criterion simulated exactly as registered (rho ≥ 0.60
AND percentile bootstrap CI over methods excludes 0; vectorised bootstrap reproduces the registered
one exactly). Probability the criterion is met:

| true rank correlation | 12 methods (ours) | 20 | 30 | 50 |
|---|---|---|---|---|
| 0.0 (false positive) | 1.7% | 0.5% | 0.0% | 0.0% |
| 0.6 (= the bar) | **43%** | 48% | 51% | 47% |
| 0.7 | **65%** | 75% | 77% | 86% |
| 0.8 | **83%** | 95% | 97% | — |

Reading: with 12 methods a true correlation of 0.7 is missed 1 time in 3; a true value exactly at
the bar is met ~half the time at ANY panel size (the criterion thresholds the estimate at the bar),
so "80% power at 0.60" is impossible by design and must never be written. 80% power at 0.70 needs
30–50 methods. Manuscript Discussion point 5 is generated from the artefact. 5 tests.

### 2.12 Re-running old analyses to check them (user, 2026-10-01 17:20: "rerunning some old stuff to check")
Lane-B slot used while ARIC is paused (`scratchpad/slot_b_now.sh`): each re-run writes to scratch,
then `scripts/compare_artefacts.py` (new, 4 tests) diffs it leaf by leaf against the archive.
- `darmanis_constraint_check`: **REPRODUCED** (128 leaves). `albiach_constraint_check`: **REPRODUCED**
  (and proves the path fix). `ish_constraint_check` (current script): **REPRODUCED**.
- I first re-ran the SUPERSEDED `fetch_ivygap_ish.py` by mistake — its output differs because it is
  the retired method, not because the archive drifted. It defaulted to the live artefact's path —
  now `results_superseded/` (OPEN_DEFECTS D25 addendum). Inventory D18 now credits the right script.
- The user renamed `GEO_GSE182109/` -> `GSE182109/` (17:15). `config.GSE182109_DIR` accepts both;
  the Abdelfattah builder goes through it (the queued seed check would otherwise have failed).
- Then, in the same slot: LinDeconSeq and RNA-Sieve Ivy GAP anatomic re-runs (pulled forward), and
  ARIC resumes only after both these and lane A's GBM authors run have finished.
- User asked about GEO GSE231746: private, embargoed until 2027-05-05 (owner/reviewer only); not
  referenced anywhere in the project; not usable.

### 2.13 GEO GSE182109 raw archive (user added `GSE182109/GSE182109_RAW/`, 2026-10-01 19:20) — provenance CLOSED
`scripts/verify_gse182109_raw.py` → `results/gse182109_raw_verification.json` (+ 4 tests). All 44 GSMs:
complete Cell Ranger triplets (132 files, 2.5 GB, sha256 each), names match GEO's GSM table, dims agree;
Cell Ranger 5.0.0 ×19 / 5.0.1 ×25 (two 36,601-gene feature lists; their union = the processed matrix's
38,224 genes). 264,951 raw barcodes → 201,986 processed cells (76.2% retained); **201,893 found in their own
sample's raw barcodes**; the 93 others (GSM5518600 ×69, GSM5518609 ×24) are absent from Cell Ranger's
filtered set — a more permissive cell call, not a mix-up (stray matches elsewhere = chance collisions).
**Exactness (first 3 samples, 16,457 cells, 34.5 M nonzeros): every recovered count EQUALS GEO's raw
count (max |diff| 7e-12), nonzero counts identical (nothing dropped), library = raw UMI total (7e-16).**
→ The independent atlas rests on GEO GSE182109; the unrecorded SCP accession no longer matters (only
`Cluster_GBM.txt`, the cluster labels, is SCP-only — and its barcodes are verified GSE182109 cells).
User's GSE231746 question: private on GEO until 2027-05-05; not in the project or its 19 papers.

### 2.14 Cleanup (user: "delete unnecessary files"), 2026-10-01 20:25 — 10.0 GB freed, nothing unique lost
`scripts/dedupe_verified.py` deletes ONLY a file proven identical (sha256) to a kept original or to its
decompressed stream; log with hashes and re-create commands: `docs/DELETED_FILES.md`.
- Deleted (identical): LGG HM450 decompressed copy (1.62 GB), Darmanis csv copy (0.18 GB), Neftel GSM3828672
  Smart-seq2 (0.98 GB; byte-identical to the SCP-named `IDHwtGBM.processed.SS2.logTPM.txt` the build reads).
- `GSE182109/.../matrix.mtx` (7.21 GB) was **NOT** a decompressed copy (my earlier notes were wrong): it was a
  TRUNCATED decompression (ends mid-line at cell 98,613 of 201,986; 7.21 of 16.56 GB) and an exact byte-prefix
  of the kept .gz's stream → deleted. The builder's `assert seen == dims[2]` would have caught a truncated read.
- NOT deleted (genuine datasets read by nothing; re-downloadable): Neftel 10x TPM (1.56 GB), GSE84465_RAW
  (0.24 GB), ABSOLUTE segtabs (0.25 GB), LGG GDC clinical TSVs. Listed for the user to decide.
- Inventory fixed on the way: D15 now lists the Neftel file the build actually reads.
Lane B: raw/X TCGA arm DONE 20:20 (LGG: FARDEEP +0.351, LinDeconSeq +0.557, ARIC +0.443, RNA-Sieve +0.637
purity rho; all B over T); ReCIDE gate running; the real ReCIDE run now uses 1 core
(`scratchpad/recide_real_1core.sh`, replaces the 2-core trigger; same log/marker).

### 2.15 The user believed ACS was an accuracy measure -- corrected; tumour factors; new methods (2026-10-02)
- **User said:** "all this time I've been thinking [ACS is] a tumor microenvironment cell composition
  accuracy measure". It is NOT: it is agreement with pre-registered anatomic orderings. Full explainer,
  worked example from tumour 703393 (SVR: 8/8 constraints, yet B>T in 90% of TCGA GBM):
  `docs/WHAT_ACS_IS.md`. Also corrected: deconvolution can stand in for single-cell / flow / IMC
  composition measurement, NOT for CT/MRI/ultrasound imaging.
- **Tumour factors** (`prespecified/tumour_acs_factors.md` → `scripts/tumour_acs_factors.py` →
  `results/tumour_acs_factors.json`; per-tumour ACS reproduced first). ACS minus chance per tumour
  +0.30..+0.57 (all tumours well above chance). No factor associated (all Holm p = 1.0; n = 9 →
  INCONCLUSIVE). The pre-declared prediction (RNA anatomy strength, F1) was NOT supported (rho -0.23).
  **Negative control finding:** sparse sampling lets even random controls beat chance (control rho with
  n structures -0.81, p 0.016) → per-tumour ACS is least reliable in sparsely sampled tumours.
- **Methods research** → `docs/METHOD_CANDIDATES.md`. Nguyen 2024 verbatim top-10-in-all-scenarios:
  DWLS, DESeq2 unmix, MuSiC, CIBERSORT, ARIC, LinDeconSeq, MIXTURE — 5/7 already covered; **DESeq2 unmix
  and MIXTURE were not**. Best reference-free: Linseed. Top-15 extras: MethylResolver, AdRoit, Scaden.
  Spatial review (NRG 2025): Cell2location, RCTD (off-label on LCM regions).
- **DESeq2 unmix added** (installed DESeq2 1.52.0): config fixed first
  (`prespecified/deseq2_unmix_config.md`: shift by automated mean-SD flatness, power 1);
  `R/run_deseq2_unmix.R`; wired into `extension.py` + `r_bridge.py`; **planted-truth gate PASSED**
  (`tests/test_deseq2_unmix_gate.py`: min per-type r >= 0.90, max |err| <= 0.10; shuffled-gene-label
  reference FAILS the bar). Real runs queued at the end of lane B.
- **Full `pytest tests/` is OWED** after the `ivygap/deconv/` change (only the extension tests were run:
  9 passed) — run it when a lane frees; never alongside two heavy jobs.
- **Lane B scheduling bug fixed:** old step 5 resumed both LGG diagnostics at once and declared lane B
  complete, which would have started the seed build alongside them. `scratchpad/lane_b_cont.sh` now
  runs ReCIDE → diagnostic 32054 → diagnostic 36099 → seeds → DESeq2 unmix, one at a time.
- External asks for the user (one at a time): MIXTURE (github elmerfer/MIXTURE), Linseed
  (github ctlab/LinSeed), InstaPrism, MethylResolver, AdRoit, spacexr.

### 2.16 LGG diagnostics complete; Figure 10 (2026-10-02 06:20)
Both paused LGG SVR fits ran (one at a time): **B removal under SVR in LGG → 99% to Tumor (465 samples)**,
replicating GBM's 92% (NNLS → T in both cohorts); **ribosomal removal under SVR in LGG → B>T 96.5% → 99.8%**
— mechanism 8 rejected in all four cohort×method cells. `Figure_lymphoid_mechanisms` rendered (4 panels;
legends/labels fixed after one look; panel B states SVR and that NNLS sends B to T) and registered as
Figure 10 (caption from artefacts; manuscript placement after the independent-atlas paragraph). Rebuild:
10 of 10 figures; checkers CLEAN. Lane B now on the seed sensitivity (started 06:17).

### 2.17 Independent-atlas seed sensitivity (Addendum 2) -- the informative cells are SEED-STABLE (2026-10-02 06:53)
Seeds 1-3 rebuilt in one streaming pass (identical cell counts per type; the cap binds only on large groups,
so B is mostly the same cells), refit: GBM SVR, **LGG NNLS (T > B; B>T 2.9-10.1% of 208-248)**, LGG SVR
SEED-STABLE; **GBM NNLS SEED-DEPENDENT** (seed 1 → T>B>NK; only 7-17 samples with signal) -- loses its
reading. MIXED stands; its informative cells hold. WHY_B_OVER_T §7i; manuscript sentence auto; test pinned.

### 2.18 DESeq2 unmix on TCGA, and the loss-scale ablation -- an identifiability finding (2026-10-02 07:10)
unmix (genuine, pre-declared shift = 1): **frozen GBM T>B>NK (B>T 23.8%), purity rho 0.780 (best in the study),
bias -0.09; frozen LGG T>B>NK (B>T 48.2%, 261/510 zero-lymphoid), rho 0.592; raw/X GBM B>NK>T (73.7%), rho 0.598;
raw/X LGG B>NK>T (83.3%), rho 0.452.** First method to meet the registered T>B criterion on the frozen signature
in both cohorts; not robust across reference builds (raw/X inverts) -- "no method robustly recovers" survives.
**Loss-scale ablation (post hoc; rule fixed first; `prespecified/unmix_loss_scale_ablation.md`;
`results/unmix_loss_ablation_gbm.json`):** shift 1, 10 → T>B; 100, 1000, 10000 → B>T (B>T 24 → 33 → 91 → 100 → 91%);
power 2 at shift 1 → B>T. Control reproduced. **Reading: LOSS SCALE IMPLICATED** (+ L1 needed). Purity rho
0.75-0.81 throughout. → **Tumour content is identified by the bulk data; the T-vs-B split follows the loss.**
WHY_B_OVER_T §7l. LGG ablation running.

### 2.19 Extension E2 -- truth-free identifiability diagnostics (user: "let's get this extension down", 2026-10-02)
Design + literature: `docs/EXTENSION_IDENTIFIABILITY.md`; binding rules: `prespecified/identifiability_diagnostics.md`
(written before computing). D1 = unmix loss-scale stability (7 settings); D2 = registered-method agreement;
truths = ABSOLUTE (Tumor), Thorsson LF (Leukocytes), EpiDISH (Lymphoid; T/B/NK secondary). H1: stability ranks
the 6 primary units as truth agreement does (SUPPORTED needs rho >= 0.83 at 6 units); H2: leukocytes identifiable
(D1 >= 0.8 and rho >= 0.4, both cohorts); S1: log-data NNLS vs linear NNLS (Avila Cobos's case).
Literature anchors (verbatim): Sturm 2019 multicollinearity/spillover; Nguyen 2024 "unstable results or multiple
solutions"; Li 2026 truth-free reproducibility benchmarking; NRG 2025 "ground truth ... often not available";
Avila Cobos 2020 "linear scale" (tension resolved: transform the LOSS, not the data).
Script `scripts/identifiability_diagnostics.py` queued after the ablation re-runs (`scratchpad/e2_queue.sh`).
**RESULT (2026-10-02 08:47; `results/identifiability_diagnostics.json`; controls reproduced exactly):**
- **H1 SUPPORTED:** D1 rho 0.886, exact p 0.0167 (6 units); D2 rho 0.657, p 0.088.
- **H1s (12 units):** D1 0.587 (p 0.024); D2 0.518 (p 0.044).
- **H2 PASS in both cohorts:** leukocyte stability 0.954 / 0.927. `unmix` vs the methylation LF:
  0.780 / 0.603, above EVERY registered method (best nu-SVR 0.688 GBM; DWLS-reimpl. 0.590 LGG,
  margin 0.013 not tested).
- **S1 expectation FAILED:** log1p-data NNLS ranks purity better than linear (0.715 vs 0.562;
  0.557 vs 0.275). It is a contrast of criteria with Avila Cobos (RMSE on simulated mixtures), not a
  refutation.
- **Units:** Tumor and Leukocytes stable (>= 0.93) and accurate; Lymphoid least stable (0.713 / 0.494)
  and inaccurate (-0.20 / -0.18). **T in GBM is stable (0.893) and does not track (-0.12; -0.006 matched)**: stability is not
  sufficient; 36% of samples are tied at zero.
- **Robustness (Addendum 1, declared before computing; `results/identifiability_robustness.json`):**
  - contrast alone has null P = 1/15;
  - near-tie swap keeps p = 0.029;
  - bootstrap B = 2000: lymphoid lowest in 100% (D1) / 98.8% (D2); H1 rho 95% interval 0.54-1.00;
  - **design flaw 1:** GBM `predeclared` = `s1`, byte-identical (pre-declared shift is 1). Over
    distinct fits H1 is unchanged (0.886, p 0.0167);
  - **design flaw 2:** D2 kept degraded Bisque and EPIC. On the immune arm's 8-method panel, D2 rho is
    0.771, p 0.051 (still n.s.).
- **Negative control** (`tests/test_identifiability.py`, 8 tests): shuffled truths gave SUPPORTED in 2 of
  20 shuffles (median p 0.67).
- **Write-up:** EXTENSION_IDENTIFIABILITY §7 and §7.1; MANUSCRIPT §4.12, Discussion beat 2d, beat 4
  qualified, Appendix D (+4 lines); Figure 11 `Figure_identifiability`; REFERENCES [45] DESeq2
  (Crossref).

### 2.20 FARDEEP's default path crashed (OPEN_DEFECTS D26) -- found by the frozen-arm regeneration (2026-10-02, fixed 09:11 by file time)
REGEN (`extension_tcga.py --methods fardeep lindeconseq aric deseq2_unmix --references frozen`) ran with
`IVYGAP_FARDEEP_CORES` unset. `run_fardeep.R` called `cut(x, 1)`, which R rejects, and both cohorts
failed. Every earlier run had exported 2 workers.
- **Fixed:** one block = joint fit.
- **Test:** `tests/test_fardeep_default_path.py` (3 genuine-package tests incl. 1-vs-2-worker identity and
  the shuffled-reference control). Shown to FAIL on the pre-fix driver.
- **Re-run queued** on the default path after REGEN: `scratchpad/fardeep_regen.sh` → `/tmp/fardeep_regen.log`,
  marker FARDEEP_REGEN.
- **REGEN exit 0 (12:23). Every regenerated summary is IDENTICAL, field by field, to the pre-REGEN
  backup** (`/tmp/tcga_*_frozen_before_regen.json`):
  - ARIC, LinDeconSeq and unmix in both cohorts;
  - **FARDEEP GBM on the fixed single-worker path: 0.6997, identical**;
  - FARDEEP LGG: re-run in progress.

  The regenerated per-sample files are therefore validated. GBM frozen holds 4 methods x 154; LGG
  frozen 3 x 510 plus FARDEEP pending. The raw/X per-sample gap (FARDEEP, LinDeconSeq, ARIC) remains
  and needs an atlas re-run.

---

### 2.21 The 2022–2026 literature, read in full; check 7; DECEPTICON's criterion tested (2026-10-02 09:15–12:00, file times)
- **Nine papers read in full text** (Europe PMC OA; `scratchpad/papers_modern/`): DREAM 2024 [47],
  omnideconv 2026 [48], DLPFC benchmark 2025 [49], DECEPTICON 2025 [46], Xu 2025, Ivich 2025, DeconX
  2026, Sutton 2022, Dai 2024. Written up in EXTENSION_IDENTIFIABILITY §2b.
  - **Four truth-free criteria in use:** agreement, reproducibility, fit/residuals, stability. None had
    been validated per compartment against orthogonal per-sample truth in tumours.
  - **omnideconv 2026 (lung, IHC):** "all methods struggled to robustly estimate tumor-infiltrating
    lymphocytes". DREAM: coarse types detected at 3–4%, CD4 subsets at 6% best case.
- **Check 7 (a third E2 design flaw):** EpiDISH (blood reference) gives shares of the immune
  compartment, and E2 compared them with tissue fractions. On matched denominators
  (`results/identifiability_denominators.json`, reproduced on re-run):
  - lymphoid non-agreement STANDS (-0.103 / -0.123);
  - lymphoid share less stable (0.552 / 0.341);
  - **H1 strengthens: D1 0.943 (p 0.008), D2 0.829 (p 0.029)**;
  - matched D2 in LGG rests on 263 of 510 samples (BayesPrism-reimpl. and NNLS return zero leukocytes
    for 201 and 156).
- **DECEPTICON's selection rule vs DNA truth** (`prespecified/agreement_selection_test.md`, written
  first; `scripts/agreement_selection_test.py`). P1 = agreement-selected ensemble vs random pairs; P2 =
  centrality vs accuracy across methods.
  - The first run was STOPPED at the 30-min background limit: it re-mapped keys inside every draw.
    Rewritten on aligned arrays; logic unchanged.
  - A 20-draw smoke run was seen before the registered 2000-draw run. P2 does not depend on draws.
  - **RESULT (2000 draws, `results/agreement_selection_test.json`):**
    - **P1 INCONCLUSIVE:** the selected ensemble beats random pairs in 4/6 units (10-method panel 5/6).
    - **P2 SUPPORTED, weakly:** centrality ~ accuracy, mean rho 0.28, p 0.035 (10-method 0.27, p 0.023).
    - **The top pair was an algorithmic-sibling pair in 9 of 10 selections**: `cibersortx` subclasses
      `svr`; the BayesPrism reimplementation iterates NNLS. Agreement = shared construction.
    - **LGG leukocytes:** ensemble 0.172 (= random) vs DWLS-reimpl. 0.590, which ranks 5th of 8 by
      centrality.
    - GBM lymphoid truths are all within noise (n = 56), so that unit is uninterpretable.
    - LGG lymphoid share: ensemble 0.395 (99th percentile, above every single method) -- weakly
      recoverable.
  - Write-up: EXTENSION_IDENTIFIABILITY §7.2; MANUSCRIPT §4.12 paragraph.
  - **Two of my sentences failed verification and were corrected:** "every secondary value within
    ±0.27" (SCDC-reimpl. reaches 0.39) and "DWLS agreed with no one" (it ranks 5th of 8).
- **References [45]–[49]** added, Crossref-verified; the list now has 49.

### 2.22 S2, the anatomy arm (Addendum 2 written 12:55, run 12:56–13:04; file times)
- **Run:** `remeasure_method.py --method deseq2_unmix --e2-s2` (new opt-in mode). The inputs are built
  once, all equivalence conditions pass, and output goes only to `results/identifiability/`. The guard
  test refuses other methods.
- **Control:** the pre-declared setting reproduces the reported row exactly (ACS 0.9077, estimates
  within 1e-16). Ivy GAP's pre-declared shift is also 1, so there are 6 distinct fits.
- **ACS by setting:** 0.908, 0.908 (s1), 0.877 (s10), 0.939 (s100), 0.939 (s1000), 0.969 (s10000),
  0.877 (p2). Every setting beats its null.
- **Reading INCONCLUSIVE:** range 0.092, between the 0.05 and 0.10 bounds. Only C6 (macrophage
  MVP > CT) and C7 (the tumour gradient) change; C1–C5 are identical.
- **Descriptive, NOT pre-declared:** the methylation-consistent fits (T>B in TCGA-GBM) score
  0.908 / 0.877, the most inverted fits 0.939 / 0.939 / 0.969; Spearman(ACS, B>T share) = 0.67 over 6
  fits. Anatomy does not move toward the right lymphoid answer.
- **Ivy GAP stability over 6 fits:** constrained types 0.944–0.983; T 0.760, B 0.546, NK 0.551 (88–93%
  exact zeros).
- **Written up:** EXTENSION_IDENTIFIABILITY §7.3 and MANUSCRIPT §4.12. **Caught before publishing:** my
  first criterion for "reproduces T>B" used a per-sample share < 0.5 and counted the power-2 fit, whose
  cohort ordering is B>T. Fixed to the registered cohort-level flag.

### 2.23 ReCIDE finished on Ivy GAP; the full test suite passes (2026-10-02 13:24–13:38)
- **ReCIDE (unbudgeted, 1 core, 14,568 s; `results/extension/recide_anatomic.json`):** ACS **0.692**
  [0.548, 0.829] on 57 pairs, null p 1e-4, no types dropped.
  - **The lowest anatomic ACS of any real method in the study** (registered comparable 0.708–0.969;
    extension 0.815–0.985). It is not broken by ACS's own criterion.
  - Li et al. 2026 [4] rank ReCIDE top on truth-free reproducibility. They are its developers ("our
    recently developed ReCIDE"). **Two truth-free criteria disagree.** Adjudication needs DNA truth:
    `prespecified/recide_tcga.md` (born 13:32:45) registers a TCGA raw/X run.
- **ReCIDE TCGA-GBM raw/X DONE 21:49** (18,461 s of R; controls 1 and 2 passed):
  - purity rho **0.522**; bias -0.069; B>T>NK at cohort level, but B>T in only 46.4% of samples;
    L1 to the prior 1.027.
  - **Adjudication (pre-declared): rank 7 of 8 GBM raw/X rows, bottom third.** Anatomy's low ranking
    agrees with DNA here; the reproducibility ranking does not.
  - Li et al.'s other pick, genuine BayesPrism, ranks 1st and 3rd of 8. One method, one cohort.
  - In MANUSCRIPT §4.10.
- **Full `pytest tests/ -q`: 417 passed, 1 skipped, 0 failed** (5 min 45 s, 13:38). The first run
  (13:30) had 2 ERRORs in the new `test_identifiability.py`: its real-data fixture read
  `config.RESULTS_DIR` at call time, after conftest's per-test isolation had re-pointed it. Fixed by
  capturing the real directory at import, as `test_anatomic_auc` does, and restoring config after a
  read-only load.
- **ReCIDE on TCGA raw/X, GBM: LAUNCHED 13:40:19** (`scripts/recide_tcga.py`; lane script
  `scratchpad/recide_tcga_lane.sh`, pid 27898; caffeinate 27899).
  - Status: `/tmp/recide_tcga.status`; log: `/tmp/recide_tcga_gbm.log`.
  - Raw ReCIDE output is saved the moment R returns: `results/extension/recide_tcga_gbm_raw.csv`.
  - **Design:** ReCIDE goes through `extension_tcga.run` unchanged via an adapter. `run_recide.py` was
    refactored into `donor_draw` / `recide_reference` / `call_recide`, shared by both runs (diff
    checked; the argument block is identical).
  - **Controls:** (2) the draw reproduces the anatomic run's 20 donors -- PASSED in the live run and
    beforehand; (1) the arm reproduces `deseq2_unmix`'s raw/X row field by field -- running.
  - **Tests:** `tests/test_recide_drivers.py` (3).
  - **Expected:** ~5 h (Ivy GAP's 122 samples took 4 h). LGG afterwards, ~17 h, if the machine allows.

### 2.24 MIXTURE: installed, gated, running (user: "go ahead", 2026-10-02 15:30)
- **Install.** MIXTURE 0.0.1, genuine, vendored at `vendor/MIXTURE/` (github.com/elmerfer/MIXTURE
  @ 6332e160, 2024-08-19; COMMIT, SHA256SUMS, PROVENANCE.md; MIT per DESCRIPTION, no LICENSE file).
  Installed with `R CMD INSTALL` from that directory.
  - Plotting/IO dependencies installed first as binaries: ComplexHeatmap 2.28.0, ggExtra 0.11.0,
    openxlsx 4.2.9.
  - Plus new packages shape, rjson, circlize, GetoptLong, clue, GlobalOptions, colourpicker, shinyjs,
    zip. No existing package changed version (`update = FALSE`).
- **Configuration** fixed before any output: `prespecified/mixture_config.md` (born 15:34:57). Package
  defaults; proportions get the central cell-size conversion.
- **Citation** verified on Crossref: Fernandez EA et al., Brief Bioinform 2021;22:bbaa317,
  doi:10.1093/bib/bbaa317. It is Nguyen 2024's ref 105.
- **Wiring:** `R/run_mixture.R` (asserts sample and cell-type names, never maps by position);
  `ExtensionSpec("mixture")`; `R_PACKAGES`; `SIGNATURE_ONLY_METHODS`. **Unbudgeted** (`R_METHOD_TIMEOUTS`).
  I added that only after launch, having NOT checked the 1-h default first.
- **Gate (`tests/test_mixture_gate.py`): PASSED.**
  - Planted min r 0.996, max error 0.027.
  - Shuffled reference: MIXTURE returns NO estimate (all NaN), which fails the bar.
  - **The first run of the gate "failed" its negative control** because NaN makes both comparisons
    False. Fixed in both gate tests with `_bar()`, which treats NaN as no recovery. The `unmix` gate had
    the same latent bug. The criterion is unchanged.
- **Tests:** bridge and extension tests 133 passed, 1 skipped. **Full suite with MIXTURE registered
  (lane B, 21:56): 422 passed, 1 skipped, 0 failed.**
- **Runs:**
  - TCGA frozen started 15:50 (`/tmp/mixture_tcga_frozen.log`). That process still carries the old
    1-h limit; re-run LGG if it times out.
  - **TCGA frozen DONE** (GBM 16:03; LGG 19:52 after the relaunch). Other methods' rows are identical
    to the backups.
    - **GBM:** purity rho 0.712, B>T>NK, B>T 100%, 21/56 with no lymphoid estimate.
    - **LGG:** purity rho 0.548, B>T>NK, B>T 100%, 164/510 with no lymphoid estimate.
    - **MIXTURE inverts the lymphoid ordering in both cohorts, like the registered panel.**
  - **Ivy GAP anatomic DONE 21:59** (`results/extension/mixture_anatomic.json`):
    - ACS **0.954** [0.891, 1.000] on 57 pairs, null p 1e-4, 0 unestimated; all 7 hard equivalence
      conditions pass.
    - Top tier on anatomy, with FARDEEP 0.985 and LinDeconSeq 0.969.
  - **TCGA raw/X DONE 22:35.**
    - GBM: purity rho 0.661, B>T>NK, but **51/56 have no lymphoid estimate**.
    - LGG: purity rho 0.471, **492/510 have no lymphoid estimate**.
    - Its raw/X "B above T 100%" rests on 5 and 18 samples.
  - **Extended agreement, n = 17** (MIXTURE added): rho 0.285, p 0.27, 95% CI -0.24 to 0.70.
    Exploratory; registered control 0.081 reproduced.
  - Queued in lane B after ReCIDE-GBM (`scratchpad/laneB2_queue.sh`, pid 32649, `/tmp/laneB2.status`):
    full pytest, then MIXTURE Ivy GAP anatomic, then MIXTURE TCGA raw/X, then ReCIDE TCGA LGG.

### 2.25 Linseed: installed, gate running (user: "continue", 2026-10-03 07:25)
- **Configuration** fixed first: `prespecified/linseed_config.md` (born 07:29:05), plus an addendum on
  the signature scale written before any output.
  - The README tutorial step for step; k = 8 a priori, as CDSeq's T; topGenes stays at the
    tutorial's 10,000.
  - CDSeq's two labellings reused; the gate bar is matched r >= 0.80.
  - Real runs wait for ReCIDE-LGG (memory).
- **Install.** linseed 0.99.3 vendored at `vendor/linseed/` (github.com/ctlab/LinSeed @ 2435a8ce,
  2025-08-22; MIT with LICENSE). Compiled from that directory; build products removed afterwards;
  sources re-verified against SHA256SUMS.
  - Imports installed first: GEOquery 2.80.0 and combinat 0.0.9.
  - **Also added as dependencies:** clipr, vroom, tzdb, selectr, rlang 1.3.0, readr, rentrez, rvest,
    httr2. **rlang may have replaced an older version (not recorded).** Every method package
    re-checked to load.
- **Citation** Crossref-verified: Zaitsev K et al., Nat Commun 2019;10:2209,
  doi:10.1038/s41467-019-09990-5.
- **Driver** `R/run_linseed.R` asserts sample names.
- **Gate (`tests/test_linseed_gate.py`): PASSED.**
  - Planted 4 types / 40 samples: matched min r **1.000**; 1,346/3,000 genes pass Linseed's filter.
  - Per-gene-shuffled control: an estimate IS returned, matched r **0.137**; 15/3,000 genes pass. It
    fails on substance, not by error.
- **Run script** `scripts/linseed_run.py`.
  - Arms: anatomic, and tcga with --cohort.
  - Reuses CDSeq's labellings unchanged.
  - The lymphoid ordering is reported only if T, NK and B are all claimed; never zero-filled.
- **Smoke test** on synthetic data: the profile labelling named exactly the 4 planted types; the marker
  labelling added one spurious NK component. Files removed.
- **Queued** after ReCIDE-LGG (memory): `scratchpad/linseed_queue.sh` (pid 67420, `/tmp/linseed.status`),
  running anatomic, then TCGA-GBM, then TCGA-LGG.

- **Real runs DONE 2026-10-03 13:21:03** (anatomic 12:57-13:03, TCGA-GBM 13:03-13:11, TCGA-LGG
  13:11-13:21; all exit 0). Descriptive, as registered:
  - **Ivy GAP anatomic (122 samples; 1,189 of 10,000 genes pass the filter):**
    - profile-named: ACS **0.677** [0.578, 0.789], p 1e-4, 57 pairs; covers Tumor, Endo, Mac, Oligo;
    - marker-named: ACS **0.659** [0.600, 0.725], p 0.0056, 41 pairs; covers Astro, Endo, Mac, Oligo;
    - CDSeq on the same test: 0.846 profile-named, 0.683 marker-named.
    - Both reference-free instruments beat their nulls; Linseed is lower on profile naming.
  - **TCGA-GBM (154 samples):** tumour component against ABSOLUTE rho **0.365** (profile-named) and
    0.089 (marker-named).
    - Marker naming claims 2 T and 2 B components but no NK, which is implausible at < 1% lymphoid.
    - No lymphoid ordering is reported, by the rule.
  - **TCGA-LGG (510 samples):** tumour rho **0.210** (profile-named) and **-0.226** (marker-named).
    No lymphoid type is claimed.
  - Manuscript: a data-driven sentence after CDSeq's paragraph in `build_manuscript.py`.
  - Anatomic and TCGA are separate fits, so "meets the constraints, tracks purity weakly" is
    descriptive, not a within-fit test.

### 2.26 Confirming the lymphoid truth with GIMiCC (user incorporated it, 2026-10-03 12:38: "spend as much time as you need confirming the truth")
- **Why.** The headline's T > B truth is EpiDISH with a BLOOD reference on brain tumour DNA (D23).
  GIMiCC (Pike et al., Acta Neuropathol Commun 2024) is glioma-specific and hierarchical:
  InfiniumPurify tumour layer, then neuronal / glial / angiogenic / immune, then lymphoid NK / B / T,
  then CD4 / CD8. Outputs are tissue fractions.
- **Software.** The user's copy `pipeline_packages / repos/GIMiCC-main` is byte-identical to GitHub @
  26cb8a15 and was installed with `R CMD INSTALL`.
  - Dependencies: ExperimentHub 3.2.2, FlowSorted.Blood.EPIC 2.16.0, InfiniumPurify 1.3.1, minfi
    1.58.0 (+32). **36 added, 0 version changes** (recorded before and after).
  - GIMiCC calls `query()` without importing it, so the driver attaches ExperimentHub, as the
    vignette does.
- **Reproduction control: PASSED.** The vignette example (EH9482, 3 GBMs) matches 54/54 printed values.
- **Rules:** `prespecified/gimicc_truth_confirmation.md` (born 12:51:36), plus Addendum 1 (12:56)
  after the first run's diagnostics and before any estimate was opened.
  - Controls: tumour layer vs ABSOLUTE >= 0.40; Immune + Microglia vs LF >= 0.40; shuffled CpGs
    |rho| < 0.20; structural invariance.
  - Q1 = cohort-level within-lymphoid T vs B (the registered criterion's computation).
  - Robustness: 20 repeats dropping 20% of probes; the direction must hold in >= 18.
- **Coverage.** All 4,022 library CpGs are present in Xena's 450K matrices, but 406 (GBM) / 450 (LGG)
  have missing values (masked probes), dropped by the complete-case rule, never imputed. Per-layer
  loss is 10-21%; **L3A (T/NK/B) loses 46/266**.
- **Code:** `R/run_gimicc.R`, `scripts/gimicc_truth.py` (resumable). Runs launched 12:57:58 (lane
  `scratchpad/gimicc_lane.sh`, `/tmp/gimicc.status`, log `/tmp/gimicc_runs.log`).
- **ReCIDE TCGA-LGG raw/X DONE 12:55:58:** purity rho 0.316, T>B>NK (T>B True; B>T in 38.6%),
  0/510 with no lymphoid estimate. **LGG rank (descriptive, outside the GBM-only rule), computed 13:29:
  9 of 9** raw/X rows against ABSOLUTE (rho 0.316), consistent with GBM's 8 of 9 (bottom third). In the
  manuscript via `_recide_lgg_rank()` in `build_manuscript.py`.

- **GIMiCC RESULT (runs 12:58-13:33, analysis 13:35:36; `results/gimicc_truth_confirmation.json`,
  `results/gimicc_secondary.json`, `results/gimicc_diagnostics.json`).** Addenda 2-6 were each
  written before the quantities they govern were computed (times from `stat`: 13:03:51, 13:09:23,
  13:17:11, 13:36:28, 13:40:27).
  - **Execution defect, caught and fixed (no rule changed).** The first analysis (13:34:00) matched
    **0 samples** in the tumour-vs-ABSOLUTE and shuffled-CpG controls: ABSOLUTE's IDs keep the vial
    letter (`-01A`), the methylation IDs do not, and the key function kept it.
    - An n = 0 control was scored as a FAIL, giving that run's INCONCLUSIVE.
    - Fix: key both sides with `k4s`. A guard now blocks any control that matches fewer than 100
      samples.
    - The defective output is kept as `results/gimicc_truth_confirmation.superseded_keybug_133400.json`.
    - **Q1's values were visible before the fix.** Q1 does not use purity, so they are identical.
  - **Controls:**

    | control | GBM | LGG |
    |---|---|---|
    | tumour layer vs ABSOLUTE (>= 0.40) | 0.935 (n 146) PASS | 0.630 (n 519) PASS |
    | immune + microglia vs LF (>= 0.40) | 0.959 (n 152) PASS | 0.921 (n 527) PASS |
    | shuffled CpGs (\|rho\| < 0.20) | -0.116 PASS | **0.648 FAIL** |
    | structural invariance | 2e-15 PASS | 1e-14 PASS |
    | probe-drop robustness (>= 18/20) | T>B 18/20 | B>T 19/20 |

  - **Registered Q1 reading: INCONCLUSIVE** (LGG negative control). The symmetric reading is the same.
    - GBM: **T > B > NK** (0.462 / 0.261 / 0.277), 148 samples; T > B in 68% of samples.
    - LGG: **B > T > NK** (0.493 / 0.405 / 0.102), 431 samples; T > B in 39%.
    - Level 3 gives the same orderings: GIMiCC zeroed level-4 T in only 5 LGG samples (plus 1 NaN)
      and 0 GBM samples.
  - **Why the LGG negative control fails (Addendum 5, D-c, plus the library):**
    - GIMiCC's IDH-mutant purity libraries are almost all tumour-hypermethylated CpGs: AST 996/1000,
      AST-HG 979, OLG 997; GBM's is balanced, 465/1000.
    - With shuffled labels the LGG tumour layer reads the sample's mean beta (rho 0.922; GBM 0.025),
      and mean beta tracks purity (0.576).
    - So the control is uninformative in LGG by construction. That explains the failure; it does not
      rescue it.
  - **Diagnostics (exploratory):**
    - **D-a: shuffling REVERSES both orderings** (GBM T>B>NK -> B>NK>T; LGG B>T>NK -> T>NK>B). Both
      are carried by which CpGs are which, not by global methylation.
    - **D-b:** the B share falls slightly with mean beta (LGG -0.139).
    - **D-d:** GIMiCC's T-dominant samples have a higher RNA T-minus-B index (GBM p 0.043, LGG p 0.019).
    - **D-e: not a tumour overflow.** The B share FALLS with purity (LGG -0.145, p 0.003; GBM -0.05).
  - **Q2:** under GIMiCC's direction the headline is unchanged in GBM (0/12 and 3/14 match). In LGG it
    is "reversed": 12/12 (frozen) and 11/13 (raw/X) methods match GIMiCC's B > T. This inherits
    GIMiCC's INCONCLUSIVE.
  - **Q3 (reproduction control exact): "FAILS".** `unmix`'s lymphoid total tracks GIMiCC's lymphoid
    tissue fraction.
    - rho 0.647 (GBM, n 56) and 0.526 (LGG, n 509); **partial on the leukocyte fraction 0.364 and
      0.316 (D-f)**.
    - The methods' median does not track: 0.077 / -0.247, partial 0.121 / -0.047.
    - unmix's T against GIMiCC's T: 0.643 / 0.491; B 0.253 / 0.296; NK 0.258 / 0.284.
    - H1 with GIMiCC's lymphoid truth: D2 rho **1.00** (p 0.0014); D1 0.771 (p 0.051).
    - E2's "the lymphoid estimate does not track methylation" must be narrowed. It holds against
      EpiDISH's share of leukocytes and for the methods' median. **It does not hold for unmix against
      GIMiCC's tissue fraction.**
  - **Q4: FAILS by rule.** GBM rho +0.185, p 0.15, n 62. LGG rho +0.100, p 0.035: weakly in the right
    direction, where EpiDISH's was significantly opposite (-0.123).
  - **What this means:**
    - **In GBM the truth is corroborated three ways:** a glioma-specific methylation instrument with
      every control passing and a CpG-specific T > B; direct counts (atlas 16/16 GBM patients; GBmap
      core, all glioblastoma, 54,257 T vs 1,250 B); flow cytometry (Klemm [53]).
    - **In LGG the methylation truth depends on the instrument.** EpiDISH says T > B; GIMiCC says B > T
      (robust, CpG-specific, not tumour overflow, but with an uninformative negative control). Direct
      counts favour T > B (atlas LGG patients 14 vs 0 and 118 vs 3; Klemm pooled), but only 2 LGG
      patients.
    - Both methylation instruments put B above the direct counts (GIMiCC several-fold): B share GIMiCC 0.28 / 0.49,
      EpiDISH 0.11 / 0.20, atlas 0.05-0.06, GBmap about 0.02.
    - **The LGG replication of the lymphoid finding therefore rests on a contested truth. The GBM
      discovery does not.**
  - **Provenance gap found and closed.** The CpG extraction (library union and Xena rows) had existed
    only as one-off commands.
    - It is now `scripts/extract_gimicc_cpgs.py`, with paths from the new `config.TCGA_HM450`.
    - `--verify` checks it reproduces the inputs byte for byte.
  - Literature: GIMiCC's lymphoid layers were validated by its authors on blood-cell mixtures only
    (GSE182379), not in tumour tissue [52].

- **Also, 13:41-13:50 (file times):**
  - GBmap per-donor counts (`scripts/gbmap_t_vs_b.py` -> `results/gbmap_t_vs_b.json`, obs only):
    all glioblastoma; T above B in **96 of 98** donors with lymphocytes; B absent in 57.
  - The data inventory registers GIMiCC as **D31**; the new scripts are listed as consumers of
    D07/D08, the atlas and GBmap.
  - Real-artefact tests (`tests/test_gimicc_truth.py`, `tests/test_atlas_t_vs_b.py`, 20 pass):
    controls n >= 100, and the reading recomputed from the controls. Shown to FIRE on the defective
    13:34 artefact.
  - The manuscript, WHY_B_OVER_T §7n, EXTENSION_IDENTIFIABILITY §7.4, OPEN_DEFECTS D23, METHODS,
    REFERENCES [52] [53] and PROJECT_ACCOUNT are updated. Every lymphoid "tracks nothing" claim now
    names its instrument.

- **D27 found by the independent recomputation (13:52-13:59, file times).**
  - `tests/test_gimicc_independent.py` (plain pandas; ABSOLUTE matched through its own `array`
    column) found more samples than the analysis: 146 against 142 (GBM), 519 against 510 (LGG).
  - Cause: 795 ABSOLUTE rows carry Broad-internal IDs in `sample` (`GBM-TCGA-06-5416-Tumor-SM-1QETM`).
  - **The registered yardstick lost 1 GBM and 9 LGG samples.**
    - Audit (`scripts/absolute_join_audit.py`): every recorded rho reproduces exactly; the lost
      samples were never deconvolved.
    - Worst-case shift of the registered agreement: +0.081 -> +0.056 (bar 0.60). Conclusion unchanged.
    - OPEN_DEFECTS D27.
  - GIMiCC now keys purity by `array`. Controls moved <= 0.013 (GBM tumour 0.935, shuffled -0.116;
    LGG 0.630, 0.648); no reading changed. Superseded output kept.
  - **Full pytest 13:50-13:59: 444 passed, 1 skipped**, plus the 24 GIMiCC/atlas tests added during
    the run.

### 2.27 Direct cell counts: T against B in an independent single-cell atlas (2026-10-03 13:21-13:27)
- **Rule** `prespecified/atlas_t_vs_b.md`, born 13:14:27.
  - Written after seeing only the C3 (T+NK) and C11 (B) cluster sizes per patient (C3 larger in
    18/18), before any T/NK split.
  - The split is the one already declared on 2026-10-01 in `abdelfattah_cluster_mapping.md`: the sign
    of T panel minus NK panel, log1p-CPM per cell.
  - Ties are counted as unresolved; the comparison uses strict T.
- **Code** `scripts/atlas_t_vs_b.py`; it reuses the reference builder's streaming and identity check.
  Output `results/atlas_t_vs_b.json`.
- **Identity check PASSED:** integer counts (max relative deviation 2.2e-15); recovered totals against
  library size, max relative error 7.9e-16; 201,986 cells.
- **Result: SUPPORTS T > B.**
  - Pooled: strict T 17,311; NK 1,302; ties 957; B 1,205; C12 unmapped 285. Within-lymphoid shares
    are 0.874 / 0.066 / 0.061.
  - **Strict T > B in 18 of 18 patients.** The smallest T/B ratio is 2.0 (ndGBM-06: 4 against 2);
    LGG-03 has 14 against 0.
  - Mean of per-patient shares: T 0.733, NK 0.213, B 0.054. That is **T > NK > B, EpiDISH's full
    ordering**; NK > B in 16/18 patients.
- **Caveats:**
  - A different cohort (18 patients: ndGBM, rGBM, LGG), so this is not a per-sample TCGA truth.
  - Dissociation biases apply, though T and B are both lymphocytes and both CD45+.
  - The NK call depends on the marker sign: cytotoxic T cells with CD3 dropout are called NK, so NK
    is if anything overstated.
- **Literature, read:** Klemm et al. 2020 [53] (flow cytometry of brain tumours) states "the lymphocyte
  compartment was mostly composed of T cells with fewer NK cells and B cells".
  - Its per-patient Table S2 is not open access through the API; it is an external item for the user.
  - References [52] GIMiCC and [53] Klemm were added and Crossref-verified; reference tests pass.

### 2.28 Flow cytometry by IDH status: Klemm et al. 2020 inspected; LGG settled on direct measurement (2026-10-06)
- **The user placed the paper and supplements** in `celldecov_reference_papers/Klemm_et_al_2020/`
  (gitignored; files dated 16:18-16:19): `PIIS0092867420305699.pdf` (paper), `mmc1.pdf`, `mmc2.pdf`.
  - **`mmc2` = Table S2** is the flow-cytometry gating definitions (T = CD3+ in four subsets; B =
    CD19/CD20+ CD3-; NK = CD3- CD56+). It has **no per-patient values**, unlike what I had expected
    on 2026-10-03.
  - **`mmc1` = Table S1** is the cohort: 17 IDH-mutant gliomas, 40 IDH-wildtype, 37 BrM, 6 non-tumour.
  - The values are in **Figure 1F** (cohort mean of each population, % of CD45+).
- **Rule** `prespecified/klemm_t_vs_b.md`, born **16:26:00** from `stat`, after reading only the
  tables and legends, before any panel was viewed.
- **Method:**
  - `pdftoppm` is absent, and nothing was installed. Page 3 was rendered with macOS `sips` only to
    locate the panel.
  - The values come from the **vector geometry**: `scripts/klemm_figure1f.py` interprets the page's
    content stream, keys colours to the legend swatches, and calibrates on the panel's gridlines
    (31.3125 pt per 25%).
  - Paths come from `config.KLEMM_2020_PDF`.
- **Controls PASSED:** melanoma-BrM CD8+ 33.011 against the printed 33.01; all-BrM lymphocytes
  46.231 against the printed 46.23 (weighted 13/16/8); every bar sums to 100%.
  - A wrong calibration fails these checks (`tests/test_klemm_t_vs_b.py`, 7 pass).
- **Result (`results/klemm_t_vs_b.json`): SUPPORTS T > B.**
  - **IDH-mutant (n 17):** T 7.28%, B 0.35%, NK 1.09% of CD45+; **21 : 1, T > NK > B**. Without DNT,
    T is 4.87%.
  - IDH-wildtype (n 40): 21.02 against 2.25 (9 : 1).
  - Non-tumour (n 6): 5.33 against 0.39.
- **Consequence, as pre-declared:**
  - GIMiCC's LGG B > T is contradicted by direct measurement, and is reported as an instrument failure
    in IDH-mutant tissue. GIMiCC's registered reading stays INCONCLUSIVE.
  - The lymphoid headline holds in both cohorts on direct measurement.
  - Updated: the manuscript (§4.4 paragraph, Reading, top caution, Appendix D), WHY_B_OVER_T §7n.4
    and §7e, OPEN_DEFECTS D23, SUBMISSION_CHECKLIST.
- **Verification:** full pytest **455 passed, 1 skipped** (7 min 31 s); both doc checkers CLEAN;
  inventory 32 datasets, 0 problems. Not committed (the user did not ask this time).
- **Noticed, not acted on:** 5 paper PDFs in `celldecov_reference_papers/` were committed before
  that folder was gitignored (publisher PDFs in a public repository). That is for the user to decide.

### 2.29 Manuscript blueprint: the main results, and why each analysis was done (2026-10-06)
- **The user asked:** commit and push; then "start working on a manuscript ... figure out the main
  results and find a way to devise why you did what you did". Committed `753256f` (Klemm), with `main`
  kept in step.
- **`docs/MANUSCRIPT_BLUEPRINT.md`:** the paper in four sentences; **four main results** (R1 detects,
  R2 does not rank, R3 lymphoid failure against a truth confirmed four ways, R4 identifiability),
  with everything else mapped to supporting or supplementary.
  - **The decision chain** has 18 rows: question, then step, then why this way, then the
    timestamped rule file. Methods should be written from it.
  - It also proposes the structure and figures (two new: a design schematic, and the truth by four
    instruments) and lists the user's decisions.
- **Errors caught in my own draft before handing it over:**
  - I credited GIMiCC's D2 rho 1.00 to stability (stability is 0.771).
  - "Every real method beats its null" (quanTIseq is not evaluable).
  - "Returns no lymphocytes at all" (it is 55/56 and 443/510).
  - A loose mechanism list.
  - "The field moved" for methods that are not new.
  - "Confirmed by DNA methylation" (GIMiCC disagrees in LGG).
  - The protocol's date (committed 09-03, not 09-01).
  - Both doc checkers are CLEAN on the corrected file.
- **`docs/PROJECT_ACCOUNT.md` is now out of date** (it predates raw/X, the extension, E2 and the truth
  tests); the blueprint says so at its top.

## 3 · Corrections made this session (do not repeat)
- **The §4.10 extension table printed B-above-T shares without denominators**, against Appendix D's
  own rule. On the raw/X arm most methods return no lymphocytes for most samples: FARDEEP's "100%"
  rests on 7 GBM samples, MIXTURE's on 5. Fixed 2026-10-02 22:40: every cell shows n, and a footnote
  says the shares are not comparable across methods. ReCIDE's frozen cells now read n/a (it needs
  cells), not "pending".
- **Clock times written into today's pre-registration addenda were estimates, and several were wrong.**
  "13:20" was given for a file last written at 12:55; "11:05" for one written at 09:32; "10:35", "09:20"
  and "08:58" were all later than the facts. They were corrected (2026-10-02, after 13:04) from
  file-system times. **The order -- every rule written before the computation it governs -- holds on
  that evidence.** The correction itself re-modified those files, so the measured times are recorded
  here (`stat`, measured before the correction, all 2026-10-02):

  | file | born | last modified |
  |---|---|---|
  | `results/identifiability_diagnostics.json` (E2 output) | 08:47:00 | 08:47:00 |
  | `results/identifiability_robustness.json` (first robustness output) | 09:00:36 | 09:14:05 |
  | `scratchpad/papers_modern/white2024_dream.txt` | 09:21:26 | 09:21:26 |
  | `scripts/identifiability_denominators.py` (check 7) | 09:24:44 | 09:30:34 |
  | `prespecified/agreement_selection_test.md` | 09:32:05 | 09:32:05 |
  | `scripts/agreement_selection_test.py` (rewrite) | 11:51:08 | 11:51:08 |
  | `prespecified/identifiability_diagnostics.md` (last write = Addendum 2) | 07:52:43 | 12:55:13 |
  | `results/identifiability/ivygap_unmix_ablation.json` (S2 output) | 13:03:57 | 13:03:57 |

  **Read `date` before writing any time into a rule file.**
- **A Crossref refresh overwrote seven verified citations with HTTP 429 errors (2026-10-03, before the paced refetch that finished 13:49; times from file mtimes).**
  `build_data_inventory.py --verify-dois` re-fetched every DOI in a burst and stored each failure over
  the good record. The cache is untracked, so git could not restore it.
  - Fixed in the script: a failed fetch never replaces a good record; 429 backs off (5/15/45 s);
    requests are paced at 1 s.
  - Then re-fetched. **Never run a bulk external refresh that can destroy what it refreshes.**
- **A sample-key mismatch made two GIMiCC controls match 0 samples, and n = 0 was scored as a FAIL
  (2026-10-03 13:34).** ABSOLUTE IDs carry the vial letter; methylation IDs do not.
  - Fixed by keying with `k4s`; a guard now blocks any control matching fewer than 100 samples.
  - **A control with no data is a defect, never a result.** Every comparison that can return
    n = 0 must assert a minimum n before its verdict is read.
- **`sed -i` destroyed a rule file's birth-time evidence (2026-10-03).** At 13:03:58 a `sed -i` fixed a
  wrongly typed time ("13:10") in `prespecified/gimicc_truth_confirmation.md` Addendum 2. BSD `sed -i`
  writes a new file and renames it over the old one, so the file's birth time became 13:03:58.
  - The original birth time, **12:51:36**, was read by `stat` at 13:03:51, just before the sed:
    `Oct  3 12:51:36 2026  Oct  3 13:03:51 2026`. That output is in this session's transcript.
  - Every other rule file still has its original birth time (checked 13:10).
  - **Never `sed -i` a file under `prespecified/`.** Append with `>>`, or rewrite in place with a
    method that keeps the inode (Python `write_text` truncates in place).
  - The "13:10" itself was typed 6 min after `date` printed 13:03:51 -- the same error as above.
    **Copy the time from the `date` output; never type it from memory.**
- **E2's registered design had two flaws, both found after the run and checked under rules written
  first (Addendum 1):**
  - the GBM pre-declared shift (1) duplicates the forced shift-1 setting;
  - D2 kept the degraded-mode methods that the immune arm excludes.

  Neither changed a reading. **Before registering any "vary setting X" design, check whether the
  pre-declared value coincides with a grid point. Reuse an existing comparable-method list instead of
  writing a new drop list.**
- **I wrote "above every registered method's median"** for `unmix`'s leukocyte rho; it is above every
  method (checked). **I attributed S1's linear-vs-registered NNLS gap to the cell-size conversion**
  without checking; the attribution was removed.
- **The user's email went into a Crossref request header once** (polite-pool `mailto`, one DOI lookup,
  2026-10-02). Never put it in a request; project scripts use a neutral User-Agent (checked: none
  contains it).
- **A paper's first author was taken from a website name in search results** ("Duan 2024"; it is Dai et al.).
  Take author lists from Europe PMC or Crossref metadata, never from a search snippet.
- **A background analysis did per-draw work that belonged outside the loop** (key re-mapping in
  8,000 null draws) and hit the 30-min limit. Align once; time a smoke run first.
- **The lid was closed at 16:42 on 2026-10-02 ('Clamshell Sleep', on battery). Every job paused
  until ~19:03.** caffeinate cannot prevent lid-close sleep. Nothing was lost to the pause itself;
  BayesPrism and ReCIDE resumed. Tell the user: lid open and plugged in while runs are going.
- **A long job launched with the Bash tool's `run_in_background` is killed at the tool's 2-h limit,
  and that limit counts sleep.** The MIXTURE TCGA frozen run lost its Python parent that way at
  ~17:50. GBM was already saved; LGG was lost, the orphaned R child was stopped, and the run was
  relaunched detached at 19:06 (`scratchpad/mixture_lgg_frozen.sh`, `/tmp/mixture.status`). **Anything
  that may exceed ~1 h goes in a nohup lane script with its own caffeinate, never in
  run_in_background.**
- **MIXTURE was launched before its driver's time limit was checked** (1-h default), against
  the standing rule. Caught at 16:01; now unbudgeted. Grep `timeout_for` / `R_METHOD_TIMEOUTS`
  before launching any new method.
- **Queue scripts that set `IVYGAP_FARDEEP_CORES=2` hid a crash on the default path (D26).** When a
  driver reads an env var with a default, test the default.
- **`extension_tcga.run` overwrote the shared per-sample file** (`results/extension/estimates_full*_extension.csv`)
  with only the methods of each call, and the unmix ablation stored override settings under the method's
  name (found 2026-10-02). Summaries (`tcga_*.json`) were merged and are intact. Fixed: merge by method;
  never write while an ablation override is set. Frozen-arm per-sample estimates regenerate after E2;
  raw/X per-sample estimates for FARDEEP/LinDeconSeq/ARIC/RNA-Sieve need an atlas re-run (gap logged).
- **The unmix ablation's control forced shift = 1** -- right for GBM (its pre-declared rule picks 1; control
  reproduced exactly) but wrong for LGG, whose rule picks a different shift: the LGG "control" (B>T 45.2%,
  rho 0.602) did not match the reported row (48.2%, 0.592). Fixed: the control now sets NO override, so the
  per-cohort rule chooses as in the reported row; LGG (and GBM, for a consistent artefact) re-run.
- **ReCIDE's first real run was killed at 06:05 on 2026-10-02 by `run_recide.py`'s 6-h default `--budget`**
  (launched without overriding it, against the user's "no budget" directive); ~9 h lost. Default is now no
  limit. Unbounded re-run queued LAST in lane B (`scratchpad/recide_rerun_tail.sh`, log
  `/tmp/recide_real2.log`, marker RECIDE_REAL2). Before launching any long job: grep its driver for
  `budget|timeout` defaults.
- **The Mac idle-sleeps overnight (Maintenance Sleep), freezing every job** (found 2026-10-02 03:21: CPU
  time grew ~4 min in 40 min of wall clock). `caffeinate -is -w <lane pid>` now holds the machine
  awake while lanes A (44794) and B (82433) run, on AC power; it ends with them; `pkill caffeinate`
  undoes it. Re-start it for any new lane script.
- **The BayesPrism `pipeline` variant does NOT free the atlas** (only `authors` does): it runs through
  `r_bridge.run_r_method` in-process, so the parent holds ~2 GB for the whole R run. Do not schedule
  another atlas-building job while a pipeline run is in its BUILD phase; after it enters R the parent
  is idle and swaps out. Fixing it means changing `ivygap/deconv/r_bridge.py` (registered path, full
  pytest) -- consider before the LGG pipeline run, which will hold a larger atlas for ~1-2 days.
- **2026-10-01 15:15 — I started the full `pytest tests/` while BayesPrism (lane A) and the raw/X
  arm (lane B) ran. Its genuine-package tests start R (`requireNamespace("BayesPrism")` alone
  ~680 MB); swap hit 7.4 / 8 GB and the BayesPrism master paged out to 7 MB. Stopped the suite;
  both lanes recovered. RULE: with two heavy jobs running, run only the Python-only tests
  for what changed; the full suite waits for a free lane.** Also: the R `RSOCKnode.R` process
  with PPID 1 is BayesPrism's own snowfall worker (lsof: connected to the master's port) —
  never kill it as an "orphan".
- REFERENCES.md must keep the exact sentence "Nothing in this list is now unverified"
  (`tests/test_references_ordered.py`).
- FARDEEP and LinDeconSeq raw/X lines were misattributed once (log order) — corrected.
- 3 of 4 extension citations wrong from memory → Crossref + DOI test.
- "Independent reference" for Neftel/Darmanis → they are GBmap constituents.
- RNA-Sieve integer check absolute → relative (float32 spacing).
- LinDeconSeq: deconSeq mangles numeric/hyphenated sample names → order-proven renaming.
- Ribosomal hypothesis was post hoc → tested with controls → rejected; SVR nulls dropped *after* the
  test failed (they guard false positives only), reason recorded in the script.
- BayesPrism probe failed once from machine overload (Edge) → keep ≤ 2 heavy jobs.

## 4 · Files (committed in 971aeb3 on 2026-10-03)
**2026-10-03 (GIMiCC truth check, atlas counts, Linseed).**
- **New:**
  - `R/run_gimicc.R`, `R/run_linseed.R`;
  - `prespecified/{gimicc_truth_confirmation, atlas_t_vs_b, linseed_config}.md`;
  - `scripts/{gimicc_truth, gimicc_secondary, gimicc_diagnostics, extract_gimicc_cpgs, atlas_t_vs_b,
    gbmap_t_vs_b, linseed_run}.py`;
  - `tests/{test_gimicc_truth, test_atlas_t_vs_b, test_linseed_gate}.py`;
  - `results/{gimicc_truth_confirmation, gimicc_secondary, gimicc_diagnostics, atlas_t_vs_b,
    gbmap_t_vs_b, linseed_*}.json`, `results/gimicc/` (gitignored; 300 MB, regenerable);
  - superseded: `results/gimicc_truth_confirmation.superseded_keybug_133400.json`.
- **Modified:**
  - `ivygap/config.py` (`TCGA_HM450`);
  - `scripts/build_data_inventory.py` (D31; consumers; 429-safe Crossref refresh);
  - `scripts/build_manuscript.py` (§4.4 truth block, Linseed, ReCIDE LGG rank, instrument-qualified
    lymphoid claims, top caution, Appendix D);
  - docs: `{METHODS, REFERENCES [52][53], WHY_B_OVER_T §7n, EXTENSION_IDENTIFIABILITY §7.4,
    OPEN_DEFECTS D23, PROJECT_ACCOUNT, METHOD_CANDIDATES, SUBMISSION_CHECKLIST, DATA_INVENTORY,
    MANUSCRIPT(+_APPENDICES)}`.

New: `R/run_{fardeep,lindeconseq,recide,bayesprism_authors}.R`, `ivygap/deconv/extension.py`,
`prespecified/abdelfattah_cluster_mapping.md`, `scripts/{atlas_cell_fractions, b_column_diagnostics,
b_profile_tissue_likeness, bayesprism_tcga, bisque_anchoring, bisque_anchoring_synthetic.R,
build_abdelfattah_reference, extension_agreement, extension_tcga, independent_atlas_test,
lymphoid_pooled_tnk, lymphoid_tracking, nk_reannotation, ribosomal_test, run_recide}.py`,
`tests/{test_extension_panel, test_independent_atlas, test_lymphoid_resolution}.py`, `vendor/ReCIDE/`,
`docs/figures/Figure_bisque_anchoring.*`, `data/reference/abdelfattah_2022/`, `results/extension/`,
many `results/*.json`.
Modified: `CLAUDE.md` (external-inputs section), `docs/{OPEN_DEFECTS (D22–D24), WHY_B_OVER_T (§7),
MANUSCRIPT(+_APPENDICES) via build_manuscript, REFERENCES ([36]–[44]), INTRODUCTION_OUTLINE,
METHOD_LIMITATIONS (§6), SUBMISSION_CHECKLIST, FIGURES, supplementary S5/S6}`, superseded-labels in 7
docs, `ivygap/{data/reference.py, deconv/r_bridge.py, bench/run_benchmark.py}`, several scripts/tests.

## 5 · Next, in order
**(2026-10-03 13:51) Current state, read this first:**
- GIMiCC truth check DONE. Registered reading INCONCLUSIVE; GBM corroborated, LGG contested (§2.26).
- Atlas and GBmap direct counts DONE (§2.27).
- Linseed DONE in all three arms. ReCIDE LGG rank DONE.
- Manuscript rebuilt. Both doc checkers CLEAN; inventory 31 datasets, 0 problems; targeted tests pass.
- Full pytest: see the end of this list.
- **Next:**
  1. ~~Klemm 2020~~ DONE 2026-10-06 (§2.28): LGG supported on direct measurement.
  2. E3 (gene-resampling stability) is still unregistered.
  3. The user decides on any abstract change (the LGG caveat; see MANUSCRIPT §4.4).
- **COMMITTED AND PUSHED 2026-10-03 at the user's request ("commit and push and sync"):** `971aeb3`
  on `anatomy-test-full-database`, `main` fast-forwarded to it, both pushed to origin.
  - 250 files, 8.8 MB, none over 5 MB.
  - GSE182109 raw data (4.2 GB) are ignored by an anchored rule; only `gsm_to_sample.tsv` is tracked.
  - `vendor/` (MIT, 5.2 MB, with provenance) is committed.
  - `results/` stays gitignored by design.

0. (2026-10-02, ~12:40) DONE:
   - REGEN and the FARDEEP re-run: every frozen summary identical to the backup;
   - manuscript (§4.12, beats 2d/3/4, Appendix D), supplementary (S9 added), figure index (11/11) rebuilt;
   - both doc checkers CLEAN; 49 doc and E2 tests plus 25 targeted tests pass.
   STILL OWED:
   - ~~the full `pytest tests/ -q`~~ PASSED 13:38 (417 passed, 1 skipped);
   - BayesPrism LGG authors → pipeline → verdict;
   - ReCIDE: anatomic DONE (§2.23); TCGA raw/X run registered (`prespecified/recide_tcga.md`), GBM then LGG;
   - the user's answer on MIXTURE, then Linseed, InstaPrism, MethylResolver;
   - ~~S2~~ DONE 13:04 (§2.22): INCONCLUSIVE by its rule; anatomy does not favour the right lymphoid answer.
1. ~~Ig-removed write-up~~ DONE. Render Figure 9 (`build_paper_figures.py`) once LGG diagnostics exist;
   add a test pinning the independent-atlas result.
2. Lane B continues automatically (raw/X arm → reruns → ReCIDE → diagnostics). Write each up:
   §4.10 table fills from artefacts; ReCIDE into §4.10; LGG diagnostics into WHY §7 / Figure 9.
3. Lane A: BayesPrism GBM authors → pipeline → LGG. Then the **pre-declared prediction verdict**
   (authors' config lowers B share?) into §4.10 and WHY_B_OVER_T.
4. Then: `extension_agreement.py` (n=16), rebuild manuscript / supplementary / figures / figure index,
   `check_doc_numbers --against-current --strict`, `check_doc_statistics --strict`, full
   `pytest tests/ -q`, update memory and this log.
5. External asks queued (ask one at a time): faster machine for BayesPrism LGG; a cohort with
   scRNA + bulk from the same tumours (Bisque paired mode); Abdelfattah authors' cluster metadata
   (SCP "Meta" file) if it exists; Ajaib 2023 per-sample IMC data.
