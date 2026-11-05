# Study guide: understand your own project in one sitting

> **THIS IS A STUDY AID, NOT REPORT TEXT.** It was written by Claude (AI) on 2026-11-02 for you to learn from.
> Do not paste any of it into your report, abstract or application (STS rule: the report must be your own
> writing). Read it, close it, then explain each section out loud in your own words. Log its existence in your
> AI prompt log.

Deadline: **Thu 2026-11-05, 8 PM ET.** Email STS about the AI-use rules first if you have not (see the end).

---

## 1. The whole project in five statements (check you can say each one without looking)

1. **The problem.** Bulk RNA sequencing mixes every cell in a tissue together. *Deconvolution* methods try to
   unmix it and estimate how much of each cell type is there. Nobody knows the true answer for real tumors, so
   nobody can easily check which method is accurate.
2. **The common shortcut.** Researchers often say "the result matches known biology, so the method works." That
   check can't fail, because you decide what counts as "matches" *after* seeing the output.
3. **Your idea (the Anatomy Test).** Use the tumor's own anatomy as a test. Write down, *before running anything*,
   which cell types must be higher in which zone of a glioblastoma. Then score each method on how many of those
   predictions its estimates satisfy. The score is the **Anatomic Concordance Score (ACS)**.
4. **The real question.** Can ACS stand in for ground truth when choosing a method? And does a method that gets
   the anatomy right also get the biology right?
5. **The answer you found.** ACS is a good *detector* of broken methods, but it is *not* a good *ranker* of
   which working method is most accurate. And getting the anatomy right does not mean getting the immune biology right.

---

## 2. Vocabulary: know these cold

| term | plain meaning |
|---|---|
| **Bulk RNA-seq** | RNA from a chunk of tissue, so a blend of all its cells |
| **Deconvolution** | statistics that estimate the cell-type mix from that blend |
| **Reference / signature matrix** | a table of what each pure cell type's RNA looks like, built from single-cell data. Here: **GBmap** |
| **Ivy GAP** | Allen Institute atlas: glioblastoma samples laser-cut from five anatomic zones and sequenced |
| **LE / IT / CT / MVP / PAN** | leading edge (mostly normal brain), infiltrating tumor, cellular tumor (dense core), microvascular proliferation (new blood vessels), pseudopalisading cells around necrosis (cells around dead, oxygen-starved tissue) |
| **Constraint** | one pre-written prediction, e.g. "endothelial cells are higher in MVP than in CT". Seven in total (C1–C7) |
| **ACS** | the (weighted) share of constraint–tumor pairs a method satisfies. Ranges 0–1 |
| **Permutation null** | shuffle zone labels *within each tumor* 10,000 times and rescore. Shows what a score looks like by luck. It comes out near **0.37**, not 0.5 |
| **Negative controls** | two fake methods (random numbers; a scrambled reference). If they scored well, the test would be meaningless |
| **Pre-registration** | the constraint file was written and hashed before any output was inspected, then registered publicly (OSF `dm2t8`, 2026-09-10). Pilot runs came before that date, so the registration is *retrospective*, and it says so. Say this plainly |
| **Ground truth / yardstick** | an independent measurement of the true composition (e.g. DNA-based tumor purity) |
| **Spearman rho** | rank correlation: do two rankings agree? 1 = identical, 0 = unrelated |
| **Underpowered / INCONCLUSIVE** | too little data to tell "no effect" from "an effect we can't see". An honest answer, not a failure |

---

## 3. How the experiment works (two rankings, then compare)

```
Ivy GAP tissue (122 samples, 10 tumors, zone labels by a pathologist)
        |
   15 deconvolution methods + 2 fake controls
        |
   Ranking A: ACS (do estimates follow the zone predictions?)
        |
        +--------- compare (Spearman rho) ---------+
        |                                          |
   Ranking B: accuracy against an independent truth (e.g. DNA purity)
```

Pre-set bar: rho >= 0.60 with a confidence interval excluding zero would mean "anatomy can choose methods."
Every possible outcome was declared publishable in advance.

Key design rules (each answers a likely judge question):
- **Only orderings are scored, never values.** Zones are still mixtures, so a target like "40% tumor" would be indefensible.
- **Comparisons are inside one tumor**, so different patients are never compared. Samples from one tumor are averaged per zone first, because blocks from the same tumor aren't independent.
- **148 other Ivy samples are never scored**: their zone labels came from expression, which would be circular. Code aborts if you try.
- **No method is chosen by outcome or by ACS**, because ACS is what's being tested.
- **T cells were excluded from the constraints before any data**, since they sit at the detection floor. So ACS *cannot see* lymphoid errors.

---

## 4. The three results, with the numbers and the catch

### Result 1: Anatomy detects broken methods
- Registered panel: real methods score ACS 0.708 (Bisque) to 0.969 (MuSiC). FARDEEP's 0.985 comes from the later extension panel. Controls score 0.400 and 0.369 and do not beat their own nulls. Margin about **+0.31**.
- Every real method beats its null at p < 1e-4.
- *Catch:* the controls only prove the constraint set isn't trivially easy. It says nothing about which working method is best.

### Result 2: Anatomy does not rank methods by accuracy
- Registered test against the synthetic-mixture truth: rho = **0.6372** (CI 0.104–0.943), meets the bar. **But** both rankings use the same atlas (GBmap), so part of that agreement is agreement about the reference, not the tissue.
- Against **DNA-measured purity** (ABSOLUTE, independent of RNA): rho = **0.081**, n = 12 methods, p = 0.80. Bar not met.
- A graded version (AUC) gave -0.06. A protein-imaging test (IMC) saw the ranking *reverse*.
- It is **underpowered**: with 12 methods a true rho of 0.7 would be detected only 65% of the time. So say "not shown to rank", never "proven it can't."
- Also: ACS can't resolve its own top. On the current leaderboard the top confidence intervals overlap (MuSiC 0.930–1.000, SVR 0.891–1.000). The README's "four-way tie at 0.9846" comes from an older run, so don't quote it.

> **The trap to avoid:** you have two numbers, 0.6372 and 0.081. They are *different tests* (synthetic-mixture truth vs DNA truth).
> Be able to say which is which, and why 0.6372 is weaker evidence (shared atlas). A judge will check.

### Result 3: Right anatomy does not mean right biology
- Pre-registered prediction (2026-09-17, about 94 minutes before the LGG methylation measurement was produced): DNA methylation will show T cells > B cells.
- Methylation puts T above B in **94.8%** of 155 GBM samples and **93.2%** of 530 LGG samples.
- **0 of 12 methods** reproduce T > B on the registered signature, in either cohort.
- *Catch 1:* the top-ACS method (MuSiC) doesn't invert it, it returns zero T/NK/B in most samples. Don't claim "the best method inverts it" (the repo corrected this on 2026-09-27). SVR and SCDC do invert it.
- *Catch 2:* GBM was the discovery cohort; LGG was the pre-registered replication. Present them that way, not as two equal confirmations.
- *Catch 3:* on a different reference build a few methods (3/14 GBM, 2/13 LGG) do get T > B, so "no method **robustly**" is the careful wording.
- *Catch 4:* EpiDISH (the methylation tool) uses a *blood* reference on brain tumors and is validated only at cohort level. The worklog (§2.26–2.28) says it was cross-checked with other methylation tools, direct single-cell counts and flow cytometry (Klemm 2020). I did not read those sections in detail: read them yourself.

### Extra (exploratory, not pre-registered; mention briefly, if at all)
Why is B above T? Nine explanations were tested and rejected; an independent atlas gave a mixed answer. Also: a "stability" check (re-fitting under different settings) tracked which compartments are determined by the data (rho 0.886, only 6 units). Label these as exploratory.

---

## 5. Numbers cheat sheet (all read from the repo's results docs)

| item | value |
|---|---|
| Ivy GAP anatomic cohort | 122 samples, 10 tumors, 9 usable (57 constraint–tumor pairs) |
| Structures | CT 30, MVP 25, PAN 24, IT 24, LE 19 |
| Methods | 15 in the registered panel (+ extension methods added after registration, labelled as such) |
| Constraint-file hash | `2d1fb47c…807a`, registered on OSF 2026-09-10 |
| Chance level of ACS | about 0.37 (not 0.50) |
| Permutations | 10,000 within-tumor |
| Control scores | 0.4000 / 0.3692 |
| Registered rho (synthetic truth) | 0.6372, p 0.0143, 14 methods (13 distinct) |
| ABSOLUTE-purity rho | 0.081, p 0.80, n = 12, CI about -0.64 to 0.70 |
| Methylation T>B | 94.8% (GBM, 155) and 93.2% (LGG, 530) |
| Methods reproducing T>B | 0 of 12, both cohorts |
| Verification | 75 of 75 re-run steps; results frozen (677 files, 0 changed) |

Caution: `quanTIseq` was scored on only 15 of 57 pairs, so it is INCONCLUSIVE, not "bad." Bisque's good T-vs-B showing is just its reference's composition.

---

## 6. The instructor's four questions: what the repo says, and what is still yours to answer

**Q1. What is "tumor anatomy" and how is it computed?**
Not computed from the data. It's the pathologist-drawn zone labels (Ivy GAP) plus histology-based expectations. Your job: explain why each of C1–C7 must hold (see `docs/WHAT_ACS_IS.md`, section "Seven predictions").

**Q2. How were the methods chosen, how do they work, do they cover the field?**
- Registered panel: published packages (MuSiC, SCDC, Bisque, EPIC, quanTIseq, DWLS, BayesPrism, CIBERSORTx) plus the project's own NNLS, SVR, elastic net and Bayesian versions. Later additions (FARDEEP, ARIC, LinDeconSeq, RNA-Sieve, ReCIDE, DESeq2 unmix, MIXTURE, Linseed) came from the top tiers of two published benchmarks (Nguyen 2024, Li 2026).
- Gaps to admit: no deep-learning methods (Scaden and DAISM-DNN need TensorFlow, which doesn't build here, and train on simulated mixtures); some methods fell back to the project's own reimplementations, or were later re-run as genuine packages. Several methods in the panel are this project's own implementations, not the published ones.
- Plain-language "how they work" (general knowledge; verify with each paper): NNLS/SVR/CIBERSORT fit non-negative mixtures of reference profiles to the bulk; DWLS reweights genes; MuSiC and SCDC use variation across donors; Bisque adapts the reference to the bulk; BayesPrism uses a Bayesian model. Write your own 1–2 sentences per method after reading the abstracts.
- Whether the panel covers the field is *your judgment call*. Decide, and say why.

**Q3. Are the data human-annotated? Is there circularity?**
- Ivy zone labels: pathologist, from tissue, no RNA involved. Constraints: from histology and in-situ hybridization, fixed before results. Good.
- Samples whose labels came from expression are excluded by code.
- Synthetic-mixture truth is built from the same atlas the methods use: **circular**, which the project says and why the DNA/methylation truths matter.
- DNA purity (ABSOLUTE) and methylation (EpiDISH) are computed, not hand-counted, but independent of RNA.
- GBmap's cell labels come from the atlas authors' annotation (`annotation_level_3`). I did not verify how those were made. **Look it up and state it.**

**Q4. Why these evaluation methods?**
- Within-tumor permutation null: because cell fractions sum to 1 and samples are nested, so a coin-flip null is wrong.
- Spearman: ranks are what a method chooser needs. Bootstrap over tumors for intervals.
- Pre-registration, negative controls, power analysis: so every outcome is reportable.
- You still need a literature precedent for rank-correlation of benchmark rankings. I have not verified one: find one.

---

## 7. Hard questions a judge may ask

| question | where the facts are |
|---|---|
| "Isn't this circular?" | section 6, Q3; the 148 excluded samples; pre-registration |
| "Why isn't your headline result positive?" | `docs/RESEARCH_QUESTION.md`: "anatomy detects, it does not rank" is a result |
| "Only 9 tumors?" | acknowledged; permutation design uses repeated measures; ACS per tumor is least reliable where few zones were sampled |
| "You scored your own implementations?" | `docs/METHOD_LIMITATIONS.md`: reimplementations are labelled, never reported as the published package |
| "What did you do yourself vs the AI?" | your prompt log and the STS disclosure. Answer truthfully |
| "What would a harder test need?" | more tumors; constraints for T/B cells; a per-sample cell-count truth |

---

## 8. Say the weaknesses first (reviewers trust that)

- Underpowered ranking test (n = 12 methods).
- Registered 0.6372 shares an atlas between both rankings.
- ACS can't see lymphoid cells at all.
- Single tumor type (glioma); small cohort (9 tumors).
- Some "methods" are the project's own reimplementations.

---

## 9. Where the evidence for 20 pages lives (planning aid only)

| report part | go read |
|---|---|
| Background and question | `README.md` sections 1–2, `Anatomy_Test.md` |
| Methods and data | `docs/WHAT_ACS_IS.md`, `docs/DATA_SOURCES.md`, `ivygap/anatomic/constraints.py` |
| Result 1 and 2 | `RESULTS.md` sections 0–3, `docs/RESEARCH_QUESTION.md`, `docs/WHAT_ACS_IS.md` |
| Result 3 | `docs/RESEARCH_QUESTION.md`, clause 2 |
| Limits | `docs/METHOD_LIMITATIONS.md` |
| Figures and captions | `docs/FIGURES.md`, `docs/figures/` |

---

## 10. A realistic schedule (sleep included)

- **Today:** email STS (if not sent). Send the CMU instructor your reply. Then work through sections 1, 2, 3 above until you can say them aloud. Talk to me whenever a term or number doesn't click. Sleep.
- **Tue:** sections 4 and 6. For each result, write 3 sentences from memory, then check them against section 5. Write your outline. Sleep.
- **Wed:** write your own draft, rough first. Then have me fact-check numbers. Disclose AI use. Sleep.
- **Thu:** application forms and uploads before 8 PM ET.

Rest of the application (essays, activities, recommendation, the AI disclosure) also needs your own words. Budget time for it.

---

## 11. What this guide did and did not read

Read in full: `STS_WORKLOG.md` (the early and late parts), `WHAT_ACS_IS.md`, `RESEARCH_QUESTION.md`, `Anatomy_Test.md`, `TWO_PROBLEMS.md`, `METHOD_LIMITATIONS.md`, the README start, the first third of `RESULTS.md`, parts of `DATA_SOURCES.md`, the close-out documents. **Not read:** worklog sections 2.26–2.37, `METHODS.md`, `MANUSCRIPT*.md`, `WHY_B_OVER_T.md`, `EXTENSION_IDENTIFIABILITY.md`, `OPEN_DEFECTS.md`, the code. Treat anything outside the sections above as unverified here.
