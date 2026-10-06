# Manuscript blueprint: the main results, and why each analysis was done

Written 2026-10-06, against commit `753256f`. A plan for writing the paper, not paper prose.
- The generated draft, with every number read from an artefact, is `docs/MANUSCRIPT.md`
  (`scripts/build_manuscript.py`). This file decides what that draft should argue, in what order, and
  why each step was taken.
- Numbers here are quoted from the current artefacts. If an artefact changes, the draft updates
  itself and this file must be re-checked.

**Out of date: `docs/PROJECT_ACCOUNT.md` and the older outlines.** They predate:
- the raw/X reference arm (D23);
- the extension panel;
- Extension E2;
- the three tests of the lymphoid truth (GIMiCC, direct cell counts, flow cytometry).

Where they disagree with this file, this file is current.

---

## 1 · The paper in four sentences

1. Choosing a deconvolution method needs ground truth, which is the measurement that settings using
   deconvolution cannot afford. So in practice estimates are trusted when they "match known biology".
   That check has never been tested.
2. We wrote down, before scoring anything, what a glioblastoma's anatomy guarantees about its cell
   composition (seven hashed, publicly registered constraints), and scored 15 estimators and 2
   deliberately broken controls against them.
3. Anatomy separates working methods from broken ones. It does not identify the accurate one, and
   the method it ranks first returns no lymphocytes in 55 of 56 GBM and 443 of 510 LGG samples.
4. No method robustly recovers the T-over-B lymphoid ordering. That ordering is confirmed by direct
   measurement (two single-cell atlases, flow cytometry) and by DNA methylation, although in LGG one
   of two methylation instruments disagrees. Matching known biology is not evidence that a
   composition estimate is correct.

**One-line thesis:** *A biological-plausibility check catches broken deconvolution, cannot choose
correct deconvolution, and passes methods that are reproducibly wrong about the immune compartment.*

---

## 2 · The main results: the paper's spine (four results, in this order)

The draft has twelve results subsections (§4.1-4.12). A paper carries four. Everything else
becomes supporting evidence inside one of them, or supplementary.

### R1 · Anatomy detects: it separates working methods from broken ones *(registered)*
- **Claim:** every real method satisfies the anatomic constraints far more often than chance, and two
  deliberately broken estimators do not.
- **Evidence:**
  - ACS of the best method (`music`) is 0.9692; of the lowest fully evaluable method (`bisque`)
    0.7077.
  - The random-proportion control scores 0.4000 and the shuffled-signature control 0.3692.
  - **Margin +0.31** on a 0-1 scale. Every fully evaluable real method beats its within-tumour
    permutation null at p < 1e-4; both controls fail theirs.
- **Figure:** `Figure_detects_not_ranks` (left panel), or an ACS leaderboard with the controls marked.
- **Caveat to state:** quanTIseq is *not evaluable* (it scores only 15 of 57 pairs, where chance
  reaches 0.50). That is different from scoring badly.

### R2 · Anatomy does not rank: it does not identify the accurate method *(registered; INCONCLUSIVE)*
- **Claim:** the ACS ordering of methods does not predict their accuracy against an instrument that
  shares nothing with the anatomic arm.
- **Evidence:**
  - Against DNA-measured tumour purity (ABSOLUTE): Spearman **rho = +0.081** (p 0.80, CI −0.64 to
    +0.70, 12 methods).
  - The registered bar (rho >= 0.60) is met only against the synthetic yardstick (+0.637, p 0.014),
    which is built from **the same atlas** the methods deconvolve against.
  - Every yardstick sharing nothing with ACS fails. That is a pattern, not one bad run.
- **Say "inconclusive" and mean it:** with 12 methods, a true rho of 0.70 meets the bar in only 65% of
  studies; 80% power needs 30-50 methods (`scripts/agreement_power.py`).
- **Figure:** `Figure_detects_not_ranks` (right panel); purity scatters as supplementary.

### R3 · Getting the anatomy right does not mean getting the biology right *(registered prediction; robustness exploratory)*
This is the finding the paper is built on.
- **Claim:** no method robustly recovers the lymphoid compartment of a glioma. The truth they miss is
  confirmed by every independent instrument that measures it directly.
- **The prediction:** "DNA methylation will show T cells > B cells", with its falsifier, was committed
  2026-09-17 23:37:27. The LGG methylation measurement was produced 95 minutes later.
- **The methods:**
  - **0 of 12 order T above B on the registered signature, in GBM and in LGG.**
  - On the donor-level raw/X rebuild, 3 of 14 and 2 of 13 do. The one consistent success (Bisque)
    returns its reference's own composition by its declared assumption, shown on planted truth.
  - Across 9 more genuine packages (the extension panel), at most 1 per cohort and reference does,
    never robustly.
  - Two failure modes: **absence** (MuSiC and NNLS return exactly zero lymphocytes in 55 of 56 GBM
    samples) and **misassignment** (the rest put B above T).
  - The method ranked first by anatomy is an absence case.
- **The truth, tested four ways:**

  | instrument | what it is | GBM / IDH-wildtype | LGG / IDH-mutant |
  |---|---|---|---|
  | EpiDISH (methylation, blood reference) | the registered truth | T > NK > B | T > NK > B |
  | GIMiCC (methylation, glioma-specific) | second instrument | T > B, every control passing | **B > T** (registered reading INCONCLUSIVE) |
  | single-cell atlases (direct counts) | Abdelfattah 2022; GBmap | T > B in 16/16 patients; 96/98 donors | T > B in 2/2 patients |
  | flow cytometry (Klemm 2020, Figure 1F) | direct measurement | T 9x B (40 tumours) | **T 21x B (17 tumours)** |

  - GBM is corroborated by every instrument.
  - In LGG one methylation instrument (GIMiCC) disagrees, and direct measurement settles the question
    against it.
  - Report that disagreement as a finding in its own right: methylation-derived lymphoid "truths" are
    instrument-dependent, and both place B above every direct count.
- **Figures:** `Figure_lymphoid_failure` and `Figure_bisque_anchoring`. **New, recommended:** "The
  lymphoid truth by four instruments", within-lymphoid T / NK / B shares from each instrument beside
  the methods' estimates. It does not exist yet.

### R4 · Why: the bulk data identify the immune compartment, not its lymphocyte types *(exploratory)*
- **Claim:** the lymphoid split is not determined by bulk RNA under this reference, so each method's
  modelling choices decide it. A truth-free check can flag this.
- **Evidence:**
  - **Stability:** re-fitting one genuine package (DESeq2 `unmix`) under seven loss scales leaves
    tumour and total-leukocyte content in place and moves the lymphoid split. Stability ranks
    compartments the way agreement with DNA does: rho 0.886 (p 0.017), and 0.943 on matched
    denominators.
  - With GIMiCC as the truth, stability gives 0.771 (p 0.051) and agreement between methods 1.00
    (p 0.0014).
  - **Mechanisms:** nine candidate mechanisms were each tested and rejected (WHY_B_OVER_T §3, §7g,
    §7h):
    - the signature cannot separate T from B;
    - macrophage spillover;
    - tumour purity;
    - per-sample model fit;
    - misfit propagation;
    - CD3+ cells in the NK column;
    - immunoglobulin, ribosomal and myelin content of the B profile.
  - Two more are partial: B absorbing tumour signal is estimator-specific (true for SVR, not NNLS),
    and an independent atlas gives a mixed result. No single cause.
- **Figure:** `Figure_identifiability`; `Figure_lymphoid_mechanisms` as supplementary.
- **Caveat to state:** stability says what *not* to trust. It certifies nothing.

### Supporting and supplementary (moved out of the main line)

| result | where it goes |
|---|---|
| Tumour-content recovery (median 33% GBM / 23% LGG of true variation); purity is the one failure factor that replicates (12/12) | inside R2, short; detail supplementary |
| Mixing-model fit (64% / 76% of marker-space variance unexplained) | supplementary; it explains *why* R2-R3 can happen |
| Reference parity; genuine package against reimplementation | Methods, plus a supplementary table |
| Extension panel (9 genuine packages, post-registration) | one paragraph in R3; full table supplementary |
| Reference-free arm (CDSeq, Linseed: they meet the constraints without a single-cell reference in the deconvolution) | supporting R1 (detection does not depend on the atlas) |
| ACS vs a per-method AUC; the IMC protein test; constraint validation against Albiach, Darmanis and ISH | supplementary |
| DECEPTICON's agreement rule tested against DNA (weak, rho 0.28) | Discussion: agreement between methods is not accuracy either |
| Spatial mesenchymal programme | drop, or supplementary |

---

## 3 · Why we did what we did: the decision chain

Each step answers a question the previous step raised. The record column is the evidence that the
reasoning came first: a rule file written, and timestamped, before the analysis it governs.

| # | the question that forced the step | what was done | why this way, not another | record |
|---|---|---|---|---|
| 1 | Deconvolution is validated against the measurement it replaces. Can it be checked without that? | Use the tumour's anatomy as the test, not the result | The "matches known biology" check is assumed everywhere and never measured. Turning anatomy from an *output* into an *input* is what makes it testable. | `Anatomy_Test.md` (protocol), committed 2026-09-03 with the constraint file |
| 2 | Where does anatomy guarantee composition? | Glioblastoma, Ivy GAP | The only public RNA-seq whose regions were labelled by a pathologist on H&E. GBM's regions have near-definitional compositions (microvascular proliferation *is* proliferating endothelium). | protocol |
| 3 | How can histology be stated without overclaiming? | Seven **ordinal, within-tumour** constraints, hashed and registered before scoring | Histology asserts direction, not size. Pairing within a tumour removes patient effects. Registering first prevents tuning. | `constraints.py` committed 2026-09-03; OSF `dm2t8` 2026-09-10 (hash `2d1fb47c...`) |
| 4 | Which samples can be scored? | Only the 122 histology-labelled samples, not the 148 ISH-cluster samples | The ISH labels were assigned using expression, so scoring them would be circular | `is_anatomic_sample` |
| 5 | What is chance? | A within-tumour permutation null; samples collapsed to (tumour, structure) | Cell fractions are compositional and correlated, so 50% is the wrong null. Samples are nested, so naive counting inflates n. | protocol; METHODS |
| 6 | Can the score fail? | Two deliberately broken estimators scored alongside | A score that rates broken and real inputs alike measures nothing, and real inputs alone cannot show it | registration |
| 7 | Should T cells be constrained? | **Excluded in advance** | The pipeline's own benchmark put T cells at the detection floor. R3 later showed they are systematically wrong. | constraint file |
| 8 | Is a method that gets the anatomy right also accurate? | Truths that share nothing with the anatomic arm: DNA copy number (ABSOLUTE) and DNA methylation | The synthetic yardstick shares the atlas, so agreement with it can be inherited, not earned (R2) | registration |
| 9 | How should the ranking result be read? | Registered bar rho >= 0.60; power simulated afterwards | A null at n = 12 is not a refutation. The simulation turns "underpowered" into a number. | `scripts/agreement_power.py` |
| 10 | GBM estimates put B above T. Real, or a method artefact? | A prediction with a falsifier, registered before the LGG methylation data existed | Lets the anomaly fail. A post hoc story could not. | `immune_failure_factors.md`, 2026-09-17 23:37 |
| 11 | Was the reference itself right? | Rebuilt from genuine counts (D16); both references reported (D23) | The atlas layer was log-transformed data treated as linear. Reporting one reference hid the arm that disagreed. | OPEN_DEFECTS D16, D23 |
| 12 | Is it the implementation, not the method? | Genuine published packages replace reimplementations; failures are labelled | A reimplementation can fail where the package does not | D18, D19, D22; `remeasurement_verdict.md` |
| 13 | Is the failure specific to this panel? | Nine more genuine packages, post-registration and labelled as such | Methods that recent benchmarks rank highly were missing from the registered panel (DESeq2, MIXTURE, Linseed, FARDEEP), one was published during the study (ReCIDE), and BayesPrism had not been run in its authors' configuration. A failure that survives them is general. | `deseq2_unmix_config.md`, `recide_tcga.md`, `mixture_config.md`, `linseed_config.md`, `bayesprism_authors_prediction.md` |
| 14 | *Why* do methods put B above T? | Nine mechanism tests, each with a rule and a control | Localise the failure: the algebra, the reference, or the gap between reference and real tissue | WHY_B_OVER_T §3-§7k |
| 15 | No single cause. Is the split even determined by the data? | Truth-free stability across loss scales, tested against DNA agreement (E2) | If the split is not identified, any method's answer is its modelling choice, and practitioners need a check that needs no truth | `identifiability_diagnostics.md`, 2026-10-02 07:52 |
| 16 | Does agreement *between methods* pick the accurate one? | DECEPTICON's rule tested against DNA | It is the field's other truth-free criterion | `agreement_selection_test.md`, 2026-10-02 09:32 |
| 17 | The headline rests on EpiDISH with a *blood* reference. Is the truth right? | Three independent tests: GIMiCC (glioma-specific methylation), direct single-cell counts (two atlases), flow cytometry | A reviewer's first objection. Direct measurement shares no assumption with methylation deconvolution. | `gimicc_truth_confirmation.md` (+6 addenda), `atlas_t_vs_b.md`, 2026-10-03; `klemm_t_vs_b.md`, 2026-10-06 |
| 18 | How do we know the numbers are right? | Every check ships with a negative control; independent recomputation; document checkers; timestamps from `stat` | Errors were caught only this way: a sign flip (D21), a sample-key mismatch, Broad-ID joins (D27) | OPEN_DEFECTS; `independent_verification.py`; tests |

**The Methods section should be written from this table.** Each design choice is stated with its
reason. The Introduction's argument is rows 1-3 and 8. The Discussion's limitations are rows 9, 13
and 17.

---

## 4 · Proposed structure

| section | content | figures |
|---|---|---|
| Introduction | The circularity (row 1); the untested plausibility check; the two-clause question | -- |
| Methods | Rows 2-9 and 11-12, each choice with its reason; registration and corrections | Fig 1: the design (anatomy as input; constraints; controls; independent truths) -- **new, schematic** |
| Results R1 | Detection | Fig 2: ACS leaderboard with controls |
| Results R2 | No ranking; the shared-atlas pattern | Fig 3: `Figure_detects_not_ranks` |
| Results R3 | Lymphoid failure; two modes; the truth by four instruments | Fig 4: `Figure_lymphoid_failure` + **new** truth panel; `Figure_bisque_anchoring` |
| Results R4 | Identifiability; mechanisms excluded | Fig 5: `Figure_identifiability` |
| Discussion | What is new; the practical reading; where it sits (reproducibility vs correctness); what would settle the ranking question | -- |
| Limitations | One tumour type; n = 12; one atlas; cohort-level truths; post-registration extensions labelled | -- |
| Supplementary | Everything in the supporting table above; defects log; verification | S-figures and S-tables |

---

## 5 · Decisions that are yours

1. **Title.**
   - The current title: "Anatomic Concordance Detects Broken Bulk-RNA Deconvolution: No Method
     Recovers the Lymphoid Compartment in Two Glioma Cohorts".
   - Since D23, the exact form is "...No Method *Robustly* Recovers...".
2. **The abstract's truth sentence.** It says "the ordering that methylation resolves". One
   methylation instrument now disagrees in LGG. Suggested: "...the T-cell-over-B-cell ordering that
   DNA methylation resolves and that direct flow-cytometric and single-cell counts confirm".
3. **Venue and length.** The plan above assumes a full paper. Tell me the page or word limit, and
   which results become supplementary follows from it.
4. **The two new figures** (the design schematic; the truth by four instruments). I can build both
   from the artefacts.
5. **The five paper PDFs** committed to `celldecov_reference_papers/` before it was gitignored:
   publisher PDFs in a public repository.
