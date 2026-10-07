# Introduction — paragraph outline with sources attached

Write from this. Each paragraph lists the claims it should make and the reference that
supports each one. Where I have read the source, the supporting line is **quoted** so you can
see the claim is really in it and phrase your own sentence around it.

**Three kinds of marking, and the difference matters:**

| mark | meaning |
|---|---|
| **[n] + quote** | I read this in the PDF. The quote is verbatim. Safe to cite for that claim. |
| **[n] — verify** | The paper is the right source for this claim, but I did not read the passage. Check it before citing. |
| **own data** | Comes from this study's artefacts. No external citation; cite your own Results or Methods. |
| **NEEDS A SOURCE** | Neither this study nor any paper in the list supports it. You must find one or drop the claim. |

---

## ¶1 — Why cell composition matters

**Claims to make:**

1. The composition and density of immune cells in the tumour microenvironment influences
   tumour progression and the success of anti-cancer therapy.
   → **[2] Sturm 2019.** Verbatim from the abstract:
   > *"The composition and density of immune cells in the tumor microenvironment (TME)
   > profoundly influence tumor progression and success of anti-cancer therapies."*

   This is the strongest single sentence you have for opening the paper, and it is a direct
   quote from a *Bioinformatics* benchmark — not a review, not a claim of mine.

2. Glioblastoma specifically — its immune compartment, why it matters clinically.
   → **Candidates found 2026-10-01: [42] Abdelfattah 2022, [43] Mathewson 2021, [44] Hara
   2021** (`docs/REFERENCES.md`). Metadata verified against Crossref; **not read here — verify
   the passage before citing.** [43] is the most direct for T cells in glioma specifically;
   [44] for immune cells shaping glioblastoma cell states.

3. Tumour purity matters in glioblastoma and is hard to measure — the rationale for the DNA
   purity yardstick. → **[41] Thomas 2025** (*Neuro-Oncology*; PDF on disk). Verbatim:
   > *"Tumor purity, the proportion of malignant cells within a tumor, is an important
   > covariate for understanding the disease, having direct clinical relevance or obscuring
   > signal of the malignant portion in molecular analyses of bulk samples. However, current
   > methods for estimating tumor purity are nonspecific and technically demanding."*

   Independently useful in Methods: GBMPurity *"was trained using simulated pseudobulk tumors
   of known purity from labeled single-cell data acquired from the GBmap resource"* — a
   published glioblastoma tool built on the same atlas this study uses as its reference.
   (Its hyphenation in the PDF reads "un- derstanding"; quote it as "understanding".)

*Keep this paragraph to two or three sentences. It is the only part a clinical reader needs
before the problem statement.*

---

## ¶2 — Why deconvolution exists: direct measurement is unavailable or unaffordable

This is the paragraph that carries your equity argument, and you have two independent sources
for it.

**Claims to make:**

1. Direct measurement of composition — flow cytometry, immunohistochemistry, single-cell
   sequencing — is frequently not available, so computational estimation is used instead.
   → **[2] Sturm 2019.** Verbatim:
   > *"Flow cytometry, immunohistochemistry staining or single-cell sequencing are often
   > unavailable such that we rely on computational methods to estimate the immune-cell
   > composition from bulk RNA-sequencing (RNA-seq) data."*

2. Single-cell experiments are **expensive**, and deconvolution is the route to cell-type
   resolution from bulk data that already exists.
   → **[3] Nguyen 2024.** Verbatim:
   > *"However, performing single-cell experiments is expensive. Although cellular
   > deconvolution cannot provide the same comprehensive information as single-cell
   > experiments, it can extract cell-type information from bulk RNA data, and therefore it
   > allows researchers to conduct studies at cell-type resolution from existing bulk
   > datasets."*

   **This resolves the costing gap.** You were going to need a per-sample dollar figure to say
   "expensive"; [3] lets you make the claim on a citation instead. If you still want a number,
   it must come from a source you cite — this project holds no costing artefact and you should
   not assert one.

3. The equity consequence: bulk RNA-seq is already routine, so a trustworthy deconvolution
   method makes cell-type analysis available wherever bulk sequencing is.
   → **inference from [2] + [3]**, and it is a fair one — both state the premise. Phrase it as
   *your* inference, not as something either paper claims.

---

## ¶3 — Many methods exist, and they disagree

**Claims to make:**

1. A large number of methods has been developed.
   → **[1] Avila Cobos 2020**, verbatim: *"Many computational methods have been developed to
   infer cell type proportions from bulk transcriptomics data."*
   → **[3] Nguyen 2024** reviews **53 methods** across 28 datasets — cite it for the scale.
   *(The 53/28 figures are from `docs/RELATED_WORK.md`; **verify** them against the paper.)*

2. Choosing between them is genuinely hard for a non-specialist.
   → **[3] Nguyen 2024**, verbatim:
   > *"The large number of methods available, the requirement of coding skills, inadequate
   > documentation, and lack of performance assessment all make it extremely difficult for
   > life scientists to ch[oose]…"*

3. Results depend strongly on choices that are not the method itself — transformation,
   normalisation, marker selection.
   → **[1] Avila Cobos 2020**, verbatim: *"the choice of normalization has a dramatic impact
   on some, but not all methods"*, and methods *"perform best when applied to data in linear
   scale."*

   **This one is worth foregrounding**, because your own OPEN_DEFECTS D16 is exactly this
   failure — the reference was log-transformed data treated as linear, and rebuilding from
   counts moved one method by 0.26. You independently reproduced a published warning. Say so
   in the Discussion, not here.

---

## ¶4 — The circularity: validation needs the measurement it replaces

This is the paragraph the whole paper turns on. Build it carefully.

**Claims to make:**

1. Benchmarks establish method performance using known composition — pseudo-bulk mixtures
   built from single-cell data, or flow cytometry on the same specimen.
   → **[1] Avila Cobos 2020**, verbatim: *"Using five single-cell RNA-sequencing (scRNA-seq)
   datasets, we generate pseudo-bulk mixtures to evaluate…"*
   → **[2] Sturm 2019** — *"We used a single-cell RNA-seq dataset of ~11,000 cells from the
   TME to simu[late]…"*

2. Performance on simulated mixtures does not transfer to real bulk.
   → **[4] Li et al. 2026.** This is the strongest citation in the paper for your central
   objection, and it is independent of you. Verbatim from their introduction:
   > *"current benchmarking studies for deconvolution methods invariably lean on pseudobulk
   > data or flow cytometry […] In real bulk RNA-expression deconvolution, such precision is a
   > mirage."*

   They measured it: on a cohort profiled by **both** scRNA-seq and bulk, pseudobulk gave
   systematically higher correlations than real bulk, and **method rankings transferred poorly
   between the two (p < 0.001)**, across 18 cohorts and 5,891 samples.

3. **The circularity, stated in your own words.** Every established route to choosing a method
   requires the measurement deconvolution is used in place of. A laboratory that cannot obtain
   single-cell data also cannot validate the method it is using instead — and by [2], that is
   the common case, not the edge case.
   → **your synthesis of [1][2][3][4].** No single paper states it this way; it follows from
   all four. This is your contribution to the framing and you should own it.

---

## ¶5 — The proposal: use anatomy the pathologist already labels

**Claims to make:**

1. Glioblastoma has anatomically distinct regions that a neuropathologist identifies on
   histology, and Ivy GAP provides RNA-seq from those regions with the labels attached.
   → **[5] Puchalski 2018**, *Science* 360:660–663 — *"An anatomic transcriptional atlas of
   human glioblastoma."* **Verify** the specific structures (LE, IT, CT, MVP, PAN) are named
   there as you describe them; I have not read this paper.

2. Some compositional facts about those regions are not in dispute and do not depend on any
   algorithm — microvascular proliferation is *defined* by proliferating endothelium; the
   leading edge is infiltrated brain and retains oligodendrocytes.
   → **own data / definitional.** These are the rationales recorded in your constraint file
   (`ivygap/anatomic/constraints.py`, C1–C7). Cite your own Methods. If a reviewer wants a
   neuropathology reference for the definitions themselves, a WHO CNS tumour classification
   citation would serve — **you would need to add it.**

3. Therefore: a method's estimates can be scored against those facts using only a slide and a
   pathologist's label, neither of which requires single-cell data.
   → **your argument.** State it as the study's premise.

4. *(Added 2026-10-07.)* Anatomy has already served as supporting evidence for deconvolution, but only
   qualitatively.
   - BayesPrism's authors deconvolved Ivy GAP and reported endothelium enriched in microvascular
     proliferation → **[19] Chu et al. 2022**.
   - Varn et al. characterised Ivy GAP regions by deconvolution and immunofluorescence → **[57] Varn et
     al. 2022**, *Cell* 185:2184.
   - Neither treats anatomic agreement as a test that could fail, nor checks it against an independent
     truth. That is the gap this study fills.
   - Robustness studies so far use simulated mixtures → **[58] Xu et al. 2025**.
   - Single-cell and single-nucleus proportions carry known dissociation biases → **[59] Slyper et al.
     2020**, the context for the CPTAC finding (§4.13).
   - All four were found by the 2026-10-07 literature check (`docs/CLAIMS_LEDGER.md`) and resolved on
     Crossref. **Read [19] and [57] for the exact wording before quoting them.**

---

## ¶6 — What this study did, and the two-clause question

**Claims to make:**

1. Constraints were fixed, hashed and publicly registered before any method was scored.
   → **own data.** OSF `dm2t8`, 2026-09-10, constraint hash `2d1fb47c…`. Cite the
   registration.

2. Fifteen estimators and two deliberately broken negative controls were ranked.
   → **own data**, Methods.

3. The ranking was then tested against two ground truths that share no input with the
   anatomic arm — DNA copy-number purity and DNA methylation.
   → **own data** for the design; **[7] Carter 2012** for ABSOLUTE and **[12] Teschendorff
   2017** for EpiDISH, both cited in Methods rather than here.

4. **State the registered question verbatim** and do not paraphrase it:
   > *Can a tumor's own anatomy stand in for ground truth when choosing a cell-type
   > deconvolution method — and does a method that gets the anatomy right also get the biology
   > right?*

   It is pre-registered, so its wording is part of the record. Quote it, mark it as the
   registered question, and let the two clauses set up Results 4.1–4.2 and 4.4.

---

## Citation order, if you write it in this sequence

`[2] · [3] · [1] · [3] · [1] · [4] · [5]` — which is close to the numbering already in
`docs/REFERENCES.md`. If your final paragraph order differs, renumber; the list is ordered by
first use and should be regenerated to match.

---

## The three gaps you must close yourself

1. **A glioma-immunology reference for ¶1.2** — candidates [42]–[44] now identified and
   metadata-verified (2026-10-01); you still need to read the passage you cite.
2. **Verify [5] Puchalski names the five anatomic structures** as you describe them. It is the
   source of your entire cohort and it is currently cited without having been read here.
3. **Verify the "53 methods across 28 datasets" figure** in [3] Nguyen before quoting it. It
   comes from this project's own related-work notes, not from the paper's abstract.

Everything else in this outline is either quoted from a source I read or comes from your own
artefacts.

---

## Quote verification

Every quotation above was checked against the PDF it is attributed to, not recalled. All nine
match verbatim:

| source | quotes checked |
|---|---|
| [1] Avila Cobos 2020 | 3 of 3 |
| [2] Sturm 2019 | 2 of 2 |
| [3] Nguyen 2024 | 2 of 2 |
| [4] Li 2026 | 2 of 2 |

*One note on transcription.* The Nguyen PDF renders its sentence as *"performing single- cell
experiments is expensive"* — a spurious space after the hyphen, introduced by the typesetting,
not by the authors. Quote it as **"performing single-cell experiments is expensive"**, which is
the sentence as written. The same artefact appears elsewhere in that PDF ("le v el" for
"level"), so check any further quotation you pull from it by eye.
