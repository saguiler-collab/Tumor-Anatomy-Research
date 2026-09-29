# Why do the methods put B cells above T cells?

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
