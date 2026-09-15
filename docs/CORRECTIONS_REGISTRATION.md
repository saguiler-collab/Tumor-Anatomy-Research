# Corrections to the OSF registration

**Registration:** <https://osf.io/dm2t8> · submitted `2026-09-10T03:32:03Z` · constraint file
`2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`

**The registration is permanent and is deliberately NOT amended.** That is the point of one.
This file is the published correction that travels beside it, and every item belongs in the
paper's methods section too. A registration with a correction a reader can check is worth more
than one nobody re-examined.

Each entry gives the registered wording verbatim, what is actually true, how that was
established, and whether the primary result moves.

---

## C1 · "Fifteen deconvolution methods"

> **Registered:** *"Fifteen deconvolution methods and two negative controls are then scored on
> how well…"* (and again: *"Fifteen cell-type deconvolution methods plus two negative controls
> are run…"*)

**Correction.** Fifteen is the plan and is what the confirmatory run executed. **Both pilot
runs the registration cites ran fourteen.** BayesPrism was added after them.

**Effect on the primary result:** none. The confirmatory run — the only one that postdates the
registration and the only one the primary result is computed from — ran fifteen.

---

## C2 · "Two published packages exceeded their wall-clock budget in one stage"

> **Registered:** *"Two published packages exceeded their wall-clock budget in one stage and
> fell back to reimplementations; the fallback is disclosed per method rather than absorbed."*

**Correction, in two parts.**

**The scope was understated.** True of `2026-09-05T2154`. In `2026-09-06T1103` and in the
confirmatory run the fallbacks occurred in **two** stages, not one.

**And both have since been measured as the genuine packages**, on inputs verified equivalent to
the leaderboard's — same 657-gene space by hash, same 88/22 donor split by name, same samples in
the same order, same normalisation:

| method | genuine R package | the reimplementation on the leaderboard | elapsed vs budget |
|---|---|---|---|
| `R:DWLS` | **0.7846** [0.657, 0.906] | 0.7385 | 2,474 s vs 2,400 s |
| `R:BayesPrism` | **0.8154** [0.710, 0.915] | 0.8769 | 2,045 s vs 2,400 s |

**BayesPrism's fallback was load-dependent, not inherent** — it completes inside the budget when
the machine is not otherwise busy. DWLS genuinely exceeds it.

The registration's own commitment — *"A Python reimplementation is never reported under a
published package's name"* — is upheld: the leaderboard rows are labelled
`python-reimplementation`, and these measurements are reported beside them, not substituted in.

**Effect on the primary result:** the archived leaderboard is unchanged. A future full re-run
would replace those two rows.

---

## C3 · "Converted from mRNA share to cell share by a per-type mRNA-content factor, applied once and centrally"

> **Registered:** *"A composition over eight cell types … Converted from mRNA share to cell
> share by a per-type mRNA-content factor, applied once and centrally."*

**Correction. The conversion is applied NOWHERE. It is the identity.**

`cell_size` on the reference every published number was computed against is **1,000,000 for all
eight cell types — spread 1.0000**. Dividing by a constant vector and renormalising does
nothing.

**Cause.** `build_from_h5ad` captures each cell's raw library size before normalising and warns
in a comment that failing to do so "would turn the cell-size correction into a no-op".
`run_benchmark` then rebuilds the reference from the already-normalised matrix **without passing
`cell_totals`**, so it falls back to `expression.sum(axis=0)`, which is 1e6 per cell by
construction.

**Established two ways**, the second requiring no code reading: the `cell_size` vector is
recorded in `results/bisque_remeasured.json`; and Bisque re-measured with the conversion
**skipped** reproduced the archived table produced with it **applied**, bit-identical — max
absolute difference 0.000000 across 122 samples and 8 cell types. A transform applied and
skipped cannot agree to the bit unless it is the identity.

**Effect on the primary result: none, and this is the important part.** Every published number
was computed with an identity transform in the chain, which is the same as without it. Nothing
needs re-running to *correct* anything. What the correction changes is the **description**: this
study reports **mRNA proportions**, not cell proportions.

**Consequence for a reader.** Where the paper compares against studies reporting cell fractions,
that difference is real and must be stated. Whether to start applying the conversion is a design
decision, recorded in `docs/OPEN_DEFECTS.md` D12 with both defensible options, and it must not
be settled by which option improves a score.

---

## C4 · The ordering is reference-dependent — not a correction to registered text, but material and unforeseen

The registration pre-specified the yardstick, the threshold and the constraint file. It did not
foresee, and did not test, whether the method ordering survives a change of **reference atlas**.
It does not.

Two independent alternative references, each built from another group's data using **that
group's own** cell-type labels, each compared against GBmap on the same sub-roster and the same
gene space:

| comparison | Spearman between ACS orderings |
|---|---|
| GBmap vs Darmanis 2017 (5 types, 4 donors) | **+0.138** |
| GBmap vs Neftel 2019 (4 types, 20 donors) | **+0.038** |
| **Darmanis vs Neftel** | **+0.710** |
| GBmap 5-type vs GBmap 4-type *(same atlas, different roster)* | **+0.908** |

The two alternatives agree with each other and neither reproduces GBmap's ordering; changing the
roster alone does not do this. A confound is not excluded — both alternatives are Smart-seq2 and
GBmap is 87% 10x — so the finding is either "the ordering depends on the atlas" or "it depends on
the atlas's platform". Both are first-order.

**This bears on the primary result.** The registered primary outcome is Spearman rho between the
ACS ranking and the accuracy ranking on held-out synthetic mixtures. **Both arms use GBmap** —
the ACS arm deconvolves Ivy GAP against it, and the mixtures in the accuracy arm are built from
its cells. So part of the reported rho = 0.750 is agreement about the reference rather than about
the tissue.

**rho = 0.750 is not wrong.** It is narrower than "anatomy tracks truth": it is agreement between
two GBmap-based rankings. The registered threshold was never designed to test reference
invariance, and the paper must say so rather than let the reader assume otherwise. See
`docs/OPEN_DEFECTS.md` D14.

---

## C5 · Cohort figures quoted from the project documentation

> **Registered / protocol:** *"270 laser-microdissected RNA-seq samples across 41 tumors"*

**Correction.** The 2014-11-25 RNA-seq release contains **37** tumours: 10 anatomic + 34
ISH-cluster, with 7 contributing to both. The 270 / 122 / 10 figures are exact. The 41 likely
counts the whole Ivy GAP project including tumours that contributed no RNA-seq.

**Effect on the primary result:** none. No constraint, threshold or cohort rule reads the number.

---

## C6 · Neftel GSE131928 described as 21 donors

**Correction.** GSE131928 is **9 patients / 21 samples**; the sample prefixes include repeated
samples of MGH105. Adult cells used for the reference built here come from **20 samples across
the adult subset**.

**Effect on the primary result:** none. Neftel is not part of the registered analysis.

---

## What is NOT corrected, and why that matters

The following were pre-registered and are unchanged, which is what makes the corrections above
readable as corrections rather than as a rewrite:

- the **constraint file** and its hash — never edited, and the registered hash still matches;
- the **primary outcome** and the **pre-registered threshold** (rho ≥ 0.6 with a bootstrap CI
  excluding zero);
- the **negative controls** and the within-tumour permutation null;
- the rule that **ACS never selects the method**, and that the frozen method is chosen by the
  outcome-blind pseudobulk benchmark;
- the **nesting rule**, the **anatomic-only scoring rule**, and the exclusion of ISH-cluster
  samples from ACS.
