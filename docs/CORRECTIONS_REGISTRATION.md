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

## THE REGISTRATION UPDATE, IN SHORT — current as of 2026-09-15

*This section is the one to post. Everything below it is the evidence.*

**The question the study registered is unchanged:** can a tumour's own anatomy substitute for
ground truth when choosing a cell-type deconvolution method — and is a method that reproduces
known biology actually more accurate?

**The registered primary outcome still stands as computed.** Spearman **rho = 0.7501** between
the ACS ranking and the ranking from real ground truth (bootstrap CI [0.328, 0.957], p = 0.0020,
14 methods). No correction below changes that number. Three change what it *means*, and they
are the reason this update exists.

**What a reader must now be told, in order of how much it matters to the registered question:**

**1. The agreement is narrower than "anatomy tracks truth." Both arms share one reference
atlas** (C4). The ACS arm deconvolves Ivy GAP against GBmap; the accuracy arm builds its
mixtures from GBmap's cells. So part of rho = 0.7501 is two rankings agreeing about the
reference rather than about the tissue. This was not foreseen in the registration and is the
single most important thing to add to it.

**2. ACS rankings are reference-dependent, and about half the ordering survives an independent
atlas** (C4, corrected). Against two independent references built from other groups' data with
those groups' own labels, Spearman is **+0.5099** (Neftel) and **+0.3655** (Darmanis). An earlier
version of this correction reported 0.038 and 0.138 and called the ordering "collapsed" — that
was confounded with expression space and is **withdrawn**. Changing sequencing platform within
one atlas costs 0.183; changing expression space costs 0.084; changing the atlas costs about
half. A leaderboard computed against a single reference must be reported with that dependence
stated.

**3. The reference was built from log-transformed values treated as linear expression** (C7).
Every method solved a linear mixing model against a log-space signature. The pseudobulk arm is
internally consistent and survives this; the **anatomic arm is not**, and that is the arm the
anatomy test rests on. The full re-run on linear data is pending. This is a limitation to
declare now, not a result to wait for.

**Two further items are declarations rather than corrections.** The mRNA-content vector for
converting RNA shares to cell abundance is now **pre-specified and committed before use** (C8),
because choosing it afterwards would be selection on an outcome. And the Ivy GAP read counts
that several methods require **do not exist publicly** (C9), which bounds the method panel by
data availability rather than by choice.

**What has NOT changed:** the constraint file and its hash, the anatomic sample definition, the
within-tumour permutation null, the exclusion of ISH cluster samples, the negative controls, and
the separation of real methods from controls — which holds under every reference tested.

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
It survives only in part — see the correction to this entry's own numbers below, which reduced
the effect from "almost none survives" to "about half does".

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
roster alone does not do this.

**Both confounds named there have since been tested, and the table above is superseded.**

*Platform* was ruled out: GBmap's own Smart-seq2 subset against its 10x subset preserves the
ordering at **+0.817**, so it is not "the atlas's platform".

*Expression space* was not foreseen at all, and it was real. GBmap was read from a
log-transformed matrix while both alternatives were linear (**C7**, `OPEN_DEFECTS.md` D16).
Rebuilding GBmap from its own counts and repeating the comparisons:

| comparison | atlas | expression space | Spearman |
|---|---|---|---|
| GBmap_log vs GBmap_linear | same | DIFFERENT | **+0.9161** |
| GBmap_linear vs Neftel | different | same (linear) | **+0.5099** |
| GBmap_linear vs Darmanis | different | same (linear) | **+0.3655** |

**The direction of C4 stands and its strength is reduced.** Expression space costs almost
nothing (0.9161), so the atlas is the driver — but with the confound removed, **about half the
ordering survives a change of atlas** rather than almost none. The wording *"it does not
[survive]"* at the head of this entry is **too strong** and is corrected to: the ordering is
**reference-dependent**, with roughly half preserved across independent atlases.

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

## C7 · The reference is built from log-transformed values, and the registration's expression-space assumption is therefore false

> **Registered:** the deconvolution model is stated as bulk expression being a non-negative
> mixture of per-type reference profiles — the linear mixing model every one of the fifteen
> methods assumes.

**Correction.** The reference profiles are **not in linear expression space.** The GBmap
`.h5ad` carries two matrices and the pipeline reads `X`, which is
`log1p(counts x one size factor per cell)`. Verified to float32 precision: the implied size
factor is constant *within* a cell to 3.67e-07 (float32 epsilon), `corr(s_i, 1/total counts)`
is 0.999847, and `corr(expm1(X), raw counts)` within a cell is 1.000000 against 0.631419 for
`X` untransformed. The genuine integer counts sit unused in `raw/X`. Full measurement in
`docs/OPEN_DEFECTS.md` D16.

So every method solved `bulk (linear FPKM) ~= profile (log space) @ w`. `log1p` compresses a
100x marker to about 4.6x, which is a violation of the mixing model, not a change of units.

**Effect on the primary result: under measurement, and it may be substantial.** It predicts
two things the project already observes — the systematic tumour under-call of 0.33–0.51 at high
purity, and ordinal results (ACS) looking healthier than magnitudes (MAE, bias), since `log1p`
is monotone and rank statistics survive it while least-squares magnitudes do not. It also puts
**C4's attribution** in question: every agreeing arm of the reference 2x2 was log-vs-log and
every disagreeing arm was log-vs-linear, so "atlas" and "expression space" were confounded.
`scripts/build_gbmap_linear_reference.py` and the `gbmap_linear` arm of
`scripts/reference_sensitivity.py` break that confound; the result will be reported here
whichever way it falls.

---

## C8 · PRE-SPECIFICATION of the mRNA-content vector for Problem 2

**This is not a correction. It is a pre-specification, committed 2026-09-15 before the vector
was used to compute any ranking — the commit date is the evidence, which is why it lives in
`prespecified/` rather than in the gitignored `results/`.** Problem 2 of the two-problem design (`docs/TWO_PROBLEMS.md`)
converts RNA contributions into cellular abundance, and that conversion needs a per-type
mRNA-content vector. Choosing it after seeing which choice moves the ranking would be selecting
on an outcome, which this project's invariants forbid.

**Specified choice: the per-type mean total UMI count of GBmap's 10x cells, computed from
`raw/X`.**

**Chosen on mechanism, stated before any effect is known.** UMIs count molecules, which is
what mRNA content means. Three alternatives are rejected for stated reasons, not measured ones:

| candidate | rejected because |
|---|---|
| pooled over GBmap's platforms | **platform, not biology.** Smart-seq2 cells carry ~100x the reads of 10x cells (Tumor 785,065 vs 8,301) and the platform mix differs sharply by type (T_cell 52,590 10x vs 101 SS2). The pooled spread is **73.9x** against 2.897x within 10x. A ~737x figure in circulation has this signature. |
| Smart-seq2 read counts | full-length read count tracks library prep, not input RNA. Measured: Tumor/T_cell = **0.926** — the ordering does not even survive. |
| sums of `X` (what the code currently passes as `cell_totals`) | sums of **log** values. Compresses Tumor/T_cell from 2.245 to 1.292, so the correction would be applied at ~43% of its magnitude (D16). |

**The specified vector**, 4,062 10x cells stratified by type, in
`results/gbmap_cell_size_by_space.json`:

| type | mean UMI / cell | relative to median |
|---|---|---|
| Tumor | 8,301 | 1.672 |
| Macrophage_Microglia | 7,005 | 1.411 |
| T_cell | 3,697 | 0.745 |
| NK_cell | 4,691 | 0.945 |
| B_cell | 3,473 | 0.700 |
| Endothelial | 4,577 | 0.922 |
| Oligodendrocyte | 5,236 | 1.055 |
| Astrocyte | 10,060 | 2.027 |

**Sanity check, not a selection criterion:** Tumor/T_cell = **2.245** here against **1.69** in
the predecessor's independent frozen `L_median` (5,897/3,483) and 1.58–1.63x in the literature
figures — same direction, same order of magnitude, from a different atlas and a different
pipeline. Convergence of that kind is reassurance; it is not why the vector was chosen, and the
vector is not being swapped if the ranking dislikes it.

**Nothing has yet been scored with it.** The conversion it feeds has never been applied
(D12), so no published number moves when it is adopted; the first ranking computed with it will
be the first, and whatever it does to the ordering stands.

**Two caveats registered with it.** Astrocyte rests on 77 cells and is a diagnostic-only
sidecar known to be unidentifiable against Tumor, so no definitional claim will rest on it.
And this vector is **not** what the published leaderboard used: that used a uniform vector, so
the conversion was the identity there (D12).

---

## C9 · The Ivy GAP read counts do not exist publicly, which bounds the method panel

> **Registered:** the method panel is specified as a fixed list of deconvolution methods, with
> no statement that any were excluded for want of a compatible input.

**Correction: one class of method is excluded by data availability, and it should be named
rather than left as a gap in the panel.** Reference-free deconvolution — CDSeq is the
representative — models the bulk as multinomial **read counts**. Ivy GAP publishes **FPKM
only**.

**Established twice, from two independent distributions of the same data.**

The 2014-11-25 Allen Institute archive contains four files — `fpkm_table.csv`,
`columns-samples.csv`, `rows-genes.csv`, `README.txt` — and the expression values are
non-integer with columns summing to exactly 1000000.0.

**GEO GSE107559** (*Ivy Glioblastoma Atlas Project (RNA-Seq)*, Puchalski et al., **Science**
2018;360(6389):660–663, PMID 29748285) carries the same three data files and states plainly:
**"Raw data not provided for this record"**, and in the series summary, *"The raw RNA-Seq and
SNP array data will be submitted to dbGaP."* Its FPKM table is **bit-identical** to the
archive's — 25,874 lines, SHA-256 `72a97ea792a81098…` on both — so the two distributions are
the same data, and neither carries counts.

**Why this matters to the registered question rather than being mere housekeeping.** The
central worry about the ACS ordering is that it may be a property of the reference atlas (C4).
A **reference-free** method is the one instrument that could test that from the other side,
because it uses no reference at all. That instrument cannot be pointed at Ivy GAP. The study
therefore cannot fully separate "this is how methods rank" from "this is how methods rank
against this atlas" on the anatomic arm, and says so.

**What remains possible, and is not a substitute.** CDSeq can be run on the pseudobulk arm,
whose mixtures rebuild from the single-cell atlas's own integer counts. That tests the
reference-free comparison on synthetic mixtures, not on tissue.

**What would unblock it.** Counts or aligned reads for the 122 anatomic RNA-seq samples, which
per the GEO record go to **dbGaP** (BioProject PRJNA420740) under controlled access requiring an
institutional signing official. That is outside this project's reach and is recorded as a
limitation, not as pending work.

**Effect on the primary result:** none. It bounds what the panel could ever have contained.

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
