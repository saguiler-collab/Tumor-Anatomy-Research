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

**Which registered sentences are affected.** Each correction below quotes the registered
wording verbatim so a reader can check it against the entry at <https://osf.io/dm2t8>:

| registered sentence | section | correction |
|---|---|---|
| *"Counts are normalised to CPM, **never log-transformed for deconvolution**, the mixing model is additive on the linear scale."* | Analysis Plan → Transformations | **C7** — true of the bulk, **false of the reference** |
| *"Converted from mRNA share to cell share by a per-type mRNA-content factor, **applied once and centrally**."* | Variables → Measured variables | **C3** — applied **nowhere**; the conversion is the identity |
| *"**Accuracy: mean absolute error against known composition**"* | Variables → Measured variables | **C12** — the two sides are in different units |
| *"T cells, because this project's own benchmark places them at the **detection floor**"* | Overview → Foreknowledge | **C11** — backwards; T cells are the largest **over**-call, 4.6–7.9x |
| *"GBmap … **is never scored and contributes no anatomic claim**."* | Sampling → Data collection | **C13** — true, and it understates GBmap's role in the result |
| *"**Fifteen** cell-type deconvolution methods…"* | Research Design → Study design | **C1** — both cited pilot runs ran fourteen |
| *"Two published packages exceeded their wall-clock budget in **one stage**"* | Other → Deviations | **C2** — two stages; both since measured as the genuine packages |
| *"A wider 270-sample cohort spanning **37 patients**"* | Research Design → Study design | **C5** — 37 **tumours**; the Sample-size section says it correctly |
| *"C7 carries double weight because a three-step monotone chain is far harder to satisfy by chance"* | Research Design → Study design | **C14** — sound, but **C6** carries the ranking |
| *"Genes are restricted to a marker subset selected from the reference alone"* | Analysis Plan → Transformations | **C14** — outcome-blind as claimed; subset **size** moves ACS by 0.31 |

Three further entries are not corrections to registered text: **C4** and **C10** are material
findings the registration did not foresee, and **C8** is a pre-specification.

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

**One further item is a declaration rather than a correction.** The mRNA-content vector for
converting RNA shares to cell abundance is now **pre-specified and committed before use** (C8),
because choosing it afterwards would be selection on an outcome.

**And one item is the project's first external validation — which did not pass.** C10: a
pre-specified prediction that ACS would reproduce a ranking established by imaging mass
cytometry failed, in the direction named in advance as most damaging, and the follow-up showed
ACS moving by 0.31 on nothing but marker-set size. Reported as INCONCLUSIVE on nine tumours. It
does not change rho = 0.7501; it changes what that number can be said to mean.

**And one item is a retraction of a correction, made the same day it was written.** C9 claimed
no public read counts exist for Ivy GAP. They do — the Allen Institute portal serves per-sample
RSEM `genes.results` files carrying `expected_count`, which GEO does not mirror. All 270 have
been fetched. This also closed the provenance chain on the published matrix for the first time:
it is exactly raw FPKM times one scale factor per sample (`corr = 1.00000000`, factor CV
1.4e-15), and that scaling is a no-op for this pipeline, which renormalises every bulk column
anyway. **Reference-free deconvolution can now be pointed at Ivy GAP itself**, which is the one
instrument that can test C4 without a reference.

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

> **Registered, under *Research Design → Study design*:** *"A wider 270-sample cohort spanning
> **37 patients**"* — and under *Sampling → Sample size*: *"The release contains 270 RNA-seq
> samples from **37 tumors**."*
>
> **Project protocol, separately:** *"270 laser-microdissected RNA-seq samples across 41 tumors"*

**The registration's *Sample size* figure is correct and its *Study design* figure is not, and
they disagree with each other.** Measured against the archive and reconciled against the live
portal: **270 samples, 37 TUMOURS** — 10 anatomic + 34 ISH-cluster, with 7 contributing to both.
So *"37 tumors"* is exact; *"37 patients"* is the wrong unit. Ivy GAP as a whole is 42 tumours
from 41 patients, which is where the protocol's 41 comes from — it counts tumours that
contributed no RNA-seq.

The 270 / 122 / 10 / 9 figures are exact, and the anatomic structure counts agree with the
portal's own table exactly (CT 30, IT 24, LE 19, MVP 25, PAN 24). Nine samples exist on the live
portal but not in the 2014-11-25 archive — a release skew, all in the ISH-cluster study, not a
parsing error.

**Effect on the primary result:** none. No constraint, threshold or cohort rule reads the number.

---

## C6 · Neftel GSE131928 described as 21 donors

**Correction.** GSE131928 is **9 patients / 21 samples**; the sample prefixes include repeated
samples of MGH105. Adult cells used for the reference built here come from **20 samples across
the adult subset**.

**Effect on the primary result:** none. Neftel is not part of the registered analysis.

---

## C7 · The reference is built from log-transformed values, and the registration's expression-space assumption is therefore false

> **Registered, under *Analysis Plan → Transformations*, verbatim:** *"Counts are normalised to
> CPM, **never log-transformed for deconvolution**, the mixing model is additive on the linear
> scale."*

**Correction. This registered sentence is false as it applies to the reference.** The intent
was right and the bulk side is exactly as described — Ivy GAP FPKM, rescaled per column, never
logged. But the **reference profiles are not in linear expression space.** The GBmap `.h5ad`
carries two matrices and the pipeline reads `X`, which is
`log1p(counts x one size factor per cell)`. Verified to float32 precision: the implied size
factor is constant *within* a cell to 3.67e-07 (float32 epsilon), `corr(s_i, 1/total counts)`
is 0.999847, and `corr(expm1(X), raw counts)` within a cell is 1.000000 against 0.631419 for
`X` untransformed. The genuine integer counts sit unused in `raw/X`. Full measurement in
`docs/OPEN_DEFECTS.md` D16.

So every method solved `bulk (linear FPKM) ~= profile (log space) @ w`. `log1p` compresses a
100x marker to about 4.6x. **The mixing model is not additive on the linear scale on both
sides**, which is precisely what the registered sentence asserts, and it is a violation of the
model rather than a change of units.

The asymmetry is worth stating because it is diagnostic: the **pseudobulk arm is internally
consistent** — its mixtures are summed from the same log matrix the profile is averaged from, so
a linear model genuinely holds there and that arm survives. The **anatomic arm is not**, because
its bulk is real linear FPKM. That is one defect explaining why the benchmark looked healthy
while accuracy on tissue was poor.

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

## C9 · RETRACTED THE SAME DAY IT WAS WRITTEN — the Ivy GAP read counts ARE public, and they are now in hand

**This entry originally claimed that no public read counts exist for Ivy GAP, and that the
count-based methods were therefore bounded by data availability. That was wrong. It is kept
here rather than deleted, because a correction file that quietly edits its own errors is worth
nothing.**

### What the retracted entry said, and what it rested on

It asserted: *"Ivy GAP publishes FPKM only"*, on two pieces of evidence — the 2014-11-25 Allen
Institute archive contains four files and none of them are counts, and **GEO GSE107559** states
**"Raw data not provided for this record"** with *"The raw RNA-Seq and SNP array data will be
submitted to dbGaP."* Both of those statements are true. The conclusion drawn from them was not.

**The error: GEO is one distribution channel and the Allen Institute portal is another.** The
portal's own download page publishes, per sample, both *"un-normalized gene-level FPKM and TPM
values"* and *"anonymized BAM files"*. Neither appears in GEO. I concluded from the most
authoritative-looking source instead of enumerating the channels.

### What is actually available

Per sample, the portal serves the **RSEM `genes.results`** file — 25,873 genes, ~1.5 MB:

```
gene_id  transcript_id(s)  length  effective_length  expected_count  TPM  FPKM
1        NM_130786_3       1766.00  1604.54          45.56           3.41  2.28
```

**`expected_count` is the count layer.** All 270 samples have been fetched
(`scripts/fetch_ivygap_counts.py`, provenance and per-file SHA-256 in
`data/raw/ivygap_counts/provenance.json`). Anonymized BAMs are also public and unauthenticated
(414 MB each, ~112 GB for 270); they are **not** needed, because RSEM has already counted, and
re-counting would substitute this pipeline's choices for the authors'.

### And it closes the provenance chain on the published matrix, which nothing else had

Mapping RSEM's Entrez ids onto Ivy GAP's internal `gene_id` through `rows-genes.csv`:

| check | result |
|---|---|
| genes mapped / shared with the published matrix | **25,873 / 25,873** |
| `corr(published FPKM, per-sample raw FPKM)` | **1.00000000** |
| published ÷ raw, across genes within a sample | **a single constant**, CV **1.4e-15** |
| the constant, three samples | 0.8596, 1.0576, 1.0647 |

So the published "normalized" matrix is exactly **raw FPKM × one scale factor per sample** —
the normalisation the authors describe as being based on genes not enriched in particular
anatomic structures. Two consequences:

1. **The per-sample RSEM files are provably the source of the matrix this study deconvolves.**
   Verified to float precision, not assumed.
2. **That normalisation is a no-op for this pipeline**, which rescales every bulk column to sum
   to 1e6 before deconvolving. A per-sample constant is removed exactly by that step. Worth
   stating because it was an untested assumption until now.

*(The alignment matters and was nearly got wrong: both id spaces are numeric, so a
position-based join silently "works" — it shares 7,741 of 25,873 ids by coincidence and gives
`corr = 0.12`. The same trap as the reference row-alignment bug this project has already been
bitten by once.)*

### Effect on the registered result and on the panel

**No published number moves.** What changes is what is *possible*:

- **CDSeq is unblocked on the anatomic arm.** That matters specifically for **C4**: the central
  worry is that the ACS ordering is a property of the reference atlas, and a reference-free
  method is the one instrument that can test it without a reference. It can now be pointed at
  Ivy GAP itself, not only at synthetic mixtures.
- **The bulk's true library sizes are known** (12.5M, 14.9M, 16.8M expected counts in the three
  samples checked), which the two-problem framing needs.
- `expected_count` is RSEM's **posterior expectation** and is fractional (45.56, not 46). It is
  not a raw integer count and is not presented as one; rounding, where a model requires
  integers, is declared at the point of use.

**Registered wording affected:** none directly. This corrects a correction.

---

## C10 · The first external validation of the registered hypothesis, and it did not pass

> **Registered:** the primary outcome is Spearman rho between the ACS ranking and the ranking
> from real ground truth on held-out synthetic mixtures. The registration did not foresee any
> test of ACS against a ground truth measured **outside** this project.

**This is new evidence, not a correction to registered text, and it is material.**

Ajaib et al. (**Neuro-Oncology** 2023;25(7):1236–1248) applied single-cell resolution **imaging
mass cytometry** — 33 antibodies, protein — to ten IDHwt GBM samples with matched bulk RNA-seq,
and scored deconvolution approaches against it. Their core contrast is one algorithm,
MCPcounter, differing only in marker set: GBM-specific **r = 0.37** (immune) / 0.43
(neoplastic), MCPcounter default 0.27, **GBmap-derived 0.06** / 0.22, CIBERSORTx 0.05 / 0.02.

A prediction was pre-specified and committed before the run: **ACS(MCP_GBM) > ACS(MCP_GBmap)**.

**It failed.** ACS gave MCP_GBmap **0.6667** against MCP_GBM's **0.5128**; paired over the same
nine tumours, delta **+0.1389**, 95% CI [+0.0444, +0.2500], **p = 0.0090**. ACS ranked the
marker set that protein ground truth places near chance **above** the one it places first, and
MCP_GBM did not beat its own permutation null (p = 0.0694). The pre-specification had named this
exact outcome in advance as the most damaging one.

**A confound named in advance turned out to be real**, and it does not rescue the prediction so
much as deepen the problem. Rebuilding GBmap's markers at Ajaib's density and changing nothing
else reverses the direction — GBmap falls to **0.3846**, below MCP_GBM, and becomes
indistinguishable from its null (p = 0.4357). But the same table shows ACS moving from
**0.3846 to 0.6923** on *nothing but the number of marker genes per cell type*, with
non-overlapping CIs, non-monotonically. **That swing of 0.31 is about three-quarters of the
range across which the published leaderboard ranks fifteen methods.**

**Effect on the registered primary result.** The number does not move: rho = 0.7501 stands as
computed. What moves is its interpretation, and this is the correction a reader most needs:

> The registered outcome measures agreement between two rankings that share a reference atlas.
> On its first test against a ground truth from outside the project, ACS did not reproduce the
> external ordering — and ACS was shown to be movable, by an amount comparable to the spread it
> is used to rank methods across, by a preprocessing choice unrelated to the tissue.

Reported as **INCONCLUSIVE** on nine evaluable tumours rather than as a refutation, which is
what the power supports. Full analysis, including what is explicitly *not* claimed, in
`docs/EXTERNAL_VALIDATION.md`. Nothing in `constraints.py` was touched; its hash is unchanged.

---

## C11 · The T-cell exclusion was registered for a reason the data contradicts

> **Registered, under *Overview → Explanation of foreknowledge*, verbatim:** *"two exclusions
> made in advance (**T cells, because this project's own benchmark places them at the detection
> floor**; and the 148 expression-labelled samples, as circular)"* — and in the frozen constraint
> file itself: *"This pipeline's own synthetic benchmark puts T cells at the detection floor, so
> a T-cell constraint would score noise."*

**Correction: the stated mechanism is not merely inaccurate, it is backwards.** T cells are not
at a detection floor. They are the **largest over-call in the entire panel.** In high-purity
mixtures, true mean T-cell content is **0.044** and the methods predict:

| MuSiC | BayesPrism | NNLS | CIBERSORTx | SCDC |
|---|---|---|---|---|
| 0.202 (**4.6x**) | 0.234 (5.3x) | 0.263 (6.0x) | 0.339 (7.7x) | 0.350 (**7.9x**) |

The methods are not failing to *detect* T cells. They are **inventing** them — T cells are where
the missing tumour mass goes when every method under-calls tumour by 0.33–0.51.

**The conclusion drawn from the false premise happens to survive, and that should be said
plainly rather than used as cover.** A constraint on a population estimated at five to eight
times its true value would indeed score noise, so excluding `T_cell` was the right call. It was
made for the wrong reason, and the registration states that reason.

**The exclusion stands and the constraint file is not edited.** Changing it would alter the
hash and void every result computed under it, which is the rule this project set itself.

**Effect on the primary result: none numerically, and it is material to interpretation.** The
exclusion is exactly what makes the anatomy test **structurally blind to the panel's largest
error mode**. A method could invent T cells without limit and still score ACS = 1.000. Any claim
that this framework can certify a method for clinical use has to carry that sentence.
`docs/CLINICAL_READINESS.md` §3.

---

## C12 · The accuracy arm scored mRNA proportions against cell-fraction truth

> **Registered, under *Variables → Measured variables*, verbatim:** *"**Accuracy: mean absolute
> error against known composition**, donor-equally aggregated, excluding the Astrocyte column."*

**Correction.** The two sides of that error are in **different units**. `pseudobulk.py` declares
its truth as `# mixtures x cell types, rows sum to 1 (CELL fractions)` and its own docstring
warns that *"scoring cell-fraction estimates against RNA-fraction truth is a classic way to get
this wrong."* Because the mRNA-to-cell conversion is the identity in this pipeline (**C3**), the
estimates being scored are **mRNA proportions**. The benchmark made the exact error its source
file warns against, and it arrived by way of the bridge that was never built.

This follows from C3 rather than being independent of it: had the conversion been applied, the
units would have matched.

**Effect on the primary result: it affects the accuracy half.** Every MAE, RMSE and bias in the
benchmark compares mRNA share with cell share. The **rank** correlation is more robust than the
magnitudes — the two quantities are monotonically related within a type — but the accuracy
*values* are not the quantity the registration names. `docs/OPEN_DEFECTS.md` D15. The fix is
implemented (`truth_mrna`, computed from pre-normalisation library sizes and gated so it is
never derived from a log matrix) and awaits the full re-run.

---

## C13 · "GBmap ... is never scored and contributes no anatomic claim" — true, and it understates GBmap's role

> **Registered, under *Sampling → Data collection procedures*, verbatim:** *"Reference dataset:
> GBmap Core, a 338,564-cell single-cell atlas of glioblastoma, used only to build the cell-type
> reference profiles and the synthetic mixtures. **It is never scored and contributes no anatomic
> claim.**"*

**Both sentences are literally true and the second is misleading**, which is why it needs a
correction rather than a retraction. GBmap is never scored, and it makes no anatomic claim. But
it **substantially determines the ordering that the study reports**:

- with expression space controlled, only **+0.5099** (Neftel) and **+0.3655** (Darmanis) of the
  ACS ordering survives a change of atlas (**C4**);
- an independent group, using **imaging mass cytometry** on matched tissue, scored GBmap-derived
  markers at **r = 0.06** for immune cells — near chance, and worst of the four approaches they
  tested (**C10**);
- and both arms of the registered primary outcome use it, so rho = 0.7501 is substantially two
  GBmap-based rankings agreeing about GBmap.

**Effect: the registered sentence should be read as scoped to anatomic claims only.** It is not
a statement that the reference is neutral with respect to the result, and a reader would
reasonably take it as one. The atlas is the single most influential choice in the study.

---

## C14 · Two registered rationales that measurement qualifies rather than refutes

Neither is an error. Both are statements the registration makes as *reasons*, which later
measurement puts a number on, and a reader checking the registration against the results will
notice both.

### "C7 carries double weight because a three-step monotone chain is far harder to satisfy by chance"

> **Registered, under *Research Design → Study design*, verbatim.**

The premise is sound — a three-step chain *is* harder to satisfy by chance. But the weighting
does not make C7 the constraint that decides anything. **C6 does.** Dropping C6 and rescoring
moves the ordering to rho = 0.678 with **13 of 14 methods changing rank**; dropping any other
single constraint moves it far less, and C3 is unanimous across methods and so separates
nothing. So the ranking rests largely on **one** pairwise constraint that carries **single**
weight, while the double weight sits on a constraint that is not load-bearing for the ordering.
`docs/OPEN_DEFECTS.md` D8.

### "Genes are restricted to a marker subset selected from the reference alone, with no reference to any outcome or structure"

> **Registered, under *Analysis Plan → Transformations*, verbatim.**

**True as written, and outcome-blindness is not the whole of the problem.** The selection is
genuinely blind to outcome and to structure. What it is not blind to is the **reference**, and
the *size* of the resulting subset turns out to move the score on its own: holding the atlas,
the algorithm, the bulk, the constraints and the tumours fixed and changing **only markers per
cell type**, ACS moved **0.3846 → 0.6923** — non-overlapping CIs, non-monotone, and from
indistinguishable-from-null (p = 0.4357) to significant (p = 0.0016). That swing of **0.31** is
about three-quarters of the range across which the leaderboard ranks fifteen methods. **C10**,
`docs/EXTERNAL_VALIDATION.md`.

So the registered sentence describes a real safeguard against one failure mode, and a different
one was left uncontrolled: not *which* genes, but *how many*.

---

## What the registration got RIGHT, including one prediction it nailed

A corrections file that lists only errors misrepresents the document it corrects. Three things
in this registration were pre-specified, load-bearing, and vindicated:

**The "no single best method" prediction was correct, and it was specific.**

> **Registered:** *"Is any single method identifiable as best? Hypothesis: no. With nine
> evaluable tumors and a weighted denominator of 65, ACS is expected to be too coarse to
> separate adjacent methods; any such claim is reported INCONCLUSIVE unless a paired bootstrap
> separates them."*

Measured: **a four-way tie at 0.9846** (NNLS, SVR, elastic net, EPIC) and a three-way tie at
0.9692 (CIBERSORTx, CIBERSORTx S-mode, Bisque). MuSiC's first place rests on **one satisfied
pair out of 57**. The registration predicted the exact failure mode, gave the mechanism, and
fixed the reporting rule in advance. That is what a registration is for.

**The negative-control commitment held.** *"If controls score highly, that is published as the
finding rather than answered by retuning the constraints."* The controls scored **0.3385**
(random) and **0.1385** (shuffled signature), neither beating its null, and the separation from
real methods holds under **every** reference atlas tested. The constraint file was never edited
— including when C5 failed its own external check, and when the pre-specified IMC prediction
failed (**C10**).

**The declared exclusions were enforced in code, not by convention.** The 148
expression-labelled samples never reached the scoring set; an assertion aborts the run if one
does.

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
