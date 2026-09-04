# Methods

Ten methods, one contract. Every one receives the same frozen `DeconvolutionInput` and
returns cell fractions over the same eight-type roster in the same order. Simplex
projection and the RNA-to-cell-size correction happen once, centrally, in
`DeconvolutionMethod.fit_predict` — so no method can quietly report a different
estimand from its neighbours.

## The estimand

Every method reports **cell fractions**, not RNA fractions.

Reference-based deconvolution against an expression signature natively recovers each
type's share of the *transcript pool*. A large, transcriptionally active tumour cell
contributes more mRNA than a small lymphocyte, so RNA share and cell share are
different numbers. Dividing by a per-type mean transcript content and renormalising
converts one to the other, and `to_cell_fractions` does it identically for all ten.

This matters for scoring as much as for reporting: pseudobulk truth is recorded as cell
fractions too. Scoring cell-fraction estimates against RNA-fraction truth manufactures a
systematic bias that looks exactly like a method difference.

## The cell-type roster

Seven primary types — Tumor, Macrophage_Microglia, T_cell, NK_cell, B_cell,
Endothelial, Oligodendrocyte — carried over unchanged from the TCGA project so the two
cohorts stay comparable.

**Astrocyte is a diagnostic-only sidecar**, never summed into the seven. In the
predecessor project GFAP alone carried 38.8% of the squared reconstruction error and
Astrocyte proved unidentifiable against Tumor under NNLS and SVR. It is retained here
rather than dropped for two reasons: Ivy GAP's leading edge is genuinely
astrocyte-rich normal brain, so the identifiability question is materially different
from TCGA's cellular-tumour samples; and Tumor/Astrocyte collinearity is precisely the
failure mode Elastic Net and the Bayesian prior are supposed to improve on. Dropping it
would delete the test case.

No definitional anatomic claim may rest on the sidecar column —
`constraints.validate()` enforces that, so a known-unresolvable column cannot carry
a constraint and cannot move the ACS.

---

## Baselines

### NNLS — `deconv/classical.py`
`min ||Sw - b||²`, `w ≥ 0`. No regularisation, so when two reference columns are
near-collinear the solution is essentially arbitrary along that direction. It is here as
the control: the improvement claimed for the regularised and Bayesian methods has to be
measured against something, and this is the solver that failed informatively on TCGA.

### SVR — `deconv/classical.py`
CIBERSORT-style ν-support-vector regression. The ε-insensitive loss ignores residuals
inside a tube around the fit, so a few wildly-off genes cannot drag the whole solution —
which is why it beat NNLS on the TCGA held-out benchmark and became that project's
frozen method. ν is swept over `{0.25, 0.5, 0.75}` and chosen by reconstruction RMSE, a
criterion that uses only expression data.

**Iteration is bounded, and that is not a detail.** libsvm defaults to `max_iter=-1`.
On synthetic data a fit takes 0.03 s; on the real frozen GBM signature it takes **5.2 s**
— 175× slower — because SMO grinds against exactly the collinear geometry (Tumor against
Astrocyte above all) that this project exists to study. An unbounded solver inside a
pipeline stage turns a benchmark into a hang with no diagnostic, so `MAX_ITER` caps it
and `n_hit_iteration_cap_` counts the fits that hit the cap. On the 122-sample Ivy GAP
cohort none did, so the bound costs nothing here; it exists so a harder reference
produces a reported number instead of a stall.

---

## The two additions

### Elastic Net — `deconv/elastic_net.py`

```
minimise  (1/2n)||Sw - b||²  +  α[ ρ||w||₁ + (1-ρ)/2·||w||² ]   subject to w ≥ 0
```

**The problem it targets.** GBM signature columns are not independent: tumour cells
adopt an astrocyte-like state, TAM-BDM and TAM-MG share most of their myeloid
programme, NK and T cells share cytotoxic machinery. Near-collinear columns make the
least-squares problem ill-conditioned — many quite different `w` reconstruct `b` almost
equally well, and noise decides which one the solver lands on. The L2 term makes the
problem strictly convex and *shares* weight smoothly between collinear columns instead
of letting the solver pick one arbitrarily. The L1 term drives genuinely absent types to
exact zero, which matters for leading-edge samples where lymphocytes really are absent.

**An honest note on L1 under a simplex constraint.** After renormalising to sum 1,
`||w||₁ ≡ 1` — so on the normalised scale the L1 term is a constant and does nothing.
It is not useless, because the fit happens on the *unnormalised* scale where `w` carries
magnitude, and the pattern of zeros survives renormalisation. But ρ behaves differently
than in ordinary regression, and pure lasso (ρ=1) mostly rescales rather than selects.
The ρ grid is therefore ridge-leaning: `{0.05, 0.2, 0.5}`.

**Tuning is outcome-blind and per sample.** α is chosen from a path anchored at
`α_max = max|Sᵀb| / (n·ρ)` — the smallest penalty that zeroes every coefficient — and
descending by a fixed ratio. A grid in *absolute* units would apply completely different
amounts of shrinkage to a sparse leading-edge sample and a dense cellular-tumour one,
i.e. silently a different method per anatomic structure. Selection is by held-out-*gene*
cross-validation on reconstruction error: no survival, no anatomic label, no known
composition. The chosen (α, ρ) per sample is retrievable via `selection_report()`.

### Bayesian — `deconv/bayesian.py`

```
w      ~ Dirichlet(α)
p_gk   = w_k·S_gk / Σ_k' w_k'·S_gk'
z_g·   ~ Multinomial(b_g, p_g·)          latent: which type each read came from
b_g    = Σ_k z_gk
```

Conditioning on the latent allocation `z` makes the composition update conjugate —
`w | z ~ Dirichlet(α + Σ_g z_·k)` — which gives an **exact** Gibbs sampler in two lines:
no variational gap, no Metropolis acceptance rate to tune. It is the same data-augmentation
scheme underneath BayesPrism's per-sample step.

**Why it is more robust, concretely:**

1. **The prior regularises.** With α < 1 a type with no real evidence settles near zero
   instead of absorbing residual noise, without L1's hard zeros.
2. **It reports uncertainty.** NNLS says "T cells: 0.8%". This says "0.8%, 95% credible
   interval 0.1–2.4%". Ivy GAP's leading-edge dissections are small and low-complexity;
   a method that cannot say when it does not know will produce confident nonsense there.
3. **Collinear types share posterior mass** rather than fighting over it. Where Tumor
   and Astrocyte are unidentifiable the posterior becomes wide and correlated — the
   correct answer to an ill-posed question, rather than NNLS's arbitrary corner.

Ivy GAP ships FPKM, and the multinomial step needs counts, so each profile is rescaled
to a fixed pseudo-count total. That total is a real modelling choice, not a detail: too
small and the prior dominates; too large and credible intervals become implausibly tight,
because FPKM is not a count of independent reads.

**`bayesian_hierarchical` pools by patient, never by structure.** Pooling across samples
that share an anatomic structure would hand the model the label it is scored on
reproducing. Pooling by patient is well motivated — a tumour's blocks share that
patient's immune baseline — and because a patient's blocks span *different* structures,
it shrinks toward the patient mean and makes between-structure differences *harder* to
detect. The bias runs conservative.

---

## The four published tools

Faithful reimplementations in `deconv/reference_based.py`; the genuine R packages run
via `deconv/r_bridge.py` when three conditions hold: `Rscript` on PATH, the package
installed, and a **cell-level** single-cell export on disk.

That third condition is binding and worth stating plainly. MuSiC, Bisque and SCDC do not
consume a signature matrix — they consume individual cells labelled with donor and cell
type, because computing statistics *across* donors is the whole method. Handing them a
collapsed signature matrix would run the code while removing the thing that makes it
that method. `release/implementation_report.json` records, per tool, which
implementation ran and why.

| Tool | What it fixes | Why it is in this roster |
|---|---|---|
| **MuSiC** | genes whose expression is inconsistent *across donors* | Ivy GAP spans ~37 tumours with strong between-patient heterogeneity — cross-subject consistency is the right currency |
| **DWLS** | a few enormously expressed genes monopolising the residual | Ivy GAP lymphocyte fractions are low single digits; without this they are numerically invisible |
| **Bisque** | the platform gap between droplet scRNA-seq and bulk | Ivy GAP bulk is 2014-era laser-capture RNA-seq — about as wide as that gap gets |
| **SCDC** | dependence on which reference you happened to pick | the only method here that can measure its own reference sensitivity |

### Deliberate deviations from the published algorithms

- **MuSiC** — the recursive tree (`music_prop.cluster`), for hierarchically related cell
  types, is not implemented. The eight-type roster is flat by construction, so the tree
  would have nothing to recurse on. Extreme gene weights are capped at the 99th
  percentile; without the cap a near-zero-variance gene takes an astronomical weight and
  the "weighting" degenerates into fitting one gene.
- **DWLS** — the published implementation selects the dampening constant by a
  multinomial-likelihood criterion over bootstrapped subsets. This uses the same
  subset-resampling idea with a solution-variance criterion. It selects the same
  constant on the test fixtures but is not guaranteed to agree in every case. Selection
  runs short chains, since the *ordering* of stability across constants is established
  long before the fits converge.
- **Bisque** — runs in `use.overlap = FALSE` mode, because Ivy GAP shares no subjects
  with any public GBM single-cell atlas. The gene-wise assay transform is estimated from
  marginal distributions rather than paired subjects. **This is Bisque's documented
  degraded mode, and any conclusion about Bisque's ranking here carries that caveat.**
- **SCDC ENSEMBLE** — ensemble weights are searched on a simplex grid minimising
  reconstruction error of the observed bulk, an unsupervised criterion. Given a single
  reference it degenerates to plain SCDC and sets `degenerate_ = True` rather than
  reporting a one-member "ensemble". Reconstructions are blended, not signatures —
  blending signatures would invent a chimeric reference no dataset supports.

---

## The two negative controls

A leaderboard without controls cannot distinguish "the constraint set detects
competence" from "the constraint set is easy." Both run through the identical pipeline
on identical inputs, scored identically.

| Control | What it is | What it tests |
|---|---|---|
| `control_random` | Dirichlet draws (α = 0.7), ignoring the expression data entirely | The floor: whatever the constraint set gives away for free |
| `control_shuffled_signature` | Real NNLS against a signature whose **gene rows** are permuted | The sharper control: structured, bulk-dependent output that means nothing |

The shuffle permutes *rows*, not entries. Permuting entries independently would also
destroy each gene's overall expression level, making the signature obviously pathological
and the control too easy to beat. Row permutation keeps the signature statistically
plausible and breaks only the gene-to-cell-type correspondence — precisely the
information deconvolution is supposed to use.

The concentration is 0.7 rather than 1.0 for the same reason. A flat Dirichlet produces
near-uniform compositions, which an ordinal constraint set rejects easily; 0.7 gives the
lumpy compositions a real-but-wrong method might produce.

**If a control scores well, that is reported as-is.** The protocol's instruction is
explicit: publish the constraint file anyway and state what a harder set would need. The
runner emits that verdict itself rather than leaving it to the reader, and it never
retunes the constraints.

## Degenerate methods are labelled, not silently renamed

Three methods can run to completion while not actually being themselves. Each detects
and reports it, because a leaderboard row saying "MuSiC" that is arithmetically NNLS is a
mislabelled result, not a result.

| Method | Degenerates when | Degraded on this data? | Reported via |
|---|---|---|---|
| MuSiC | the reference carries no cross-donor variance — its gene weighting collapses to a constant | **yes**, the frozen signature is collapsed | `degenerate_`, `degeneracy_reason_` |
| SCDC ENSEMBLE | only one reference is supplied | **yes** | `degenerate_`, `degeneracy_reason_`, `ensemble_weights_` |
| Bisque | (a) no subjects assayed both ways; (b) the reference carries fewer than `min_donors` donor profiles | **yes, both** | `degenerate_`, `degeneracy_reason_` |

The vendored `reference_frozen/` signature is a collapsed matrix, so **MuSiC is
degenerate on it**. Restoring a genuine MuSiC run needs the cell-level single-cell
export — see `reference_frozen/PROVENANCE.json`.

**Two gaps found by test, 2026-09-03.**

*SCDC ENSEMBLE* set `degenerate_ = True` but never set `degeneracy_reason_`. Every
report surfaces degeneracy by printing the *reason*, so the flag alone disclosed
nothing: `scdc_ensemble` would have appeared in the leaderboard as though an ensemble
had happened.

*Bisque* set **neither**, while its own class docstring said "it is running in its
degraded mode, and any conclusion about Bisque's ranking here carries that caveat." The
caveat had no way to travel: every report reads `degenerate_` / `degeneracy_reason_`,
and a docstring reaches no artefact. `implementation_report.json` therefore recorded
`degenerate: false` for a method that is degraded twice over on this data —
a direct contradiction of the invariant *"Bisque without overlapping subjects is in its
degraded mode. Each reports this."*

The two degradations are independent and are now reported separately, because claiming
both unconditionally would mislead as badly as claiming neither:

1. **No-overlap mode.** Unconditional here: Ivy GAP overlaps no public GBM single-cell
   atlas, so the assay transform is estimated from marginal distributions rather than
   from paired subjects. This is Bisque's documented `use.overlap = FALSE` path and is
   weaker than the published method.
2. **The donor-profile fallback.** The marginal the transform matches is specified as
   being taken across reference *donors*. A collapsed signature matrix has none, so the
   spread across *cell types* stands in — a different quantity that happens to have the
   right shape. Fires only when the reference carries fewer than `min_donors` (3) donor
   profiles, which the vendored frozen reference does (it carries zero).

### The disclosure now reaches the headline artefact

`acs_leaderboard.csv` used to carry no degeneracy information at all, so a reader
looking at the table — rather than at a sidecar JSON beside it — saw **ten independent
real methods where there are eight**: MuSiC is NNLS on a collapsed reference, SCDC
ENSEMBLE is SCDC with one reference, and Bisque is running two degraded modes at once.
The leaderboard now carries `degenerate`, `degeneracy_reason`, `implementation` and
`acs_tie_group`, the last naming every method sharing an ACS value so a reader counting
points for a rank correlation counts the right number. `test_implementation_disclosure_names_every_python_reimplementation`
now requires a reason wherever the flag is set.

### The disclosure is written by every run that produces a leaderboard

`implementation_report.json` used to be written only by the benchmark stage — which does
not run without a cell-level reference. On this machine that stage is skipped, so a real
run produced a leaderboard with no record of which methods had fallen back to a Python
reimplementation and which were degenerate. `run_anatomic` now writes the disclosure
itself, into whichever output directory that run owns, containing per method:
`implementation`, `fallback_reason`, `degenerate`, `degeneracy_reason`.

Two invariants depend on this file existing:

- a Python reimplementation is never reported as the published package;
- a method running degenerately says so under its own name.

**Current state on this machine:** R is installed but none of MuSiC, DWLS, BisqueRNA or
SCDC are, so all four published tools run as Python reimplementations. Every artefact
says so. They are not reported as the published packages anywhere.

## The gene space every method receives

Deconvolution runs on the informative subset of the signature, not on every gene the
reference and the bulk share. `config.py` states both reasons and both are load-bearing:
least squares is dominated by the largest residuals, which belong to highly expressed
housekeeping genes that are similar in every cell type and carry no information about
composition; and nu-SVR is superlinear in the number of rows.

`reference.frozen_gene_space()` applies it: top `SIGNATURE_GENES_PER_TYPE` (200) genes
per cell type by a one-vs-rest expression ratio, plus every marker gene, ranked **within
the shared space** so no selected gene is chosen and then discarded for being absent from
the bulk. On this cohort that is **1,614 of 16,007** shared genes. The choice is made
once, for every method, and written to `benchmark/gene_space.json`.

### The defect this replaced, found 2026-09-03

`config.USE_SIGNATURE_GENE_SUBSET = True` was declared, documented as load-bearing, and
**never read**. `select_signature_genes` was called only inside `run_benchmark`, which
does not run without a cell-level single-cell reference — which is exactly the
configuration this project ships in. On the frozen-signature path `run_all.py` passed the
full 16,007-gene intersection to every method.

Measured on this machine, one NuSVR fit:

| gene space | nu=0.25 | nu=0.5 | nu=0.75 | hit the 200,000-iteration cap |
|---|---|---|---|---|
| 16,007 shared | 21.1 s | 37.2 s | 43.6 s | **all three** |
| 1,614 subset | 0.26 s | 0.68 s | 1.28 s | none |

At 366 fits per cohort that is 3.5 hours for the 122-sample run and 7.7 for the
270-sample one, against about 15 minutes for both together.

**The cost was not only time.** Every fit on the full intersection hit the iteration cap,
so libsvm returned its best iterate rather than a converged solution — for every sample,
at every nu. SVR's entire column would have been non-converged estimates, and
`n_hit_iteration_cap_` would have reported 100%. A leaderboard row for "SVR" built from
that is not a result about SVR.

Four tests now guard it, including the negative control that the full intersection *does*
come back when the switch is turned off — otherwise the positive test would pass for
reasons unrelated to the flag.

**One difference from the previous run.** `RESULTS.md` recorded 1,475 genes; this path
selects 1,614. The earlier number came from ranking over the whole frozen profile and
intersecting with the bulk afterwards, which discards top-scoring genes the bulk does not
carry. Ranking inside the shared space keeps every gene it selects. The two runs are
therefore not on an identical gene set and their ACS values are not directly comparable —
which is why the superseded numbers are labelled as superseded rather than compared.

## Running the published packages, rather than reimplementations of them

For most of this project's life, MuSiC, DWLS, Bisque and SCDC were this project's Python
reimplementations wearing the published tools' names. Every artefact said so, but the
leaderboard was not a comparison of the published methods. Making them actually run took
more than installing them, and each obstacle failed in the same direction: silently,
into the fallback.

### The packages

| package | obstacle | resolution |
|---|---|---|
| BisqueRNA | **archived from CRAN**, though `install_deps.R` said "this one IS on CRAN" | authors' GitHub |
| DWLS | on CRAN, but depends on **MAST** (Bioconductor), which `install.packages()` cannot see — the error reads `dependency 'MAST' is not available`, like a network fault | install MAST from Bioconductor first |
| SCDC | needs xbioc, which needs **pkgmaker**, also archived from CRAN — the chain fails at its first link | both from the maintainer's GitHub |
| all | ~half of a first pass lost to transient `cannot download any files` | every step retries |

A single attempt reports a package as unavailable when it is merely unlucky, and the
pipeline then runs a reimplementation under the published tool's name. Retrying is not
politeness; it is the difference between a true and a false statement in the artefacts.

### MuSiC changed its interface

`music_prop()` in MuSiC 1.x takes a bulk **matrix** and a `SingleCellExperiment`; 0.x
took two `ExpressionSet`s. The old call fails with `argument "bulk.mtx" is missing`.
`run_music.R` detects which interface is present rather than pinning a version.

### The cell-level export

All four drivers consume cells, not a signature matrix — `run_music.R` states it
directly: passing pre-averaged profiles "would run the function while removing the thing
that makes it MuSiC". So the export is what the R side actually reads, and it is where
most of the remaining problems lived.

- **Per gene set, not per reference.** Exporting all 16,758 genes x 15,311 cells takes
  ~49 minutes and 1.1 GB gzipped, measured — and R re-reads it for every method on every
  cohort. No method needs those genes: each receives `data.bulk`, already restricted to
  the run's gene space, and MuSiC, Bisque and SCDC all intersect internally. Exports are
  written per gene set and cached.
- **Keyed on the cells as well as the genes.** Two references can share a gene set and
  hold entirely different cells.
- **Metadata always rewritten.** A stale `sc_meta` from a different build surfaces in R
  as `no shared cell ids between counts and metadata`, naming neither the stale file nor
  the cause.
- **Written to `.part` and renamed.** An interrupted run left a truncated 11.8 MB counts
  file where 1.1 GB was expected, and `check()` treats existence as availability — the
  next run would have handed R a silently truncated reference.
- **Cells with zero expression in the gene set are dropped.** BisqueRNA rejects the run
  outright (`Zero expression in selected genes for 16 cells`, from `CountsToCPM`). Those
  cells carry no compositional information in the space being deconvolved. 16 of 11,755
  on the real atlas.

### DWLS was asked to select genes twice

`run_dwls.R` stopped with *"only 31 genes shared between the DWLS signature and the bulk
matrix"* — its own guard, and the diagnosis only became visible once fallback reasons
started carrying R's error line.

The cause is this pipeline, not DWLS. The published workflow selects signature genes
from the whole transcriptome; equal footing requires every method here to receive the
same **pre-filtered** gene space (top-N per cell type plus markers, 676 genes on this
run). DWLS's MAST step then selects again from an input already reduced to the
informative genes, and at the published cutoffs the second pass left 31 genes for 8 cell
types — about four each.

Deconvolving on 31 genes would have made DWLS look bad for a reason that has nothing to
do with DWLS. The driver now relaxes `diff.cutoff` and `pval.cutoff` progressively —
starting at the published defaults — until enough genes survive, printing which
thresholds were used. This is a deviation from the published defaults, made because of
the pre-filtering this pipeline imposes and recorded here; it was decided from the gene
count alone, before any DWLS score existed.

A second deviation is for tractability. `buildSignatureMatrixMAST` fits a hurdle model
per gene across every cell — 758 s on 11,739 cells — and the pipeline calls DWLS in three
stages. The signature does not depend on the bulk, so it is cached beside the export and
keyed on that export's checksum; a different reference or donor split gets its own.

The cost is stubborn. Capping cells to 2,000 still took ~35 minutes, and DWLS's other
signature builder (`buildSignatureMatrixUsingSeurat`) took the same — because the
dominant cost is the condition-number search over gene counts that both builders share,
not the differential-expression step. MAST is kept as the one verified end to end here.
A run therefore spends roughly half an hour per gene set on this one method's
preprocessing, and that is a property of DWLS on this hardware, not something the
pipeline can tune away.

Cells are capped **per cell type** (250), not as a share of the whole. A signature is a
per-cell-type summary, so type proportions are irrelevant to it — what matters is how
many cells estimate each type. A proportional subsample keeps abundant types abundant
and leaves rare ones thin, which is backwards: Astrocyte, at ~1% of cells, would have
been estimated from a handful. Capping per type spends the budget where it buys accuracy,
keeps every cell of any type below the cap, and costs less. The realised counts are
printed into the run log.

A third: a sample that fails to solve becomes NA rather than aborting the method. DWLS
built its signature successfully and then failed on particular mixtures with
`NA/NaN argument`; one such sample was losing the whole cohort to the Python fallback.
NA is already what this pipeline means by a failed sample, and the count is reported.

### Donor leakage, in two places

The R methods do not read the `ReferenceBundle`. They read the export, written from
whatever cell source is registered — so the donor-held-out design has to be enforced
there too, and it was not.

1. **Inside the benchmark.** `run_benchmark` builds its reference from training donors
   only, but the caller had registered the full atlas. MuSiC, Bisque and SCDC would have
   seen the held-out donors' cells — the exact cells the test mixtures are pooled from —
   while NNLS and SVR saw only the training signature. The R methods would have been
   scored against their own reference, in the stage whose entire job is to select a
   method.
2. **In the ACS leaderboard.** `run_benchmark` restores the caller's full source when it
   finishes, but stage 4 deconvolves with the training-only reference. The leaderboard
   would have given the R methods 110 donors and the Python methods an 88-donor
   signature.

Neither is detectable by the equal-footing certificate, which hashes the
`DeconvolutionInput`: the difference lives outside it, in module state the bridge reads.
Both are now enforced explicitly, with tests, and the leakage assertion was seen to fail
against the unguarded code before the guard existed.

### Degradation is a property of the data, not of the implementation

`run_bisque.R` passes `use.overlap = FALSE`, so the genuine BisqueRNA package is degraded
too — but the flag was set inside the Python solver, and when R succeeded that solver
never ran. `implementation_report.json` would have recorded `degenerate: false` for
`R:BisqueRNA`: the same invariant violation, reappearing on the other path, and appearing
precisely when the real package started working.

`degradation_for(data)` now answers it without solving, and `RMethod` sets the flags on
both the R-success and fallback paths. MuSiC's answer flips correctly: degenerate on the
collapsed frozen signature, **not** degenerate on the cell-level atlas.

### A fallback must say why

The reason recorded was `str(exc).split("\n")[0]` — `run_music.R exited 1`, the exit code
and nothing else. Every "this ran as a Python reimplementation" disclosure therefore
omitted the one detail that makes it actionable. R writes `Error in <call> : <message>`,
so that line is now lifted out. The difference is between "music fell back" and "music
fell back because `music_prop()` wanted `bulk.mtx`".

## Reading a real single-cell atlas

`build_from_h5ad` called `X.toarray()` unconditionally. Core GBmap is 338,564 cells x
27,632 genes; dense float32 is tens of TB against 9 GB of RAM. Not a slow path, an
impossible one — and it made the atlas unusable, which in turn kept the four published
tools, method selection and two of three yardsticks blocked.

- `anndata.read_h5ad(backed="r")` does not help. It leaves only `X` on disk and eagerly
  loads everything else, and this file carries a `layers['scaled']` matrix and a `raw/X`
  matrix, each the size of `X`. Reading `obs`, `var` and `X` through `anndata.io`
  primitives skips both.
- **Cells are subsampled donor-balanced within cell type**, seeded, and recorded to
  `reference_sampling_<name>.json`. This follows the predecessor's own documented choice
  — `reference_frozen/PROVENANCE.json` records the frozen signature as "100 donors,
  donor-balanced" — rather than inventing one. Balance matters beyond memory: an
  unbalanced panel lets deeply-sequenced donors dominate every cell-type mean and skews
  MuSiC's cross-donor variance.
- The total is a **target, not a cap**: the trim is proportional within cell type with a
  floor of 25 cells per type, so a rare population cannot be thinned out of existence.
  The overshoot is reported rather than hidden.
- **CELLxGENE indexes `var` by Ensembl id**; this project uses HGNC symbols. Matching
  them directly yields an empty intersection, so symbols come from `feature_name`, and a
  disjoint gene space raises rather than building on whatever survived.

Measured on the real file: 338,564 cells -> 314,700 after roster mapping -> 15,311 kept
across 110 donors and 16,758 genes, in 237 s, with `has_cross_donor_variance` True.

## The negative controls are distributions, not single draws

Both controls are specified as one draw: `control_shuffled_signature` permutes the
signature's gene rows once from `config.RANDOM_SEED`, `control_random` draws one
Dirichlet sample per row from the same seed. A single draw is one sample from a random
variable, and reporting it as *the* control value treats a point estimate as a
parameter.

On real Ivy GAP data that variable is wide. Measured over 200 permutations of the same
signature on the same cohort, the shuffled-signature control has **mean 0.361, sd
0.156**, ranging 0.000–0.708 — against a total spread of 0.29 across every real method
on the leaderboard. The pre-specified draw landed at the **91st percentile**, scoring
0.585.

That single number, read as "the control", ties a real method and beats its own
permutation null, and invites the conclusion that the constraint set is too permissive.
Read as one draw from a distribution centred near chance, it says the opposite. A
superseded run of the same control on a slightly different gene set recorded 0.215 —
also a legitimate draw, and the reason those two runs appeared to disagree.

**So a single-draw control cannot support "the controls behave" in either direction**,
and whichever way it lands, the reading is an artefact of the seed.

`anatomic/control_calibration.py` therefore scores every control over `--control-draws`
independent draws (default 200) as part of each run. What it changes and does not
change:

- The **constraint file is untouched.**
- The **controls' definitions are untouched**, and the **leaderboard rows are untouched** —
  the pre-specified single draw is still reported, because that is what was registered.
- What is added is each control's own sampling distribution, so the leaderboard row can
  be located in it and so the question a reader actually has — *is this method
  distinguishable from a meaningless composition?* — is answered against a percentile.

Methods are judged against the **hardest** control, not whichever happens to be
convenient. Both controls are calibrated even though `control_random` never came near a
real method on this cohort: "it did not matter this time" is not a reason to leave one
of two controls uncalibrated while arguing that single-draw controls are the defect.

Note the direction of travel. Replacing a point estimate with a distribution makes the
control **harder** to beat and the test more conservative, never less. That is the only
direction it is safe to move an instrument after seeing what it scored, and it is why
this is not the retuning the protocol forbids.

### The verdict wording had the same problem

`control_verdict` compared the best control against the **median** real method. A
control could therefore tie the *worst* real method **and** beat its own permutation
null while the run announced "CONTROLS BEHAVE" — which is exactly what happened on the
first full-database run. There are now three outcomes, and the reassuring one requires
all three facts to hold: below the median, ties no real method, and does not beat its
own null.

## The ACS scoring set is enforced, not remembered

`run_anatomic.run(deconvolve_all=True)` solves every sample in the archive — all 270 —
while still scoring ACS on the 122 H&E-selected anatomic samples alone. The separation
matters because the two sets answer different questions:

- **Composition** is well defined on all 270 samples, and the wider set is what takes the
  prognostic cohort from 7 tumours to 29.
- **ACS** may only ever be scored on the anatomic subset. Ivy GAP's other 148 samples had
  their structure assigned by in-situ hybridization — from expression — so an anatomic
  constraint scored on them is circular by construction.

The function raises if any sample lacking `is_anatomic_study` reaches the scoring set.
That is a guard, not a comment: the failure it prevents produces a *higher* ACS on a
larger apparent cohort, which is exactly the kind of error nobody notices.

There is a second, less obvious reason to run both. Several methods use cross-sample
statistics — Bisque's per-gene transformation, SCDC's weighting — so the same sample does
not receive the same estimate depending on what it was solved alongside. Comparing the
two runs measures how much of the leaderboard is a property of the method rather than of
the cohort it was handed; `acs_cohort_sensitivity.csv` reports it per method with the
rank change.

The two runs write to separate directories and cannot overwrite each other. This is the
same isolation rule the synthetic/real split already enforces, applied one level down.

## Per-tumour ACS

The protocol asks for per-tumour results explicitly, and until now only the pooled score
and the per-constraint breakdown were reported. `acs_per_tumor.csv` gives, for every
method and every tumour, which constraints that tumour could evaluate and which it
satisfied.

It is a decomposition of the pooled score, not a second statistic:
`test_per_tumor_table_reweights_exactly_to_the_pooled_acs` requires the weighted sum of
the parts to equal the whole to 1e-12. With nine evaluable tumours a pooled ACS can be
carried by one or two of them, and the pooled number alone cannot show that.

### The ceiling is a property of the instrument

Maximum evaluable tumours per constraint, from Ivy GAP's H&E study alone — before any
method runs:

| C1 | C2 | C3 | C4 | C5 | C6 | C7 |
|---|---|---|---|---|---|---|
| 8 | 8 | 9 | 9 | **6** | 9 | 8 |

Maximum weighted denominator: **65.0**, across 9 distinct tumours. One tumour supplied
only CT and PAN; two supplied no PAN. No deconvolution method, hyperparameter, gene set
or reimplementation can raise any of these numbers, and no wider cohort can either —
Ivy GAP has 10 H&E-annotated tumours and that is all it has.

This bounds what "improving ACS" can mean. A method can satisfy more constraints; the
evidence base it is judged on cannot grow.

## Uncertainty calibration — a measured limitation

The Bayesian methods report 95% credible intervals. An interval is a claim, and this
project checks the claim rather than trusting it. `ivygap/bench/calibration.py` measures
realised coverage against known-truth pseudobulk on every run and writes
`uncertainty_calibration.csv` into the release bundle.

**It caught a real defect during development.** The sampler's original default rescaled
each FPKM profile to 100,000 pseudo-counts, treating it as that many independent
multinomial reads. FPKM is nothing of the sort — it is already depth- and
length-normalised, and neighbouring genes are strongly correlated — so the likelihood
was far too informative and the intervals came out about five times too narrow:

| pseudocount | realised coverage | MAE | mean CI width |
|---|---|---|---|
| 1e3 | 0.86 | 0.0083 | 0.038 |
| 1e4 | 0.60 | 0.0064 | 0.013 |
| 1e5 | **0.23** | 0.0063 | 0.004 |
| 1e6 | 0.08 | 0.0063 | 0.001 |

Point accuracy is flat above 1e4, so the large values bought nothing but overconfidence.
The default is now 1e3.

**A limitation that remains, and cannot be fixed by widening intervals.** Scored against
mixtures pooled from *held-out donors*, the dominant error is not scatter but systematic
per-type bias — on the test fixture about −0.13 for Tumor and +0.15 for
Macrophage_Microglia, against a mean interval half-width of 0.024. The posterior is over
composition *given the reference*; it cannot represent the possibility that the reference
is wrong for this patient. No posterior conditional on the reference can cover a bias in
the reference, and tuning the intervals until it did would be fitting the interval to the
answer.

`propagate_reference_uncertainty` (default on) draws the signature matrix itself each
sweep from its cross-donor standard error, which admits the *estimable* part of that
uncertainty. The effect is real but modest — averaging over hundreds of genes cancels
most per-gene reference noise.

**How to read these intervals.** Comparatively, not absolutely. A wide interval reliably
means "this sample is harder", which is genuinely useful for Ivy GAP's small
leading-edge dissections. The absolute width understates total error whenever the
reference does not match the patient, and `uncertainty_calibration.csv` states by how
much on the run in front of you.

## What the benchmark measures

Per method, per cell type, on donor-held-out mixtures with known composition:

- **MAE / RMSE / bias** — absolute accuracy and systematic over/under-estimation.
- **Pearson / Spearman** — reported *alongside* error, never instead. A method reporting
  every macrophage fraction as exactly twice the truth correlates at r = 1.0 and is
  wrong about every sample.
- **false_positive_fraction** — where a type is *truly absent*, the mean fraction the
  method invented. Ivy GAP's leading edge really does lack lymphocytes; a method that
  cannot represent absence manufactures immune infiltrate at the tumour margin.
- **limit_of_detection** — the smallest true fraction still recovered above noise. This
  is the number that decides whether a 1–2% T-cell population is measurable at all.
- **per-niche breakdown** — the same metrics restricted to mixtures shaped after each
  anatomic structure. A method can win overall and still fail exactly where this project
  operates.

The headline selection criterion is `mae_primary`, which **excludes the Astrocyte
sidecar** so an acknowledged-unresolvable column cannot move the ranking.
