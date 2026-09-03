> **SUPERSEDED — preserved verbatim, do not cite.**
>
> This is `RESULTS.md` as it stood on 2026-09-03 before the re-run. Its ACS numbers were
> computed on a **1,475-gene** space; the current run uses **1,614** genes, selected by
> the same rule applied inside the shared space rather than before intersecting with it.
> Different gene space, different estimates — the two are not comparable and are not
> compared.
>
> It is kept because it is the only surviving record of that run: the artefacts it
> describes were destroyed by the synthetic-fixture overwrite documented in
> `SUPERSEDED.md`, and these hand-transcribed numbers were what remained. That is also
> why `scripts/summarize_results.py` now exists — so no future write-up depends on
> anyone copying a number correctly.
>
> Replaced by `RESULTS.md` at the repository root.

# Results — first real run

**Run date:** 2026-09-03 · **Cohort:** Ivy GAP anatomic study, 122 samples / 10 tumors
**Constraint freeze hash:** `2d1fb47c98832adfae20e5b79a97b731ac3cced25fa14c1dd7bb02da7895807a`

Read [Anatomy_Test.md](Anatomy_Test.md) first. This file reports what the pipeline
produced and, more importantly, what it does not yet support.

---

## 1. What ran, and what did not

| Protocol step | Status |
|---|---|
| 1. Restore frozen inputs | **done** — `reference_frozen/`, hashes recorded |
| 2. Write + timestamp the constraint file | **done in repo** — `constraints.py`, hash above. *Not yet publicly registered (OSF).* |
| 3. Ingest Ivy GAP | **done** — 122/10 anatomic, 148/34 ISH-cluster, verified against the Allen README |
| 4. Build + prove the scorer | **done** — known-answer fixtures pass; permutation null centred at chance |
| 5. Score every method + controls | **done** — 10 methods, 2 controls, 10,000 permutations each |
| 6. The agreement test | **partially** — all three protocol yardsticks UNAVAILABLE; run on a simulated substitute |
| 7. Prove it transfers | not started |

---

## 2. The ACS leaderboard

122 samples, 10 tumors, 1,475 signature genes, 10,000 within-tumor label permutations
per method. 9 tumors contribute at least one evaluable constraint (W26-1-1 supplied only
CT and PAN, so no constraint is evaluable for it).

| method | ACS | 95% CI | null mean | null p |
|---|---|---|---|---|
| scdc / scdc_ensemble | **0.877** | 0.761 – 1.000 | 0.373 | <0.001 |
| nnls / music | 0.815 | 0.746 – 0.875 | 0.379 | <0.001 |
| elastic_net | 0.769 | 0.672 – 0.855 | 0.380 | <0.001 |
| bisque | 0.769 | 0.646 – 0.875 | 0.382 | <0.001 |
| svr | 0.754 | 0.672 – 0.847 | 0.367 | <0.001 |
| bayesian | 0.662 | 0.561 – 0.779 | 0.382 | <0.001 |
| dwls | 0.631 | 0.535 – 0.755 | 0.369 | <0.001 |
| bayesian_hierarchical | 0.631 | 0.545 – 0.734 | 0.382 | 0.001 |
| **control_random** | **0.400** | 0.288 – 0.530 | 0.384 | **0.437** |
| **control_shuffled_signature** | **0.215** | 0.157 – 0.294 | 0.372 | **0.996** |

### The controls behave

This is the result that makes the rest readable. Both controls sit below every real
method, and **neither beats its own permutation null** (p = 0.437 and p = 0.996). The
constraint set is therefore discriminating, not permissive — the outcome the protocol
says would otherwise force reporting it as too easy and publishing the constraint file
anyway.

`control_shuffled_signature` scoring *below* the null (0.215 vs 0.372) is expected: a
permuted signature does not produce random orderings, it produces systematically wrong
ones.

### Two degeneracies confirmed exactly

`nnls` and `music` returned **identical** scores, and so did `scdc` and
`scdc_ensemble`. Both were predicted in advance by the degeneracy flags: MuSiC's entire
contribution is weighting genes by cross-donor variance, and the vendored signature is a
collapsed matrix that carries none, so MuSiC *is* NNLS here; SCDC ENSEMBLE with a single
reference *is* SCDC. Neither is presented as a distinct method.

### Per constraint (best method, scdc)

| ID | claim | tumors evaluable | satisfied |
|---|---|---|---|
| C1 | Tumor: CT > LE | 8 | 8 (1.00) |
| C2 | Oligodendrocyte: LE > CT | 8 | 8 (1.00) |
| C3 | Endothelial: MVP > CT | 9 | 9 (1.00) |
| C4 | Endothelial: MVP is the maximum | 9 | 9 (1.00) |
| C5 | Macrophage: PAN > LE | 6 | 6 (1.00) |
| C6 | Macrophage: MVP > CT | 9 | 7 (0.78) |
| C7 | Tumor: LE < IT < CT (weight 2) | 8 | 5 (0.63) |

C7 is hardest, as designed — a three-step monotone chain is a conjunction and much
harder to satisfy by chance than any single step. Note the evaluable counts: **no
constraint rests on more than 9 tumors, and C5 on only 6.**

---

## 3. The agreement test

**The protocol's three yardsticks are all UNAVAILABLE**, and this is a finding about the
available data, not a failure of the run:

| yardstick | why |
|---|---|
| `synthetic_mixtures` | **real and vendored**, but scored only NNLS and SVR — 2 methods, Spearman needs ≥6. See the pairwise comparison below. |
| `absolute_purity` | the source project records this as **BLOCKED — no ABSOLUTE purity file was ever located** |
| `sc_pseudobulk` | requires paired single-cell data, not present |

So the study's primary result **is not yet computable**, and `agreement_report.json`
says so in `primary_result_status`.

### Preliminary reading on a simulated substitute

`simulated_donor_mismatch` generates 250 mixtures from donor-perturbed copies of the
frozen signature while methods solve against the unperturbed one, so the dominant error
is reference mismatch rather than sampling noise. It is **not** the protocol's yardstick
1 and every artefact says so.

> **ρ = 0.494**, bootstrap CI **[−0.293, 1.000]**, 10 methods / **8 distinct**
> → **NULL RESULT** — below the pre-registered 0.60, and the CI includes zero.

| method | ACS rank | truth rank | gap |
|---|---|---|---|
| scdc / scdc_ensemble | 1 | 1 | 0 |
| bayesian | 8 | 8 | 0 |
| bayesian_hierarchical | 9 | 9 | 0 |
| nnls / music | 3 | 5 | 2 |
| elastic_net | 5 | 7 | 2 |
| svr | 7 | 4 | 3 |
| **bisque** | 5 | **10** | **5** |
| **dwls** | 9 | **3** | **6** |

**DWLS is the case that carries the finding.** It is essentially tied for most accurate
(MAE 0.0072 against SCDC's 0.0070) and simultaneously tied for *worst* on anatomic
concordance (0.631). A method can be numerically excellent and still fail the constraints
a pathologist would predict. Bisque is the mirror image: middling ACS, worst accuracy.

If this held on a real yardstick it would be the protocol's second branch — *"it matches
known biology" is not evidence of correctness* — which the protocol calls "arguably more
useful". **It does not hold yet.** See the caveats before treating it as anything.

### One comparison against REAL ground truth

The protocol's yardstick 1 does exist — 500 held-out pseudobulk mixtures with known
composition, vendored in `reference_frozen/tcga_benchmark/`. It covers only **NNLS and
SVR**, so it cannot support a rank correlation. But a *pairwise* comparison on two
methods is still meaningful, and it is reported rather than discarded:

| method | ACS (real Ivy GAP tissue) | MAE (500 real mixtures) |
|---|---|---|
| nnls | **0.815** | 0.1039 |
| svr | 0.754 | **0.0970** |

> **ACS prefers NNLS. Real ground truth prefers SVR. They disagree.**

This is one comparison on two methods — it is not a result. It is notable for two
reasons: it uses genuine held-out mixtures rather than a simulation, and it points the
same way as the simulated yardstick, which also ranks SVR (MAE 0.0083) well above NNLS
(0.0380). Two independent accuracy measures agree with each other and both disagree with
anatomic concordance.

If a full yardstick reproduced this over ≥6 methods, it would be the protocol's second
branch. On two methods it is a reason to prioritise obtaining that yardstick, nothing
more.

### Caveats that must travel with that number

1. **The yardstick is simulated.** Mixtures come from perturbed copies of the same
   signature the methods deconvolve with. Real Ivy GAP tissue has a far larger reference
   mismatch (a modern droplet scRNA-seq reference against 2014 laser-capture bulk), so
   the true accuracy ranking may differ.
2. **8 distinct points.** The CI spans −0.29 to 1.00. This design cannot separate
   ρ = 0.49 from ρ = 0.0 or from ρ = 0.9.
3. **ACS and accuracy were measured on different data** — real tissue and simulated
   mixtures. That is the intended design (ACS exists for tissue without ground truth),
   but it means a method could rank differently for reasons of cohort, not method.
4. **9 tumors, and C5 on 6.** The ACS side is itself thinly supported.

---

## 4. A negative result on the project's opening premise

Elastic Net (0.769) and both Bayesian variants (0.662, 0.631) did **not** beat plain
NNLS (0.815) on ACS, and on the simulated yardstick Elastic Net (MAE 0.0437) was worse
than NNLS (0.0380). The regularisation and probabilistic modelling that motivated this
rebuild did not improve either measure here.

That is reported because it is what happened. It is one cohort, one collapsed reference,
and a simulated accuracy yardstick — not a general claim about the methods.

---

## 5. Defects found and fixed during this run

Recorded because each would have produced plausible-looking wrong numbers — or, in the
last case, silently lost correct ones.

| defect | consequence had it shipped |
|---|---|
| ISH-cluster samples matched as anatomic (`CT-CD44` starts with `CT`) | would have tripled the cohort with samples whose structure labels came *from expression* — the exact circularity the study exists to eliminate |
| observed ACS and its null computed by different estimators | every p-value would compare two different statistics the moment a method emitted NaN |
| vectorised null propagated NaN where pandas skipped it | null mean 0.324 vs 0.386 on affected data |
| simulated yardstick divided columns by their totals | ~1.33× systematic error in *every* method's accuracy; leaderboard would still have looked plausible |
| `tumor_details.csv` served as HTML with HTTP 200 | a `<style>` column parsed as clinical data |
| unbounded `NuSVR` iteration | pipeline hang with no diagnostic on the real collinear signature |
| `NICHE_TEMPLATES` kept pre-rename `CTmvp`/`CTpan` | per-niche tables labelled with structures that exist nowhere else |
| degenerate duplicates counted as independent points | Spearman judged on 10 points that were really 8 |
| benchmark declared a winner on a 0.0006 MAE gap with 2 donor groups | a coin flip reported as a method selection |
| **synthetic and real runs shared output paths** | **a fixture run silently destroyed a completed 30-minute real run.** The manifest still said SYNTHETIC FIXTURE, so nothing was mislabelled — but labelling is not isolation, and the real numbers survived only because they had already been transcribed here by hand. Synthetic runs now write to `results_synthetic/` and `release_synthetic/`. |

---

## 6. What would move this forward, in order

1. **Publicly register the constraint file** (OSF). The hash is above; the timestamp must
   precede any result anyone cites. This is cheap and it is what makes the design
   credible to a reviewer.
2. **Obtain a real yardstick covering ≥6 methods.** A GBM single-cell atlas at
   `data/reference/gbmap_core.h5ad` enables both `synthetic_mixtures` (real
   donor-held-out mixtures) and `sc_pseudobulk`, and turns the primary result from
   "not computable" into an answer.
3. **Locate ABSOLUTE purity for TCGA-GBM.** The source project never found it. It is the
   one yardstick that never touches RNA.
4. **Widen the ACS cohort.** 9 evaluable tumors is the binding limit on every CI here.
