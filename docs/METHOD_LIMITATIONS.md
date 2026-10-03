# What this study could not measure, and why

Every limit below is a fact about **this machine or this method's coverage**, not a finding
about deconvolution. They are collected here so a reader does not have to infer them from a
leaderboard, and so that none of them is quietly reported as a result.

The distinction that runs through the whole page: **"the method produced a bad estimate" and
"the method could not be evaluated" are different statements.** Collapsing them would make a
method look worse than the evidence supports.

---

## 1. quanTIseq is INCONCLUSIVE on the anatomic arm — not bad

quanTIseq appears on the leaderboard at **ACS 0.600**, at the bottom of the real
methods. That number must not be read as a rank against the others.

| | quanTIseq | every other real method |
|---|---|---|
| constraint-tumour pairs scored | **15** | 57 |
| its own permutation null (mean) | 0.502 | ~0.375 |
| p against that null | **0.309** | < 1e-4 |
| beats its own null | **no** | yes |
| `comparable` flag in the leaderboard | **false** | true |

**Why the coverage is low.** quanTIseq ships its own TIL10 signature and ignores the
reference it is handed. TIL10 covers four immune populations; most of this study's anatomic
constraints name cell types quanTIseq does not estimate at all, so only
15 of 57 constraint-tumour pairs can be scored.

**Why that raises its null rather than lowering its score.** A smaller constraint set is
easier to satisfy by chance, so quanTIseq's null sits at 0.502 where everyone
else's sits near 0.375. With 15 pairs there is not enough
power to clear it. The honest verdict is **INCONCLUSIVE — underpowered**, which the protocol
treats as an answer, not a failure.

**It has no purity estimate either.** quanTIseq returns no usable tumour-content estimate in
**either** cohort, so it is absent from the accuracy yardstick in GBM (n=154) and LGG (n=510)
alike. It is one of 3 methods evaluated on ACS but not on purity.

**What this does to the headline.** The control-separation claim was stated as *worst real
method 0.600 vs best control 0.400*, a margin of +0.200 — and quanTIseq is that
"worst real method". Excluding it as not-evaluable, the worst **fully-scored** method is
`bisque` at **0.708**, so the real margin is **+0.308**.
Reporting quanTIseq as the floor understated the separation. Both figures are now given
wherever the claim appears.

---

## 2. Two methods lost their genuine R package on this laptop

This is the one place where **this machine's speed changed a published number.** Both
packages run correctly here and both produce a *higher* ACS than the Python reimplementation
the leaderboard actually carries. For DWLS the cause is a plain timeout. For BayesPrism it is
not, and that is stated rather than smoothed over.

> ⚠ SUPERSEDED (OPEN_DEFECTS D22, 2026-09-30): the genuine-package re-measurements quoted below ran on GBmap's log `X` layer while the leaderboard ran on `raw/X`, so these figures are cross-matrix and are withdrawn. Current values: MANUSCRIPT §4.7 and `results/*_remeasured.json` (DWLS: no estimate for 73 of 122 samples on raw/X; Bisque and EPIC reproduce their leaderboard rows).

| method | budget | genuine package took | outcome | ACS reimplementation | ACS genuine |
|---|---|---|---|---|---|
| **DWLS** | 2400 s | **2593 s** (+8%) | fell back | 0.7231 | **0.7846** |
| **BayesPrism** | 2400 s | 2045 s with 3 workers; **4076 s** at `n.cores=1` | fell back — *cause not a clean timeout, see below* | 0.8000 | **0.8154** |

> ⚠ SUPERSEDED (OPEN_DEFECTS D22, 2026-09-30): the genuine-package re-measurements quoted below ran on GBmap's log `X` layer while the leaderboard ran on `raw/X`, so these figures are cross-matrix and are withdrawn. Current values: MANUSCRIPT §4.7 and `results/*_remeasured.json` (DWLS: no estimate for 73 of 122 samples on raw/X; Bisque and EPIC reproduce their leaderboard rows).

DWLS misses its budget by 193 seconds —
about 8%. **On a faster
machine the leaderboard would carry 0.7846, not 0.7231.** That is
a hardware artefact, and it is disclosed rather than absorbed.

Both were re-measured separately under a 14400 s budget, on the
**leaderboard's own gene space** (657 genes, hash-matched) with
training donors, held-out donors and the reference composition all verified identical
(`input_equivalence` in each artefact). An earlier version of these re-measurements used a
different gene space and their conclusions were withdrawn — see OPEN_DEFECTS **D10**, now
closed.

**BayesPrism's fallback is not explained by the budget.** Run on its own it finishes in
2,045 s, comfortably *inside* the 2,400 s budget, yet it fell back in the confirmatory run.
The artefact records that it fell back and not why, so no cause is claimed here. The likely
neighbour is OPEN_DEFECTS **D18** — BayesPrism's R socket cluster hangs on this machine
without the timeout firing — but that is an association, not a demonstrated mechanism, and
D18's own history is a record of three confident diagnoses that were wrong. At `n.cores=1`
it takes 4,076 s and *would* exceed the budget honestly.

> ⚠ SUPERSEDED (OPEN_DEFECTS D22, 2026-09-30): the genuine-package re-measurements quoted below ran on GBmap's log `X` layer while the leaderboard ran on `raw/X`, so these figures are cross-matrix and are withdrawn. Current values: MANUSCRIPT §4.7 and `results/*_remeasured.json` (DWLS: no estimate for 73 of 122 samples on raw/X; Bisque and EPIC reproduce their leaderboard rows).

Its score is **unchanged by core count** (0.8154 at three workers and at `n.cores=1`), so the
declared `n.cores` deviation in `docs/METHODS.md` is result-neutral. Only the runtime moves,
and it roughly doubles.

---

## 3. Two methods could not be run at all — toolchain, not CPU

Worth separating, because it would be wrong to blame the processor for either.

| method | blocker | nature |
|---|---|---|
| **CDSeq** | `ld: library 'emutls_w' not found` — R is configured for `/opt/gfortran/bin/gfortran`, which is not installed | missing **Fortran toolchain**; one installer fixes it |
| **Scaden** | `tensorflow>=2.0` publishes no wheel for Python 3.14.5 | missing **dependency build**; unrelated to hardware |

CDSeq was later unblocked and **did run** on the anatomic arm: ACS **0.6829** reference-free
(no atlas at all, p = 0.0009) and **0.7561** with the atlas used only for labelling. It cuts
both ways — the anatomy is recoverable without an atlas, and the atlas still earns its place.

Scaden is **recommended for dropping** rather than chased with a second Python environment;
the reasoning is in `docs/ROAD_TO_PAPER.md`.

> **Scaden, re-verified 2026-09-24.** The blocker is unchanged and was re-tested rather than
> re-quoted. The source is present at `~/Downloads/scaden-master`; its `setup.py` declares
> `install_requires=[... "tensorflow>=2.0" ...]`; the project runs **Python 3.14.5**; and
> `pip install --dry-run --no-deps "tensorflow>=2.0"` returns **`ERROR: No matching
> distribution found`** today. TensorFlow publishes no wheel for this interpreter.
>
> **The practical blocker is the weaker of the two reasons, and it is the one that can
> expire.** If TensorFlow ships a 3.14 wheel tomorrow, Scaden becomes installable and the
> substantive objection still stands: Scaden is *trained on simulated bulk mixtures*, which is
> precisely the construction this study already discloses as a limitation of its own accuracy
> arm. Adding a method that shares the simulator's assumptions to a benchmark already
> criticised for leaning on simulated mixtures compounds the weakness instead of testing
> against it. CDSeq was pursued instead because, being reference-free, it is independent of
> both the simulator and the signature matrix — and it ran.
>
> So Scaden is a **deliberate exclusion with a stated rationale**, not a method that failed.
> The paper should say which of the two reasons it is resting on.

---

## 4. Three methods are on ACS but not on the accuracy yardstick

`bayesian_hierarchical, cibersortx_smode, quantiseq` return no usable tumour-content estimate, in **both** cohorts —
GBM (n=154) and LGG (n=510). So the ACS-vs-accuracy correlation is computed over 12 methods,
not 15. This is coverage, not failure, and it is why the agreement test reports `n_methods`
beside every rho.

- `cibersortx_smode` needs the single cells the signature was built from, and no cell source
  is registered on the frozen path.
- `bayesian_hierarchical` requires a `patient_id` the frozen arms predate (OPEN_DEFECTS D20).
- `quantiseq` — section 1.

---

## 5. The end-to-end pipeline check does not complete on this machine

`python scripts/run_all.py --synthetic` is named in `CLAUDE.md` as the pre-flight validation.
It hangs in BayesPrism at **11x its own 2,400 s budget** — reproduced 2026-09-23 at 7 h 25 m,
writing nothing. The cause is **unknown** and is deliberately not guessed at; the workers were
sampled directly and are computing rather than blocked on a dead peer, which rules out the
leading theory. See OPEN_DEFECTS **D18**.

Consequence: the validation that actually runs is `pytest tests/` (298 tests, ~4-6 min),
`scripts/independent_verification.py` and `scripts/check_doc_numbers.py --against-current`.

---

## 6. The post-registration extension panel (added 2026-09-30 / 10-01)

Five top-tier methods from Nguyen et al. 2024 and Li et al. 2026 were pursued after
registration. Each runs as the **genuine published package**, with no fallback, in a labelled
panel that enters no registered statistic (`ivygap/deconv/extension.py`). What it took, method
by method — every workaround is a disclosed compatibility shim, never a reimplementation, and
each package passed a planted-truth gate before touching real data:

| method | route | what it needed | planted-truth gate |
|---|---|---|---|
| **FARDEEP** | CRAN 1.0.1 | nothing; `permn = 0` (p-values only — estimates verified identical); blocked over samples for 2 cores (verified identical) | — (package defaults; parallel = joint, max diff 0) |
| **LinDeconSeq** | GitHub, pinned commit `20f1aec` | DESeq2 from Bioconductor; `deconSeq()` with this project's signature, since glioma has no pure bulk samples for `findMarkers()` | max err 0.023, r ≥ 0.998 |
| **ARIC** | PyPI 1.0.1 | nothing | max err 0.013 |
| **RNA-Sieve** | PyPI 0.1.4, installed without its stale pins | **shim:** flatten `x0` for SciPy ≥ 1.11 (what older SciPy did silently); **serial path** instead of 10 spawned workers (8 GB, 2 cores; same per-gene problems); exact counts recovered from per-cell CPM × library size | max err 0.011, r ≥ 0.999 |
| **ReCIDE** | GitHub, vendored at commit `31bbd6b` (MIT) | **shim:** `summarize_each` (defunct in dplyr ≥ 1.2; its pinned dplyr 1.1.x no longer compiles on R 4.6); sparse input so Seurat 5 builds the v3 assays its code reads; **> 500 cells per donor**, so the capped reference (≤ 400/donor) cannot feed it — per its README, a donor subsample with all their cells | *in progress* |

**Still blocked, for stated reasons:**

| method | reason | kind |
|---|---|---|
| Scaden, DAISM-DNN | no TensorFlow build for Python 3.14; both train on simulated bulk (§3) | toolchain + design |
| CDSeq | a multinomial read-count model; the Ivy GAP bulk is FPKM | input-type mismatch |
| AutoGeneS | pins scipy 1.3 / matplotlib 3.0; primarily a GA gene-selection method needing the full gene space | toolchain + scope |
| MIXTURE | published as an application repository of scripts (last updated 2020), not a package | packaging |
| Linseed | reference-free: components are unlabelled, so placing it on a labelled roster needs a component-to-type mapping that would itself be a major deviation | design |

**One reading note.** ReCIDE is the top method in Li et al. 2026 [4], which comes from ReCIDE's
own laboratory. An independent evaluation of it is worth more than the original benchmark,
and should be framed that way.

## What none of this changes

The study's conclusions do not rest on any method whose evaluation was limited here. The
lymphoid failure replicates in both cohorts against orthogonal methylation truth; the control
separation survives on fully-scored methods with a **wider** margin than reported; and the
ACS-vs-accuracy null holds on the 12 methods that carry both scores.
