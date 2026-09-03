# Project instructions

## Mission

Test whether a tumor's own anatomy can substitute for ground truth when choosing a
cell-type deconvolution method — and whether a method that reproduces known biology is
actually more accurate. See `Anatomy_Test.md`; it is the protocol, not a summary.

The predecessor TCGA-GBM project lives under `~/Downloads/cancer_judging_machine-master 22/`
and is **read-only from here**. It supplies the frozen signature matrix vendored in
`reference_frozen/` and, later, the ABSOLUTE-purity yardstick. Do not modify it.

## Sources of truth

1. `Anatomy_Test.md` — the protocol.
2. `ivygap/anatomic/constraints.py` — the pre-registered constraint file. Hashed.
3. `docs/DATA_SOURCES.md` — what the archive actually contains, verified against the
   2014-11-25 release.
4. `docs/METHODS.md` — the methods and every deviation from published algorithms.
5. `ivygap/config.py` — the control panel. Never type a path or cell-type name as a raw
   string anywhere else.

## Scientific invariants

Breaking one invalidates results rather than merely making them untidy.

- **Never use an outcome to select anything.** Not the method, not a hyperparameter, not
  a gene set. `run_survival.assert_selection_frozen` aborts if the recorded criterion
  mentions an outcome.
- **Never let ACS select the method.** ACS is the object under test. Selecting on it
  before the agreement test has answered whether that is safe assumes the conclusion.
- **Never edit `constraints.py` because a method or a control scored badly.** If controls
  score high, the protocol's instruction is explicit: report it, publish the constraint
  file anyway, state what a harder set would need. Do not retune and re-report.
- **Never score the ISH cluster samples as anatomic.** Their structure labels were
  assigned using expression; scoring anatomic constraints on them is circular.
  `is_anatomic_sample` requires the `-reference-histology` suffix.
- **Never aggregate over samples where tumors are nested.** Collapse to (tumor,
  structure) first. Block counts run 6 to 15.
- **Never use a 50% coin-flip null.** Cell fractions are compositional and correlated;
  the null is within-tumor permutation of structure labels.
- **Never report a degenerate method under its own name.** MuSiC without cross-donor
  variance is NNLS; SCDC ENSEMBLE with one reference is SCDC; Bisque without overlapping
  subjects is in its degraded mode. Each reports this.
- **Never report a Python reimplementation as the published package.**
- **Never impute a missing input.** Missing clinical data means BLOCKED.
- **Never validate counts against hard-coded literals.** Measure and record.
- **Never let a definitional claim rest on the Astrocyte sidecar**, which is known to be
  unidentifiable against Tumor.

## Evidence required before calling something done

An exit code is not evidence. As applicable:

- File presence, schema, dimensions, hashes, provenance.
- Sample, patient and event counts, plus the reason for every exclusion.
- Non-negativity, row sums, finite values, category mappings, missingness.
- Fold coverage, and that no fold split a patient.
- Negative controls, not just positive ones. `test_inverted_signal_is_rejected` matters
  more than `test_planted_signal_is_recovered`: a scorer that reports high agreement on
  both is measuring nothing.
- No failed folds, no silently dropped samples, no method quietly missing from a ranking.

## Working practices

- Run `pytest tests/ -q` before and after any change to `deconv/`, `bench/` or
  `anatomic/`. It is ~40 seconds.
- Validate end-to-end with `python scripts/run_all.py --synthetic` before touching the
  real-data path.
- Add new work under explicit versioned paths. Label superseded work as superseded
  through additive documentation rather than deleting it.
- Do not commit, push, open a pull request, publish or deploy unless explicitly asked.
- Never place matrices, model binaries or complete large logs into context. Use `rg` and
  focused file ranges.

## When something fails

- Never retry the same failure blindly. After three occurrences of one failure
  signature, stop and form a new evidence-based hypothesis.
- A blocker report must contain: the failure, the evidence, what was attempted, the
  preserved state, the smallest required human action, and the exact resume command.
- Report `PASSED`, `IN_PROGRESS` or `BLOCKED` truthfully. An underpowered result is
  reported as INCONCLUSIVE, which is an answer, not a failure.
