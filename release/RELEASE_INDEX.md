# Release bundle — Ivy GAP GBM anatomic deconvolution

Built 2026-09-04 · run type: **real data**

Selected method: **music**
Selected on: mean absolute error on primary cell types, donor-held-out pseudobulk, aggregated donor-equally

## Files

- `acs_leaderboard.csv` — THE LEADERBOARD. Anatomic Concordance Score per method, with a bootstrap CI over tumours and a p-value against the within-tumour label-permutation null. The two rows flagged is_control are negative controls and are what make the rest interpretable - read the control verdict in anatomic_report.json first.
- `agreement_test.csv` — THE PRIMARY RESULT. Spearman rank correlation between the ACS ranking and the ranking from real ground truth, per yardstick, against a threshold fixed before it was computed. A yardstick marked UNAVAILABLE produced no scores in this run - that is not the same as disagreeing.
- `agreement_report.json` — The primary result with its pre-registered threshold and the requirement it was judged against (rho >= 0.60 AND a bootstrap CI excluding zero).
- `agreement_method_ranks.csv` — Each method's ACS rank beside its ground-truth rank, with the gap. This is where the agreement test is actually legible: a large gap means a method the anatomy likes and the truth does not, or vice versa.
- `acs_per_constraint.csv` — Every constraint x every method: how many tumours could evaluate it and how many satisfied it. A constraint no method satisfies is evidence about the constraint, or about the reference - read it before reading the leaderboard.
- `constraint_file.json` — The frozen, pre-registered constraint file: C1-C7 with directions, weights, histological citations, the declared T-cell exclusion, and the facts this signature is structurally unable to test. Its SHA-256 is the registration hash.
- `anatomic_report.json` — Machine-readable summary, including the control verdict and the note that ACS was never used to select a method.
- `composition.csv` — Composition per (tumour, anatomic structure) from the best-scoring real method. One row per tumour per structure - already collapsed, so uneven block counts cannot skew anything downstream.
- `benchmark_summary.csv` — Per-method accuracy on donor-held-out pseudobulk with known composition. This is the ground-truth ranking the agreement test correlates ACS against, and the only evidence any method is selected on.
- `benchmark_per_type.csv` — Per method, per cell type: MAE, RMSE, bias, correlation, the fraction invented where the type is truly absent, and the limit of detection.
- `benchmark_per_niche.csv` — The same accuracy broken down by which anatomic niche each synthetic mixture was shaped after. A method can win overall and still fail at the leading edge.
- `uncertainty_calibration.csv` — For every method reporting credible intervals: how often the interval actually contained the truth, against a nominal 0.95. Read before quoting any interval.
- `method_selection_decision.json` — Which method was selected, on what criterion, with the full ranking - and the explicit record that survival and the anatomic labels were not used.
- `equal_footing_certificate.json` — Content hashes proving every method and control was given identical inputs, identical tumours and identical folds.
- `implementation_report.json` — For each of MuSiC, DWLS, Bisque and SCDC: whether the genuine R package ran, and if not, exactly why the Python reimplementation was used instead.
- `registration_status.json` — Whether the constraint file is actually pre-registered. A hash proves the constraints have not changed; it does not prove they were written before the results, which is what the pre-registration claim rests on. UNREGISTERED is the honest default and is stated as such.
- `method_configs.json` — Every method's configuration as run, with anything differing from its own defaults marked, and every declared departure from a tool's PUBLISHED defaults listed with the reason and when it was decided. The fairness rule is defaults with no per-method tuning, so a departure that is invisible is not defensible.
- `reference_coverage.json` — What the cell-type roster drops from the atlas, counted. Dropped populations do not leave the tissue, only the model, so their expression is absorbed by the nearest retained type - and which constraints that lands on is reported.
- `control_calibration.json` — Each negative control scored over many independent draws, not once. The leaderboard's control rows are single draws; on real data the shuffled-signature control has a standard deviation of about 0.16 against a real-method spread of 0.29, so one draw cannot establish that the controls behave in either direction. Read this before reading a control row.
- `acs_per_tumor.csv` — Every method x every tumour: which constraints that tumour could evaluate and which it satisfied. The protocol asks for per-tumour results explicitly, because with 9 evaluable tumours a pooled ACS can be carried by one or two of them and the pooled number alone cannot show that.
- `acs_cohort_sensitivity.csv` — ACS for each method when only the 122 anatomic samples are deconvolved, beside ACS when all 270 archive samples are, with the rank change. Several methods use cross-sample statistics, so this measures how much of the leaderboard is a property of the method rather than of the cohort it was handed.
- `data_reconciliation.json` — Protocol step 3's gate: our parse of the frozen 2014-11-25 archive checked sample-by-sample and structure-by-structure against the Allen Institute portal's own live metadata export, which is produced independently of this pipeline.
- `clinical_missingness_audit.json` — Whether a blank survival_days is missing at random. It is not: the blanks are strongly enriched for MGMT methylation and younger age, the two strongest favourable prognostic factors in GBM. Read before any survival number.
- `full_database_acs_leaderboard.csv` — The leaderboard recomputed when every archive sample is deconvolved rather than only the anatomic ones. ACS is still scored on the anatomic subset alone - the ISH-cluster samples had their structure assigned from expression, so scoring them would be circular.
- `full_database_anatomic_report.json` — The full-database run's machine-readable summary, including how many samples were deconvolved versus how many were scored.
- `survival_power.json` — The canonical survival verdict. Ivy GAP publishes survival_days and no vital-status column, so who was censored is unknowable and this is BLOCKED. Nothing is imputed.
- `survival_power_declared_policy.json` — Present only if the run was explicitly asked to read recorded survival times as observed deaths. Carries the assumption, the missingness audit, and what the event count can actually detect - computed before any model was fitted.
- `survival_metrics_declared_policy.csv` — Out-of-fold prognostic evaluation per method under that declared assumption. NOT the canonical result. Read survival_power_declared_policy.json first; on this many events every row is INCONCLUSIVE by design.

## Not in this build

These stages had not produced output when the bundle was built:

- `benchmark/gene_space.json`

## Survival

BLOCKED: tumor_details.csv is present and joins correctly, but it publishes survival_days and NO vital-status column. Who was censored is not knowable from the release, so a C-index cannot be computed. Nothing is imputed and no prognostic claim is made.

## The primary result

AVAILABLE

## How the published tools were run

| tool | genuine R package | reason |
|---|---|---|
| music | yes | MuSiC available with cell-level reference |
| dwls | yes | DWLS available with cell-level reference |
| bisque | yes | BisqueRNA available with cell-level reference |
| scdc | yes | SCDC available with cell-level reference |
| scdc_ensemble | yes | SCDC available with cell-level reference |
