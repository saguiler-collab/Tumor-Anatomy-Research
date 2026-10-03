# MIXTURE -- configuration fixed before any output

**Written 2026-10-02.** The time is this file's birth time (`stat`); no time is typed here. At writing,
MIXTURE had been installed and loaded but had produced no estimate on any data. Post-registration
extension panel, the same standing as FARDEEP, LinDeconSeq, ARIC, RNA-Sieve, ReCIDE and DESeq2
`unmix`: reported under its own name, never in a registered statistic.

## Why MIXTURE

Nguyen et al. 2024 [3] list it among the seven methods in the top 10 of reference-based methods in
all three of their scenarios. With `unmix` added it is the last of the seven not yet run here. The
user approved the install on 2026-10-02.

## The package

- MIXTURE 0.0.1, genuine and unmodified, vendored from github.com/elmerfer/MIXTURE @ 6332e160
  (`vendor/MIXTURE/`: COMMIT, SHA256SUMS, PROVENANCE.md) and installed from that directory.
- **Algorithm (from its source):** nu-SVR with a linear kernel, nu tuned over {0.25, 0.5, 0.75};
  negative coefficients clipped; normalised to sum to 1; cell types under 0.1% dropped recursively, up
  to 6 iterations.
- The signature is z-scored globally, CIBERSORT's convention, so inputs are linear-scale.

## The call (package defaults throughout; nothing tuned, nothing chosen on any outcome)

```
MIXTURE(expressionMatrix = bulk, signatureMatrix = reference profiles,
        functionMixture = nu.svm.robust.RFE (default), useCores = 1,
        nullDist = "none" (default), verbose = FALSE)        # no fileSave
GetMixture(result, type = "proportion")
```

- `bulk` and the reference profiles are the arm's own linear inputs on its shared gene space, exactly
  as every other signature-only extension method receives them.
- **The output is a proportion of the signature's profiles,** an mRNA proportion on the signature's
  scale, so the central cell-size conversion applies (NOT in `R_RETURNS_CELL_FRACTIONS`), as for
  CIBERSORT and DESeq2 `unmix`.
- A sample MIXTURE cannot estimate stays NaN, never zero-filled.

## Gate before any real data (`tests/test_mixture_gate.py`)

- **Planted:** mixtures built from the frozen signature with Dirichlet proportions must be recovered:
  minimum per-type r >= 0.90 and maximum absolute error <= 0.10. This is the bar `unmix` met.
- **Negative control:** the same mixtures against a reference whose gene labels are shuffled must
  FAIL that bar.
- **Reading:** if the planted set fails, MIXTURE does not proceed to real data, and that is reported.

## Real runs (after the gate)

- **Ivy GAP anatomic:** `scripts/remeasure_method.py --method mixture`, on the leaderboard's own
  inputs with all equivalence conditions.
- **TCGA:** `scripts/extension_tcga.py --methods mixture` on the frozen and raw/X arms, in both cohorts.
- Reported as every extension method is: anatomic ACS with its null; purity rho; lymphoid ordering
  against methylation.
