"""
OPEN_DEFECTS D31: the genuine DWLS package's first step (solveOLSInternal) hands quadprog an unscaled
t(S) %*% S, and on counts-per-million inputs quadprog reports "constraints are inconsistent, no solution!"
for a feasible, well-conditioned problem. R/run_dwls.R divides the signature and the bulk by one constant.

Pinned on a real case (tests/fixtures/dwls_conditioning_ivygap.csv: the raw/X Ivy GAP DWLS signature and two
Ivy GAP samples, one that failed in the recorded run and one that solved):
- the failure is real: unscaled, the failing sample errors with quadprog's message;
- the fix removes it: rescaled, the same sample solves;
- the fix changes nothing else: where the unscaled problem solves, the rescaled solution is identical, for the
  OLS step and for a dampened step.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FIXTURE = ROOT / "tests" / "fixtures" / "dwls_conditioning_ivygap.csv"


# The two DWLS steps are copied verbatim from the installed package (DWLS 0.1.0, namespace functions
# solveOLSInternal and solveDampenedWLSj) so that the test needs only quadprog: loading DWLS itself pulls
# in Seurat and MAST and costs over a minute.
R_CODE = r"""
if (!requireNamespace("quadprog", quietly = TRUE)) { cat("SKIP: quadprog\n"); quit(status = 0) }
library(quadprog)
solveOLSInternal <- function(S, B) {
  D <- t(S) %%*%% S; d <- t(S) %%*%% B; A <- cbind(diag(dim(S)[2])); bzero <- c(rep(0, dim(S)[2]))
  solution <- solve.QP(D, d, A, bzero)$solution; names(solution) <- colnames(S); solution }
solveDampenedWLSj <- function(S, B, goldStandard, j) {
  multiplier <- 1 * 2^(j - 1); sol <- goldStandard; ws <- as.vector((1/(S %%*%% sol))^2)
  wsScaled <- ws/min(ws); wsDampened <- wsScaled; wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened); D <- t(S) %%*%% W %%*%% S; d <- t(S) %%*%% W %%*%% B
  A <- cbind(diag(dim(S)[2])); bzero <- c(rep(0, dim(S)[2])); sc <- norm(D, "2")
  solution <- solve.QP(D/sc, d/sc, A, bzero)$solution; names(solution) <- colnames(S); solution }
fx <- read.csv("%s", row.names = 1, check.names = FALSE)
S <- as.matrix(fx[, setdiff(colnames(fx), c("bulk_fails_unscaled", "bulk_solves_unscaled"))])
bf <- fx$bulk_fails_unscaled; bs <- fx$bulk_solves_unscaled
k <- max(S)
try_ols <- function(S, b) tryCatch(solveOLSInternal(S, b), error = function(e) conditionMessage(e))
a <- try_ols(S, bf); b <- try_ols(S / k, bf / k)
cat("FAIL_UNSCALED:", if (is.character(a)) a else "solved", "\n")
cat("FAIL_RESCALED:", if (is.character(b)) b else "solved", "\n")
o1 <- try_ols(S, bs); o2 <- try_ols(S / k, bs / k)
cat("OLS_DIFF:", max(abs(o1 / sum(o1) - o2 / sum(o2))), "\n")
d1 <- solveDampenedWLSj(S, bs, o1, 3); d2 <- solveDampenedWLSj(S / k, bs / k, o2, 3)
cat("DAMPENED_DIFF:", max(abs(d1 / sum(d1) - d2 / sum(d2))), "\n")
"""


def _run() -> dict[str, str]:
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    r = subprocess.run(["Rscript", "-e", R_CODE % FIXTURE], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    if r.stdout.startswith("SKIP"):
        pytest.skip("R package quadprog not installed")
    return {ln.split(":", 1)[0]: ln.split(":", 1)[1].strip() for ln in r.stdout.splitlines() if ":" in ln}


def test_quadprog_failure_is_real_and_rescaling_removes_it_without_changing_any_solution():
    out = _run()
    assert "constraints are inconsistent" in out["FAIL_UNSCALED"]      # the defect, reproduced
    assert out["FAIL_RESCALED"] == "solved"                             # the fix
    assert float(out["OLS_DIFF"]) < 1e-10                               # nothing else moves
    assert float(out["DAMPENED_DIFF"]) < 1e-10


def test_the_bridge_applies_the_rescaling():
    src = (ROOT / "R" / "run_dwls.R").read_text()
    assert "signature <- signature / dwls_scale" in src and "bulk_mat <- bulk_mat / dwls_scale" in src
    assert 'Sys.getenv("IVYGAP_REPAIRED", "0"), "1"' in src           # opt-in until the verification finishes


def test_repaired_reimplementation_reproduces_the_package_given_its_dampening_constant(monkeypatch):
    """OPEN_DEFECTS D32: with the package's j, the published-dampening reimplementation must return the
    package's proportions on real inputs; the legacy reimplementation (weights capped relative to the
    LARGEST weight) is the negative control and must not."""
    import json

    import numpy as np
    import pandas as pd

    from ivygap.deconv.reference_based import DWLSDeconvolution
    fx = pd.read_csv(FIXTURE, index_col=0)
    ref = json.loads((ROOT / "tests" / "fixtures" / "dwls_package_outputs.json").read_text())
    S = fx.drop(columns=["bulk_fails_unscaled", "bulk_solves_unscaled"])
    k = S.to_numpy().max()
    m = DWLSDeconvolution()
    rng = np.random.default_rng(0)
    for sample, col in ref["column_map"].items():
        want = np.array([ref["samples"][sample]["proportions"][t] for t in S.columns])
        b = fx[col].to_numpy(float) / k
        got = m._solve_published(S.to_numpy(float) / k, b, rng, j=ref["samples"][sample]["j"])
        assert np.max(np.abs(got - want)) < 1e-8
        legacy = m._solve_dampened(S.to_numpy(float) / k, b, power=8)
        assert np.max(np.abs(legacy - want)) > 0.01          # the legacy dampening is a different algorithm
