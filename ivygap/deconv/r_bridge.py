"""
r_bridge.py — run the genuine MuSiC / DWLS / Bisque / SCDC R packages.

WHEN THIS IS USED
-----------------
The registry prefers the real R package over the Python reimplementation whenever
both are usable, and records which one ran. "Usable" means three things are all true:

  1. `Rscript` is on PATH,
  2. the specific package is installed,
  3. a CELL-LEVEL single-cell reference export exists on disk.

The third condition is the binding one and is worth being explicit about. MuSiC,
Bisque and SCDC do not consume a signature matrix — they consume individual cells,
each labelled with a donor and a cell type, because their whole contribution is
computing statistics *across* donors and cells that a collapsed signature matrix has
already averaged away. Handing them a signature matrix would run the code but not the
method.

So the reference builder writes a cell-level export alongside the signature matrix
when it has the underlying single-cell data. Without that export this module reports
itself unavailable and the registry falls back to Python, recording the reason. It
never silently substitutes one for the other.

WHY SUBPROCESS AND NOT rpy2
---------------------------
rpy2 pins tightly to both the Python and R ABI and is a common cause of environments
that work on one machine and not the next. A subprocess boundary with CSV in and CSV
out is slower and completely robust, and this is not an inner loop.
"""

from __future__ import annotations

import gzip
import os
import signal
import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ivygap import config
from .base import DeconvolutionInput, DeconvolutionMethod

R_DIR = config.PROJECT_ROOT / "R"

#: R package required by each method
R_PACKAGES = {
    "music": "MuSiC",
    "dwls": "DWLS",
    "bisque": "BisqueRNA",
    "scdc": "SCDC",
    "scdc_ensemble": "SCDC",
    "epic": "EPIC",
    "quantiseq": "quantiseqr",
    "bayesprism": "BayesPrism",
}


def sc_export_paths(ref_name: str) -> tuple[Path, Path]:
    """(counts, metadata) paths for a cell-level reference export."""
    return (config.REFERENCE_DIR / f"sc_counts_{ref_name}.csv.gz",
            config.REFERENCE_DIR / f"sc_meta_{ref_name}.csv")


# -----------------------------------------------------------------------------
# The cell-level export, written per gene set instead of once for every gene
# -----------------------------------------------------------------------------
# Exporting the whole reference is not viable. Measured on the real atlas subsample,
# 16,758 genes x 15,311 cells takes ~49 minutes to write and 1.1 GB gzipped — and then
# R re-reads it for every method, on every cohort.
#
# The methods never need those genes. Each receives `data.bulk`, already restricted to
# the run's gene space (~1,600 on the frozen path, ~640 in the benchmark), and MuSiC,
# Bisque and SCDC all intersect against it internally anyway. Exporting the same gene
# set cuts the write by an order of magnitude and the R-side read with it.
#
# The source lives in memory because a ReferenceBundle deliberately does not carry
# cells — it carries summaries of them. Whoever built the reference registers it here.
_CELL_SOURCE: dict[str, tuple] = {}


def set_cell_source(ref_name: str, expression, meta) -> None:
    """Register the cell-level reference so per-gene-set exports can be written."""
    _CELL_SOURCE[ref_name] = (expression, meta)


def has_cell_source(ref_name: str) -> bool:
    return ref_name in _CELL_SOURCE


def get_cell_source(ref_name: str):
    """The registered (expression, meta) for `ref_name`, or None."""
    return _CELL_SOURCE.get(ref_name)


def clear_cell_source(ref_name: str) -> None:
    _CELL_SOURCE.pop(ref_name, None)


def export_for_genes(ref_name: str, genes) -> tuple[Path, Path]:
    """
    Write (or reuse) a cell-level export restricted to `genes`.

    Cached on a hash of the gene list, so the several methods in one run share a single
    write rather than repeating it.
    """
    import hashlib

    expression, meta = _CELL_SOURCE[ref_name]

    # Align metadata to the cells that actually have expression. The pipeline registers
    # matched pairs, but a caller need not, and a mismatch surfaces far away as a
    # KeyError naming cell barcodes — which says nothing about the cause.
    shared_cells = [c for c in expression.columns if c in meta.index]
    if not shared_cells:
        raise RBridgeError(
            f"no cell ids shared between the expression and metadata registered for "
            f"{ref_name!r} ({expression.shape[1]} cells vs {len(meta)} metadata rows)")
    if len(shared_cells) < expression.shape[1]:
        expression = expression[shared_cells]
    meta = meta.loc[shared_cells]

    genes = [g for g in genes if g in expression.index]
    if not genes:
        raise RBridgeError(
            f"none of the requested genes are in the cell-level reference {ref_name!r}")

    # Key on the gene list, the cell ids, AND THE VALUES.
    #
    # The values used to be absent from the key, which made this a cache that ignored the
    # content it was caching. Re-registering the same cells and genes at a different
    # scaling -- raw library sizes instead of 1e6 CPM, say -- produced an identical key and
    # silently handed R the stale file. That is not a hypothetical: it is the mechanism the
    # cell-size probe had to vary in order to measure anything, and it is the same shape as
    # D10, where a re-measurement ran on a gene space nobody had recorded.
    #
    # Hashed in column chunks rather than through one `to_numpy().tobytes()`, because the
    # production block is roughly 657 genes x 80,000 cells and a single materialised copy
    # of that is several hundred megabytes on a machine this pipeline already pushes into
    # swap. The hash is exact either way -- chunking changes only the peak memory.
    digest = hashlib.sha256()
    digest.update("\n".join(map(str, genes)).encode())
    digest.update(b"\x00")
    digest.update("\n".join(map(str, expression.columns)).encode())
    digest.update(b"\x00")
    _block_for_key = expression.loc[genes]
    digest.update(str(_block_for_key.shape).encode())
    digest.update(str(_block_for_key.to_numpy(dtype="float64", copy=False).dtype).encode())
    _step = max(1, 4096)
    for _i in range(0, _block_for_key.shape[1], _step):
        chunk = _block_for_key.iloc[:, _i:_i + _step].to_numpy(dtype="float64")
        digest.update(np.ascontiguousarray(chunk).tobytes())
    del _block_for_key
    key = digest.hexdigest()[:16]
    counts_path = config.REFERENCE_DIR / f"sc_counts_{ref_name}_{key}.csv.gz"
    meta_path = config.REFERENCE_DIR / f"sc_meta_{ref_name}.csv"

    # Restricting to a signature gene set leaves some cells with zero counts across
    # every gene kept. They carry no information about composition, and BisqueRNA
    # refuses outright — "Zero expression in selected genes for N cells" from
    # CountsToCPM — which sends the genuine package to the Python fallback for a reason
    # that has nothing to do with Bisque.
    #
    # Dropped here rather than inside one driver, so every R method sees the same cells.
    block = expression.loc[genes]
    nonzero = block.sum(axis=0) > 0
    n_dropped = int((~nonzero).sum())
    if n_dropped:
        block = block.loc[:, nonzero]
        meta = meta.loc[block.columns]
        print(f"  {n_dropped} cell(s) have zero expression across the {len(genes):,} "
              f"exported genes and are dropped from the R export "
              f"({block.shape[1]:,} remain)", flush=True)

    if not counts_path.exists():
        config.ensure_dirs()
        print(f"  writing cell-level export for {len(genes):,} genes "
              f"x {block.shape[1]:,} cells ...", flush=True)
        # Write-then-rename. An interrupted run previously left a partial
        # sc_counts_*.csv.gz on disk, and `check()` treats existence as availability —
        # so the next run would have handed R a silently truncated reference.
        tmp_path = counts_path.with_suffix(counts_path.suffix + ".part")
        with gzip.open(tmp_path, "wt") as fh:
            block.to_csv(fh)
        tmp_path.replace(counts_path)
    # ALWAYS rewrite the metadata, never reuse it. The counts file is cached on a gene
    # hash, but the metadata belongs to the CELL SOURCE — and a stale sc_meta left by a
    # different reference build lists different cell ids, which surfaces in R as
    # "no shared cell ids between counts and metadata" with no hint of the cause.
    meta.to_csv(meta_path)
    return counts_path, meta_path


def rscript_available() -> bool:
    return shutil.which("Rscript") is not None


#: Result of every package probe made in this process: package -> (available, why).
#: Whether an R package is installed cannot change while the process runs, so probing
#: once is not an optimisation with a correctness cost — it removes one.
_PACKAGE_PROBE: dict[str, tuple[bool, str]] = {}


def probe_r_package(package: str) -> tuple[bool, str]:
    """
    Ask R itself whether the package loads, once per process, and say how it knows.

    Two things this fixes, both real.

    **The probe was uncached.** Every `check()` spawned a fresh R process with a 120 s
    timeout. A 15-method run makes that call many times over, and the test suite more; a
    single pytest run was observed spending over twenty minutes almost entirely in
    repeated `requireNamespace` probes while the machine was busy.

    **A failed probe was indistinguishable from an absent package.** The old code caught
    every exception — `TimeoutExpired` included — and returned False, which `check()` then
    reported as "R package X is not installed". Under load the probe can time out while
    the package is installed and perfectly usable, and the method would quietly run its
    Python reimplementation instead. The fallback is recorded, so nothing is hidden, but
    the recorded *reason* would be false and the leaderboard would depend on how busy the
    machine was. The three states are now distinct, and a probe that fails for any reason
    other than the package being absent says so.
    """
    if package in _PACKAGE_PROBE:
        return _PACKAGE_PROBE[package]
    if not rscript_available():
        result = (False, "Rscript is not on PATH")
    else:
        try:
            proc = subprocess.run(
                ["Rscript", "-e",
                 f'cat(if (requireNamespace("{package}", quietly=TRUE)) "yes" else "no")'],
                capture_output=True, text=True, timeout=120,
            )
            if proc.stdout.strip().endswith("yes"):
                result = (True, f"R package {package} loads")
            elif proc.stdout.strip().endswith("no"):
                result = (False, f"R package {package} is not installed")
            else:
                result = (False,
                          f"the probe for R package {package} returned neither yes nor no "
                          f"(exit {proc.returncode}); treating it as unavailable, but this "
                          f"is a PROBE FAILURE, not evidence the package is absent")
        except subprocess.TimeoutExpired:
            result = (False,
                      f"the probe for R package {package} timed out after 120 s. This is a "
                      f"PROBE FAILURE, not evidence the package is absent — most likely the "
                      f"machine was loaded. A method falling back for this reason is running "
                      f"its Python reimplementation for a reason unrelated to the method.")
        except Exception as exc:                          # noqa: BLE001
            result = (False,
                      f"the probe for R package {package} raised "
                      f"{type(exc).__name__}; PROBE FAILURE, not evidence of absence")
    _PACKAGE_PROBE[package] = result
    return result


def r_package_available(package: str) -> bool:
    """Ask R itself whether the package loads — presence on disk is not enough."""
    return probe_r_package(package)[0]


def clear_package_probe_cache() -> None:
    """For tests, and for a process that installs a package after probing for it."""
    _PACKAGE_PROBE.clear()


@dataclass
class Availability:
    """Why a method can or cannot run through R. Recorded in the results manifest."""
    method: str
    available: bool
    reason: str

    def as_dict(self) -> dict:
        return {"method": self.method, "r_available": self.available, "reason": self.reason}


def check(method_name: str, ref_name: str) -> Availability:
    pkg = R_PACKAGES.get(method_name)
    if pkg is None:
        return Availability(method_name, False, "no R implementation is wired for this method")
    if not rscript_available():
        return Availability(method_name, False, "Rscript is not on PATH")
    ok, why = probe_r_package(pkg)
    if not ok:
        # `why` distinguishes "not installed" from "the probe failed". Recording the
        # difference matters: a fallback caused by a loaded machine is not a fact about
        # the method, and must not read like one in implementation_report.json.
        return Availability(method_name, False, why)
    # Either a registered in-memory cell source (from which a per-gene-set export is
    # written on demand) or a pre-written full export will do.
    # EPIC and quanTIseq consume a signature matrix, which every run already has, so
    # they do not require a cell-level export the way MuSiC/Bisque/SCDC/BayesPrism do.
    if method_name in SIGNATURE_ONLY_METHODS:
        return Availability(method_name, True, f"{pkg} available (signature-based)")

    counts, meta = sc_export_paths(ref_name)
    if not has_cell_source(ref_name) and not (counts.exists() and meta.exists()):
        return Availability(
            method_name, False,
            f"cell-level reference export missing ({counts.name}) and no cell source "
            f"registered. These methods need individual cells with donor and cell-type "
            f"labels, not a signature matrix.",
        )
    return Availability(method_name, True, f"{pkg} available with cell-level reference")


#: Wall-clock budget per method, in seconds. A method that cannot finish inside its
#: budget falls back to its Python reimplementation with the timeout recorded as the
#: reason — the same disclosure path as any other R failure.
#:
#: This exists because DWLS can consume the entire run. buildSignatureMatrixMAST's cost
#: is dominated by a condition-number search over gene counts, not by the cells or the
#: differential-expression method, so neither subsampling nor swapping the DE step
#: bounds it: one observed call ran for over two hours at 100% CPU without finishing.
#:
#: A benchmark that never finishes produces no result at all, which is strictly worse
#: than a disclosed fallback. The budget is generous enough that the methods measured to
#: complete here (MuSiC ~95 s, SCDC ~90 s, Bisque ~30 s) are nowhere near it.
#: Diagnostic files an R script wrote beside its output on the most recent call, per
#: method. Populated by `run_r_method`; read by the run reporters so a sidecar is an
#: artefact rather than something that existed for the lifetime of a temp directory.
LAST_SIDECARS: dict[str, list] = {}

#: What gene space each R method actually received on its most recent call. Read by the
#: run reporters so "which genes did this method see" is an artefact, not a claim.
LAST_GENE_SPACE: dict[str, dict] = {}

#: R packages MEASURED to return CELL fractions rather than mRNA shares, so this project's
#: central conversion would be the second one and must be skipped.
#:
#: Measured by `scripts/verify_cell_size_semantics.py --mixture all_types`: two cell types
#: differing 3x in mRNA, mixed at equal cell counts, so 0.750 is an mRNA share and 0.500 a
#: cell share. MuSiC 0.7502, SCDC 0.7503, EPIC 0.7502, BayesPrism 0.7502 -- all mRNA share.
#: Bisque 0.5000, and its off-pair mass of 0.750 against the others' 0.600 confirms it
#: independently. Artefact: results/cell_size_semantics_normalized_all_types.json.
#: See docs/OPEN_DEFECTS.md D1.
#:
#: Keyed on the R method, NOT on the method name, because the convention belongs to the
#: implementation that actually ran. This project's Python reimplementation of Bisque has
#: NOT been measured on that probe, so the fallback path is deliberately left unflagged:
#: assuming it shares the R package's convention would be exactly the unmeasured inference
#: D1 has already been corrected for twice.
R_RETURNS_CELL_FRACTIONS = frozenset({"bisque"})

#: Methods that consume a reference PROFILE rather than individual cells.
SIGNATURE_ONLY_METHODS = frozenset({"epic", "quantiseq"})

#: Methods that bring their OWN published signature and ignore the one they are handed.
#:
#: EPIC is signature-only but reads `args$signature` — this project's reference — so the
#: project's gene space is the right one for it and equal footing is untouched. quanTIseq
#: reads nothing but its built-in TIL10, so genes ranked against a different reference
#: are simply the wrong genes for it. These methods get `bulk_full` when the caller
#: supplied it; the choice is recorded per run rather than assumed.
OWN_SIGNATURE_METHODS = frozenset({"quantiseq"})

#: Roster types each method structurally models, for the methods that model only some.
#: Mirrors `DeconvolutionMethod.models_cell_types`; kept here so the all-NaN guard below
#: can tell "this method does not do tumour cells" from "this method returned nothing".
OWN_COVERAGE: dict[str, frozenset[str]] = {
    "quantiseq": frozenset({"T_cell", "NK_cell", "B_cell", "Macrophage_Microglia"}),
}

R_METHOD_TIMEOUTS: dict[str, int] = {
    # 40 minutes. Measured: DWLS's signature build completed in 1,625 s at best on this
    # data, so the budget clears a successful run with room, and bounds the case where
    # the condition-number search does not converge.
    "dwls": 2400,
    # BayesPrism's Gibbs sampler scales with cells x genes x types.
    "bayesprism": 2400,
}
DEFAULT_R_TIMEOUT = 3600


def timeout_for(method_name: str) -> int:
    return R_METHOD_TIMEOUTS.get(method_name, DEFAULT_R_TIMEOUT)


def _summarise_r_failure(exc: Exception) -> str:
    """
    A fallback reason that names the CAUSE, not just the exit code.

    The first line of an RBridgeError is "run_music.R exited 1", which records that R
    failed and nothing about why. Every artefact that discloses "this ran as a Python
    reimplementation" then says so without the one detail that would let anyone fix it —
    and the disclosure exists precisely so a reimplementation is never mistaken for the
    published package.

    R writes its diagnostics as "Error in <call> : <message>", so that line is lifted
    out when present and appended to the exit line.
    """
    text = str(exc)
    head = text.split("\n")[0]

    lines = [ln.strip() for ln in text.splitlines()]
    for i, ln in enumerate(lines):
        if ln.startswith("Error"):
            # R wraps long messages onto the following line.
            detail = ln
            if i + 1 < len(lines) and lines[i + 1] and not lines[i + 1].startswith(
                    ("Error", "Calls:", "Execution halted", "---")):
                detail = f"{detail} {lines[i + 1]}"
            return f"{head} | {detail[:400]}"
    return head


def _run_bounded(cmd: list[str], budget: int):
    """
    Run `cmd` with a timeout that actually fires.

    `subprocess.run(..., timeout=N)` is not sufficient here, and the failure is silent.
    On timeout it kills the DIRECT child and then calls `communicate()` a second time,
    with no timeout, to reap the output. If that child spawned grandchildren holding the
    inherited stdout pipe, the second call blocks forever: the TimeoutExpired never
    propagates and the "budget" is unbounded in practice.

    That is exactly what happened to DWLS. MAST's workers keep the pipe open, so a
    2,400 s budget was still running at 5,280 s — verified in isolation that
    `timeout_for('dwls')` returned 2400, that the signature default was None, and that a
    plain sleeping R child DOES time out correctly. The mechanism was fine; the pipe was
    the problem.

    So: put the child in its own process group, and on timeout kill the GROUP, which
    reaches the grandchildren holding the pipe. Then reap with a short bounded wait.
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, start_new_session=True)
    try:
        out, err = proc.communicate(timeout=budget)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):       # pragma: no cover
            proc.kill()
        try:
            proc.communicate(timeout=30)
        except subprocess.TimeoutExpired:                   # pragma: no cover
            pass
        raise
    return subprocess.CompletedProcess(cmd, proc.returncode, out, err)


class RBridgeError(RuntimeError):
    pass


def run_r_method(method_name: str, data: DeconvolutionInput,
                 ref_name: str | None = None,
                 timeout: int | None = None) -> pd.DataFrame:
    """
    Execute one R deconvolution method and return samples x cell_types RNA proportions.

    Raises RBridgeError on any failure. Callers decide whether to fall back; this
    function never quietly returns a Python-computed result under an R method's name.
    """
    ref_name = ref_name or data.primary.name
    avail = check(method_name, ref_name)
    if not avail.available:
        raise RBridgeError(f"{method_name} unavailable via R: {avail.reason}")

    script = R_DIR / f"run_{method_name.replace('_ensemble', '')}.R"
    if not script.exists():
        raise RBridgeError(f"missing driver script {script}")

    # Restrict the export to the genes the method is actually given. See the note on
    # export_for_genes: the full-reference export is ~49 minutes and 1.1 GB, and none of
    # it beyond `data.bulk`'s gene space is ever used.
    if method_name in SIGNATURE_ONLY_METHODS:
        counts_path = meta_path = Path("")          # unused by these drivers
    elif has_cell_source(ref_name):
        counts_path, meta_path = export_for_genes(ref_name, list(data.bulk.index))
    else:
        counts_path, meta_path = sc_export_paths(ref_name)

    with tempfile.TemporaryDirectory(prefix=f"ivygap_{method_name}_") as tmp:
        tmp = Path(tmp)
        bulk_path = tmp / "bulk.csv"
        out_path = tmp / "proportions.csv"

        bulk_for_r = data.bulk
        gene_space_note = "shared marker subset"
        if method_name in OWN_SIGNATURE_METHODS:
            if data.bulk_full is not None:
                bulk_for_r = data.bulk_full
                gene_space_note = "full shared space (method supplies its own signature)"
            else:
                gene_space_note = (
                    "shared marker subset — NO bulk_full was supplied, so this method "
                    "ran on genes selected against a signature it does not use"
                )
        bulk_for_r.to_csv(bulk_path)
        LAST_GENE_SPACE[method_name] = {
            "n_genes": int(bulk_for_r.shape[0]),
            "gene_space": gene_space_note,
        }

        # Methods that consume a signature matrix rather than cells (EPIC, quanTIseq)
        # get the same reference profile every Python method solves against, written to
        # the same temp dir. Written for every method so the payload is one shape.
        sig_path = tmp / "signature.csv"
        data.primary.profile.to_csv(sig_path)

        # PER-GENE VARIABILITY, as a STANDARD DEVIATION.
        #
        # EPIC weights genes by `rowSums(refProfiles / (refProfiles.var + 1e-12))`, a
        # signal-to-noise ratio that is only dimensionally sensible if refProfiles.var is on
        # the same scale as the profile. EPIC's own bundled TRef confirms it: its
        # refProfiles.var has median 30.9 against profile median 9.2 and max 78,469 against
        # 71,786 -- the same order of magnitude, not a squared one.
        #
        # `ReferenceBundle.sigma` is a VARIANCE, so it is square-rooted here. Passing it raw
        # would divide a mean by a squared spread and systematically over-weight
        # low-expression genes, which is a silent mis-weighting rather than an error.
        #
        # Written for every method; only the ones that ask for it read it. Until 2026-09-14
        # it was written for none, so EPIC ran with `using identical weights for all genes`
        # -- its published contribution switched off. See docs/OPEN_DEFECTS.md D13.
        var_path = tmp / "signature_var.csv"
        (data.primary.sigma.clip(lower=0.0) ** 0.5).to_csv(var_path)

        payload = {
            "bulk": str(bulk_path),
            "signature": str(sig_path),
            "signature_var": str(var_path),
            "sc_counts": str(counts_path),
            "sc_meta": str(meta_path),
            "out": str(out_path),
            "cell_types": list(data.cell_types),
            "seed": config.RANDOM_SEED,
            "ensemble": method_name.endswith("_ensemble"),
        }
        cfg_path = tmp / "args.json"
        cfg_path.write_text(json.dumps(payload, indent=2))

        budget = timeout if timeout is not None else timeout_for(method_name)
        proc = _run_bounded(["Rscript", str(script), str(cfg_path)], budget)
        if proc.returncode != 0:
            raise RBridgeError(
                f"{script.name} exited {proc.returncode}\n"
                f"--- stdout ---\n{proc.stdout[-4000:]}\n"
                f"--- stderr ---\n{proc.stderr[-4000:]}"
            )
        if not out_path.exists():
            raise RBridgeError(f"{script.name} exited 0 but wrote no output file")

        result = pd.read_csv(out_path, index_col=0)

        # SIDECARS. Some scripts write diagnostics beside their output -- EPIC's
        # `otherCells` fraction and its per-sample convergence codes are the ones that
        # matter. They used to be written into this temporary directory and deleted with
        # it, so the numbers needed to judge EPIC's estimand were produced on every run
        # and never survived one. Anything matching the output's stem is copied out.
        for extra in sorted(tmp.glob(f"{out_path.stem}_*.csv")):
            dest = config.DIAGNOSTICS_DIR / f"{method_name}_{extra.name.split('_', 1)[1]}"
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(extra, dest)
            LAST_SIDECARS.setdefault(method_name, []).append(
                str(dest.relative_to(config.PROJECT_ROOT)))

    # R reorders and renames freely; force the project's roster and sample order back
    # on, and fail loudly if something is genuinely absent rather than filling zeros.
    missing = set(data.cell_types) - set(result.columns)
    if missing:
        raise RBridgeError(f"{method_name} returned no column for: {sorted(missing)}")
    result.index = result.index.astype(str)
    missing_samples = set(data.samples) - set(result.index)
    if missing_samples:
        raise RBridgeError(
            f"{method_name} returned no estimate for {len(missing_samples)} samples"
        )
    result = result.loc[data.samples, list(data.cell_types)].astype("float64")
    _reject_unusable(method_name, result, data, script.name)
    return result


def _reject_unusable(method_name: str, result: "pd.DataFrame", data, script_name: str):
    """
    Refuse a result that carries no estimate, however well-formed it looks.

    "Exited 0" is not evidence that a method produced anything. quantiseqr exited 0,
    wrote a well-formed 122x8 CSV with the right index and the right column names, and
    every cell in it was NaN — and run 2026-09-05T2154 recorded that as
    `implementation: R:quantiseqr, fallback_reason: null`. A table with no numbers in it
    is a failure of the R path, and the caller's fallback machinery is exactly what
    should handle it; silently ranking it is what must not happen.

    Checked over the types the method actually models, so a partial-coverage method's
    legitimately all-NaN columns do not trip it.
    """
    declared = OWN_COVERAGE.get(method_name)
    modelled = ([c for c in data.cell_types if c in declared] if declared
                else list(data.cell_types))
    if not modelled:
        raise RBridgeError(
            f"{method_name} declares coverage {sorted(declared or [])}, none of which "
            f"is in the roster {list(data.cell_types)}"
        )

    values = result[modelled].to_numpy(dtype="float64")
    if not np.isfinite(values).any():
        raise RBridgeError(
            f"{script_name} exited 0 and wrote a {result.shape[0]}x{result.shape[1]} "
            f"table, but every value across the {len(modelled)} cell type(s) it models "
            f"is NaN or infinite. The R path produced no estimate."
        )

    # Not fatal — a method may legitimately fail particular mixtures (DWLS does) — but
    # it is recorded, because "some samples have no estimate" must never be invisible.
    n_dead = int((~np.isfinite(values)).all(axis=1).sum())
    if n_dead:
        LAST_GENE_SPACE.setdefault(method_name, {})["n_samples_without_estimate"] = n_dead


class RMethod(DeconvolutionMethod):
    """
    Wrapper presenting an R package as a normal project method.

    `fallback` is the Python reimplementation to use when R is unavailable. Which one
    actually ran is recorded in `implementation_` and travels into the manifest — the
    two are never conflated in reporting.
    """

    requires_r = True

    def __init__(self, r_method: str, fallback: DeconvolutionMethod,
                 allow_fallback: bool = True, **params):
        super().__init__(r_method=r_method, allow_fallback=allow_fallback, **params)
        self.r_method = r_method
        self.fallback = fallback
        self.allow_fallback = allow_fallback
        self.name = fallback.name
        self.family = fallback.family
        self.uses_multiple_references = fallback.uses_multiple_references
        # Structural coverage is a property of the METHOD, not of which implementation
        # runs it: quanTIseq has no tumour population whether quantiseqr or the
        # reimplementation solves it. Carried over so the base class applies the same
        # partial-coverage handling on both paths.
        self.models_cell_types = fallback.models_cell_types
        self.implementation_: str = "unknown"
        self.fallback_reason_: str | None = None
        # Degradation is a property of the DATA, not of which implementation ran. Bisque
        # is in no-overlap mode whether the genuine R package or the reimplementation
        # solves it, so these are set on both paths — otherwise the disclosure vanishes
        # the moment the real package starts working.
        self.degenerate_: bool = False
        self.degeneracy_reason_: str | None = None

    @property
    def returns_cell_fractions(self) -> bool:
        """
        True only when the GENUINE R package ran AND that package is measured to convert
        cell size itself.

        Dynamic rather than a class attribute because the two implementations behind this
        wrapper can disagree: R:BisqueRNA is measured to return cell fractions, and the
        Python reimplementation has not been measured at all. `_solve_all` runs before the
        base class reads this, so `implementation_` is already known.
        """
        return (self.r_method in R_RETURNS_CELL_FRACTIONS
                and str(self.implementation_).startswith("R:"))

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        try:
            result = run_r_method(self.r_method, data)
            self.implementation_ = f"R:{R_PACKAGES[self.r_method]}"
            self.fallback_reason_ = None
            self.degenerate_, self.degeneracy_reason_ = self.fallback.degradation_for(data)
            return result.to_numpy(dtype="float64")
        except (RBridgeError, subprocess.TimeoutExpired) as exc:
            if not self.allow_fallback:
                raise
            self.implementation_ = "python-reimplementation"
            if isinstance(exc, subprocess.TimeoutExpired):
                self.fallback_reason_ = (
                    f"exceeded its {timeout_for(self.r_method)}s budget for the genuine "
                    f"R package and was stopped; this is a wall-clock limit, not a "
                    f"failure of the method")
            else:
                self.fallback_reason_ = _summarise_r_failure(exc)
            out = np.asarray(self.fallback._solve_all(data), dtype="float64")
            self.degenerate_ = bool(getattr(self.fallback, "degenerate_", False))
            self.degeneracy_reason_ = getattr(self.fallback, "degeneracy_reason_", None)
            return out
