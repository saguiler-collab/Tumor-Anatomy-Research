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

    # Key on the cells as well as the genes: two references can share a gene set and
    # hold entirely different cells, and reusing one's counts for the other is exactly
    # the mismatch the metadata comment above describes.
    digest = hashlib.sha256()
    digest.update("\n".join(map(str, genes)).encode())
    digest.update(b"\x00")
    digest.update("\n".join(map(str, expression.columns)).encode())
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


def r_package_available(package: str) -> bool:
    """Ask R itself whether the package loads — presence on disk is not enough."""
    if not rscript_available():
        return False
    try:
        proc = subprocess.run(
            ["Rscript", "-e",
             f'cat(if (requireNamespace("{package}", quietly=TRUE)) "yes" else "no")'],
            capture_output=True, text=True, timeout=120,
        )
        return proc.stdout.strip().endswith("yes")
    except Exception:                                    # noqa: BLE001
        return False


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
    if not r_package_available(pkg):
        return Availability(method_name, False, f"R package {pkg} is not installed")
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
#: Methods that consume a reference PROFILE rather than individual cells.
SIGNATURE_ONLY_METHODS = frozenset({"epic", "quantiseq"})

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
        data.bulk.to_csv(bulk_path)

        # Methods that consume a signature matrix rather than cells (EPIC, quanTIseq)
        # get the same reference profile every Python method solves against, written to
        # the same temp dir. Written for every method so the payload is one shape.
        sig_path = tmp / "signature.csv"
        data.primary.profile.to_csv(sig_path)

        payload = {
            "bulk": str(bulk_path),
            "signature": str(sig_path),
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
    return result.loc[data.samples, list(data.cell_types)].astype("float64")


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
        self.implementation_: str = "unknown"
        self.fallback_reason_: str | None = None
        # Degradation is a property of the DATA, not of which implementation ran. Bisque
        # is in no-overlap mode whether the genuine R package or the reimplementation
        # solves it, so these are set on both paths — otherwise the disclosure vanishes
        # the moment the real package starts working.
        self.degenerate_: bool = False
        self.degeneracy_reason_: str | None = None

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
