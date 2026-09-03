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
}


def sc_export_paths(ref_name: str) -> tuple[Path, Path]:
    """(counts, metadata) paths for a cell-level reference export."""
    return (config.REFERENCE_DIR / f"sc_counts_{ref_name}.csv.gz",
            config.REFERENCE_DIR / f"sc_meta_{ref_name}.csv")


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
    counts, meta = sc_export_paths(ref_name)
    if not (counts.exists() and meta.exists()):
        return Availability(
            method_name, False,
            f"cell-level reference export missing ({counts.name}). These methods need "
            f"individual cells with donor and cell-type labels, not a signature matrix.",
        )
    return Availability(method_name, True, f"{pkg} available with cell-level reference")


class RBridgeError(RuntimeError):
    pass


def run_r_method(method_name: str, data: DeconvolutionInput,
                 ref_name: str | None = None, timeout: int = 7200) -> pd.DataFrame:
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

    counts_path, meta_path = sc_export_paths(ref_name)

    with tempfile.TemporaryDirectory(prefix=f"ivygap_{method_name}_") as tmp:
        tmp = Path(tmp)
        bulk_path = tmp / "bulk.csv"
        out_path = tmp / "proportions.csv"
        data.bulk.to_csv(bulk_path)

        payload = {
            "bulk": str(bulk_path),
            "sc_counts": str(counts_path),
            "sc_meta": str(meta_path),
            "out": str(out_path),
            "cell_types": list(data.cell_types),
            "seed": config.RANDOM_SEED,
            "ensemble": method_name.endswith("_ensemble"),
        }
        cfg_path = tmp / "args.json"
        cfg_path.write_text(json.dumps(payload, indent=2))

        proc = subprocess.run(["Rscript", str(script), str(cfg_path)],
                              capture_output=True, text=True, timeout=timeout)
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

    def _solve_all(self, data: DeconvolutionInput) -> np.ndarray:
        try:
            result = run_r_method(self.r_method, data)
            self.implementation_ = f"R:{R_PACKAGES[self.r_method]}"
            self.fallback_reason_ = None
            return result.to_numpy(dtype="float64")
        except (RBridgeError, subprocess.TimeoutExpired) as exc:
            if not self.allow_fallback:
                raise
            self.implementation_ = "python-reimplementation"
            self.fallback_reason_ = str(exc).split("\n")[0]
            return np.asarray(self.fallback._solve_all(data), dtype="float64")
