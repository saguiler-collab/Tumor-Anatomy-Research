"""
real_yardsticks.py — load the protocol's ground-truth yardsticks from actual data.

WHY THIS MODULE EXISTS SEPARATELY
---------------------------------
There is a large difference between "we could not obtain ground truth" and "ground truth
exists and covers 2 of the 10 methods". The first is an unknown; the second is a measured
limitation with a specific remedy. Reporting them the same way — an empty dict and an
`UNAVAILABLE` row — throws away the distinction and makes the project look less finished
than it is.

So the yardsticks are wired to real files, and the agreement test receives whatever
coverage actually exists. When that coverage is below what a rank correlation needs,
`agreement.py` says so with the number, rather than saying nothing was found.

WHAT IS ACTUALLY AVAILABLE, AND WHY
-----------------------------------
`synthetic_mixtures` — the protocol's yardstick 1. Real: 500 held-out pseudobulk
    mixtures with known composition, vendored in `reference_frozen/tcga_benchmark/`.
    Covers **NNLS and SVR only**, because those were the two methods the predecessor
    project ran against it. Adding the other eight requires the cell-level single-cell
    data the mixtures were pooled from, which is not in the repository.

`absolute_purity` — yardstick 2. The predecessor recorded a Tumor-vs-ABSOLUTE
    correlation of rho = 0.76 (n = 120) for one frozen method, and separately recorded
    the ABSOLUTE purity file itself as BLOCKED — never located. One method cannot
    contribute to a ranking and the correlation cannot be recomputed for others, so this
    yardstick returns nothing and says why.

`sc_pseudobulk` — yardstick 3. Requires paired single-cell data. Not present.

NAME MAPPING IS EXPLICIT
------------------------
The predecessor labelled its methods `NNLS` and `SVR`; this project uses `nnls` and
`svr`. The mapping is written out rather than inferred by lowercasing, because a silent
name-matching rule is how a method's score ends up attached to a different method.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ivygap import config

TCGA_BENCHMARK_DIR = config.FROZEN_REFERENCE_DIR / "tcga_benchmark"

#: predecessor label -> this project's method name. Explicit by design.
METHOD_NAME_MAP = {
    "NNLS": "nnls",
    "SVR": "svr",
}


def load_synthetic_mixtures() -> tuple[dict[str, float], dict]:
    """
    Yardstick 1: mean MAE over the seven primary cell types, per method, on 500 real
    held-out mixtures with known composition.

    Returns (scores, provenance). An empty `scores` means the file is absent; a small
    one means the file is present but thin, and the provenance says which.
    """
    path = TCGA_BENCHMARK_DIR / "benchmark_summary.csv"
    if not path.exists():
        return {}, {"available": False,
                    "reason": f"{path.name} not vendored in reference_frozen/"}

    df = pd.read_csv(path)
    primary = [c for c in config.PRIMARY_CELL_TYPES]
    df = df[df["cell_type"].isin(primary)]

    scores, unmapped = {}, []
    for label, block in df.groupby("method"):
        name = METHOD_NAME_MAP.get(str(label))
        if name is None:
            unmapped.append(str(label))
            continue
        # The source reports MAE in PERCENT; this project works in fractions. Scaling
        # is harmless for a rank correlation but wrong for any absolute comparison, so
        # convert rather than leave two unit systems in one table.
        scores[name] = float(block["mae"].mean()) / 100.0

    return scores, {
        "available": True,
        "source": str(path.relative_to(config.PROJECT_ROOT)),
        "sha256": config.sha256_file(path),
        "n_mixtures": 500,
        "is_real_ground_truth": True,
        "methods_covered": sorted(scores),
        "methods_not_mapped": unmapped,
        "units": "MAE as a fraction (source reports percent; divided by 100 here)",
        "limitation": (
            f"covers {len(scores)} method(s). A rank correlation needs at least "
            f"{6}. Extending it requires the cell-level single-cell data the mixtures "
            f"were pooled from."
        ),
    }


def load_absolute_purity() -> tuple[dict[str, float], dict]:
    """
    Yardstick 2. Returns nothing, with the reason — see the module docstring.

    Kept as a real function rather than an omission so the pipeline reports a specific
    blocker instead of an absence, and so wiring it later is a change in one place.
    """
    path = TCGA_BENCHMARK_DIR / "orthogonal_validation.csv"
    prov = {
        "available": False,
        "is_real_ground_truth": True,
        "reason": "the predecessor project recorded ABSOLUTE purity as BLOCKED — no "
                  "purity file was ever located — so the correlation exists for one "
                  "frozen method only and cannot be recomputed for the others.",
    }
    if path.exists():
        df = pd.read_csv(path)
        row = df[df["name"].str.contains("ABSOLUTE", case=False, na=False)]
        if not row.empty:
            prov["recorded_for_one_method"] = {
                "method": str(row.iloc[0]["method"]),
                "spearman_rho": float(row.iloc[0]["rho"]),
                "n": int(row.iloc[0]["n"]),
            }
    return {}, prov


def load_sc_pseudobulk() -> tuple[dict[str, float], dict]:
    """Yardstick 3. Requires paired single-cell data, which is not present."""
    return {}, {
        "available": False,
        "is_real_ground_truth": True,
        "reason": "requires paired single-cell GBM data where composition is counted "
                  "rather than estimated. Not present in this repository.",
    }


def load_all() -> tuple[dict[str, dict[str, float]], dict]:
    """
    Every protocol yardstick, plus a provenance block naming exactly what each one
    contributed and why.
    """
    scores, prov = {}, {}
    for name, loader in (("synthetic_mixtures", load_synthetic_mixtures),
                         ("absolute_purity", load_absolute_purity),
                         ("sc_pseudobulk", load_sc_pseudobulk)):
        s, p = loader()
        scores[name] = s
        prov[name] = p
    return scores, prov


def write_provenance(path: Path | None = None) -> dict:
    _, prov = load_all()
    path = path or (config.BENCH_DIR / "yardstick_provenance.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(prov, indent=2))
    return prov
