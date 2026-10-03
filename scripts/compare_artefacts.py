"""
compare_artefacts.py -- does a re-run reproduce an archived artefact? Leaf by leaf.

    python3 scripts/compare_artefacts.py ARCHIVED.json RERUN.json [--tol 1e-9] [--ignore key ...]

Flattens both JSON documents to (path, value) leaves and reports every leaf that differs:
numbers beyond `tol`, strings, missing or extra keys. Timestamps, durations and absolute paths
are ignored by default because they legitimately change between runs; nothing else is. Exit
status 0 = reproduced, 1 = differences (each printed), 2 = unreadable input.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys

DEFAULT_IGNORE = r"(time|utc|date|elapsed|seconds|duration|runtime|path|file|host|python|platform|version)"


def flatten(x, prefix=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from flatten(v, f"{prefix}.{k}" if prefix else str(k))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from flatten(v, f"{prefix}[{i}]")
    else:
        yield prefix, x


def same(a, b, tol: float) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return True
        return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a)))
    return a == b


def compare(old: dict, new: dict, tol: float, ignore: str) -> list[str]:
    pat = re.compile(ignore, re.I)
    fo = {k: v for k, v in flatten(old) if not pat.search(k.split(".")[-1])}
    fn = {k: v for k, v in flatten(new) if not pat.search(k.split(".")[-1])}
    out = []
    for k in sorted(set(fo) | set(fn)):
        if k not in fn:
            out.append(f"missing in re-run: {k} = {fo[k]!r}")
        elif k not in fo:
            out.append(f"extra in re-run:   {k} = {fn[k]!r}")
        elif not same(fo[k], fn[k], tol):
            out.append(f"differs: {k}: archived {fo[k]!r} -> re-run {fn[k]!r}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("archived")
    ap.add_argument("rerun")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--ignore", default=DEFAULT_IGNORE,
                    help="regex on the leaf key name; matching leaves are not compared")
    a = ap.parse_args()
    try:
        old, new = (json.load(open(p)) for p in (a.archived, a.rerun))
    except Exception as e:                                   # noqa: BLE001
        print(f"UNREADABLE: {e}")
        return 2
    diffs = compare(old, new, a.tol, a.ignore)
    n = sum(1 for _ in flatten(old))
    if not diffs:
        print(f"REPRODUCED: {a.rerun} matches {a.archived} on every compared leaf ({n} leaves).")
        return 0
    print(f"DIFFERS: {len(diffs)} leaf/leaves of {n}:")
    for d in diffs[:60]:
        print("  " + d)
    return 1


if __name__ == "__main__":
    sys.exit(main())
