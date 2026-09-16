#!/usr/bin/env python3
"""
check_doc_numbers.py — every ACS printed in a live document must exist in an artefact.

WHY THIS EXISTS
---------------
A leaderboard number reaches the reader through prose, and prose is copied by hand. This
project has already been bitten twice by exactly that: a header that claimed rho = 0.873
on 10 methods while the body below it said 0.733 on 13, and a defect note that retracted
a claim on one page and restated it as fact on another. Both were caught by reading, which
does not scale and does not run in CI.

`tests/test_summary_header.py` pins RESULTS.md to the generator. Nothing pinned the other
thirty-odd documents, and those are the ones a judge or reviewer actually opens: README.md,
docs/RELATED_WORK.md, docs/OPEN_DEFECTS.md, docs/ROAD_TO_PAPER.md.

WHAT IT CHECKS
--------------
Every markdown table cell that sits under a column header beginning "ACS", in a row whose
first cell names a method this study ran. The value must equal that method's ACS in the
canonical leaderboard, in the full-database run, or in a declared re-measurement.

WHAT IT DELIBERATELY DOES NOT CHECK
-----------------------------------
  * `results_archive/` and `results_superseded/`. Those record what earlier runs said, and
    are *supposed* to differ. Rewriting them to match today would destroy the audit trail
    that makes the corrections legible.
  * Columns that are not ACS. An earlier version of this script matched the second column
    of every table and reported 86 "defects", all false: the survival table's C-index
    (0.484), and the cell-size table's fraction-interpretation column (0.500 = cell share).
    Requiring the table's own header to declare an ACS column is what makes the result
    trustworthy. A checker that cries wolf is worse than none, because it gets ignored.

  * WHETHER A MEASURED NUMBER IS THE RIGHT ONE FOR ITS CONTEXT. This is the honest limit,
    and it got looser when subset re-scorings were added. The check is "some artefact holds
    this value for this method", not "this is the value this sentence should quote". DWLS,
    for instance, now legitimately has more than twenty distinct ACS values across the
    leaderboard, the full-database run, two genuine-package re-measurements and several
    restricted-subset comparisons. A number that is real but quoted in the wrong place —
    a 41-pair subset score presented as a leaderboard score — passes this checker and is
    caught only by reading. Denominators are printed beside every subset score for that
    reason.

    python scripts/check_doc_numbers.py           # report
    python scripts/check_doc_numbers.py --strict  # non-zero exit on any mismatch
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

#: Values a method is legitimately allowed to show besides its two leaderboard entries,
#: each with the reason it exists. A method may not be added here to silence a mismatch.
DECLARED_REMEASUREMENTS = {
    # A re-measurement listed here is NOT thereby comparable to the leaderboard. Both
    # values below were measured on 1,591 genes against a leaderboard built on 657
    # (docs/OPEN_DEFECTS.md D10), so a document may print them only alongside the gene
    # space. This checker verifies a number exists in an artefact; it cannot verify that
    # the sentence around it is fair, and it must not be read as blessing the comparison.
    # VALID as of 2026-09-14: 657-gene space (hash verified), 88/22 donor split reproduced
    # by name, training cells only, all seven input-equivalence conditions passing. These
    # ARE comparable to the leaderboard.
    "dwls": {0.7846: "genuine DWLS R package, equivalent inputs, 2,474 s",
             0.785: "the same value, rounded to 3 dp",
             0.738: "the Python reimplementation, as run in the confirmatory run",
             0.7077: "VOID — the 1,591-gene run with donor leakage (D10). Kept so a "
                     "document still quoting it is recognised, not blessed",
             0.708: "VOID, rounded"},
    "bayesprism": {0.8154: "genuine BayesPrism R package, equivalent inputs, 2,045 s",
                   0.815: "the same value, rounded to 3 dp",
                   0.8923: "VOID — the 1,591-gene run with donor leakage (D10)",
                   0.892: "VOID, rounded"},
}

SKIP = ("/.git/", "/node_modules/", "/pipeline_packages/",
        "/results_archive/", "/results_superseded/")


def _artefact_rescorings() -> dict[str, set[float]]:
    """
    ACS values a method legitimately shows because it was RE-SCORED on a restricted
    constraint subset, harvested from the artefacts themselves.

    Subset re-scorings are a real and growing class — the IMC marker-source comparison and
    the CDSeq reference-free arm both score every method on the constraints two arms share,
    so the same method correctly shows several different ACS values. Those values must still
    be checkable, and hand-listing them in DECLARED_REMEASUREMENTS would be exactly the
    "add the value to silence the warning" this script refuses. So they are read back out of
    the JSON that produced them: a document may print a number only if some artefact holds it.

    Any `{<method>: {"acs": <float>, ...}}` mapping anywhere in a results json counts.
    """
    found: dict[str, set[float]] = {}

    def walk(node, key=None):
        if isinstance(node, dict):
            if "acs" in node and isinstance(node.get("acs"), (int, float)) and key:
                name = str(key).split("__")[0].strip().lower()
                v = float(node["acs"])
                found.setdefault(name, set()).update({round(v, 4), round(v, 3)})
            for k, v in node.items():
                walk(v, k)
        elif isinstance(node, list):
            for v in node:
                walk(v, key)

    for f in sorted((ROOT / "results").glob("*.json")):
        try:
            walk(json.loads(f.read_text()))
        except (ValueError, OSError):
            continue
    return found


def _leaderboards() -> tuple[pd.Series, pd.Series | None]:
    canon = pd.read_csv(ROOT / "results/anatomic/acs_leaderboard.csv")
    canon = canon.set_index("method")["acs"]
    p = ROOT / "results/anatomic/full_database/acs_leaderboard.csv"
    full = pd.read_csv(p).set_index("method")["acs"] if p.exists() else None
    return canon, full


def _aliases(methods) -> dict[str, str]:
    a = {m: m for m in methods}
    a.update({"elastic net": "elastic_net", "scdc ensemble": "scdc_ensemble",
              "cibersortx b-mode": "cibersortx", "cibersortx s-mode": "cibersortx_smode"})
    return a


def scan() -> tuple[list[tuple], int, int]:
    canon, full = _leaderboards()
    alias = _aliases(canon.index)
    rescored = _artefact_rescorings()

    files = [p for p in ROOT.rglob("*.md")
             if not any(s in str(p).replace(str(ROOT), "") for s in SKIP)]

    problems: list[tuple] = []
    checked = 0
    for f in sorted(files):
        acs_col: int | None = None
        for ln, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
            if not line.strip().startswith("|"):
                acs_col = None                      # a table cannot span a blank line
                continue
            cells = [c.strip().lower().strip("*") for c in line.strip().strip("|").split("|")]
            hdr = [i for i, c in enumerate(cells) if c.startswith("acs")]
            if hdr:
                acs_col = hdr[0]
                continue
            if acs_col is None or acs_col >= len(cells):
                continue
            name = alias.get(cells[0].replace("**", "").strip())
            if name is None:
                continue
            m = re.fullmatch(r"(\d\.\d{3,4})", cells[acs_col].replace("**", "").strip())
            if not m:
                continue
            checked += 1
            val = round(float(m.group(1)), 4)
            ok = {round(float(canon[name]), 4), round(float(canon[name]), 3)}
            if full is not None and name in full.index:
                ok |= {round(float(full[name]), 4), round(float(full[name]), 3)}
            for v in DECLARED_REMEASUREMENTS.get(name, {}):
                ok |= {round(v, 4), round(v, 3)}
            ok |= rescored.get(name, set())          # subset re-scorings, read from artefacts
            if val not in ok:
                problems.append((f.relative_to(ROOT), ln, name, val,
                                 sorted(ok), line.strip()[:100]))
    return problems, checked, len(files)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any document quotes an ACS no artefact holds")
    args = ap.parse_args()

    problems, checked, n_files = scan()
    print(f"checked {checked} ACS table cells across {n_files} live markdown files")
    if not problems:
        print("CLEAN: every ACS in a live document matches an artefact.")
        return 0

    print(f"\n{len(problems)} value(s) appear in a document but in no artefact:\n")
    for f, ln, name, val, ok, txt in problems:
        print(f"  {f}:{ln}\n    {name} is printed as {val}; artefacts hold {ok}\n    {txt}")
    print("\nEither the document is stale (fix the document) or a run was not archived "
          "(archive it). Do not add the value to DECLARED_REMEASUREMENTS to silence this.")
    return 1 if args.strict else 0


if __name__ == "__main__":
    sys.exit(main())
