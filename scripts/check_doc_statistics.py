"""Verify the paper's headline statistics against the specific artefact field each comes from.

WHY THE FIRST VERSION OF THIS FILE WAS WORTHLESS, recorded because the failure is instructive.
It harvested every float in the results tree -- 755,490 distinct values -- and asked whether a
quoted number appeared anywhere in that set. At that density essentially every 3- and 4-decimal
number matches something, so the check reported CLEAN when a rho was replaced with an invented
0.4242 AND when a sign was flipped. It passed because it could not fail. Its own docstring had
warned that harvesting every float "would make every number matchable and the check worthless",
and it did it anyway.

WHAT THIS VERSION DOES INSTEAD. Each headline statistic is pinned to the exact artefact and
field that produces it. A document may quote that statistic only with the current value. The
superseded value of the same statistic is ALSO recorded, and quoting it outside a region
explicitly marked as superseded is reported -- that is what staleness actually looks like, and
it is the shape OPEN_DEFECTS D21 and the log-`X` carryover both had.

This trades coverage for teeth. It checks fewer numbers than the harvesting version pretended
to, and unlike that version it fails when something is wrong.

    python3 scripts/check_doc_statistics.py
    python3 scripts/check_doc_statistics.py --strict
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

DOCS = ["README.md", "RESULTS.md", "docs/ENDPOINT.md", "docs/FIGURES.md",
        "docs/MANUSCRIPT.md", "docs/SUPPLEMENTARY.md", "docs/METHOD_LIMITATIONS.md"]

#: Regions a document declares as another build's or as superseded. Inside them a superseded
#: value is correct and the current value would be wrong, so the polarity flips.
SUPERSEDED_MARKERS = (
    "Reference build for this entire section: GBmap `X`",
    "SUPERSEDED TABLE — values below are from the log-`X` arm.",
    "⚠ SUPERSEDED",
    "log `X` (superseded)",
)


def _j(rel: str):
    p = ROOT / "results" / rel
    return json.loads(p.read_text()) if p.exists() else None


def _dig(obj, path: str):
    cur = obj
    for key in path.split("."):
        if cur is None:
            return None
        if isinstance(cur, list):
            try:
                cur = cur[int(key)]
                continue
            except (ValueError, IndexError):
                return None
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def _yardstick(name: str, field: str = "rho"):
    d = _j("anatomic/agreement_report.json") or {}
    for row in d.get("by_yardstick", []):
        if row.get("yardstick") == name:
            return row.get(field)
    return None


#: (label, current, superseded-or-None, decimals, required-context regex)
#:
#: THE CONTEXT REGEX IS LOAD-BEARING. Without it a bare digit search fires on any statistic
#: that happens to share digits: the model-fit FLOOR is -0.0071 and the true-parity rho is
#: +0.0071, so two lines about the floor were reported as stale parity values. A number alone
#: does not identify which quantity a sentence is about.
def canonical() -> list[tuple[str, float | None, float | None, int, str]]:
    fz, h5 = _j("yardstick_agreement.json") or {}, _j("yardstick_agreement_h5ad.json") or {}
    mx = _j("matrix_arm_comparison.json") or {}
    CORR = r"(?:rho|ρ|Spearman|correlation|tau)"
    return [
        ("pseudobulk yardstick rho",
         _dig(fz, "arm_1_pseudobulk.rho"), 0.7501, 4, CORR),
        ("ACS vs purity, frozen, all methods",
         _dig(fz, "arm_2_absolute_purity.all_methods.spearman"), 0.0698, 4,
         rf"(?:{CORR}|purity|ABSOLUTE)"),
        ("ACS vs purity, frozen, excl degenerate",
         _dig(fz, "arm_2_absolute_purity.excluding_degenerate.spearman"), -0.1044, 4,
         rf"(?:{CORR}|degenerate|purity)"),
        ("ACS vs purity, true parity, all methods",
         _dig(h5, "arm_2_absolute_purity.all_methods.spearman"), -0.0156, 4,
         rf"(?:{CORR}|parity|h5ad|purity)"),
        ("ACS vs purity, true parity, excl degenerate",
         _dig(h5, "arm_2_absolute_purity.excluding_degenerate.spearman"), 0.0071, 4,
         rf"(?:{CORR}|parity|h5ad)(?!.*floor)"),
        ("absolute_purity rho in the agreement report",
         _yardstick("absolute_purity"), -0.0810, 4, rf"(?:{CORR}|purity)"),
        ("synthetic_mixtures rho in the agreement report",
         _yardstick("synthetic_mixtures"), None, 4, CORR),
        ("Kendall tau, log X vs raw/X leaderboards",
         mx.get("kendall_tau"), None, 4, r"(?:tau|Kendall)"),
        ("Spearman, log X vs raw/X leaderboards",
         mx.get("spearman"), None, 4, r"Spearman"),
    ]


def _regions(text: str) -> list[bool]:
    """Per-line: is this line inside a superseded/other-build region?"""
    out, flag = [], False
    for line in text.split("\n"):
        if any(m in line for m in SUPERSEDED_MARKERS):
            flag = True
        elif line.startswith("## "):
            flag = False
        out.append(flag)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true", help="exit non-zero on any problem")
    args = ap.parse_args()

    stats = canonical()
    print("CANONICAL HEADLINE STATISTICS, read from the artefact field each comes from\n")
    missing = []
    for label, cur, sup, nd, _ctx in stats:
        if cur is None:
            missing.append(label)
            print(f"  {label:<46} ARTEFACT FIELD ABSENT")
        else:
            s = f"(superseded {sup:+.4f})" if sup is not None else ""
            print(f"  {label:<46} {cur:+.{nd}f}  {s}")
    print()

    problems = []
    for rel in DOCS:
        p = ROOT / rel
        if not p.exists():
            continue
        text = p.read_text(errors="replace")
        lines = text.split("\n")
        sup_region = _regions(text)
        for label, cur, sup, nd, ctx in stats:
            if sup is None or cur is None:
                continue
            ctx_re = re.compile(ctx, re.I)
            # SIGN-ONLY DIFFERENCE gets its own matcher. When the superseded and current
            # values share a magnitude and differ only in sign -- exactly OPEN_DEFECTS D21,
            # where -0.0810 was published for +0.0810 -- a magnitude search cannot tell them
            # apart, and the "current value is alongside it" skip below matched the FLIPPED
            # value and waved it through. Caught by the negative control, not by reading.
            sign_only = abs(abs(cur) - abs(sup)) < 10 ** -nd and (cur < 0) != (sup < 0)
            if sign_only:
                # require the explicit wrong sign, immediately before the digits
                # NOTE THE GROUP. Written as "-|−" the alternation binds looser than the
                # concatenation, so the pattern meant "any hyphen OR the number" and matched
                # every line containing a dash. Found by the negative control again.
                pat = re.compile(r"(?<![\d.])" + (r"(?:-|−)" if sup < 0 else r"\+")
                                 + re.escape(f"{abs(sup):.{nd}f}"))
            else:
                pat = re.compile(re.escape(f"{abs(sup):.{nd}f}"))
            for i, line in enumerate(lines):
                if sup_region[i] or not pat.search(line):
                    continue
                if not ctx_re.search(line):
                    continue          # the line shares digits but is about something else
                # INLINE LABELLING IS CORRECT, not stale. A sentence that gives the current
                # value AND names the old one as superseded is exactly what a corrected
                # document should look like; flagging it would push the author to delete the
                # record instead of labelling it.
                if not sign_only and re.search(r"supersed|log-`X`|log `X`", line, re.I):
                    continue
                if not sign_only and re.search(re.escape(f"{abs(cur):.{nd}f}"), line):
                    continue          # the current value is right here alongside it
                # the two may round to the same printed string; then there is nothing to say
                if not sign_only and abs(abs(cur) - abs(sup)) < 10 ** -nd:
                    continue
                problems.append((rel, i + 1, label, sup, cur, line.strip()[:110]))

    if problems:
        print(f"{len(problems)} superseded value(s) quoted OUTSIDE a marked region:\n")
        for rel, ln, label, sup, cur, txt in problems:
            print(f"  {rel}:{ln}")
            print(f"    {label}: document prints {abs(sup):.4f}; current is {cur:+.4f}")
            print(f"    {txt}\n")

    if missing:
        print(f"{len(missing)} canonical statistic(s) have no artefact field — the check "
              f"cannot verify them: {', '.join(missing)}\n")

    if not problems and not missing:
        print("CLEAN: no document quotes a superseded headline statistic outside a marked "
              "region, and every canonical statistic resolves to an artefact field.")
        return 0
    return 1 if args.strict else 0


if __name__ == "__main__":
    raise SystemExit(main())
