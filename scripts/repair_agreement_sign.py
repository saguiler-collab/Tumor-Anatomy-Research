"""Regenerate `agreement_report.json` after the yardstick-direction defect.

THE DEFECT. `run_anatomic` passed `higher_is_better={k: False for k in yardsticks}` to the
agreement test, under the comment "every yardstick here is an error metric". That is true of
`synthetic_mixtures`, whose score is mean absolute error, and false of `absolute_purity`,
whose score is `spearman_vs_purity` -- a correlation between estimated tumour content and
DNA-measured purity, where HIGHER IS BETTER. `test_agreement` negates the truth vector when
told lower is better, so the orthogonal yardstick's rho was written with its sign inverted.

WHAT IT DID AND DID NOT AFFECT. |rho|, the p-value, the method count and the pre-registered
threshold test (rho >= 0.60, failed either way) are unchanged. The SIGN and the bootstrap
confidence interval are not. `results/yardstick_agreement.json`, written by a different
script that does its own Spearman, was correct throughout: +0.0810. The two artefacts
disagreed in sign on the same statistic, which is how this was found.

WHY A SCRIPT AND NOT A RE-RUN. Only the agreement STAGE consumed the wrong flag. Nothing
upstream of it -- deconvolution, the ACS leaderboard, the per-constraint scores -- reads it,
so every input to this stage is unchanged and a full re-run would recompute hours of
identical numbers. This calls the same production functions on the same inputs:
`_load_yardsticks`, `agreement.run_all_yardsticks` and `agreement.write_report`. It does not
hand-edit the artefact, and it recomputes every yardstick rather than only the one that was
wrong, so a second wrong direction could not survive it.

Usage:  python3 scripts/repair_agreement_sign.py [--check]
        --check reports what would change and writes nothing.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ivygap import config                                      # noqa: E402
from ivygap.anatomic import agreement                          # noqa: E402
from ivygap.anatomic.run_anatomic import _load_yardsticks      # noqa: E402
from ivygap.bench.real_yardsticks import HIGHER_IS_BETTER      # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report the change without writing")
    args = ap.parse_args()

    anat = config.PROJECT_ROOT / "results" / "anatomic"
    lb_path, rep_path = anat / "acs_leaderboard.csv", anat / "agreement_report.json"
    if not lb_path.exists():
        print(f"ABORT: {lb_path} absent; nothing to repair.")
        return 1

    lb = pd.read_csv(lb_path).set_index("method")
    if "comparable" not in lb.columns:
        print("ABORT: leaderboard has no `comparable` column; cannot reproduce the "
              "method set run_anatomic used.")
        return 1

    # EXACTLY the restriction run_anatomic applies: a partial-coverage method's ACS comes
    # from a different denominator over a different subset of constraints, so it may not
    # enter the rank correlation.
    acs_scores = {m: float(r["acs"]) for m, r in lb.iterrows() if bool(r["comparable"])}
    yardsticks, _ = _load_yardsticks()

    unknown = sorted(set(yardsticks) - set(HIGHER_IS_BETTER))
    if unknown:
        print(f"ABORT: no declared direction for {unknown}.")
        return 1

    table = agreement.run_all_yardsticks(
        acs_scores, yardsticks,
        higher_is_better={k: HIGHER_IS_BETTER[k] for k in yardsticks})

    before = {}
    if rep_path.exists():
        for y in json.loads(rep_path.read_text()).get("by_yardstick", []):
            before[y.get("yardstick")] = y

    print(f"{'yardstick':<26} {'was':>10} {'now':>10}   {'n':>3}  change")
    print("-" * 74)
    changed = []
    for y, row in table.iterrows():
        old = before.get(y, {}).get("rho")
        new = row["rho"]
        o = "n/a" if old is None or not np.isfinite(old) else f"{old:+.4f}"
        n = "n/a" if not np.isfinite(new) else f"{new:+.4f}"
        note = ""
        if old is not None and np.isfinite(old) and np.isfinite(new):
            if abs(old - new) < 1e-9:
                note = "unchanged"
            elif abs(old + new) < 1e-9:
                note = "SIGN INVERTED"
                changed.append(y)
            else:
                note = "DIFFERS"
                changed.append(y)
        print(f"{y:<26} {o:>10} {n:>10}   {int(row['n_methods']):>3}  {note}")

    print()
    print(f"direction used: {json.dumps({k: HIGHER_IS_BETTER[k] for k in yardsticks})}")

    if args.check:
        print("\n--check: nothing written.")
        return 0

    if rep_path.exists():
        # Keep the wrong artefact rather than destroying the record of the defect.
        bak = rep_path.with_suffix(".json.pre_sign_fix")
        if not bak.exists():
            shutil.copy2(rep_path, bak)
            print(f"preserved the pre-fix artefact at {bak.relative_to(config.PROJECT_ROOT)}")

    report = agreement.write_report(table, rep_path)
    # Say in the artefact itself that it was regenerated, and why.
    report["regenerated"] = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "by": "scripts/repair_agreement_sign.py",
        "reason": "the agreement stage was given higher_is_better=False for every "
                  "yardstick, which is wrong for `absolute_purity` (a Spearman "
                  "correlation, higher is better) and inverted its rho's sign. Only the "
                  "agreement stage read that flag; all inputs are unchanged.",
        "yardsticks_whose_value_changed": changed,
        "direction_used": {k: HIGHER_IS_BETTER[k] for k in yardsticks},
    }
    rep_path.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {rep_path.relative_to(config.PROJECT_ROOT)}")

    # THE ORDERING GATE. `registration_status.json` is written LAST by `run_anatomic`,
    # deliberately: it is the only point at which every result file carries this run's
    # mtime. Rewriting any artefact after it leaves it judging an earlier state of the
    # tree, which `tests/test_registration_ordering_gate.py` detects and which is how a
    # stale "every result file postdates the registration" verdict gets reprinted as
    # though it were current. So it is RECOMPUTED here -- by the same authoritative
    # function the pipeline calls, never edited.
    from ivygap.anatomic import registration
    reg = registration.status()
    (anat / "registration_status.json").write_text(reg.to_json())
    print(f"recomputed registration_status.json -> {reg.state}")
    if getattr(reg, "results_predating_registration", None):
        print(f"  WARNING: {len(reg.results_predating_registration)} file(s) predate the "
              f"registration; investigate rather than regenerate.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
