#!/usr/bin/env python3
"""
archive_run.py — freeze a completed run so it cannot be lost or quietly altered.

WHY THIS EXISTS
---------------
This project has lost a completed run twice. Once when a `--synthetic` fixture wrote to
the same paths as a real run; once when the test suite called `run_benchmark.run()`
against the real results tree while a pipeline was mid-flight. Both were fixed at the
cause — separate synthetic paths, a run lock, autouse output isolation in conftest — and
both fixes prevent a repeat of *that* mechanism. Neither makes a finished run durable.

An archive does. Once a run is archived it is a read-only snapshot with a hash for every
file, so a later run cannot overwrite it, and any drift is detectable rather than
assumed. The paper's numbers should cite an archive id, not `results/`, because
`results/` is by definition the thing the next run replaces.

    python scripts/archive_run.py --label "14-method run"
    python scripts/archive_run.py --verify 2026-09-05T1412
    python scripts/archive_run.py --list
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ivygap import config

ARCHIVE_ROOT = config.PROJECT_ROOT / "results_archive"

#: Files that make a run interpretable. A snapshot without these is not an archive of a
#: result, it is a pile of numbers.
REQUIRED = (
    "anatomic/acs_leaderboard.csv",
    "anatomic/anatomic_report.json",
    "anatomic/agreement_report.json",
    "anatomic/constraint_file.json",
)


def _hash_tree(root: Path) -> dict[str, str]:
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            out[str(p.relative_to(root))] = config.sha256_file(p)
    return out


def archive(label: str = "", results: Path | None = None) -> Path:
    src = results or config.RESULTS_DIR
    if not src.exists():
        raise SystemExit(f"no results tree at {src}")

    missing = [f for f in REQUIRED if not (src / f).exists()]
    if missing:
        raise SystemExit(
            f"refusing to archive an incomplete run — missing {missing}. "
            f"An archive is a citable object; a partial one is worse than none.")

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M")
    dest = ARCHIVE_ROOT / stamp
    if dest.exists():
        dest = ARCHIVE_ROOT / f"{stamp}-{datetime.now().strftime('%S')}"
    dest.mkdir(parents=True)

    shutil.copytree(src, dest / "results", dirs_exist_ok=True)
    for extra in ("RESULTS.md", "REGISTRATION.json"):
        p = config.PROJECT_ROOT / extra
        if p.exists():
            shutil.copy2(p, dest / extra)

    # provenance the archive can be judged on without opening a single CSV
    rep = json.loads((src / "anatomic/anatomic_report.json").read_text())
    agree = json.loads((src / "anatomic/agreement_report.json").read_text())
    primary = next((y for y in agree.get("by_yardstick", []) if y.get("n_methods")), {})
    lb = (src / "anatomic/acs_leaderboard.csv").read_text().splitlines()

    manifest = {
        "archived_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "label": label,
        "constraint_freeze_hash": rep.get("constraint_freeze_hash"),
        "n_samples_scored": rep.get("n_samples"),
        "n_tumors_scored": rep.get("n_tumors"),
        "n_samples_deconvolved": rep.get("n_samples_deconvolved"),
        "best_acs_method": rep.get("best_acs_method"),
        "best_acs": rep.get("best_acs"),
        "control_verdict": rep.get("control_verdict"),
        "primary_result": {k: primary.get(k) for k in
                           ("yardstick", "rho", "ci_low", "ci_high", "p_value",
                            "n_methods", "n_distinct", "meets_threshold",
                            "ci_excludes_zero")},
        "registration_state": (rep.get("registration") or {}).get("state"),
        "n_methods_on_leaderboard": max(len(lb) - 1, 0),
        "git_commit": _git_commit(),
        "files": _hash_tree(dest),
    }
    (dest / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str))

    # make it read-only so an accident has to be deliberate
    for p in dest.rglob("*"):
        if p.is_file():
            p.chmod(0o444)

    return dest


def _git_commit() -> str | None:
    import subprocess
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                           text=True, cwd=config.PROJECT_ROOT, timeout=10)
        return r.stdout.strip() or None
    except Exception:                                       # noqa: BLE001
        return None


def verify(stamp: str) -> int:
    dest = ARCHIVE_ROOT / stamp
    mpath = dest / "MANIFEST.json"
    if not mpath.exists():
        print(f"no archive at {dest}")
        return 1
    manifest = json.loads(mpath.read_text())
    recorded = manifest.get("files", {})
    now = _hash_tree(dest)
    now.pop("MANIFEST.json", None)

    changed = [f for f, h in recorded.items()
               if f != "MANIFEST.json" and now.get(f) != h]
    added = [f for f in now if f not in recorded]
    gone = [f for f in recorded if f != "MANIFEST.json" and f not in now]

    if not (changed or added or gone):
        print(f"{stamp}: INTACT — {len(recorded)} files match their recorded hashes")
        return 0
    print(f"{stamp}: ALTERED")
    for f in changed[:10]:
        print(f"  changed  {f}")
    for f in gone[:10]:
        print(f"  missing  {f}")
    for f in added[:10]:
        print(f"  added    {f}")
    return 1


def listing() -> int:
    if not ARCHIVE_ROOT.exists():
        print("no archives yet")
        return 0
    rows = sorted(ARCHIVE_ROOT.iterdir())
    for d in rows:
        m = d / "MANIFEST.json"
        if not m.is_file():
            continue
        man = json.loads(m.read_text())
        p = man.get("primary_result") or {}
        rho = p.get("rho")
        print(f"{d.name}  {man.get('n_methods_on_leaderboard','?')} methods  "
              f"best={man.get('best_acs_method','?')} {man.get('best_acs') or float('nan'):.3f}  "
              f"rho={'—' if rho is None else format(rho, '.3f')}  "
              f"{man.get('label','')}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", default="")
    ap.add_argument("--verify", metavar="STAMP")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()

    if a.list:
        return listing()
    if a.verify:
        return verify(a.verify)

    dest = archive(a.label)
    man = json.loads((dest / "MANIFEST.json").read_text())
    print(f"archived {len(man['files'])} files to {dest.relative_to(config.PROJECT_ROOT)}")
    print(f"  best ACS   : {man['best_acs_method']} {man['best_acs']}")
    print(f"  primary    : rho = {man['primary_result'].get('rho')}")
    print(f"  registration: {man['registration_state']}")
    print(f"  files are read-only; verify with --verify {dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
