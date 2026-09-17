#!/usr/bin/env python3
"""
preflight_b1.py — GO / NO-GO for the one expensive run, checked BEFORE it starts.

WHY THIS EXISTS
---------------
B1 is a single ~10 h unattended run and there is no time to do it twice. Everything it needs
is checkable in about a minute, and every check below corresponds to something that has
actually gone wrong in this project at least once: a reference built from the wrong matrix, a
gene space that could not be reproduced, an R package that was absent so a reimplementation ran
under a published package's name, a provenance sidecar silently overwritten, results replaced
without the previous run being archived.

It changes nothing. It reads, reports, and exits non-zero if the run would be unsafe.

    python scripts/preflight_b1.py                     # check the default (log) build
    python scripts/preflight_b1.py --matrix raw/X      # check the build B1 should use
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

OK, WARN, BAD = "  OK  ", " WARN ", " FAIL "
issues: list[str] = []
warns: list[str] = []


def say(state: str, what: str, detail: str = "") -> None:
    print(f"[{state}] {what}" + (f"\n         {detail}" if detail else ""))
    if state == BAD:
        issues.append(what)
    elif state == WARN:
        warns.append(what)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", choices=["X", "raw/X"], default="raw/X")
    args = ap.parse_args()
    print(f"PRE-FLIGHT for run_all.py --matrix {args.matrix}\n" + "=" * 68)

    # --- 1. inputs ------------------------------------------------------------------
    h5 = config.REFERENCE_DIR / "gbmap_core.h5ad"
    if h5.exists():
        say(OK, f"atlas present ({h5.stat().st_size / 1e9:.1f} GB)")
        try:
            import h5py
            with h5py.File(h5, "r") as f:
                have = [k for k in ("X", "raw/X") if k in f]
                if args.matrix in have:
                    say(OK, f"atlas carries {args.matrix!r}", f"present: {have}")
                else:
                    say(BAD, f"atlas has no {args.matrix!r}", f"present: {have}")
        except Exception as exc:                                   # noqa: BLE001
            say(WARN, "could not open the atlas to check its matrices", str(exc)[:120])
    else:
        say(BAD, "atlas missing", str(h5))

    if config.IVYGAP_FPKM_PATH.exists():
        say(OK, "Ivy GAP bulk present")
    else:
        say(BAD, "Ivy GAP bulk missing", str(config.IVYGAP_FPKM_PATH))

    # --- 2. the R packages, so no reimplementation runs under a published name -------
    try:
        from ivygap.deconv import r_bridge
        probes = {}
        for m, pkg in sorted(r_bridge.R_PACKAGES.items()):
            try:
                probes[m] = bool(r_bridge.probe_r_package(pkg))
            except Exception:                                      # noqa: BLE001
                probes[m] = False
        missing = sorted(k for k, v in probes.items() if not v)
        if not missing:
            say(OK, f"all {len(probes)} R packages load")
        else:
            say(WARN, f"{len(missing)} R package(s) will fall back to Python",
                f"{missing} — these run as REIMPLEMENTATIONS and are labelled as such, "
                f"which is correct behaviour but weakens the run")
    except Exception as exc:                                       # noqa: BLE001
        say(WARN, "could not probe the R bridge", str(exc)[:140])

    # --- 3. would this run destroy the current results? ------------------------------
    lb = config.RESULTS_DIR / "anatomic" / "acs_leaderboard.csv"
    if lb.exists():
        arch = ROOT / "results_archive"
        n_arch = len(list(arch.glob("*"))) if arch.exists() else 0
        say(WARN, "a results tree already exists and run_all OVERWRITES it",
            f"{lb.relative_to(ROOT)} holds the leaderboard every document cites, and "
            f"`check_doc_numbers.py` validates against it. A run with a DIFFERENT "
            f"reference produces a DIFFERENT leaderboard, so those documents would all "
            f"become stale at once. Stage 7 archives into results_archive/ "
            f"({n_arch} snapshot(s) there now) — but archive by hand first if you want "
            f"the current numbers trivially recoverable:\n"
            f"         cp -R results results_PRE_B1_$(date +%Y%m%dT%H%M)")
    else:
        say(OK, "no existing results tree to overwrite")

    # --- 4. the D16 chain: does raw/X actually reach Problem 1? ----------------------
    if args.matrix == "raw/X":
        lin = ROOT / "data/reference/gbmap_linear/profile.csv"
        say(OK if lin.exists() else WARN,
            "linear reference already built" if lin.exists()
            else "linear reference not pre-built (run_all builds its own; this is fine)")
        src = ROOT / "ivygap/bench/run_benchmark.py"
        gate = '"raw" in src or src == "counts"' in src.read_text()
        say(OK if gate else BAD,
            "Problem 1's truth is gated on a counts matrix" if gate
            else "the Problem 1 gate is missing from run_benchmark.py",
            "`truth_mrna` is populated only when library sizes come from counts; with "
            "--matrix raw/X they do, so Problem 1 becomes scoreable (D16)." if gate else "")

    # --- 5. disk and time -------------------------------------------------------------
    free = shutil.disk_usage(ROOT).free / 1e9
    say(OK if free > 20 else (WARN if free > 8 else BAD),
        f"{free:.0f} GB free",
        "a run writes results/, release/ and a results_archive/ snapshot; 20 GB is "
        "comfortable" if free <= 20 else "")

    # --- 6. the test suite is the cheapest proxy for "the code works" -----------------
    print("\nrunning the test suite (~5 min) — the cheapest proof the pipeline is sound ...")
    r = subprocess.run([sys.executable, "-m", "pytest", "tests/", "-q"],
                       capture_output=True, text=True, cwd=ROOT)
    tail = (r.stdout.strip().splitlines() or ["(no output)"])[-1]
    say(OK if r.returncode == 0 else BAD, f"test suite: {tail}")

    # --- verdict ----------------------------------------------------------------------
    print("\n" + "=" * 68)
    if issues:
        print(f"NO-GO — {len(issues)} blocking issue(s):")
        for i in issues:
            print(f"  * {i}")
        return 1
    print("GO.")
    if warns:
        print(f"\n{len(warns)} thing(s) to read before starting, none blocking:")
        for w in warns:
            print(f"  * {w}")
    print(f"\n  nohup python3 scripts/run_all.py --matrix {args.matrix} > b1.log 2>&1 &\n"
          f"\nThen leave it. Check with: tail -f b1.log")
    return 0


if __name__ == "__main__":
    sys.exit(main())
