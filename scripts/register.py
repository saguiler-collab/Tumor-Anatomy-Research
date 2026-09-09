#!/usr/bin/env python3
"""
register.py — record the public registration receipt, and check it.

WHY THIS IS A SCRIPT AND NOT AN INSTRUCTION TO EDIT A FILE
----------------------------------------------------------
`REGISTRATION.json` has to carry the freeze hash of the constraint file exactly. Typing
it by hand is how it ends up off by a character, and a wrong hash reads as
HASH_MISMATCH, which is the state the protocol treats as fatal: it says the constraints
changed after registration and every result under them is void. So the hash is read from
the constraint file rather than accepted as input, and cannot be passed in.

    python scripts/register.py --status
    python scripts/register.py --url https://osf.io/xxxxx --registry OSF
    python scripts/register.py --url https://osf.io/xxxxx --registry OSF \
        --registered-utc 2026-09-09T14:03:00Z

`--registered-utc` defaults to now. Pass the timestamp OSF shows on the registration if
it differs, because the ordering check compares it against every result file's mtime.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap.anatomic import constraints as K            # noqa: E402
from ivygap.anatomic import registration as R           # noqa: E402


def _print_status() -> int:
    s = R.status()
    print(f"state              : {s.state}")
    print(f"live freeze hash   : {s.freeze_hash}")
    if s.recorded_hash:
        print(f"recorded hash      : {s.recorded_hash}")
    if s.registry or s.url:
        print(f"registry / url     : {s.registry} {s.url}")
    if s.registered_utc:
        print(f"registered (UTC)   : {s.registered_utc}")
    if s.results_predating_registration:
        n = len(s.results_predating_registration)
        print(f"results predating  : {n} file(s) — these were produced BEFORE the "
              f"registration")
        for f in s.results_predating_registration[:10]:
            print(f"    {f}")
        if n > 10:
            print(f"    ... and {n - 10} more")
    print()
    print(s.verdict)
    return 0 if s.state == "REGISTERED" else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", action="store_true",
                    help="report the current state and exit")
    ap.add_argument("--url", help="the public registration URL (e.g. https://osf.io/xxxxx)")
    ap.add_argument("--registry", default="OSF",
                    help="registry name, default OSF")
    ap.add_argument("--registered-utc",
                    help="ISO-8601 UTC timestamp of the registration; defaults to now")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing REGISTRATION.json")
    args = ap.parse_args()

    if args.status or not args.url:
        if not args.url and not args.status:
            ap.print_help()
            print()
        return _print_status()

    path = R.REGISTRATION_PATH
    if path.exists() and not args.force:
        print(f"{path.name} already exists. Re-run with --force to overwrite it, but "
              f"read it first — overwriting a genuine receipt loses the timestamp that "
              f"the whole pre-registration claim depends on.")
        return 1

    if not args.url.startswith(("http://", "https://")):
        print(f"--url must be a public URL a reader can open; got {args.url!r}")
        return 1

    stamp = args.registered_utc or dt.datetime.now(dt.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")

    receipt = {
        "registry": args.registry,
        "url": args.url,
        "registered_utc": stamp,
        # Read from the constraint file, never accepted as input.
        "constraint_freeze_hash": K.freeze_hash(),
        "what_this_certifies": (
            "The seven anatomic constraints hashed here were registered publicly at the "
            "URL above on the date above. The hash covers the canonical payload of the "
            "claims — ids, cell types, structures, directions, kinds and weights — not "
            "the source file, so reformatting constraints.py does not change it while "
            "altering a claim does."),
        "what_this_does_not_certify": (
            "That the analyses postdate the registration. Runs archived before this date "
            "are pilot work and are labelled as such; the ordering check reports any "
            "result file older than the registration."),
    }
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path.name}\n")
    return _print_status()


if __name__ == "__main__":
    raise SystemExit(main())
