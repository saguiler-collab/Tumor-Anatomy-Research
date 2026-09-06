#!/usr/bin/env python3
"""
fetch_allen.py — pull the Allen Human Brain Atlas down to disk, and query its API.

WHY THIS EXISTS
---------------
`Anatomy_Test.md` step 7 says one tissue is a case study and two is a method. Ivy GAP
gives nine to ten evaluable tumours and a weighted denominator of 65, and no analysis
inside Ivy GAP widens that. The Allen Human Brain Atlas is the protocol's own suggested
second tissue: bulk expression, annotated to a structure ontology, across donors.

WHAT IS ACTUALLY THERE (measured 2026-09-06, no API key needed)
---------------------------------------------------------------
  microarray   6 donors, ~1.6 GB total, ~400-1000 sampling sites per brain
  RNA-seq      2 donors, ~44 MB total, 240 samples over matched structures
  ontology     1,327 structures

READ THIS BEFORE PLANNING A SECOND-TISSUE RUN
---------------------------------------------
Downloading the data is the easy part and it is not the blocker. Three things are:

  1. THE CONSTRAINTS DO NOT TRANSFER. `constraints.py` encodes claims about GBM
     anatomic structures — CT, IT, LE, MVP, PAN. Normal brain has none of them. A
     second tissue needs its OWN pre-registered constraint file about which cell types
     are enriched in which brain regions, written and hashed before any Allen
     deconvolution is inspected. Reusing the GBM file is not an option, and writing a
     new one after seeing Allen output would destroy the pre-registration that the
     whole design rests on.

  2. THE ROSTER DOES NOT TRANSFER. Normal cortex is mostly neurons. This project's
     roster has no neuron column, and `reference_coverage.json` records why one could
     not simply be added: the GBmap atlas contains 22 neurons. A normal-brain run needs
     a normal-brain single-cell reference, not this one.

  3. PLATFORM. The 6-donor microarray set is the large one, but this pipeline is built
     for RNA-seq — CPM throughout, and an RNA-seq-derived frozen signature. Cross-
     platform deconvolution is a known hard problem, not a detail. The 2-donor RNA-seq
     set is platform-compatible and small; the 6-donor microarray set is large and is
     not. That trade is a study design decision, so it is stated here rather than made
     silently by whichever file this script downloads first.

USAGE
    python scripts/fetch_allen.py --list
    python scripts/fetch_allen.py --rnaseq              # 44 MB, platform-compatible
    python scripts/fetch_allen.py --microarray          # 1.6 GB, 6 donors
    python scripts/fetch_allen.py --query structures    # ontology via the RMA API
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config                                            # noqa: E402

API = "http://api.brain-map.org/api/v2/data"
DOWNLOAD = "https://human.brain-map.org/api/v2/well_known_file_download"
ALLEN_DIR = config.RAW_DIR / "allen_human_brain"

#: donor -> well-known-file id. Verified live 2026-09-06; sizes measured, not quoted.
MICROARRAY = {
    "H0351.2001": (178238387, 406.2),
    "H0351.2002": (178238373, 382.3),
    "H0351.1009": (178238359, 158.5),
    "H0351.1012": (178238316, 230.1),
    "H0351.1015": (178238266, 206.0),
    "H0351.1016": (178236545, 219.9),
}
RNASEQ = {
    "H0351.2001": (278447594, 22.0),
    "H0351.2002": (278448166, 22.4),
}


def _get_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def rma(model: str, criteria: str = "", num_rows: int = 50, **kw) -> list[dict]:
    """
    One RMA query. Returns `msg`, and raises if the service reports failure rather than
    handing back an empty list that reads like 'no results'.
    """
    q = f"{API}/{model}/query.json?num_rows={num_rows}"
    if criteria:
        q += "&criteria=" + urllib.parse.quote(criteria, safe="[]$,'@=")
    for k, v in kw.items():
        q += f"&{k}={urllib.parse.quote(str(v), safe='[]$,@=')}"
    out = _get_json(q)
    if not out.get("success"):
        raise RuntimeError(f"Allen API refused the query: {out.get('msg')}\n  {q}")
    return out["msg"]


def download(name: str, wkf_id: int, expect_mb: float, dest_dir: Path) -> Path:
    """
    Fetch one well-known file, verifying it is the size the catalogue claims.

    The size check is not decoration. A dead Allen endpoint earlier in this project
    returned HTTP 200 with an HTML error body under a text/csv content type, and a
    downloader that only checks the status code would have written that to disk and
    called it data.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{name}.zip"
    if dest.exists():
        print(f"  {name}: already present ({dest.stat().st_size / 1048576:.1f} MB)")
        return dest

    url = f"{DOWNLOAD}/{wkf_id}"
    print(f"  {name}: fetching {expect_mb:.1f} MB from {url}")
    tmp = dest.with_suffix(".partial")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        digest = hashlib.sha256()
        while chunk := r.read(1 << 20):
            fh.write(chunk)
            digest.update(chunk)

    got_mb = tmp.stat().st_size / 1048576
    if abs(got_mb - expect_mb) > max(2.0, 0.05 * expect_mb):
        tmp.unlink()
        raise RuntimeError(
            f"{name}: expected ~{expect_mb:.1f} MB, got {got_mb:.1f} MB. Refusing to "
            f"keep it — a wrong-sized file here is a truncated download or an error "
            f"page, and either one is worse than no file."
        )
    with open(tmp, "rb") as fh:
        head = fh.read(4)
    if head[:2] != b"PK":
        tmp.unlink()
        raise RuntimeError(f"{name}: not a zip (magic {head!r}); the endpoint returned "
                           f"something other than the archive.")
    tmp.rename(dest)
    print(f"    saved {dest.relative_to(config.PROJECT_ROOT)}  "
          f"sha256 {digest.hexdigest()[:16]}...")
    return dest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show what is available")
    ap.add_argument("--rnaseq", action="store_true", help="2 donors, ~44 MB")
    ap.add_argument("--microarray", action="store_true", help="6 donors, ~1.6 GB")
    ap.add_argument("--query", choices=["donors", "structures", "probes"],
                    help="hit the RMA API instead of downloading")
    ap.add_argument("--gene", default="GFAP", help="gene acronym for --query probes")
    args = ap.parse_args()

    if args.list:
        print(f"Allen Human Brain Atlas — no API key required\n")
        print(f"RNA-seq    (platform-compatible with this pipeline)")
        for n, (i, mb) in RNASEQ.items():
            print(f"   {n:12s} {mb:6.1f} MB   {DOWNLOAD}/{i}")
        print(f"\nMicroarray (larger, but a different platform — see module docstring)")
        for n, (i, mb) in MICROARRAY.items():
            print(f"   {n:12s} {mb:6.1f} MB   {DOWNLOAD}/{i}")
        print(f"\ntotal microarray {sum(v[1] for v in MICROARRAY.values())/1024:.2f} GB, "
              f"total RNA-seq {sum(v[1] for v in RNASEQ.values()):.0f} MB")
        return 0

    if args.query:
        if args.query == "donors":
            rows = rma("Donor", "products[id$eq2]", num_rows=20)
            print(f"{len(rows)} donors in the Human Brain Microarray product:")
            for d in rows:
                print(f"  {d.get('name')}  age_id={d.get('age_id')}  "
                      f"sex={d.get('sex_full_name') or d.get('gender_id')}")
        elif args.query == "structures":
            rows = rma("Structure", "[ontology_id$eq1]", num_rows=10)
            print("first 10 of the human structure ontology:")
            for s in rows:
                print(f"  {s.get('acronym'):12s} depth={s.get('depth'):<3} "
                      f"{s.get('name')}")
        else:
            rows = rma("Probe",
                       f"[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],"
                       f"gene[acronym$eq'{args.gene}']", num_rows=10)
            print(f"{len(rows)} probes for {args.gene}:")
            for p in rows:
                print(f"  {p.get('name')}  id={p.get('id')}")
        return 0

    targets = {}
    if args.rnaseq:
        targets.update({f"rnaseq_{k}": v for k, v in RNASEQ.items()})
    if args.microarray:
        targets.update({f"microarray_{k}": v for k, v in MICROARRAY.items()})
    if not targets:
        ap.print_help()
        return 1

    print(f"downloading {len(targets)} file(s) into "
          f"{ALLEN_DIR.relative_to(config.PROJECT_ROOT)}")
    for name, (wkf, mb) in targets.items():
        download(name, wkf, mb, ALLEN_DIR)
    print("\ndone. NOTE: having the data does not unblock a second-tissue run — the "
          "constraint file, the roster and the platform all have to be dealt with "
          "first. See this module's docstring.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
