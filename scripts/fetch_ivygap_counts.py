#!/usr/bin/env python3
"""
fetch_ivygap_counts.py — get Ivy GAP READ COUNTS, which unblocks the count-based methods.

Why this exists
---------------
The published Ivy GAP expression matrix is normalised FPKM, and every count-based method
(CDSeq above all, but the assumption runs through MuSiC, SCDC and BayesPrism too) wants
counts. `docs/CORRECTIONS_REGISTRATION.md` C9 asserted no public counts exist, on the strength
of GEO GSE107559 stating "Raw data not provided for this record". **That was wrong**, and the
correction is recorded in C9 itself: GEO is one distribution channel, and the Allen Institute
portal is another. The portal publishes, per sample, the RSEM `genes.results` file:

    gene_id  transcript_id(s)  length  effective_length  expected_count  TPM  FPKM

`expected_count` is what CDSeq needs. 25,874 genes per sample, ~1.5 MB each, 270 samples.

It also publishes anonymized BAMs (414 MB each, ~112 GB for all 270). Those are NOT needed:
RSEM has already done the counting, and re-counting from BAMs would substitute my pipeline for
the authors' without adding information.

What `expected_count` IS, stated because it matters
--------------------------------------------------
It is RSEM's posterior expected read count per gene, so it is FRACTIONAL (45.56, not 46). It is
not a raw integer count and is not presented as one. Any model needing integers gets an
explicitly declared rounding step, applied here once and recorded, rather than each method
rounding differently in private.

    python scripts/fetch_ivygap_counts.py            # all 270
    python scripts/fetch_ivygap_counts.py --limit 5  # smoke test

Writes per-sample files under `data/raw/ivygap_counts/` and a provenance record beside them.
Resumable: a sample already on disk with the right gene count is skipped.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config                                          # noqa: E402

MANIFEST = ROOT / "FPKM" / "raw_fpkm.csv"
OUT = config.RAW_DIR / "ivygap_counts"
EXPECTED_GENES = 25873          # 25,874 lines minus the header; MEASURED, and checked not
                                # assumed: a short file is a truncated download, not a sample
                                # with fewer genes.


def fetch(url: str, dest: Path, tries: int = 3) -> bytes:
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                body = r.read()
            if body.count(b"\n") < EXPECTED_GENES:
                raise ValueError(f"only {body.count(chr(10).encode())} lines")
            dest.write_bytes(body)
            return body
        except Exception as exc:                                   # noqa: BLE001
            if k == tries - 1:
                raise
            print(f"      retry {k + 1}: {type(exc).__name__} {exc}")
            time.sleep(2 * (k + 1))
    raise AssertionError("unreachable")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=0, help="fetch only the first N (smoke test)")
    ap.add_argument("--pause", type=float, default=0.3,
                    help="seconds between requests; this is someone else's server")
    args = ap.parse_args()

    if not MANIFEST.exists():
        print(f"BLOCKED: {MANIFEST} not found. It is the portal's "
              f"'un-normalized gene-level FPKM and TPM values for each sample' csv.")
        return 2
    rows = list(csv.DictReader(MANIFEST.open()))
    if args.limit:
        rows = rows[:args.limit]
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(rows)} samples from {MANIFEST.relative_to(ROOT)} -> "
          f"{OUT.relative_to(ROOT)}/")

    recs, skipped, fetched = [], 0, 0
    for i, r in enumerate(rows, 1):
        well = r["rna_well"]
        dest = OUT / f"{well}.genes.results"
        if dest.exists() and dest.read_bytes().count(b"\n") >= EXPECTED_GENES:
            body = dest.read_bytes()
            skipped += 1
        else:
            print(f"  [{i:3d}/{len(rows)}] {well}  {r['structure_acronym']}")
            body = fetch(r["file_download_link"], dest)
            fetched += 1
            time.sleep(args.pause)
        recs.append({
            "rna_well": well,
            "specimen_name": r["specimen_name"],
            "tumor_name": r["tumor_name"],
            "structure_acronym": r["structure_acronym"],
            "rnaseq_total_reads": int(r["rnaseq_total_reads"]),
            "percent_aligned_mrna": float(r["rnaseq_percent_reads_aligned_to_mrna"]),
            "url": r["file_download_link"],
            "sha256": hashlib.sha256(body).hexdigest(),
            "n_lines": body.count(b"\n"),
        })

    prov = {
        "what_this_is": "Per-sample RSEM gene-level results for Ivy GAP, carrying "
                        "expected_count, effective_length, TPM and FPKM.",
        "why": "Unblocks the count-based methods, CDSeq in particular. Supersedes the claim "
               "in CORRECTIONS_REGISTRATION C9 that no public counts exist.",
        "source": "Allen Institute Ivy GAP portal, per-sample well_known_file_download URLs "
                  "listed in the portal's un-normalized FPKM/TPM manifest.",
        "expected_count_is_fractional": "RSEM posterior expected counts, not raw integers. "
                                        "Rounding, where a model needs integers, is declared "
                                        "at the point of use rather than done silently here.",
        "bams_not_used": "The portal also publishes anonymized BAMs (414 MB each, ~112 GB for "
                         "270). Not fetched: RSEM already counted, and re-counting would "
                         "substitute this pipeline's choices for the authors'.",
        "n_samples": len(recs), "n_fetched_now": fetched, "n_already_present": skipped,
        "samples": recs,
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"\n{len(recs)} samples on disk ({fetched} fetched now, {skipped} already there)")
    print(f"wrote {(OUT / 'provenance.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
