"""
dedupe_verified.py -- delete a duplicate ONLY after proving it is byte-identical to a kept file.

User instruction (2026-10-01): "once the project's files have become too overwhelming you should
delete unnecessary files". The only files this script will delete are second copies: a file whose
content (sha256) equals that of a retained original, or of the original's decompressed stream when
the original is the .gz. Nothing that is the only copy of anything is touched. Each deletion is
appended to docs/DELETED_FILES.md with both hashes, the bytes freed and the reason, so it can be
audited and the file re-created exactly (gunzip -k, or a copy) if ever needed.

    python3 scripts/dedupe_verified.py            # dry run: hashes and reports, deletes nothing
    python3 scripts/dedupe_verified.py --delete   # deletes the verified duplicates
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ivygap import config  # noqa: E402

LOG = config.PROJECT_ROOT / "docs" / "DELETED_FILES.md"

#: (duplicate, kept original, original is gzip of the duplicate?, why it is unnecessary)
CANDIDATES = [
    (config.GSE182109_DIR / "Processed_matrix.mtx.gz.download" / "matrix.mtx",
     config.GSE182109_DIR / "Processed_matrix.mtx.gz.download" / "Processed_matrix.mtx.gz", True,
     "decompressed copy of the Abdelfattah processed matrix; every script reads the .gz"),
    (config.TCGA_DOWNLOADS_DIR / "project" / "TCGA.LGG.sampleMap_HumanMethylation450",
     config.TCGA_DOWNLOADS_DIR / "TCGA.LGG.sampleMap_HumanMethylation450.gz", True,
     "decompressed copy of the LGG HM450 matrix; the .gz equals the Xena server file byte for byte"),
    (config.VENDORED_DIR / "extra_data" / "GSE84465_GBM_All_data.csv",
     config.RAW_DIR / "darmanis_2017" / "GSE84465_GBM_All_data.csv.gz", True,
     "decompressed copy of the Darmanis counts; the pipeline reads data/raw/darmanis_2017/*.csv.gz"),
    (config.VENDORED_REPOS_DIR / "SCDC" / "neftel_data" / "GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",
     config.VENDORED_REPOS_DIR / "SCDC" / "neftel_data" / "IDHwtGBM.processed.SS2.logTPM.txt", False,
     "the Smart-seq2 matrix under its GEO file name; build_neftel_reference.py reads the SCP-named copy"),
]


def sha256_stream(fh) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    for b in iter(lambda: fh.read(8 << 20), b""):
        h.update(b)
        n += len(b)
    return h.hexdigest(), n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--delete", action="store_true")
    a = ap.parse_args()
    rows, freed = [], 0
    for dup, orig, gz, why in CANDIDATES:
        rel = lambda p: str(p.relative_to(config.PROJECT_ROOT))     # noqa: E731
        if not dup.exists():
            print(f"absent (already removed?): {rel(dup)}"); continue
        if not orig.exists():
            print(f"KEEP {rel(dup)}: its original {rel(orig)} is missing -- it may be the only copy")
            continue
        with open(dup, "rb") as fh:
            h_dup, n_dup = sha256_stream(fh)
        with (gzip.open(orig, "rb") if gz else open(orig, "rb")) as fh:
            h_orig, n_orig = sha256_stream(fh)
        same = h_dup == h_orig and n_dup == n_orig
        print(f"{'IDENTICAL' if same else 'DIFFERENT'}  {rel(dup)}  ({n_dup / 1e9:.2f} GB)  vs  "
              f"{rel(orig)}{' (decompressed)' if gz else ''}", flush=True)
        if same and a.delete:
            dup.unlink()
            freed += n_dup
            rows.append(f"| {dt.datetime.now():%Y-%m-%d %H:%M} | `{rel(dup)}` | {n_dup:,} | `{h_dup}` | "
                        f"`{rel(orig)}`{' (gzip of it)' if gz else ''} | {why} |")
    if rows:
        new = not LOG.exists()
        with open(LOG, "a") as fh:
            if new:
                fh.write("# Deleted files\n\nEvery file removed from this project, with proof it was a "
                         "second copy: its sha256 equalled that of the kept original (or of the "
                         "original's decompressed stream). Re-create any of them exactly with "
                         "`gunzip -k` or a copy of the original. Written by "
                         "`scripts/dedupe_verified.py`.\n\n"
                         "| when | deleted | bytes | sha256 (deleted = original content) | kept original | why unnecessary |\n"
                         "|---|---|---|---|---|---|\n")
            fh.write("\n".join(rows) + "\n")
        print(f"deleted {len(rows)} verified duplicate(s), freed {freed / 1e9:.2f} GB; logged in "
              f"{LOG.relative_to(config.PROJECT_ROOT)}")
    elif a.delete:
        print("nothing deleted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
