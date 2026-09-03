"""
download_ivygap.py — fetch the Ivy GAP public release from the Allen Institute.

Ivy GAP (the Ivy Glioblastoma Atlas Project) is openly downloadable; no application,
credential or data-use agreement is required, which is a real practical advantage
over TCGA's controlled tiers.

WHAT GETS DOWNLOADED
--------------------
1. `gene_expression_matrix_2014-11-25.zip` (~50 MB) containing
     fpkm_table.csv       genes x RNA-seq samples, FPKM
     rows-genes.csv       gene_id -> gene_symbol / entrez_id / chromosome
     columns-samples.csv  rna_well_id -> tumour, block, anatomic structure
2. `tumor_details.csv`, the per-tumour clinical table that carries survival.

The second file is the fragile one: the Allen Institute has moved it between
well-known-file ids over the years. This module tries the known locations and, if all
of them fail, says so plainly and lets the pipeline continue *without* survival rather
than inventing a clinical table. The survival module then reports BLOCKED. Silently
degrading to fabricated outcomes would be far worse than an honest missing input.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import requests

from ivygap import config

USER_AGENT = "ivygap-anatomic-deconvolution/1.0 (research pipeline)"

# The expression matrix. Both hosts serve the same archive; the api.brain-map.org
# well-known-file id is the more stable of the two.
EXPRESSION_URLS = [
    "http://api.brain-map.org/api/v2/well_known_file_download/305873915",
    "https://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873915",
    "https://glioblastoma.alleninstitute.org/static/download/gene_expression_matrix_2014-11-25.zip",
]

# The per-tumour clinical/survival table.
#
# HOW THE WORKING URL WAS FOUND (2026-09-03)
# -------------------------------------------
# All three URLs previously listed here are dead, and two of them fail in the worst
# possible way:
#
#   .../well_known_file_download/305873922   -> HTTP 404 (honest)
#   .../static/download/tumor_details.csv    -> HTTP 200, Content-Type: text/csv,
#                                               body is an HTML landing page
#
# The second is the defect `_assert_looks_like_csv` exists to catch. The remedy was to
# stop guessing well-known-file ids and read the portal's own download page
# (https://glioblastoma.alleninstitute.org/static/download.html), which links to a
# `/api/v2/gbm/` namespace that is not documented anywhere else. Those endpoints serve
# real CSV. The retired URLs are kept below the live one so that if the Allen Institute
# moves the namespace again, the failure report names everything that was tried.
TUMOR_DETAILS_URLS = [
    "https://glioblastoma.alleninstitute.org/api/v2/gbm/tumor_details.csv",
    "https://glioblastoma.alleninstitute.org/api/v2/well_known_file_download/305873922",
    "http://api.brain-map.org/api/v2/well_known_file_download/305873922",
    "https://glioblastoma.alleninstitute.org/static/download/tumor_details.csv",
]

# The portal's per-RNA-seq-sample metadata table. This is the authoritative source for
# protocol step 3's gate ("sample counts per structure reconcile against the portal's
# own documentation") and it is also the only file that carries the donor_id <-> tumor_id
# correspondence, without which the clinical table cannot be joined to the expression
# matrix at all: `tumor_details.csv` is keyed on donor_id, `columns-samples.csv` on
# tumor_id, and neither file contains both.
RNA_SEQ_SAMPLE_DETAILS_URLS = [
    "https://glioblastoma.alleninstitute.org/api/v2/gbm/rna_seq_samples_details.csv",
]

EXPECTED_MEMBERS = {"fpkm_table.csv", "rows-genes.csv", "columns-samples.csv"}


class DownloadFailed(RuntimeError):
    """Raised when every candidate URL for a required file failed."""


class NotACSV(RuntimeError):
    """The server returned 200 with a body that is not a CSV table."""


def _assert_looks_like_csv(body: bytes) -> None:
    """
    Reject an HTML page served with a 200.

    The Allen Institute serves a styled HTML landing page — not a 404 — for several of
    the tumour-details URLs that used to serve CSV. It downloads cleanly, saves cleanly,
    and then fails deep inside the clinical parser with a baffling error about a
    `<style>` column. A tolerant parser might even extract something from it. Check the
    content here, where the diagnosis is obvious and the remedy is a one-line message.
    """
    head = body[:2048].lstrip().lower()
    for marker in (b"<!doctype", b"<html", b"<style", b"<head", b"<body"):
        if head.startswith(marker) or marker in head[:512]:
            raise NotACSV(
                "the server returned an HTML page rather than a CSV table "
                "(the Allen Institute serves a landing page for retired download URLs)"
            )
    first = body.split(b"\n", 1)[0]
    if b"," not in first and b"\t" not in first:
        raise NotACSV(f"first line has no delimiter; not a table: {first[:120]!r}")


def _get(url: str, timeout: int = 300) -> bytes:
    resp = requests.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})
    resp.raise_for_status()
    if not resp.content:
        raise RuntimeError(f"empty body from {url}")
    return resp.content


def _try_urls(urls: list[str], what: str) -> bytes:
    errors = []
    for url in urls:
        try:
            print(f"  trying {url}")
            body = _get(url)
            print(f"  ok — {len(body):,} bytes")
            return body
        except Exception as exc:                     # noqa: BLE001 — report, then try next
            errors.append(f"    {url}\n      -> {type(exc).__name__}: {exc}")
    raise DownloadFailed(
        f"could not download {what}. Tried:\n" + "\n".join(errors)
    )


def download_expression(force: bool = False) -> Path:
    """
    Download and extract the Ivy GAP expression archive.

    Returns the directory holding the extracted CSVs. Idempotent: an existing, valid
    extraction is reused unless `force=True`.
    """
    config.ensure_dirs()
    out_dir = config.IVYGAP_FPKM_PATH.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    already = all((out_dir / m).exists() for m in EXPECTED_MEMBERS)
    if already and not force:
        print(f"expression archive already extracted at {out_dir} — skipping")
        return out_dir

    print("Downloading Ivy GAP expression matrix...")
    body = _try_urls(EXPRESSION_URLS, "the Ivy GAP expression matrix")

    config.IVYGAP_ZIP_PATH.write_bytes(body)
    print(f"  saved {config.IVYGAP_ZIP_PATH} "
          f"(sha256 {config.sha256_file(config.IVYGAP_ZIP_PATH)[:16]}...)")

    with zipfile.ZipFile(io.BytesIO(body)) as zf:
        names = {Path(n).name for n in zf.namelist()}
        missing = EXPECTED_MEMBERS - names
        if missing:
            raise DownloadFailed(
                f"archive downloaded but is missing expected members: {sorted(missing)}. "
                f"Archive actually contains: {sorted(names)}"
            )
        for info in zf.infolist():
            if info.is_dir():
                continue
            name = Path(info.filename).name
            # Flatten: the archive nests everything one level deep, and nothing
            # downstream benefits from reproducing that.
            with zf.open(info) as src, open(out_dir / name, "wb") as dst:
                dst.write(src.read())
    print(f"  extracted to {out_dir}")
    return out_dir


def download_tumor_details(force: bool = False) -> Path | None:
    """
    Download the per-tumour clinical table (carries survival).

    Returns the path, or None if every candidate URL failed. A None return is a
    legitimate, reportable state — it is not an error to be swallowed. The survival
    module checks for the file and reports BLOCKED rather than proceeding.
    """
    config.ensure_dirs()
    path = config.IVYGAP_TUMOR_DETAILS_PATH
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists() and not force:
        print(f"tumour details already present at {path} — skipping")
        return path

    print("Downloading Ivy GAP tumour details (clinical/survival)...")
    try:
        body = _try_urls(TUMOR_DETAILS_URLS, "the Ivy GAP tumour details table")
        _assert_looks_like_csv(body)
    except (DownloadFailed, NotACSV) as exc:
        print(f"\nWARNING: {exc}\n")
        print(
            "  Survival analysis will be reported as BLOCKED (missing required input).\n"
            "  To supply it manually: open https://glioblastoma.alleninstitute.org/,\n"
            "  use the Download menu to fetch the tumour details / clinical table, and\n"
            f"  save it as {path}\n"
            "  Nothing is imputed and no prognostic claim is made without it.\n"
        )
        return None

    path.write_bytes(body)
    print(f"  saved {path} (sha256 {config.sha256_file(path)[:16]}...)")
    return path


def download_rna_seq_sample_details(force: bool = False) -> Path | None:
    """
    Download the portal's per-RNA-seq-sample metadata table.

    Returns the path, or None if it could not be fetched. Like the clinical table this
    is allowed to be absent — the pipeline then cannot close protocol step 3's
    reconciliation gate and cannot join clinical data, and it says so rather than
    proceeding on an internally-derived count.
    """
    config.ensure_dirs()
    path = config.IVYGAP_RNA_SEQ_DETAILS_PATH
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists() and not force:
        print(f"rna-seq sample details already present at {path} — skipping")
        return path

    print("Downloading Ivy GAP rna_seq_samples_details (portal reconciliation)...")
    try:
        body = _try_urls(RNA_SEQ_SAMPLE_DETAILS_URLS,
                         "the Ivy GAP rna-seq sample details table")
        _assert_looks_like_csv(body)
    except (DownloadFailed, NotACSV) as exc:
        print(f"\nWARNING: {exc}\n")
        print("  Step 3's portal reconciliation and the donor_id <-> tumor_id join\n"
              "  both become unavailable. Counts will be reported as INTERNAL ONLY.\n")
        return None

    path.write_bytes(body)
    print(f"  saved {path} (sha256 {config.sha256_file(path)[:16]}...)")
    return path


def main(force: bool = False) -> int:
    download_expression(force=force)
    details = download_tumor_details(force=force)
    samples = download_rna_seq_sample_details(force=force)
    print("\nDownload step complete.")
    print(f"  expression      : OK")
    print(f"  clinical        : {'OK' if details else 'MISSING (survival will be BLOCKED)'}")
    print(f"  sample details  : {'OK' if samples else 'MISSING (portal reconciliation unavailable)'}")
    return 0


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download even if present")
    raise SystemExit(main(force=ap.parse_args().force))
