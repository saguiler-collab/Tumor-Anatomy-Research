"""
fetch_cptac_gbm.py -- download the CPTAC glioblastoma tumours that have BOTH bulk RNA-seq and
single-nucleus RNA-seq from the same cryopulverised material (Wang et al., Cancer Cell 2021,
doi:10.1016/j.ccell.2021.01.006), from the NCI Genomic Data Commons (project CPTAC-3; open access).

Fetched, per case (only what the analysis needs; the 2.6 GB of loom files are not):
  snrna/      CellRanger "10x Filtered Counts" matrices (MEX, tar.gz)
  snrna_gdc/  GDC's own Seurat cluster table (seurat.analysis.tsv)
  bulk/       STAR "Counts" gene-expression quantification (the same format as the TCGA bulk here)
  methylation/ methylation beta values
Every file is checked against the md5 GDC records, and data/external/cptac_gbm/MANIFEST.json lists
case, file id, name, size and md5. Re-running skips files already present with a matching md5.

    python3 scripts/fetch_cptac_gbm.py            # list, then download
    python3 scripts/fetch_cptac_gbm.py --list     # list only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

try:
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:  # pragma: no cover
    CTX = None
UA = {"User-Agent": "TumorAnatomyResearch/1.0 (CPTAC GBM fetch)"}
API = "https://api.gdc.cancer.gov"
OUT = config.CPTAC_GBM_DIR
KINDS = {  # folder: (experimental_strategy, data_type, workflow_type or file-name filter)
    "snrna": ("scRNA-Seq", "Gene Expression Quantification", "CellRanger - 10x Filtered Counts"),
    "snrna_gdc": ("scRNA-Seq", "Single Cell Analysis", "seurat.analysis.tsv"),
    "bulk": ("RNA-Seq", "Gene Expression Quantification", "STAR - Counts"),
    "methylation": ("Methylation Array", "Methylation Beta Value", None),
}


def api(path: str, params: dict) -> dict:
    url = f"{API}/{path}?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120, context=CTX) as r:
        return json.load(r)


def snrna_cases() -> list[str]:
    f = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "CPTAC-3"}},
        {"op": "=", "content": {"field": "cases.primary_site", "value": "Brain"}},
        {"op": "=", "content": {"field": "experimental_strategy", "value": "scRNA-Seq"}}]}
    hits = api("files", {"filters": json.dumps(f), "fields": "cases.submitter_id", "size": 500, "format": "json"})
    return sorted({c["submitter_id"] for h in hits["data"]["hits"] for c in h.get("cases", [])})


def listing(cases: list[str]) -> list[dict]:
    rows = []
    for kind, (es, dt, wf) in KINDS.items():
        f = {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.submitter_id", "value": cases}},
            {"op": "=", "content": {"field": "experimental_strategy", "value": es}},
            {"op": "=", "content": {"field": "data_type", "value": dt}}]}
        hits = api("files", {"filters": json.dumps(f), "size": 2000, "format": "json",
                             "fields": "file_id,file_name,file_size,md5sum,access,analysis.workflow_type,"
                                       "cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type"})
        for h in hits["data"]["hits"]:
            w = (h.get("analysis") or {}).get("workflow_type")
            if kind == "snrna_gdc":
                if h["file_name"] != wf:
                    continue
            elif wf and w != wf:
                continue
            c = h["cases"][0]
            rows.append({"kind": kind, "case": c["submitter_id"], "file_id": h["file_id"], "file_name": h["file_name"],
                         "size": h["file_size"], "md5": h["md5sum"], "access": h["access"], "workflow": w,
                         "samples": [s.get("submitter_id") for s in c.get("samples", [])],
                         "sample_types": sorted({s.get("sample_type") for s in c.get("samples", [])})})
    return rows


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def download(row: dict) -> Path:
    dst = OUT / row["kind"] / row["case"] / f"{row['file_id']}__{row['file_name']}"
    if dst.exists() and md5(dst) == row["md5"]:
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".part")
    req = urllib.request.Request(f"{API}/data/{row['file_id']}", headers=UA)
    with urllib.request.urlopen(req, timeout=600, context=CTX) as r, open(tmp, "wb") as fo:
        for b in iter(lambda: r.read(1 << 20), b""):
            fo.write(b)
    got = md5(tmp)
    if got != row["md5"]:
        tmp.unlink()
        raise SystemExit(f"BLOCKED: md5 mismatch for {row['file_name']} ({got} != {row['md5']})")
    tmp.rename(dst)
    return dst


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    cases = snrna_cases()
    rows = listing(cases)
    print(f"{len(cases)} cases with snRNA-seq; files: " + ", ".join(
        f"{k} {sum(r['kind'] == k for r in rows)} ({sum(r['size'] for r in rows if r['kind'] == k) / 1e9:.2f} GB)"
        for k in KINDS))
    if any(r["access"] != "open" for r in rows):
        print("BLOCKED: a listed file is not open access"); return 2
    if a.list:
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    for i, r in enumerate(rows, 1):
        p = download(r)
        r["local"] = str(p.relative_to(config.PROJECT_ROOT))
        print(f"  [{i}/{len(rows)}] {r['kind']:11s} {r['case']} {r['file_name'][:60]} md5 ok", flush=True)
    (OUT / "MANIFEST.json").write_text(json.dumps({"source": "NCI GDC, project CPTAC-3 (open access)",
                                                   "paper": "Wang LB et al., Cancer Cell 39:509 (2021), doi:10.1016/j.ccell.2021.01.006",
                                                   "cases": cases, "files": rows}, indent=2))
    print(f"wrote {(OUT / 'MANIFEST.json').relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
