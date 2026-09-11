#!/usr/bin/env python3
"""
SUPERSEDED — do not use, and do not cite its output.

Replaced by scripts/ish_constraint_check.py. This version pulled ISH values one
experiment at a time from the Allen API and was wrong three ways: it used 72 of 632
experiments (11.4%); it pooled expression energy across DIFFERENT sub-block
specimens, conflating tumour-to-tumour variation with structure differences; and it
pooled markers that contradict each other. The correct input was already on disk —
gene_expression_details.csv has one row per (gene x sub-block) with all five
structures as columns, so the comparison is paired inside one block by construction.
"""

"""
fetch_ivygap_ish.py — check the pre-registered constraints against in-situ hybridization.

WHAT THIS IS FOR
----------------
The seven constraints in `ivygap/anatomic/constraints.py` are justified by neuropathology.
Several of them already *cite* Ivy GAP's ISH as evidence — C5's evidence line reads "Ivy
GAP ISH of myeloid markers by structure" — but nothing in this project ever measured it.

Ivy GAP quantified ISH signal per anatomic structure and serves it through the same API
the RNA-seq comes from. So the constraints can be checked against an independent
measurement on the same tissue, made by the people who took the images, without
downloading 130 GB of pictures: the quantification is a few MB of JSON.

If ISH agrees with a constraint, that constraint rests on measurement rather than on
expert expectation. If it disagrees, that is a finding about the constraint and must be
reported as one.

WHAT THIS IS NOT
----------------
**Not a protocol analysis, and not a way to justify the constraint file.** The constraints
are frozen and hashed (`2d1fb47c…`) and registered at <https://osf.io/dm2t8>. They do not
move because ISH agrees or disagrees. This is exploratory, belongs under the
registration's "Other planned analysis", and must be labelled that way wherever it is
reported.

**Not a second deconvolution.** ISH measures transcript abundance in a tissue neighbourhood,
not cell-type composition. A marker gene being higher in one structure is consistent with
more cells of that type, more transcript per cell, or both. It supports an ordinal claim;
it does not measure a fraction.

THE CIRCULARITY THAT MUST BE AVOIDED
------------------------------------
ISH is how Ivy GAP assigned structure labels to the 148 "cluster" samples this project
excludes from scoring as circular. Using ISH here is legitimate because the constraints
are scored on the **H&E** samples, whose labels come from histology, not from expression —
different samples, different labelling basis. Any write-up must say so explicitly, or a
careful reader will assume the circularity was overlooked rather than avoided.

MARKER GENES ARE DECLARED BEFORE THE QUERY RUNS
-----------------------------------------------
The obvious way to fake this analysis is to pick markers until the constraints pass. So
`MARKERS` below is fixed in this file, chosen as canonical textbook markers for each cell
type, and every gene queried is reported whether or not it agrees. There is no
"representative marker" selection step and no dropping of inconvenient genes.

    python scripts/fetch_ivygap_ish.py --list-markers
    python scripts/fetch_ivygap_ish.py --out results/ish_constraint_check.json
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.parse
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

API = "http://api.brain-map.org/api/v2/data/query.json"

#: Ivy GAP ISH structures (graph 15) mapped onto this project's constraint vocabulary.
#: Verified against the ontology 2026-09-10. The ISH survey uses bare abbreviations; the
#: RNA-seq survey uses the "-reference-histology" suffixed forms, and they are different
#: structure records for the same anatomy.
ISH_STRUCTURE = {
    "LE":    {"id": 122767057, "name": "Leading Edge"},
    "IT":    {"id": 122767063, "name": "Infiltrating Tumor"},
    "CT":    {"id": 122767073, "name": "Cellular Tumor"},
    "MVP":   {"id": 122767079, "name": "Microvascular proliferation"},   # CTmvp
    "PAN":   {"id": 122767097, "name": "Pseudopalisading cells around necrosis"},  # CTpan
}
ID_TO_REGION = {v["id"]: k for k, v in ISH_STRUCTURE.items()}

#: Canonical markers per roster cell type. DECLARED BEFORE ANY QUERY. Chosen because they
#: are the textbook markers for each population, not because of anything they do here.
#: Every one is reported, agreeing or not.
MARKERS: dict[str, list[str]] = {
    # Endothelium. ESM1 is the endothelial marker Ivy GAP's own API documentation uses as
    # its worked example; CD34 and KDR(VEGFR2) are canonical endothelial markers; CAV1 is
    # endothelial-enriched but less specific and is reported separately rather than dropped.
    "Endothelial":          ["ESM1", "CD34", "KDR", "CAV1"],
    # Myeloid. CD163 is the M2-like macrophage marker. LAPTM5 is pan-leukocyte rather than
    # myeloid-specific, so it is included but must be read as weaker evidence.
    "Macrophage_Microglia": ["CD163", "LAPTM5"],
    # Tumour cell density. Every available marker here is a glioma-STATE marker rather
    # than a density marker, which is a real limitation stated in the report.
    "Tumor":                ["SOX2", "PTPRZ1", "EGFR", "CD44", "BIRC5", "TOP2A"],
}

#: Cell types whose constraints CANNOT be checked against this panel, recorded so the
#: silence is visible rather than mistaken for agreement.
NOT_CHECKABLE: dict[str, str] = {
    "Oligodendrocyte": (
        "C2 (Oligodendrocyte: LE > CT) cannot be checked. The Ivy GAP ISH panel contains "
        "no myelin or oligodendrocyte-lineage marker — MBP, PLP1, MOG, MAG, CNP, SOX10, "
        "MOBP and CLDN11 all return zero experiments. OLIG2 is present and was "
        "DELIBERATELY NOT USED: in glioma it is expressed by the tumour cells themselves "
        "and by OPCs, so scoring C2 on OLIG2 would measure tumour content and report it "
        "as oligodendrocyte content, which would corrupt the very constraint it claims "
        "to check."),
}

#: Which constraints each cell type speaks to, so the report can be read against the
#: constraint file rather than requiring the reader to hold both in their head.
CONSTRAINT_FOR = {
    "Endothelial":          ["C3 (MVP > CT)", "C4 (MVP is the maximum)"],
    "Macrophage_Microglia": ["C5 (PAN > LE)", "C6 (MVP > CT)"],
    "Tumor":                ["C1 (CT > LE)", "C7 (LE < IT < CT)"],
}


def _get(criteria: str, num_rows: int = 200, tries: int = 3) -> list[dict]:
    """One RMA query. curl rather than urllib: this machine's Python has no trust store."""
    # Two encoding rules, both learned by getting them wrong, and both failing the same
    # dangerous way — an EMPTY BODY rather than an error, so every gene silently looked
    # like it had no ISH experiment and all seven constraints reported "not evaluable".
    # A clean null result caused entirely by URL escaping.
    #
    #   * the apostrophe must stay literal: %27ESM1%27 is not matched;
    #   * the brackets must be ENCODED: a literal [ ] returns an empty body.
    #
    # Verified against the service: with brackets encoded, ESM1 returns 32 experiments.
    url = (f"{API}?criteria={urllib.parse.quote(criteria, safe=chr(39) + '$,:@=')}"
           f"&num_rows={num_rows}")
    for attempt in range(tries):
        out = subprocess.run(["curl", "-s", "--max-time", "60", url],
                             capture_output=True, text=True).stdout
        try:
            d = json.loads(out)
        except json.JSONDecodeError:
            time.sleep(1.5 * (attempt + 1))
            continue
        if not d.get("success"):
            raise RuntimeError(f"Allen API refused: {str(d.get('msg'))[:200]}\n  {url}")
        return d.get("msg", [])
    raise RuntimeError(f"no parseable response after {tries} tries: {url}")


def section_data_sets(gene: str) -> list[int]:
    """Every ISH experiment for one gene in the Glioblastoma product (id 17)."""
    rows = _get(f"model::SectionDataSet,rma::criteria,products[id$eq17],"
                f"genes[acronym$eq'{gene}']")
    return [r["id"] for r in rows if r.get("id")]


def structure_energy(sds_id: int) -> dict[str, float]:
    """expression_energy per constraint region for one ISH experiment."""
    rows = _get(f"model::StructureUnionize,rma::criteria,"
                f"section_data_set[id$eq{sds_id}],rma::include,structure", num_rows=60)
    out = {}
    for r in rows:
        sid = (r.get("structure") or {}).get("id") or r.get("structure_id")
        region = ID_TO_REGION.get(sid)
        if region and r.get("expression_energy") is not None:
            out[region] = float(r["expression_energy"])
    return out


def check_ordering(values: dict[str, float], higher: str, lower: str):
    """Does ISH put `higher` above `lower`? None when either region is missing."""
    if higher not in values or lower not in values:
        return None
    return values[higher] > values[lower]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/ish_constraint_check.json")
    ap.add_argument("--list-markers", action="store_true")
    ap.add_argument("--max-experiments", type=int, default=6,
                    help="ISH experiments per gene to average over")
    args = ap.parse_args()

    if args.list_markers:
        print("Markers declared in this file, before any query:\n")
        for ct, genes in MARKERS.items():
            print(f"  {ct:22s} {', '.join(genes)}")
            print(f"  {'':22s} speaks to: {', '.join(CONSTRAINT_FOR[ct])}\n")
        return 0

    report = {
        "what_this_is": (
            "Ivy GAP ISH expression energy per anatomic structure, checked against the "
            "pre-registered constraints. EXPLORATORY — the constraints are frozen and "
            "hashed and do not move on this result."),
        "constraint_freeze_hash": None,
        "markers_declared_before_query": MARKERS,
        "structures": {k: v["name"] for k, v in ISH_STRUCTURE.items()},
        "per_gene": [],
        "per_cell_type": {},
    }
    try:
        from ivygap.anatomic import constraints as K
        report["constraint_freeze_hash"] = K.freeze_hash()
    except Exception:
        pass

    for cell_type, genes in MARKERS.items():
        agg: dict[str, list[float]] = defaultdict(list)
        for gene in genes:
            try:
                sds = section_data_sets(gene)
            except RuntimeError as exc:
                report["per_gene"].append({"gene": gene, "cell_type": cell_type,
                                           "error": str(exc)[:180]})
                continue
            if not sds:
                report["per_gene"].append({"gene": gene, "cell_type": cell_type,
                                           "n_experiments": 0,
                                           "note": "no ISH experiment in Ivy GAP"})
                print(f"  {gene:8s} ({cell_type}): no ISH experiment")
                continue

            per_exp = []
            for sid in sds[:args.max_experiments]:
                try:
                    e = structure_energy(sid)
                except RuntimeError:
                    continue
                if e:
                    per_exp.append(e)
            if not per_exp:
                continue

            mean = {r: sum(d[r] for d in per_exp if r in d) / max(1, sum(r in d for d in per_exp))
                    for r in ISH_STRUCTURE if any(r in d for d in per_exp)}
            for r, v in mean.items():
                agg[r].append(v)
            report["per_gene"].append({
                "gene": gene, "cell_type": cell_type,
                "n_experiments_found": len(sds),
                "n_experiments_used": len(per_exp),
                "mean_expression_energy": {k: round(v, 4) for k, v in mean.items()},
            })
            print(f"  {gene:8s} ({cell_type:22s}) {len(per_exp)}/{len(sds)} exp  "
                  + "  ".join(f"{k}={v:.2f}" for k, v in sorted(mean.items())))

        pooled = {r: sum(v) / len(v) for r, v in agg.items() if v}
        report["per_cell_type"][cell_type] = {
            "pooled_mean_expression_energy": {k: round(v, 4) for k, v in pooled.items()},
            "n_genes_contributing": len({g["gene"] for g in report["per_gene"]
                                         if g.get("cell_type") == cell_type
                                         and g.get("mean_expression_energy")}),
            "speaks_to": CONSTRAINT_FOR[cell_type],
        }

    # Constraint-by-constraint verdict, on the pooled marker signal.
    endo = report["per_cell_type"].get("Endothelial", {}).get("pooled_mean_expression_energy", {})
    mac = report["per_cell_type"].get("Macrophage_Microglia", {}).get("pooled_mean_expression_energy", {})
    tum = report["per_cell_type"].get("Tumor", {}).get("pooled_mean_expression_energy", {})

    verdicts = {
        "C1": {"claim": "Tumor: CT > LE", "ish_agrees": check_ordering(tum, "CT", "LE")},
        "C3": {"claim": "Endothelial: MVP > CT", "ish_agrees": check_ordering(endo, "MVP", "CT")},
        "C4": {"claim": "Endothelial: MVP is the maximum",
               "ish_agrees": (max(endo, key=endo.get) == "MVP") if endo else None},
        "C5": {"claim": "Macrophage_Microglia: PAN > LE", "ish_agrees": check_ordering(mac, "PAN", "LE")},
        "C6": {"claim": "Macrophage_Microglia: MVP > CT", "ish_agrees": check_ordering(mac, "MVP", "CT")},
        "C7": {"claim": "Tumor: LE < IT < CT",
               "ish_agrees": (tum.get("LE", 0) < tum.get("IT", 0) < tum.get("CT", 0))
                             if {"LE", "IT", "CT"} <= set(tum) else None},
    }
    verdicts["C2"] = {"claim": "Oligodendrocyte: LE > CT",
                      "ish_agrees": None,
                      "why_not_evaluable": NOT_CHECKABLE["Oligodendrocyte"]}
    report["constraint_verdicts"] = verdicts
    report["not_checkable_against_this_panel"] = NOT_CHECKABLE

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    print(f"\nwrote {out}\n")
    print("Constraint vs ISH — EXPLORATORY, the constraints do not move on this:")
    for cid, v in verdicts.items():
        mark = {True: "agrees", False: "DISAGREES", None: "not evaluable"}[v["ish_agrees"]]
        print(f"  {cid}  {v['claim']:38s} ISH {mark}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
