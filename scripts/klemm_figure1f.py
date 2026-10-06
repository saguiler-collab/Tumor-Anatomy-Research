"""
klemm_figure1f.py -- prespecified/klemm_t_vs_b.md: T against B in IDH-mutant and IDH-wildtype gliomas,
measured by flow cytometry (Klemm et al. 2020, Cell 181:1643, Figure 1F).

EXPLORATORY. Figure 1F gives, per group, the mean of each immune population as % of CD45+ cells, as
horizontal stacked bars (myeloid left of 0%, lymphoid right). The figure is VECTOR graphics, so each
bar segment is a filled path with exact coordinates. This script:
  1. interprets page 3's content stream (CTM, fill colours, path bounding boxes) -- no rendering;
  2. takes the colour of every population from the figure's own legend swatches (row-major order as
     printed: MG, MDM, Neutrophils, CD14low/CD16+ Mono, CD14+/CD16+ Mono / CD16- Gran., iMC, DC,
     CD4+ T, Treg / CD8+ T, DNT, B, NK);
  3. calibrates against the panel's own gridlines (every 25%) and its 0% (where myeloid meets lymphoid);
  4. checks itself against numbers the paper PRINTS in its text -- melanoma-BrM CD8+ T = 33.01% of
     CD45+, all-BrM lymphocytes = 46.23% (n 13 / 16 / 8 by primary, Table S1) -- and that every bar sums
     to 100%; any failure blocks the reading;
  5. applies the pre-declared reading: T = CD4+ + Treg + CD8+ + DNT against B, per glioma group.

    python3 scripts/klemm_figure1f.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from pypdf import PdfReader
from pypdf.generic import ContentStream

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

OUT = config.RESULTS_DIR / "klemm_t_vs_b.json"
PAGE = 3                                       # Figure 1 (its legend is on page 4)
PANEL_F = (255.0, 552.0, 80.0, 225.0)          # x0, x1, y0, y1 in page points (panel F's box)
BAR_HEIGHT = 8.75
LEGEND = ["MG", "MDM", "Neutrophils", "CD14low/CD16+ Mono", "CD14+/CD16+ Mono",
          "CD16- Gran.", "iMC", "DC", "CD4+ T", "Treg",
          "CD8+ T", "DNT", "B", "NK"]
MYELOID = LEGEND[:8]
LYMPHOID = LEGEND[8:]
T_POPS = ["CD4+ T", "Treg", "CD8+ T", "DNT"]
ROWS = ["Non-tumor", "Glioma IDH mut", "Glioma IDH wt", "BrM breast", "BrM lung", "BrM melanoma"]
N_BRM = {"BrM breast": 13, "BrM lung": 16, "BrM melanoma": 8}          # Table S1
PRINTED = {"melanoma_CD8_pct_CD45": 33.01, "brm_lymphocytes_pct_CD45": 46.23}   # paper's own text
RESOLUTION = 0.01                                                      # % of CD45+ (see controls)


def _mul(a, b):
    return [a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3], a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3],
            a[4]*b[0] + a[5]*b[2] + b[4], a[4]*b[1] + a[5]*b[3] + b[5]]


def shapes(pdf: Path, pno: int) -> list[dict]:
    """Every filled or stroked path on a page: bounding box in page points, and its colour."""
    r = PdfReader(str(pdf))
    ops = ContentStream(r.pages[pno - 1].get_contents(), r).operations
    ctm, fill, stroke, stack, pts, out = [1, 0, 0, 1, 0, 0], None, None, [], [], []
    num = ("cm", "re", "m", "l", "c", "v", "y", "rg", "g", "k", "sc", "scn", "RG", "G", "K", "SC", "SCN")
    for operands, op in ops:
        op = op.decode() if isinstance(op, bytes) else str(op)
        v = [float(x) for x in operands if hasattr(x, "as_numeric") or isinstance(x, (int, float))] if op in num else None
        if op == "q":
            stack.append((list(ctm), fill, stroke))
        elif op == "Q":
            ctm, fill, stroke = stack.pop()
        elif op == "cm":
            ctm = _mul(v, ctm)
        elif op in ("rg", "g", "k", "sc", "scn"):
            fill = tuple(round(x, 4) for x in v)
        elif op in ("RG", "G", "K", "SC", "SCN"):
            stroke = tuple(round(x, 4) for x in v)
        elif op == "re":
            x, y, w, h = v
            pts += [(ctm[0]*a + ctm[2]*b + ctm[4], ctm[1]*a + ctm[3]*b + ctm[5])
                    for a, b in ((x, y), (x + w, y), (x, y + h), (x + w, y + h))]
        elif op in ("m", "l", "c", "v", "y"):
            pts += [(ctm[0]*v[i] + ctm[2]*v[i+1] + ctm[4], ctm[1]*v[i] + ctm[3]*v[i+1] + ctm[5])
                    for i in range(0, len(v), 2)]
        elif op in ("f", "F", "f*", "B", "B*", "b", "b*", "S", "s"):
            if pts:
                xs, ys = [p[0] for p in pts], [p[1] for p in pts]
                kind = "stroke" if op in ("S", "s") else "fill"
                out.append({"kind": kind, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys),
                            "colour": stroke if kind == "stroke" else fill})
            pts = []
        elif op == "n":
            pts = []
    return out


def measure(S: list[dict], scale_divisor: float = 25.0) -> dict:
    """Per-row % of CD45+ for every population. `scale_divisor` is the % between gridlines (25); it is
    a parameter only so the tests can show that a wrong calibration fails the printed-number controls."""
    x0, x1, y0, y1 = PANEL_F
    inF = [s for s in S if x0 <= s["x0"] and s["x1"] <= x1 and y0 <= s["y0"] and s["y1"] <= y1]
    grid = sorted({round(s["x0"], 3) for s in inF if s["kind"] == "stroke" and abs(s["x1"] - s["x0"]) < 0.3
                   and (s["y1"] - s["y0"]) > 5 and s["colour"] and min(s["colour"]) > 0.5})
    gaps = [b - a for a, b in zip(grid, grid[1:])]
    if len(grid) < 5 or max(gaps) - min(gaps) > 0.01:
        raise SystemExit(f"BLOCKED: panel F gridlines not found as an even grid: {grid}")
    pt_per_pct = (sum(gaps) / len(gaps)) / scale_divisor
    swatch = [s for s in inF if s["kind"] == "fill" and 88 <= s["y0"] <= 116 and 5.5 <= s["x1"] - s["x0"] <= 6.0
              and 5.5 <= s["y1"] - s["y0"] <= 6.0 and s["colour"] and max(s["colour"]) > 0.05]
    swatch.sort(key=lambda s: (-round(s["y0"]), s["x0"]))
    if len(swatch) != len(LEGEND):
        raise SystemExit(f"BLOCKED: expected {len(LEGEND)} legend swatches, found {len(swatch)}")
    key = {tuple(s["colour"]): name for s, name in zip(swatch, LEGEND)}
    if len(key) != len(LEGEND):
        raise SystemExit("BLOCKED: two legend entries share a colour")
    bars = [s for s in inF if s["kind"] == "fill" and abs((s["y1"] - s["y0"]) - BAR_HEIGHT) < 0.3]
    centres = sorted({round((s["y0"] + s["y1"]) / 2, 1) for s in bars}, reverse=True)
    if len(centres) != len(ROWS):
        raise SystemExit(f"BLOCKED: expected {len(ROWS)} bar rows, found {len(centres)}")
    rows = {}
    for name, yc in zip(ROWS, centres):
        segs = sorted([s for s in bars if abs((s["y0"] + s["y1"]) / 2 - yc) < 0.5], key=lambda s: s["x0"])
        labels = [key.get(tuple(s["colour"])) for s in segs]
        if labels != LEGEND:
            raise SystemExit(f"BLOCKED: row {name} segments are not in legend order: {labels}")
        zero = segs[len(MYELOID)]["x0"]                       # where lymphoid segments begin
        rows[name] = {lab: round((s["x1"] - s["x0"]) / pt_per_pct, 4) for lab, s in zip(labels, segs)}
        rows[name]["_myeloid_total"] = round((zero - segs[0]["x0"]) / pt_per_pct, 4)
        rows[name]["_lymphoid_total"] = round((segs[-1]["x1"] - zero) / pt_per_pct, 4)
    return {"gridlines_pt": grid, "pt_per_percent": pt_per_pct, "rows": rows}


def controls(m: dict) -> dict:
    rows = m["rows"]
    mel_cd8 = rows["BrM melanoma"]["CD8+ T"]
    brm_lym = sum(rows[k]["_lymphoid_total"] * n for k, n in N_BRM.items()) / sum(N_BRM.values())
    totals = {k: round(v["_myeloid_total"] + v["_lymphoid_total"], 3) for k, v in rows.items()}
    c = {"melanoma_CD8": {"measured": round(mel_cd8, 3), "printed": PRINTED["melanoma_CD8_pct_CD45"],
                          "pass": abs(mel_cd8 - PRINTED["melanoma_CD8_pct_CD45"]) < 0.01},
         "brm_lymphocytes_weighted": {"measured": round(brm_lym, 3), "printed": PRINTED["brm_lymphocytes_pct_CD45"],
                                      "pass": abs(brm_lym - PRINTED["brm_lymphocytes_pct_CD45"]) < 0.01},
         "bar_totals_100": {"totals": totals, "pass": all(abs(t - 100) < 0.05 for t in totals.values())}}
    c["all_pass"] = all(v["pass"] for v in c.values() if isinstance(v, dict))
    return c


def reading(row: dict) -> dict:
    t = sum(row[p] for p in T_POPS)
    b, nk = row["B"], row["NK"]
    lym = t + nk + b
    # the rule's "B >= T where both are printed" applies to printed values; Figure 1F prints none, so
    # a measured difference within the reading resolution is UNRESOLVED in either direction
    verdict = ("SUPPORTS T > B" if t - b > RESOLUTION else
               "CONTRADICTS T > B" if b - t > RESOLUTION else "UNRESOLVED")
    order = ">".join(k for k, _ in sorted({"T": t, "NK": nk, "B": b}.items(), key=lambda kv: -kv[1]))
    return {"T_pct_CD45": round(t, 3), "B_pct_CD45": round(b, 3), "NK_pct_CD45": round(nk, 3),
            "T_without_DNT_pct_CD45": round(t - row["DNT"], 3),
            "T_to_B": round(t / b, 1) if b else None, "ordering": order,
            "within_lymphoid_share": {"T": round(t / lym, 4), "NK": round(nk / lym, 4), "B": round(b / lym, 4)},
            "reading": verdict}


def main() -> int:
    pdf = config.KLEMM_2020_PDF
    if not pdf.exists():
        print(f"BLOCKED: {pdf} not found"); return 2
    m = measure(shapes(pdf, PAGE))
    c = controls(m)
    out = {"rule": "prespecified/klemm_t_vs_b.md", "source": "Klemm F, et al. Cell 181:1643-1660 (2020), Figure 1F",
           "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "page": PAGE,
           "method": "vector geometry of Figure 1F (exact path coordinates), calibrated on the panel's gridlines",
           "resolution_pct_CD45": RESOLUTION, "pt_per_percent": round(m["pt_per_percent"], 5),
           "gridlines_pt": m["gridlines_pt"], "controls": c, "rows": m["rows"]}
    if not c["all_pass"]:
        out["reading"] = "BLOCKED: a printed-number control failed"
        OUT.write_text(json.dumps(out, indent=2)); print(json.dumps(c, indent=1)); return 1
    out["IDH_mutant_glioma_n17"] = reading(m["rows"]["Glioma IDH mut"])
    out["IDH_wildtype_glioma_n40"] = reading(m["rows"]["Glioma IDH wt"])
    out["non_tumor_n6"] = reading(m["rows"]["Non-tumor"])
    OUT.write_text(json.dumps(out, indent=2))
    print("controls:", json.dumps({k: (v["pass"] if isinstance(v, dict) else v) for k, v in c.items()}))
    print(f"  melanoma CD8+ {c['melanoma_CD8']['measured']} vs printed {c['melanoma_CD8']['printed']}; "
          f"BrM lymphocytes {c['brm_lymphocytes_weighted']['measured']} vs printed {c['brm_lymphocytes_weighted']['printed']}")
    for k in ("IDH_mutant_glioma_n17", "IDH_wildtype_glioma_n40", "non_tumor_n6"):
        print(k, out[k])
    print(f"wrote {OUT.relative_to(config.PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
