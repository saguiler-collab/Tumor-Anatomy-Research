#!/usr/bin/env python3
"""
build_plates.py — render the figure plates from a results tree.

The plates carry a `const D = {...}` blob that every SVG and every number in the prose is
drawn from. That blob was assembled by hand once and immediately went stale: the published
page still showed ten real methods, DWLS at 0.877 as `R:DWLS`, and 200 test mixtures, long
after the run had thirteen comparable methods, DWLS as a reimplementation, and 500
mixtures. A figure that disagrees with the artefacts is worse than no figure, because it
gets quoted with confidence.

So the page is generated. `figures/plates_template.html` holds the markup and the drawing
code with `__FIGURE_DATA__` where the data goes; this script fills it from
`build_figure_data.build()` and writes a publishable file.

    python scripts/build_plates.py -o figures/plates.html
    python scripts/build_plates.py --results-dir results_archive/<stamp>/results

Then publish `figures/plates.html` as the artifact.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "build_figure_data", ROOT / "scripts" / "build_figure_data.py")
_bfd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bfd)

TEMPLATE = ROOT / "figures" / "plates_template.html"
PLACEHOLDER = "__FIGURE_DATA__"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("-o", "--out", default="figures/plates.html")
    ap.add_argument("--template", default=str(TEMPLATE))
    args = ap.parse_args()

    tpl_path = Path(args.template)
    if not tpl_path.exists():
        print(f"no template at {tpl_path}")
        return 1
    tpl = tpl_path.read_text()
    if PLACEHOLDER not in tpl:
        print(f"{tpl_path} has no {PLACEHOLDER} placeholder — is it already rendered?")
        return 1

    data = _bfd.build(Path(args.results_dir))
    blob = json.dumps(data, separators=(",", ":"), sort_keys=False)
    # A lone "</script>" inside the JSON would close the block early.
    blob = blob.replace("</", "<\\/")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(tpl.replace(PLACEHOLDER, blob))

    print(f"wrote {out}  ({out.stat().st_size / 1024:.1f} KB)")
    print(f"  {data['n_real_methods']} real methods "
          f"({data['n_comparable_methods']} comparable), {data['n_controls']} controls")
    reg = (data.get("registration") or {}).get("state")
    print(f"  registration state in the artefacts: {reg}")
    by = (data.get("agreement") or {}).get("by_yardstick") or []
    row = next((r for r in by if r.get("yardstick") == "synthetic_mixtures"), None)
    if row:
        print(f"  primary rho = {row['rho']:.4f}  "
              f"[{row['ci_low']:.3f}, {row['ci_high']:.3f}]  n={int(row['n_methods'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
