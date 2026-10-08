"""
build_software_inventory.py -- every R and Python package this study's code loads, with the version installed on the
machine that produced the results. Each version is read from the installed package, not recalled.

For each package it records the files that load it, where it was installed from (CRAN, Bioconductor, or a GitHub
repository and commit), and, for the R deconvolution bridges, the method each runs. Vendored code is listed with its
provenance file.

It writes no citations. Each package's authors state how it should be cited: in R, `citation("<package>")`; for
Python, the package's documentation. The STS 2027 rules (Appendix 2, p.31) require the reference list to be built
without AI. This file is the record of what to look up, not the list itself.

    python3 scripts/build_software_inventory.py      # writes docs/SOFTWARE_INVENTORY.md
"""
from __future__ import annotations

import ast
import importlib.metadata as md
import platform
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "SOFTWARE_INVENTORY.md"
R_FILES = sorted((ROOT / "R").glob("*.R"))
PY_FILES = sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "ivygap").rglob("*.py"))
R_LOAD_R = re.compile(r"(?:library|require)\(\s*[\"']?([A-Za-z][A-Za-z0-9.]*)[\"']?\s*[,)]"
                      r"|requireNamespace\(\s*[\"']([A-Za-z][A-Za-z0-9.]*)[\"']|\b([A-Za-z][A-Za-z0-9.]*)::")
R_LOAD_PY = re.compile(r"(?:library|require)\(\s*[\"']?([A-Za-z][A-Za-z0-9.]*)[\"']?\s*\)"
                       r"|requireNamespace\(\s*[\\\"']+([A-Za-z][A-Za-z0-9.]*)")
R_BASE = {"base", "stats", "utils", "methods", "parallel", "graphics", "grDevices", "tools", "datasets", "grid",
          "splines", "stats4", "tcltk", "compiler"}


def superseded(p: Path) -> bool:
    head = p.read_text(errors="ignore")[:400]
    return "SUPERSEDED" in head.split("\n\n")[0] if '"""' in head else False


def r_packages() -> dict[str, set[str]]:
    used: dict[str, set[str]] = defaultdict(set)
    for f in R_FILES:
        for m in R_LOAD_R.finditer(f.read_text(errors="ignore")):
            used[next(g for g in m.groups() if g)].add(str(f.relative_to(ROOT)))
    for f in PY_FILES:
        if superseded(f):
            continue
        for m in R_LOAD_PY.finditer(f.read_text(errors="ignore")):
            used[next(g for g in m.groups() if g)].add(str(f.relative_to(ROOT)))
    return {k: v for k, v in used.items() if k not in R_BASE}


def r_versions(pkgs: list[str]) -> tuple[str, dict[str, dict]]:
    code = ("pk <- commandArgs(TRUE); for (p in pk) { d <- tryCatch(utils::packageDescription(p), warning=function(w) NULL,"
            " error=function(e) NULL); if (is.null(d) || !is.list(d)) { cat(p, 'NA', '', '', sep='\\t'); cat('\\n'); next };"
            " src <- if (!is.null(d$RemoteType) && d$RemoteType == 'github') paste0('GitHub ', d$RemoteUsername, '/',"
            " d$RemoteRepo, '@', substr(d$RemoteSha, 1, 12)) else if (!is.null(d$Repository)) d$Repository else"
            " if (!is.null(d$biocViews)) 'Bioconductor' else 'local install';"
            " cat(p, d$Version, src, sep='\\t'); cat('\\n') };"
            " cat('R', paste(R.version$major, R.version$minor, sep='.'), '', sep='\\t'); cat('\\n')")
    out = subprocess.run(["Rscript", "-e", code, *pkgs], capture_output=True, text=True).stdout.splitlines()
    rows = {}
    r_version = "?"
    for ln in out:
        parts = (ln.split("\t") + ["", "", ""])[:3]
        if parts[0] == "R":
            r_version = parts[1]
        else:
            rows[parts[0]] = {"version": parts[1] if parts[1] != "NA" else "NOT INSTALLED", "source": parts[2]}
    return r_version, rows


def py_packages() -> dict[str, set[str]]:
    local = ({"ivygap", "scripts"} | {p.stem for p in (ROOT / "scripts").glob("*.py")}
             | {p.name for p in (ROOT / "vendor").glob("*")})
    std = set(sys.stdlib_module_names) | {"__future__"}
    used: dict[str, set[str]] = defaultdict(set)
    for f in PY_FILES:
        if superseded(f):
            continue
        try:
            tree = ast.parse(f.read_text(errors="ignore"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import) else
                     [node.module] if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module else [])
            for n in names:
                top = n.split(".")[0]
                if top not in std and top not in local:
                    used[top].add(str(f.relative_to(ROOT)))
    return used


def main() -> int:
    rp = r_packages()
    r_version, rv = r_versions(sorted(rp))
    pp = py_packages()
    dists = md.packages_distributions()
    when = datetime.now().astimezone()
    L = ["# Software inventory", "",
         f"*Generated {when:%Y-%m-%d %H:%M %Z} by `scripts/build_software_inventory.py`, on {platform.platform()} with R "
         f"{r_version} and Python {platform.python_version()}. Every version is read from the installed package. Files "
         f"marked SUPERSEDED in their docstring are not counted.*", "",
         "**What this is for.** It records which software produced the results, at which version. **It is not a "
         "reference list.** The authors of each package say how it should be cited: in R, run `citation(\"<package>\")`; "
         "for a Python package, see its documentation. Read the source and build the reference list from it (STS 2027 "
         "rules, Appendix 2: \"Entrants should generate their reference lists without the use of AI\").", "",
         "## R packages", "",
         "| package | version | installed from | loaded by | runs |", "|---|---|---|---|---|"]
    for p in sorted(rp, key=str.lower):
        files = sorted(rp[p])
        runs = sorted({Path(f).stem[len("run_"):] for f in files if Path(f).name.startswith("run_")})
        via = sorted(Path(f).stem for f in files if f.endswith(".py"))
        role = ", ".join(runs + ([f"called from Python: {', '.join(via)}"] if via else [])) or "support"
        info = rv.get(p, {"version": "?", "source": ""})
        L.append(f"| `{p}` | {info['version']} | {info['source']} | {', '.join(f'`{f}`' for f in files)} | {role} |")
    L += ["", "*'runs': the method name of each `R/run_<method>.R` bridge that loads the package, and the Python scripts "
          "that run it through R; 'support' means it is loaded only for installation or data handling.*", "",
          "## Python packages", "",
          "| import | distribution | version | loaded by (files) |", "|---|---|---|---|"]
    for top in sorted(pp, key=str.lower):
        dist = (dists.get(top) or [top])[0]
        try:
            ver = md.version(dist)
        except md.PackageNotFoundError:
            ver = "NOT INSTALLED (optional import)"
        files = sorted(pp[top])
        shown = ", ".join(f"`{f}`" for f in files[:6]) + (f" and {len(files) - 6} more" if len(files) > 6 else "")
        L.append(f"| `{top}` | {dist} | {ver} | {shown} |")
    L += ["", "## Vendored code (copied into this repository, with its provenance)", "",
          "| folder | provenance file |", "|---|---|"]
    for d in sorted(p for p in (ROOT / "vendor").glob("*") if p.is_dir()):
        prov = sorted(str(f.relative_to(ROOT)) for f in d.iterdir() if re.match(r"(PROVENANCE|README|LICENSE|DESCRIPTION)",
                                                                               f.name, re.I))
        L.append(f"| `{d.relative_to(ROOT)}/` | {', '.join(f'`{f}`' for f in prov) or 'none found'} |")
    L += ["", "## Data", "", "Every public dataset, with its accession, its source check and the scripts that read it, "
          "is in `docs/DATA_INVENTORY.md` (generated by `scripts/build_data_inventory.py`)."]
    OUT.write_text("\n".join(L) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rp)} R packages (R {r_version}), {len(pp)} Python imports")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
