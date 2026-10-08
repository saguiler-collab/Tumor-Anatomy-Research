"""
freeze_results.py -- record the results freeze (docs/CORE_CLOSEOUT.md): the SHA-256 and size of every file under
results/ (the verification folder excluded), docs/figures/ and docs/supplementary/, with the git commit and the time.

Run once, after the verification report and after every figure and table has been regenerated from the frozen run.
`scripts/check_core_closeout.py --freeze` re-hashes against it; any later run that touches a frozen file shows there.

    python3 scripts/freeze_results.py
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ivygap import config  # noqa: E402

MANIFEST = ROOT / "docs" / "results_freeze_manifest.tsv"
DOC = ROOT / "docs" / "RESULTS_FREEZE.md"
TREES = (config.RESULTS_DIR, ROOT / "docs" / "figures", ROOT / "docs" / "supplementary")
EXCLUDE = (config.RESULTS_DIR / "verification",)         # logs and state of the re-run harness, which keep changing


#: R console logs (results/diagnostics/<method>_r_stdout.log) are rewritten by every call of that R method, tests
#: included: a record of the last call, not a result.
CONSOLE_LOG = "_r_stdout.log"


def frozen_files(trees=TREES, exclude=EXCLUDE) -> list[Path]:
    return [p for t in trees for p in sorted(Path(t).rglob("*"))
            if p.is_file() and not p.name.startswith(".") and not p.name.endswith(CONSOLE_LOG)
            and not any(p.is_relative_to(e) for e in exclude)]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(path: Path, files: list[Path], root: Path) -> list[tuple[str, int, str]]:
    rows = [(str(p.relative_to(root)), p.stat().st_size, sha256(p)) for p in files]
    path.write_text("path\tbytes\tsha256\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in rows))
    return rows


def main() -> int:
    rows = write_manifest(MANIFEST, frozen_files(), ROOT)
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    head, dirty = git("rev-parse", "--short", "HEAD"), bool(git("status", "--porcelain"))
    when = datetime.now().astimezone()
    by_tree = {}
    for p, b, _ in rows:
        top = "/".join(p.split("/")[:2]) if p.startswith("docs/") else p.split("/")[0]
        n, s = by_tree.get(top, (0, 0))
        by_tree[top] = (n + 1, s + b)
    L = ["# Results freeze", "",
         f"*Written by `scripts/freeze_results.py` on {when:%Y-%m-%d %H:%M:%S %Z}, from git `{head}`"
         + (" plus the uncommitted close-out changes, which are committed together with this manifest" if dirty else "")
         + ". Checked by `python3 scripts/check_core_closeout.py --freeze`.*", "",
         "From this point every number, figure and table used for the write-up comes from this run. A file that changes "
         "after the freeze shows up in the check. The definition is in `docs/CORE_CLOSEOUT.md`.", "",
         "| tree | files | size |", "|---|---|---|"]
    L += [f"| `{t}/` | {n} | {s / 1e6:,.1f} MB |" for t, (n, s) in sorted(by_tree.items())]
    L += [f"| **total** | **{len(rows)}** | **{sum(b for _, b, _ in rows) / 1e6:,.1f} MB** |", "",
          f"Every file's SHA-256 is in `{MANIFEST.relative_to(ROOT)}`. Excluded: `results/verification/` (the "
          "re-run harness's own logs and state), R console logs (`*_r_stdout.log`, rewritten by every R call) and hidden "
          "files.", "",
          "**Limit.** `results/` is not in git. The manifest records those files, but on a fresh clone they cannot be "
          "re-checked. The figures and tables in `docs/` are in git and can be."]
    DOC.write_text("\n".join(L) + "\n")
    print(f"froze {len(rows)} files ({sum(b for _, b, _ in rows) / 1e6:,.1f} MB) -> {MANIFEST.relative_to(ROOT)}, "
          f"{DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
