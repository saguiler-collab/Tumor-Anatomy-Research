"""
check_core_closeout.py -- is the core close-out done? docs/CORE_CLOSEOUT.md defines the eight items.

Reads evidence only and changes nothing. Prints PASS or FAIL for each item with what it read, and exits 1 if any item
fails: an exit code is not evidence on its own (CLAUDE.md).

    python3 scripts/check_core_closeout.py
    python3 scripts/check_core_closeout.py --freeze     # re-hash every file docs/results_freeze_manifest.tsv records
"""
from __future__ import annotations

import argparse
import filecmp
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
PY = sys.executable
R = ROOT / "results"
SNAP = ROOT / "results_snapshot_20261006"
DOCS = ROOT / "docs"
PYTEST_LOG = R / "verification" / "pytest_core.log"     # first line "START <epoch>", written by the run that counts
CODE_TREES = ("ivygap", "scripts", "R", "tests")
# what scripts/repaired_methods.py and the genuine-DWLS re-run write (docs/METHOD_REPAIRS.md)
REPAIRED_OUTPUTS = ("tcga_gbm_frozen.json", "tcga_gbm_h5ad.json", "tcga_lgg_frozen.json", "tcga_lgg_h5ad.json",
                    "cptac_frozen.json", "cptac_h5ad.json", "ivygap_dwls_genuine_rescaled.json")
D36_FILES = ("estimates_full_extension.csv", "estimates_full_lgg_extension.csv")
D36_TEST = "def test_a_gene_subset_fit_never_writes_the_shared_per_sample_file"


def mtime(p: Path) -> float:
    return p.stat().st_mtime if p.exists() else 0.0


def newest_code_mtime() -> tuple[float, str]:
    files = [p for t in CODE_TREES for p in (ROOT / t).rglob("*") if p.suffix in (".py", ".R") and p.is_file()]
    p = max(files, key=lambda q: q.stat().st_mtime)
    return p.stat().st_mtime, str(p.relative_to(ROOT))


# ------------------------------------------------------------------------------------------------ the eight items
def probe(value: str | None) -> dict:
    """config.REPAIRED_METHODS and the three budgets, read in a fresh process with IVYGAP_REPAIRED unset or set."""
    env = {k: v for k, v in os.environ.items() if k != "IVYGAP_REPAIRED"}
    if value is not None:
        env["IVYGAP_REPAIRED"] = value
    code = ("import json; from ivygap import config; from ivygap.deconv import r_bridge as b; "
            "M = ('dwls', 'bayesprism', 'quantiseq'); print(json.dumps({'repaired': config.REPAIRED_METHODS, "
            "'budgets': {m: b.timeout_for(m) for m in M}, 'unbudgeted': b.UNBUDGETED, "
            "'registered': {m: b.REGISTERED_R_METHOD_TIMEOUTS.get(m, b.DEFAULT_R_TIMEOUT) for m in M}}))")
    out = subprocess.run([PY, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def item_repaired() -> tuple[bool, str]:
    d, r = probe(None), probe("0")
    on = d["repaired"] and all(v == d["unbudgeted"] for v in d["budgets"].values())
    off = (not r["repaired"]) and r["budgets"] == r["registered"]
    r_default = 'Sys.getenv("IVYGAP_REPAIRED", "1")' in (ROOT / "R" / "run_dwls.R").read_text()
    missing = [f for f in REPAIRED_OUTPUTS if not (R / "repaired" / f).exists()]
    return (on and off and r_default and not missing,
            f"unset: repaired={d['repaired']}, budgets {d['budgets']}; IVYGAP_REPAIRED=0: repaired={r['repaired']}, "
            f"budgets {r['budgets']}; R/run_dwls.R defaults on: {r_default}; results/repaired: "
            f"{len(REPAIRED_OUTPUTS) - len(missing)}/{len(REPAIRED_OUTPUTS)} outputs" + (f", missing {missing}" if missing else ""))


def item_verification() -> tuple[bool, str]:
    import verify_rerun as vr                           # noqa: PLC0415  (reads state only; runs nothing)
    state = vr.load_state()
    core = [s for s in vr.STEPS if s["id"] not in vr.SUPPLEMENTARY]
    missing = [s["id"] for s in core if s["id"] not in state]
    failed, problems, n_diff = [], [], 0
    for s in core:
        rec = state.get(s["id"])
        if rec is None:
            continue
        if rec["exit"] != 0 or rec.get("verdict") in ("FAILED", "DISAGREES"):
            failed.append(s["id"])
        outs = {k: vr._reclassify(vr._recompare_moved(s["id"], k, v)) for k, v in rec["outputs"].items()}
        bad = [k for k, v in outs.items() if v.get("status") in ("DIFFERS", "NOT WRITTEN", "UNREADABLE")]
        n_diff += len(bad)
        if bad and s["id"] not in vr.EXPLAINED:
            problems.append(f"{s['id']}: {bad} with no stated cause")
        elif bad:
            chk = vr.explanation_check(s["id"], bad)
            if not chk.startswith("checked"):
                problems.append(f"{s['id']}: {chk}")
    report = DOCS / "VERIFICATION_RERUN.md"
    fresh = mtime(report) >= mtime(vr.STATE) - 1
    matrix = DOCS / "EVALUATION_MATRIX.md"
    repaired_beside = matrix.exists() and "repaired" in matrix.read_text().lower()
    ok = not missing and not failed and not problems and fresh and repaired_beside
    return ok, (f"{len(core)} core steps, {len(core) - len(missing)} re-run" + (f" (missing {missing})" if missing else "")
                + f"; failed {failed or 'none'}; {n_diff} differing outputs, unexplained: {problems or 'none'}; "
                f"report newer than state: {fresh}; evaluation matrix reports repaired results: {repaired_beside}")


def item_d36() -> tuple[bool, str]:
    same = {f: (R / "extension" / f).exists() and filecmp.cmp(R / "extension" / f, SNAP / "extension" / f, shallow=False)
            for f in D36_FILES}
    test = D36_TEST in (ROOT / "tests" / "test_extension_panel.py").read_text()
    return all(same.values()) and test, f"byte-identical to the snapshot: {same}; guard test present: {test}"


def pytest_evidence(log: Path, code_mtime: float) -> tuple[bool, str]:
    """The run counts if it started after the newest code change and its summary has passes and no failure or error."""
    if not log.exists():
        return False, f"{log.relative_to(ROOT) if log.is_relative_to(ROOT) else log} not found"
    text = log.read_text(errors="ignore")
    m = re.match(r"START (\d+(?:\.\d+)?)", text)
    summary = [ln for ln in text.splitlines() if re.search(r"\d+ passed", ln)]
    if not m or not summary:
        return False, "no START line or no pytest summary in the log"
    last = summary[-1].strip("= ").strip()
    started = float(m.group(1))
    clean = not re.search(r"\d+ (failed|error|errors)\b", last)
    after = started > code_mtime
    return clean and after, (f"summary '{last}'; started {datetime.fromtimestamp(started):%Y-%m-%d %H:%M:%S}, "
                             f"after the last code change: {after}")


def item_pytest() -> tuple[bool, str]:
    t, f = newest_code_mtime()
    ok, detail = pytest_evidence(PYTEST_LOG, t)
    return ok, detail + f" (newest code file {f}, {datetime.fromtimestamp(t):%Y-%m-%d %H:%M:%S})"


def item_ledger() -> tuple[bool, str]:
    led = DOCS / "CLAIMS_LEDGER.md"
    statuses = re.findall(r"^\| \*\*(L\d+)\*\*.*\| ([^|]+) \|\s*$", led.read_text(), re.M) if led.exists() else []
    bad = [f"{c}: {s.strip()}" for c, s in statuses if "NOT VERIFIED" in s or "pending" in s.lower()]
    fresh = mtime(led) >= mtime(DOCS / "VERIFICATION_RERUN.md")
    return bool(statuses) and not bad and fresh, (f"{len(statuses)} claims; not verified or pending: {bad or 'none'}; "
                                                  f"rebuilt after the verification report: {fresh}")


def item_audit() -> tuple[bool, str]:
    aud = DOCS / "REGISTRATION_AUDIT.md"
    verdicts = re.findall(r"\*\*(PRECEDES|AFTER|CANNOT CHECK|NOT PROVABLE HERE \(see note\))\*\*",
                          aud.read_text()) if aud.exists() else []
    bad = [v for v in verdicts if v in ("AFTER", "CANNOT CHECK")]
    fresh = mtime(aud) >= mtime(DOCS / "VERIFICATION_RERUN.md")
    n_prec, n_np = verdicts.count("PRECEDES"), sum(v.startswith("NOT PROVABLE") for v in verdicts)
    return bool(verdicts) and not bad and fresh, (f"{n_prec} PRECEDES, {n_np} not provable here (note recorded), "
                                                  f"{len(bad)} AFTER or CANNOT CHECK; rebuilt after the report: {fresh}")


def item_docs() -> tuple[bool, str]:
    runs = {name: subprocess.run([PY, f"scripts/{name}", *args], cwd=ROOT, capture_output=True, text=True)
            for name, args in (("check_doc_numbers.py", ["--against-current"]), ("check_doc_statistics.py", ["--strict"]))}
    clean = {n: r.returncode == 0 and "CLEAN" in r.stdout for n, r in runs.items()}
    index = (DOCS / "FIGURES.md").read_text()
    figs = sorted(p.stem for p in (DOCS / "figures").glob("*.png"))
    unindexed = [f for f in figs if f not in index]
    return all(clean.values()) and not unindexed, (f"checkers CLEAN: {clean}; {len(figs)} figures, not in "
                                                   f"docs/FIGURES.md: {unindexed or 'none'}")


def item_ai_log() -> tuple[bool, str]:
    import export_ai_log as eal                         # noqa: PLC0415
    log = Path.home() / "STS_AI_LOG" / "PROMPT_LOG.md"
    live = [f for d in eal.PROJECT_DIRS for f in (eal.PROJECTS / d).glob("*.jsonl")]
    if not log.exists() or not live:
        return False, f"log present: {log.exists()}; live transcripts: {len(live)}"
    newest = max(live, key=lambda f: f.stat().st_mtime)
    sess = eal.read_session(newest)
    words = sess["prompts"] + sess["feedback"]                  # everything the student wrote to the AI
    last = max(datetime.fromisoformat(t.replace("Z", "+00:00")) for t, _ in words) if words else None
    exported = datetime.fromtimestamp(log.stat().st_mtime, tz=timezone.utc)
    ok = last is not None and exported > last
    return ok, (f"latest prompt or feedback {last.astimezone():%Y-%m-%d %H:%M:%S} in session {newest.stem[:8]}; "
                f"log exported {exported.astimezone():%Y-%m-%d %H:%M:%S}" if last else "no typed prompt found")


ITEMS = [("1 repaired DWLS / BayesPrism / quanTIseq", item_repaired),
         ("2 registered-vs-repaired verification", item_verification),
         ("3 D36 restoration", item_d36),
         ("4 pytest", item_pytest),
         ("5 claims ledger", item_ledger),
         ("6 registration audit", item_audit),
         ("7 documentation audit", item_docs),
         ("8 AI-use log", item_ai_log)]


# ------------------------------------------------------------------------------------------------------ the freeze
def compare_manifest(manifest: Path, root: Path, current: list[Path]) -> dict:
    """Changed, missing and added files against a freeze manifest (path, bytes, sha256)."""
    from freeze_results import sha256                   # noqa: PLC0415
    rows = [ln.split("\t") for ln in manifest.read_text().splitlines()[1:] if ln.strip()]
    recorded = {p: (int(b), h) for p, b, h in rows}
    out = {"recorded": len(recorded), "changed": [], "missing": [], "added": []}
    for p, (b, h) in recorded.items():
        f = root / p
        if not f.exists():
            out["missing"].append(p)
        elif f.stat().st_size != b or sha256(f) != h:
            out["changed"].append(p)
    out["added"] = sorted(str(f.relative_to(root)) for f in current if str(f.relative_to(root)) not in recorded)
    return out


def check_freeze() -> int:
    import freeze_results as fr                         # noqa: PLC0415
    if not fr.MANIFEST.exists():
        print(f"FAIL  no freeze recorded ({fr.MANIFEST.relative_to(ROOT)} missing)")
        return 1
    res = compare_manifest(fr.MANIFEST, ROOT, fr.frozen_files())
    bad = res["changed"] or res["missing"] or res["added"]
    print(f"{'FAIL' if bad else 'PASS'}  freeze: {res['recorded']} files recorded; changed {len(res['changed'])}, "
          f"missing {len(res['missing'])}, added {len(res['added'])}")
    for k in ("changed", "missing", "added"):
        for p in res[k][:20]:
            print(f"      {k}: {p}")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--freeze", action="store_true", help="re-hash the files the results freeze recorded")
    if ap.parse_args().freeze:
        return check_freeze()
    n_fail = 0
    for name, fn in ITEMS:
        try:
            ok, detail = fn()
        except Exception as exc:                        # a check that cannot read its evidence fails, loudly
            ok, detail = False, f"could not read the evidence: {type(exc).__name__}: {exc}"
        n_fail += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name}\n      {detail}")
    print(f"\n{len(ITEMS) - n_fail}/{len(ITEMS)} core items pass")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
