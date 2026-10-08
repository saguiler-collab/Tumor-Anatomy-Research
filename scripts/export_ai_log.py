"""
export_ai_log.py -- the AI-use record that the Regeneron STS 2027 rules ask for (Appendix 4, p.34): a log of every
prompt given to the AI, and a statement of which code the AI generated.

Reads the Claude Code session transcripts for this project and its predecessor (~/.claude/projects/), keeps a raw copy
of each (Claude Code deletes transcripts older than `cleanupPeriodDays`, 30 by default), and writes:

- PROMPT_LOG.md: every prompt the student typed, and every answer the student gave to a question from the AI, in time
  order (US Eastern), by session. Nothing is paraphrased. IDE context and system text are removed; slash commands are
  kept.
- AI_CONTRIBUTION_INVENTORY.md: per session, the AI's tool calls by kind and every file it created or edited through
  its file tools; from git, which commits carry an AI co-author trailer and which tracked code files were first added
  by one.

Output goes outside the repository by default, so that nothing personal is committed or pushed by accident.

    python3 scripts/export_ai_log.py                  # writes ~/STS_AI_LOG/
    python3 scripts/export_ai_log.py --out <dir>
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parent.parent
PROJECTS = Path.home() / ".claude" / "projects"
# this project, then the predecessor TCGA-GBM project it was rebuilt from (CLAUDE.md); Claude Code names a project's
# transcript folder after its path, with every character other than a letter, digit or '-' replaced by '-'
_PREDECESSOR = Path.home() / "Downloads" / "cancer_judging_machine-master 22"
PROJECT_DIRS = [re.sub(r"[^A-Za-z0-9-]", "-", str(p)) for p in
                (ROOT, _PREDECESSOR / "GBM_cell_atlas",
                 Path.home() / "Downloads" / "cancer_judging_machine-master 16" / "cancer_judging_machine")]
TZ = ZoneInfo("America/New_York")
TAG = re.compile(r"<(system-reminder|ide_selection|ide_opened_file|ide_diagnostics)>.*?</\1>", re.S)
OPENED = re.compile(r"<ide_opened_file>The user opened the file (\S+) in the IDE")
NOT_TYPED = ("<task-notification>", "<local-command-stdout>", "Caveat: The messages below",
             "This session is being continued from a previous conversation", "[Request interrupted")
FILE_TOOLS = {"Write", "Edit", "NotebookEdit", "MultiEdit"}
CODE_DIRS = ("ivygap/", "scripts/", "R/", "tests/")


def local(ts: str) -> datetime:
    return datetime.fromisoformat(ts.replace("Z", "+00:00")).astimezone(TZ)


def text_of(content) -> str:
    if isinstance(content, str):
        return content
    return "\n".join(x.get("text", "") for x in content if x.get("type") == "text")


def typed_prompt(d: dict) -> str | None:
    """The student's own words in a user entry, or None if the entry was written by the system or a tool."""
    c = d["message"]["content"]
    if d.get("toolUseResult") is not None or (isinstance(c, list) and any(x.get("type") == "tool_result" for x in c)):
        return None
    origin = d.get("origin") if isinstance(d.get("origin"), dict) else {}
    raw = text_of(c)
    if origin:
        if origin.get("kind") != "human":
            return None
    elif d.get("isMeta") or d.get("isCompactSummary") or raw.lstrip().startswith(NOT_TYPED):
        return None
    m = re.search(r"<command-name>(.*?)</command-name>.*?<command-args>(.*?)</command-args>", raw, re.S)
    if m:
        return f"[slash command] {m.group(1).strip()} {m.group(2).strip()}".strip()
    opened = OPENED.search(raw)
    body = TAG.sub("", raw).strip()
    if not body:
        return None
    return body + (f"\n\n*(file open in the editor: `{Path(opened.group(1)).name}`)*" if opened else "")


DECLINED = "The user doesn't want to proceed with this tool use."
REJECTION = "The user provided the following reason for the rejection:"


def rejection_feedback(d: dict) -> str | None:
    """The student's own words when declining an action the AI proposed (a plan, an edit), or None.

    Only a result that IS the rejection notice counts: it starts with it. A command whose output merely quotes the
    notice (the AI searching the transcript, for one) is not the student's feedback."""
    c = d["message"]["content"]
    for x in c if isinstance(c, list) else []:
        if x.get("type") == "tool_result":
            t = x["content"] if isinstance(x["content"], str) else " ".join(
                y.get("text", "") for y in x["content"] if isinstance(y, dict))
            if t.lstrip().startswith(DECLINED) and REJECTION in t:
                return t.split(REJECTION, 1)[1].strip() or None
    return None


def read_session(path: Path) -> dict:
    s = {"id": path.stem, "file": path, "prompts": [], "answers": [], "feedback": [], "tools": collections.Counter(),
         "files": collections.Counter(), "first": None, "last": None}
    for line in path.open():
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        ts = d.get("timestamp")
        if ts and d.get("type") in ("user", "assistant"):
            s["first"] = s["first"] or ts
            s["last"] = ts
        if d.get("type") == "user":
            tur = d.get("toolUseResult")
            if isinstance(tur, dict) and isinstance(tur.get("answers"), dict):
                s["answers"].append((ts, tur["answers"]))
            p = typed_prompt(d)
            if p:
                s["prompts"].append((ts, p))
            fb = rejection_feedback(d)
            if fb:
                s["feedback"].append((ts, fb))
        elif d.get("type") == "assistant" and not d.get("isSidechain"):
            for x in d["message"].get("content", []) if isinstance(d["message"].get("content"), list) else []:
                if x.get("type") == "tool_use":
                    s["tools"][x["name"]] += 1
                    fp = (x.get("input") or {}).get("file_path") or (x.get("input") or {}).get("notebook_path")
                    if x["name"] in FILE_TOOLS and fp:
                        s["files"][fp] += 1
    return s


def backup(out: Path) -> list[Path]:
    """Copy each transcript into out/raw/<project>/, replacing an older or smaller copy; return the copies."""
    kept = []
    for proj in PROJECT_DIRS:
        for src in sorted((PROJECTS / proj).glob("*.jsonl")):
            dst = out / "raw" / proj / src.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists() or src.stat().st_size >= dst.stat().st_size:
                shutil.copy2(src, dst)
        kept += sorted((out / "raw" / proj).glob("*.jsonl")) if (out / "raw" / proj).exists() else []
    return kept


def rel(fp: str) -> str:
    p = Path(fp)
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p).replace(str(Path.home()), "~")


def git_inventory() -> list[str]:
    def git(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout
    shas = git("log", "--format=%H").split()
    ai = {h for h in shas if re.search(r"Co-Authored-By: Claude", git("show", "-s", "--format=%B", h), re.I)}
    L = ["## From git", "", f"- Commits on this branch: {len(shas)}; with an AI co-author trailer "
         f"(`Co-Authored-By: Claude ...`): {len(ai)}; without: {len(shas) - len(ai)}.", ""]
    for h in [h for h in shas if h not in ai]:
        L.append(f"  - without a trailer: `{h[:7]}` {git('show', '-s', '--format=%ad %s', '--date=short', h).strip()}")
    files = [f for f in git("ls-files").splitlines() if f.startswith(CODE_DIRS)]
    by_dir: dict[str, list[str]] = collections.defaultdict(list)
    for f in files:
        first = git("log", "--diff-filter=A", "--format=%H", "--", f).split()
        tag = "AI co-authored" if first and first[-1] in ai else "no AI trailer"
        by_dir[f.split("/")[0] + "/"].append(f"`{f}` ({tag})")
    L += ["", f"### Tracked code files ({len(files)}), by the commit that first added each", ""]
    for d in CODE_DIRS:
        n_ai = sum("AI co-authored" in x for x in by_dir[d])
        L.append(f"- **{d}** {len(by_dir[d])} files, {n_ai} first added in an AI co-authored commit")
    L += ["", "<details><summary>Every file</summary>", ""] + [f"- {x}" for d in CODE_DIRS for x in by_dir[d]]
    return L + ["", "</details>"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path.home() / "STS_AI_LOG")
    out = ap.parse_args().out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    sessions = sorted((read_session(p) for p in backup(out)), key=lambda s: s["first"] or "")
    sessions = [s for s in sessions if s["first"]]
    made = f"{datetime.now(TZ):%Y-%m-%d %H:%M %Z}"

    P = ["# Prompt log: the AI-assisted parts of this project", "",
         f"*Exported {made} by `scripts/export_ai_log.py` from the Claude Code session transcripts (raw copies in "
         f"`raw/`). Each entry is the student's prompt exactly as typed, the student's answer to a question the AI asked, "
         f"or the student's feedback when declining an action the AI proposed. AI replies are not included; they are in "
         f"the raw transcripts. Times are US Eastern.*", ""]
    total = 0
    for s in sessions:
        proj = s["file"].parent.name.split("-Repos-")[-1].split("-Downloads-")[-1]
        P += ["", f"## Session `{s['id'][:8]}` ({proj}): {local(s['first']):%Y-%m-%d %H:%M} to "
              f"{local(s['last']):%Y-%m-%d %H:%M}", "",
              f"{len(s['prompts'])} prompts, {len(s['answers'])} answered questions, {len(s['feedback'])} feedback "
              f"notes.", ""]
        events = ([(t, "prompt", p) for t, p in s["prompts"]] + [(t, "answer", a) for t, a in s["answers"]]
                  + [(t, "feedback", f) for t, f in s["feedback"]])
        for t, kind, x in sorted(events, key=lambda e: e[0]):
            if kind == "prompt":
                total += 1
                body = x if len(x) <= 6000 else x[:6000] + "\n\n*(truncated at 6,000 characters; full text in raw/)*"
                P += [f"**{local(t):%Y-%m-%d %H:%M}** -- prompt {total}", ""] + [f"> {ln}" if ln else ">"
                                                                                for ln in body.splitlines()] + [""]
            elif kind == "feedback":
                P += [f"**{local(t):%Y-%m-%d %H:%M}** -- feedback, given when declining an action the AI proposed", ""]
                P += [f"> {ln}" if ln else ">" for ln in x.splitlines()] + [""]
            else:
                P += [f"**{local(t):%Y-%m-%d %H:%M}** -- answers to the AI's questions", ""]
                P += [f"- Q: {q}\n  - **A: {a}**" for q, a in x.items()] + [""]
    (out / "PROMPT_LOG.md").write_text("\n".join(P) + "\n")

    I = ["# AI contribution inventory", "",
         f"*Exported {made} by `scripts/export_ai_log.py`. Facts only: counts of what the AI did, read from the "
         f"session transcripts and from git. Edits the AI made through shell commands (sed, heredocs, scripts) are "
         f"counted under Bash, not itemised by file.*", "", "## From the session transcripts", "",
         "| session | project | from | to | student prompts | answered questions | AI tool calls | files written or edited |",
         "|---|---|---|---|---|---|---|---|"]
    for s in sessions:
        proj = s["file"].parent.name.split("-Repos-")[-1].split("-Downloads-")[-1]
        I.append(f"| `{s['id'][:8]}` | {proj} | {local(s['first']):%Y-%m-%d} | {local(s['last']):%Y-%m-%d} | "
                 f"{len(s['prompts'])} | {len(s['answers'])} | {sum(s['tools'].values()):,} | {len(s['files'])} |")
    for s in sessions:
        I += ["", f"### Session `{s['id'][:8]}`", "",
              "Tool calls: " + ", ".join(f"{k} {v:,}" for k, v in s["tools"].most_common()), "",
              "<details><summary>Files created or edited through the AI's file tools</summary>", ""]
        I += [f"- `{rel(f)}` ({n})" for f, n in sorted(s["files"].items(), key=lambda kv: rel(kv[0]))]
        I += ["", "</details>"]
    I += [""] + git_inventory()
    (out / "AI_CONTRIBUTION_INVENTORY.md").write_text("\n".join(I) + "\n")
    print(f"{len(sessions)} sessions, {total} prompts, {sum(len(s['answers']) for s in sessions)} answered questions, "
          f"{sum(len(s['feedback']) for s in sessions)} feedback notes -> {out}/PROMPT_LOG.md, "
          f"{out}/AI_CONTRIBUTION_INVENTORY.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
