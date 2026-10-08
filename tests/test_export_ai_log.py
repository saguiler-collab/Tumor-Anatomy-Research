"""The AI-use log must hold the student's words and nothing else (scripts/export_ai_log.py).

A log that let system notifications or tool output pass as prompts, or that dropped the student's feedback on a
proposed action, would misstate who wrote what.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import export_ai_log as eal  # noqa: E402


def _user(content, **extra):
    return {"type": "user", "message": {"role": "user", "content": content}, **extra}


def test_typed_prompts_are_kept_and_system_text_is_not():
    typed = _user([{"type": "text", "text": "<ide_opened_file>The user opened the file /x/Anatomy_Test.md in the IDE."
                                            "</ide_opened_file>Keep going"}], origin={"kind": "human"})
    assert eal.typed_prompt(typed).startswith("Keep going")
    assert "Anatomy_Test.md" in eal.typed_prompt(typed)                     # noted as context, not as words
    assert eal.typed_prompt(_user("<task-notification>done</task-notification>",
                                  origin={"kind": "task-notification"})) is None
    assert eal.typed_prompt(_user("This session is being continued from a previous conversation ...",
                                  isCompactSummary=True)) is None
    assert eal.typed_prompt(_user([{"type": "tool_result", "content": "exit 0"}], toolUseResult={})) is None


def test_feedback_on_a_declined_action_is_kept_and_other_tool_results_are_not():
    declined = _user([{"type": "tool_result", "content": "The user doesn't want to proceed with this tool use. "
                       "The user provided the following reason for the rejection:  Do the commit after the tests."}])
    assert eal.rejection_feedback(declined) == "Do the commit after the tests."
    assert eal.rejection_feedback(_user([{"type": "tool_result", "content": "wrote docs/FIGURES.md"}])) is None
    quoted = _user([{"type": "tool_result", "content": "2026-10-07 search hit:  The user doesn't want to proceed with "
                     "this tool use. The user provided the following reason for the rejection:  Step 2 ..."}])
    assert eal.rejection_feedback(quoted) is None          # output that quotes the notice is not the student's words
    assert eal.rejection_feedback(_user("plain text")) is None
