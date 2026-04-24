#!/usr/bin/env python3
"""
stangene PostToolUse hook — family handoff to stancounts (conditional) or QC.

Triggers after a Bash tool call. If the command was a successful invocation
of ``stangene`` (entry point, or ``python -m stangene``) and exited 0,
inject a handoff hint. The next step in the stan* family pipeline splits on
matrix type:

  - X is log1p-normalized → stancounts (reverse log1p → recover counts)
  - X is already counts    → QC stage directly

Design rules:
- Never block. This is PostToolUse — the tool already executed.
- No-op on any non-stangene command or non-zero exit, silently.
- No-op on any parse error or missing field (hook must never break
  a working session).
- Exit 0 in all cases. Feedback is delivered via JSON stdout
  (``hookSpecificOutput.additionalContext``), with stderr as a
  best-effort fallback for older runtimes.
"""

from __future__ import annotations

import json
import re
import sys
from typing import Optional

# Match:
#   - `stangene [harmonize] ...` (entry point)
#   - `python[3] [-m] [<path-or-module>]stangene[.py|.sub] ...`
# Dotted module paths (e.g. `python -m scripts.stangene`) are supported.
# Not matched: `stangene-report.md`, `cat path/stangene.txt`.
# Package-manager commands (`pip install stangene`, etc.) are filtered out
# by `_PKG_MGR_CMD` below.
_STANGENE_CMD = re.compile(
    r"(?:^|[\s;|&\n])"
    r"(?:"
    r"stangene(?:\s|$)"
    r"|python\d?\s+(?:-m\s+)?[\w./-]*\bstangene\b"
    r")"
)

# Silence false positives from package-manager invocations.
_PKG_MGR_CMD = re.compile(
    r"\b(?:pip|pip3|conda|mamba|micromamba|uv|poetry|pipx)"
    r"\s+(?:install|uninstall|add|remove|show|list|search|info|update|upgrade|sync)\b"
)

_HANDOFF = (
    "✅ stangene 执行成功。家族流水线的下一步视数据 matrix type 而定:\n"
    "\n"
    "  • X 是 **log1p normalized** → 用 **stancounts** 反向恢复 counts\n"
    "    (eca-curation pipeline 对应 02_counts_recovery 阶段)\n"
    "  • X 已是 **整数 counts** → 直接进入 **QC 阶段**\n"
    "    (eca-curation pipeline 对应 03_qc,由 qc-iterator subagent 驱动)\n"
    "\n"
    "怎么判断 matrix type:\n"
    "  • stanobj 的 standardize 报告(JSON)里一般含 matrix type 分类;\n"
    "  • 或手动检查:`adata.X.max()` 通常远大于 20 且为整数 → counts;\n"
    "    若是非整数浮点且最大值 ~ log1p(library_size) 量级 → 已 normalized。\n"
    "\n"
    "建议的下一步:\n"
    "  • 在 eca-curation pipeline session 里:运行\n"
    "      /eca-run <dataset>\n"
    "    pipeline 会按 stage 顺序推进,conditional 逻辑由 02_counts_recovery "
    "的 stage script 内部处理;\n"
    "  • 否则主对话根据 matrix type 显式触发 stancounts skill 或 QC。"
)


def _read_payload() -> dict:
    try:
        raw = sys.stdin.read()
    except OSError:
        return {}
    if not raw:
        return {}
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return {}


def _get_bash_command(payload: dict) -> Optional[str]:
    tool_name = (
        payload.get("tool_name")
        or (payload.get("tool_use") or {}).get("name")
        or (payload.get("toolUse") or {}).get("name")
    )
    if tool_name != "Bash":
        return None
    tool_input = (
        payload.get("tool_input")
        or (payload.get("tool_use") or {}).get("input")
        or (payload.get("toolUse") or {}).get("input")
        or {}
    )
    cmd = tool_input.get("command")
    return cmd if isinstance(cmd, str) else None


def _exit_code(payload: dict) -> int:
    res = (
        payload.get("tool_response")
        or payload.get("tool_result")
        or payload.get("toolResult")
        or {}
    )
    if not isinstance(res, dict):
        return 0
    for key in ("exit_code", "exitCode", "returncode", "returnCode"):
        code = res.get(key)
        if code is not None:
            try:
                return int(code)
            except (TypeError, ValueError):
                return 1
    if res.get("is_error") or res.get("isError"):
        return 1
    return 0


def main() -> int:
    payload = _read_payload()
    cmd = _get_bash_command(payload)
    if not cmd:
        return 0
    if _exit_code(payload) != 0:
        return 0

    # Split by shell operators so chained commands like
    # ``pip install stangene && stangene harmonize ...`` still trigger on
    # the run segment while the install segment is skipped.
    matched = False
    for seg in re.split(r"(?:&&|\|\||;|\|)", cmd):
        if _PKG_MGR_CMD.search(seg):
            continue
        if _STANGENE_CMD.search(seg):
            matched = True
            break
    if not matched:
        return 0

    out = {
        "hookSpecificOutput": {
            "hookEventName": "PostToolUse",
            "additionalContext": _HANDOFF,
        }
    }
    try:
        print(json.dumps(out))
    except Exception:
        pass
    try:
        print(_HANDOFF, file=sys.stderr)
    except Exception:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
