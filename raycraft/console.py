"""控制台编码兜底。

Windows 的 ``sys.stdout`` 在输出被重定向到管道/文件时会退回本地代码页
（简体中文系统是 cp936，英文系统是 cp1252），此时打印 ``─``、``═``、``∞``
或中文会抛 ``UnicodeEncodeError``。CI 与 ``python run.py ... > out.txt``
都属于这种情况，所以在所有入口处主动把标准流切到 UTF-8。

对真正的终端没有副作用：Python 3.6+ 在 Windows 控制台上本来就用 UTF-16
写入，``encoding`` 已经是 ``utf-8``，函数会直接跳过。
"""

from __future__ import annotations

import sys

_UTF8_ALIASES = {"utf8", "utf8mb4", "cp65001", "utf_8"}


def force_utf8() -> None:
    """把 stdout/stderr 切换成 UTF-8，失败时静默放弃。"""
    for stream in (sys.stdout, sys.stderr):
        encoding = (getattr(stream, "encoding", "") or "").lower().replace("-", "")
        if encoding in _UTF8_ALIASES:
            continue
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except (OSError, ValueError):
            pass
