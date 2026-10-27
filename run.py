"""启动引导：把项目根目录加入 sys.path 后调用 CLI 主入口。

本机 Python 为嵌入式发行版（``python310._pth`` 限制），``python -m`` 不会
自动把当前目录加入模块搜索路径，因此统一用本脚本启动。
"""

import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from raycraft.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
