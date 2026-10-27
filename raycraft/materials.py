"""光学材料库：Sellmeier 色散模型、阿贝数、标准谱线。

色散公式（波长单位 µm）::

    n²(λ) = 1 + Σ_i  B_i λ² / (λ² - C_i)

该式即 Schott 玻璃目录给出的 Sellmeier 系数形式，在 0.3–2.5 µm
范围内精度优于 1e-5。所有波长对外统一使用纳米，内部换算为微米。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

# 标准夫琅禾费谱线（nm）
F_LINE = 486.1327   # 氢 F 线，蓝
D_LINE = 587.5618   # 氦 d 线，黄（光学设计主波长）
C_LINE = 656.2725   # 氢 C 线，红

# 常用的三波长组
SPECTRUM_FDC: Tuple[float, float, float] = (F_LINE, D_LINE, C_LINE)


@dataclass(frozen=True)
class Glass:
    """一种各向同性光学介质。

    Attributes:
        name: 目录名，例如 ``N-BK7``。
        nd: 厂商标称的 d 线折射率（仅用于核对，计算一律走 Sellmeier）。
        vd: 厂商标称的阿贝数。
        sellmeier: ``(B1, B2, B3, C1, C2, C3)``；``None`` 表示空气。
        note: 简短说明。
    """

    name: str
    nd: float
    vd: float
    sellmeier: Optional[Tuple[float, float, float, float, float, float]] = None
    note: str = ""

    @property
    def is_air(self) -> bool:
        return self.sellmeier is None

    def index(self, wavelength_nm: float) -> float:
        """返回给定波长下的折射率。"""
        if self.sellmeier is None:
            return 1.0
        lam = wavelength_nm * 1e-3
        l2 = lam * lam
        b1, b2, b3, c1, c2, c3 = self.sellmeier
        n2 = (
            1.0
            + b1 * l2 / (l2 - c1)
            + b2 * l2 / (l2 - c2)
            + b3 * l2 / (l2 - c3)
        )
        return math.sqrt(n2)

    def abbe(self) -> float:
        """由色散曲线反算的阿贝数 Vd = (nd-1)/(nF-nC)。"""
        if self.is_air:
            return 0.0
        nd = self.index(D_LINE)
        return (nd - 1.0) / (self.index(F_LINE) - self.index(C_LINE))

    def dispersion(self) -> float:
        """部分色散 nF - nC。"""
        if self.is_air:
            return 0.0
        return self.index(F_LINE) - self.index(C_LINE)


AIR = Glass(name="air", nd=1.0, vd=0.0, sellmeier=None, note="空气，视为无色散")

# --------------------------------------------------------------------------
# Schott 玻璃目录（Sellmeier 系数，波长单位 µm）
# --------------------------------------------------------------------------
_CATALOG: Tuple[Glass, ...] = (
    Glass(
        "N-BK7", 1.51680, 64.17,
        (1.03961212, 0.231792344, 1.01046945,
         0.00600069867, 0.0200179144, 103.560653),
        "硼硅酸盐冕牌玻璃，最通用的正透镜材料",
    ),
    Glass(
        "N-K5", 1.52249, 59.48,
        (1.08511833, 0.199562005, 0.930511663,
         0.00661099503, 0.024110866, 111.982777),
        "轻冕玻璃，低色散",
    ),
    Glass(
        "N-SK16", 1.62041, 60.32,
        (1.34317774, 0.241144399, 0.994317969,
         0.00704687339, 0.0229005, 92.7508526),
        "重冕玻璃，高折射率低色散，正透镜常用",
    ),
    Glass(
        "N-BAF10", 1.67003, 47.11,
        (1.5851495, 0.143559385, 1.08521269,
         0.00926681282, 0.0424489805, 105.613573),
        "钡火石玻璃",
    ),
    Glass(
        "N-LAK22", 1.65113, 55.89,
        (1.14229781, 0.535138441, 1.04088385,
         0.00585778594, 0.0198546147, 100.834017),
        "镧冕玻璃，消色差双胶合常用正片",
    ),
    Glass(
        "N-LASF9", 1.85025, 32.17,
        (2.00029547, 0.298926886, 1.80691843,
         0.0121426017, 0.0538736236, 156.530829),
        "重镧火石玻璃，高折射率，广角镜头常用",
    ),
    Glass(
        "F2", 1.62004, 36.37,
        (1.34533359, 0.209073176, 0.937357162,
         0.00997743871, 0.0470450767, 111.886764),
        "火石玻璃，消色差双胶合的负片",
    ),
    Glass(
        "N-F2", 1.62005, 36.43,
        (1.39757037, 0.159201403, 1.2686543,
         0.00995906143, 0.0546931752, 119.248346),
        "无铅火石玻璃",
    ),
    Glass(
        "SF2", 1.64769, 33.82,
        (1.40301821, 0.231767504, 0.939056586,
         0.0105795466, 0.0493226978, 112.405955),
        "重火石玻璃",
    ),
    Glass(
        "SF5", 1.67270, 32.25,
        (1.52481889, 0.187085527, 1.42729015,
         0.011254756, 0.0588995392, 129.141675),
        "重火石玻璃，高色散",
    ),
    Glass(
        "N-SF6", 1.80518, 25.36,
        (1.77931763, 0.338149866, 2.08734474,
         0.0133714182, 0.0617533621, 174.01759),
        "致密火石玻璃，极高色散",
    ),
    Glass(
        "N-SF11", 1.78472, 25.68,
        (1.73759695, 0.313747346, 1.89878101,
         0.013188707, 0.0623068142, 155.23629),
        "致密火石玻璃，短波吸收边靠近可见光",
    ),
)

CATALOG: Dict[str, Glass] = {g.name.upper(): g for g in _CATALOG}
CATALOG["AIR"] = AIR


def glass(name: str) -> Glass:
    """按名称取玻璃，大小写不敏感；空气可用 ``air`` / ``AIR`` / ``""``。"""
    key = (name or "air").strip().upper()
    if key in ("", "AIR", "VACUUM"):
        return AIR
    try:
        return CATALOG[key]
    except KeyError:
        raise KeyError(
            f"未知玻璃 {name!r}，可用材料：{', '.join(sorted(CATALOG))}"
        ) from None


def index_of(name: str, wavelength_nm: float) -> float:
    """便捷函数：材料名 + 波长 → 折射率。"""
    return glass(name).index(wavelength_nm)


def list_glasses() -> Tuple[Glass, ...]:
    """返回全部目录玻璃（不含空气）。"""
    return _CATALOG


def abbe_pair(crown: str, flint: str) -> float:
    """两种玻璃的色散比 (V_crown / V_flint)，消色差条件要求其等于
    负片的相对光焦度比。"""
    return glass(crown).abbe() / glass(flint).abbe()
