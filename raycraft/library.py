"""镜头库：单透镜、消色差双胶合、柯克三片式、双面非球面单透镜。

处方初值来自一阶（薄透镜/矩阵）光学解析估算，再由本项目的 DLS 优化器
精修，最终数值直接写在这里，可直接作为设计起点使用。面型数据格式::

    (radius, thickness, material, semi_diameter, is_stop, label)

``radius = inf`` 表示平面，``thickness`` 为到下一面的距离，``material``
为该面**之后**的介质。

各镜头的像质指标（RMS 点列图半径，五环采样，F/d/C 三波长）记录在
``REFERENCE`` 中，``tests/test_raycraft.py`` 会用它做回归保护——
任何把像质改坏的改动都会被 CI 拦下。
"""

from __future__ import annotations

import copy
import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import materials
from .surface import Surface
from .system import OpticalSystem

__all__ = ["SYSTEMS", "get", "names", "catalog", "set_f_number", "REFERENCE"]

Row = Tuple[float, float, str, float, bool, str]
FDC = (materials.F_LINE, materials.D_LINE, materials.C_LINE)


def _build(
    name: str,
    rows: Sequence[Row],
    fields: Sequence[float],
    description: str,
    wavelengths: Sequence[float] = FDC,
    primary: float = materials.D_LINE,
    object_distance: float = math.inf,
    f_number: Optional[float] = None,
    aspheres: Optional[Dict[int, Tuple[float, Tuple[float, ...]]]] = None,
    image_shift: float = 0.0,
) -> OpticalSystem:
    surfaces = [
        Surface(
            radius=r,
            thickness=t,
            material=m,
            semi_diameter=sd,
            is_stop=stop,
            label=lab,
        )
        for (r, t, m, sd, stop, lab) in rows
    ]
    if aspheres:
        for idx, (conic, coeffs) in aspheres.items():
            surfaces[idx].conic = conic
            surfaces[idx].coefficients = tuple(coeffs)
    sys_ = OpticalSystem(
        name,
        surfaces,
        fields=fields,
        wavelengths=wavelengths,
        primary_wavelength=primary,
        object_distance=object_distance,
        description=description,
    )
    if f_number is not None:
        set_f_number(sys_, f_number)
    sys_.autosize()
    if f_number is not None:
        set_f_number(sys_, f_number)
    sys_.image_shift = image_shift
    return sys_


def set_f_number(system: OpticalSystem, f_number: float) -> None:
    """缩放光阑口径，使近轴 F/# 等于给定值（一次线性修正即可精确）。"""
    cur = system.f_number
    if not math.isfinite(cur) or cur <= 0 or f_number <= 0:
        return
    k = system.stop_index
    system.surfaces[k].semi_diameter *= cur / f_number
    system._solve()


# --------------------------------------------------------------------------
def singlet() -> OpticalSystem:
    rows: List[Row] = [
        (51.68, 6.0, "N-BK7", 12.5, True, "L1 前"),
        (math.inf, 0.0, "air", 12.5, False, "L1 后"),
    ]
    return _build(
        "singlet",
        rows,
        fields=(0.0, 3.0, 6.0),
        description="f=100 mm 平凸单透镜，F/4。球差与轴向色差均未校正，作为对照基准："
        "像面取在近轴焦点上，轴上光斑已被球差弥散到 0.17 mm，边缘视场再叠加彗差与像散。",
        image_shift=0.0,
    )


# --------------------------------------------------------------------------
def achromat() -> OpticalSystem:
    rows: List[Row] = [
        (63.2064, 6.0, "N-BK7", 8.309, True, "冕牌 前"),
        (-36.4608, 3.0, "F2", 8.786, False, "胶合面"),
        (-129.2754, 0.0, "air", 8.804, False, "火石 后"),
    ]
    return _build(
        "achromat",
        rows,
        fields=(0.0, 1.5, 3.0),
        description="f=100 mm 双胶合消色差物镜，F/6。冕牌正片 + 火石负片，"
        "F 光与 C 光焦点重合，轴上色差被压到衍射量级。",
        f_number=6.0,
        image_shift=-0.0264,
    )


# --------------------------------------------------------------------------
def cooke_triplet() -> OpticalSystem:
    rows: List[Row] = [
        (13.9790, 3.2590, "N-SK16", 6.652, False, "正片 1 前"),
        (math.inf, 3.3375, "air", 5.964, False, "正片 1 后"),
        (-35.9837, 1.0000, "F2", 3.808, False, "负片 2 前"),
        (12.8607, 3.0184, "air", 3.576, False, "负片 2 后"),
        (math.inf, 4.5000, "air", 3.321, True, "光阑"),
        (46.2218, 2.9524, "N-SK16", 7.185, False, "正片 3 前"),
        (-33.6696, 0.0000, "air", 7.534, False, "正片 3 后"),
    ]
    return _build(
        "cooke_triplet",
        rows,
        fields=(0.0, 8.0, 16.0),
        description="f=50 mm 柯克三片式，F/5.6，全视场 ±16°（32°）。"
        "正-负-正三片把球差、彗差、像散与场曲同时压到几微米，是三阶像差平衡的教科书结构。",
        f_number=5.6,
        image_shift=-0.1222,
    )


# --------------------------------------------------------------------------
def asphere() -> OpticalSystem:
    rows: List[Row] = [
        (32.7463, 6.0, "N-SK16", 5.750, False, "非球面 前"),
        (-47.2228, 2.500, "air", 5.968, False, "非球面 后"),
    ]
    aspheres = {
        0: (-9.49824, (-1.61692e-05, -4.44989e-07, -1.32950e-09)),
        1: (5.00000, (-4.33514e-05, -3.29825e-07)),
    }
    return _build(
        "asphere",
        rows,
        fields=(0.0, 1.0, 2.0),
        description="f=32 mm 双面非球面单色物镜（633 nm），F/2.8。两个面都带二次曲面与偶次项，"
        "轴上与轴外同时逼近衍射极限（艾里半径 2.2 µm），说明「一片非球面顶一组球面」的真实含义。",
        wavelengths=(633.0,),
        primary=633.0,
        f_number=2.8,
        aspheres=aspheres,
        image_shift=-0.0187,
    )


SYSTEMS: Dict[str, Callable[[], OpticalSystem]] = {
    "singlet": singlet,
    "achromat": achromat,
    "cooke_triplet": cooke_triplet,
    "asphere": asphere,
}

_TITLES = {
    "singlet": "单透镜 · 未校正基准",
    "achromat": "消色差双胶合 · 100 mm",
    "cooke_triplet": "柯克三片式 · 50 mm",
    "asphere": "非球面单透镜 · 32 mm F/2.8",
}

#: 回归基准：各视场 RMS 光斑半径（µm），供测试比对
REFERENCE: Dict[str, List[float]] = {
    "singlet": [172.0, 205.3, 312.9],
    "achromat": [9.0, 4.9, 16.9],
    "cooke_triplet": [2.8, 7.1, 9.8],
    "asphere": [2.7, 1.3, 4.4],
}

_CACHE: Dict[str, OpticalSystem] = {}


def names() -> List[str]:
    return list(SYSTEMS)


def get(name: str) -> OpticalSystem:
    key = name.strip().lower()
    if key not in SYSTEMS:
        raise KeyError(f"未知镜头 {name!r}，可用：{', '.join(SYSTEMS)}")
    if key not in _CACHE:
        _CACHE[key] = SYSTEMS[key]()
    return copy.deepcopy(_CACHE[key])


def catalog() -> List[dict]:
    out = []
    for key in SYSTEMS:
        s = get(key)
        out.append(
            {
                "name": key,
                "title": _TITLES.get(key, key),
                "description": s.description,
                "efl": s.efl,
                "f_number": s.f_number,
                "fields": s.fields,
                "surfaces": len(s.surfaces),
            }
        )
    return out
