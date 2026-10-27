"""像质分析：点列图、光线扇形图、场曲/畸变、几何 MTF、布局数据。

关于几何 MTF
------------
本模块用「特征函数」估计几何 MTF：设像面上采样到的 N 个落点为
``r_i``，则

    MTF(f) = | (1/N) Σ exp(-2πi f·r_i) |

它是像面辐照度分布特征函数的无偏估计，等价于对点列图密度做傅里叶
变换，但不需要分箱，噪声仅随 1/√N 下降，因此用几千条光线即可得到
平滑曲线。

几何 MTF 忽略衍射，在像差主导（大像差、短焦距）时与实测吻合；衍射
受限时它恒接近 1。为此本模块同时给出由 F/# 与波长解析计算的衍射
极限 MTF（圆孔自相关），两者并列显示才能正确判断系统处于哪种极限。
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from . import materials
from .system import OpticalSystem, TraceResult, normalize

__all__ = [
    "pupil_hexapolar",
    "pupil_disk",
    "spot_diagram",
    "ray_fan",
    "field_analysis",
    "distortion",
    "geometric_mtf",
    "diffraction_mtf",
    "mtf_analysis",
    "encircled_energy",
    "layout",
    "surface_profiles",
]


# --------------------------------------------------------------------------
# 光瞳采样
# --------------------------------------------------------------------------
def pupil_hexapolar(rings: int = 6) -> List[Tuple[float, float]]:
    """六角极坐标采样：每环 6k 点，点列图常用，同心环便于看球差。"""
    pts: List[Tuple[float, float]] = [(0.0, 0.0)]
    for k in range(1, rings + 1):
        rho = k / rings
        n = 6 * k
        for j in range(n):
            a = 2.0 * math.pi * j / n
            pts.append((rho * math.cos(a), rho * math.sin(a)))
    return pts


def pupil_disk(count: int = 4000) -> List[Tuple[float, float]]:
    """圆域均匀面积采样（斐波那契螺旋），用于 MTF 与能量集中度。"""
    pts: List[Tuple[float, float]] = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(count):
        rho = math.sqrt((i + 0.5) / count)
        a = i * golden
        pts.append((rho * math.cos(a), rho * math.sin(a)))
    return pts


def _hit(
    system: OpticalSystem, field: float, px: float, py: float, wl: float, defocus: float
) -> Optional[Tuple[float, float]]:
    res = system.trace_to_image(field, px, py, wl, defocus=defocus)
    if not res.ok or res.image_point is None:
        return None
    return res.image_point[0], res.image_point[1]


# --------------------------------------------------------------------------
# 点列图
# --------------------------------------------------------------------------
def spot_diagram(
    system: OpticalSystem,
    field: float,
    rings: int = 6,
    defocus: float = 0.0,
    wavelengths: Optional[Sequence[float]] = None,
) -> dict:
    """计算单一视场的点列图。"""
    wls = list(wavelengths or system.wavelengths)
    samples = pupil_hexapolar(rings)
    groups = []
    total = 0
    failed = 0
    for wl in wls:
        xs: List[float] = []
        ys: List[float] = []
        for px, py in samples:
            total += 1
            hit = _hit(system, field, px, py, wl, defocus)
            if hit is None:
                failed += 1
                continue
            xs.append(hit[0])
            ys.append(hit[1])
        groups.append({"wavelength": wl, "x": xs, "y": ys})

    # 以主波长像斑质心为参考点
    primary = min(groups, key=lambda g: abs(g["wavelength"] - system.primary_wavelength))
    if primary["x"]:
        cx = sum(primary["x"]) / len(primary["x"])
        cy = sum(primary["y"]) / len(primary["y"])
    else:
        cx = cy = 0.0

    for g in groups:
        g["x"] = [x - cx for x in g["x"]]
        g["y"] = [y - cy for y in g["y"]]

    all_r2 = [x * x + y * y for g in groups for x, y in zip(g["x"], g["y"])]
    rms = math.sqrt(sum(all_r2) / len(all_r2)) if all_r2 else 0.0
    geo = math.sqrt(max(all_r2)) if all_r2 else 0.0
    # 主波长单独统计（色差不计入）
    p_r2 = [x * x + y * y for x, y in zip(primary["x"], primary["y"])]
    rms_primary = math.sqrt(sum(p_r2) / len(p_r2)) if p_r2 else 0.0

    return {
        "field": field,
        "reference": (cx, cy),
        "groups": groups,
        "rms": rms,
        "rms_primary": rms_primary,
        "geo_radius": geo,
        "rays": total - failed,
        "vignetted": failed,
    }


# --------------------------------------------------------------------------
# 光线扇形图
# --------------------------------------------------------------------------
def ray_fan(
    system: OpticalSystem,
    field: float,
    n: int = 41,
    defocus: float = 0.0,
    wavelengths: Optional[Sequence[float]] = None,
) -> dict:
    """子午（tangential）与弧矢（sagittal）光线扇形图。

    横轴为归一化光瞳坐标，纵轴为像面横向偏离量 (µm)。
    """
    wls = list(wavelengths or system.wavelengths)
    chief = _hit(system, field, 0.0, 0.0, system.primary_wavelength, defocus)
    ref_y = chief[1] if chief else 0.0
    ref_x = chief[0] if chief else 0.0

    tang = []
    for wl in wls:
        xs, ys = [], []
        for k in range(n):
            py = -1.0 + 2.0 * k / (n - 1)
            hit = _hit(system, field, 0.0, py, wl, defocus)
            if hit is None:
                continue
            xs.append(py)
            ys.append((hit[1] - ref_y) * 1000.0)
        tang.append({"wavelength": wl, "pupil": xs, "value": ys})

    sag = []
    for wl in wls:
        xs, ys = [], []
        for k in range(n):
            px = -1.0 + 2.0 * k / (n - 1)
            hit = _hit(system, field, px, 0.0, wl, defocus)
            if hit is None:
                continue
            xs.append(px)
            ys.append((hit[0] - ref_x) * 1000.0)
        sag.append({"wavelength": wl, "pupil": xs, "value": ys})

    return {"field": field, "tangential": tang, "sagittal": sag}


# --------------------------------------------------------------------------
# 场曲与畸变
# --------------------------------------------------------------------------
def _fan_spread(
    system: OpticalSystem, field: float, axis: str, defocus: float, n: int = 21
) -> float:
    """扇形图相对主光线的 RMS 弥散（µm 量级不影响比较）。"""
    wl = system.primary_wavelength
    chief = _hit(system, field, 0.0, 0.0, wl, defocus)
    if chief is None:
        return float("inf")
    ref = chief[1] if axis == "y" else chief[0]
    total = 0.0
    count = 0
    for k in range(n):
        p = -1.0 + 2.0 * k / (n - 1)
        px, py = (0.0, p) if axis == "y" else (p, 0.0)
        hit = _hit(system, field, px, py, wl, defocus)
        if hit is None:
            continue
        v = (hit[1] if axis == "y" else hit[0]) - ref
        total += v * v
        count += 1
    return math.sqrt(total / count) if count else float("inf")


def _best_focus(
    system: OpticalSystem, field: float, axis: str, span: float = 4.0, coarse: int = 41
) -> float:
    """在 ±span mm 内搜索使扇形图 RMS 最小的离焦量（粗扫 + 抛物线拟合）。"""
    best_d = 0.0
    best_v = float("inf")
    step = 2.0 * span / (coarse - 1)
    samples = []
    for k in range(coarse):
        d = -span + k * step
        v = _fan_spread(system, field, axis, d)
        samples.append((d, v))
        if v < best_v:
            best_v, best_d = v, d
    idx = [i for i, (d, v) in enumerate(samples) if d == best_d][0]
    if 0 < idx < len(samples) - 1:
        (d0, v0), (d1, v1), (d2, v2) = samples[idx - 1], samples[idx], samples[idx + 1]
        denom = v0 - 2.0 * v1 + v2
        if abs(denom) > 1e-15:
            offset = 0.5 * (v0 - v2) / denom * step
            if abs(offset) <= step:
                best_d = d1 + offset
    return best_d


def distortion(
    system: OpticalSystem, fields: Optional[Sequence[float]] = None
) -> List[dict]:
    """相对畸变：真实主光线像高与近轴像高之比。"""
    out = []
    wl = system.primary_wavelength
    theta_unit = 1.0
    # 近轴主光线在首面的起始高度（主光线通过光阑中心）
    m = system.system_matrix_at(wl, upto=system.stop_index)
    m00, m01 = m[0][0], m[0][1]
    for f in fields or system.fields:
        th = math.radians(f)
        u0 = math.tan(th) if system.is_infinite_object else 0.0
        y0 = -m01 * u0 / m00 if abs(m00) > 1e-15 else 0.0
        y_par = system.paraxial_height(y0, u0, wl)
        chief = system.trace_to_image(f, 0.0, 0.0, wl)
        y_real = chief.image_point[1] if chief.ok else float("nan")
        dist = (y_real - y_par) / y_par * 100.0 if abs(y_par) > 1e-9 else 0.0
        out.append(
            {
                "field": f,
                "paraxial_height": y_par,
                "real_height": y_real,
                "distortion": dist,
            }
        )
    return out


def field_analysis(system: OpticalSystem, fields: Optional[Sequence[float]] = None) -> dict:
    """场曲（子午/弧矢最佳焦面）与畸变随视场的变化。"""
    fl = list(fields or system.fields)
    curves = []
    for f in fl:
        ft = _best_focus(system, f, "y")
        fs = _best_focus(system, f, "x")
        curves.append({"field": f, "tangential_focus": ft, "sagittal_focus": fs})
    return {"curves": curves, "distortion": distortion(system, fl)}


# --------------------------------------------------------------------------
# MTF
# --------------------------------------------------------------------------
def geometric_mtf(
    system: OpticalSystem,
    field: float,
    freqs: Sequence[float],
    n_rays: int = 3000,
    defocus: float = 0.0,
    wavelength: Optional[float] = None,
) -> dict:
    """视场的几何 MTF（子午 T / 弧矢 S）。"""
    wl = system.primary_wavelength if wavelength is None else wavelength
    samples = pupil_disk(n_rays)
    # 先用通过光阑中心的光线定参考点
    ref = _hit(system, field, 0.0, 0.0, wl, defocus)
    if ref is None:
        return {"field": field, "freqs": list(freqs), "tangential": [0.0] * len(freqs),
                "sagittal": [0.0] * len(freqs)}
    rx, ry = ref
    xs: List[float] = []
    ys: List[float] = []
    for px, py in samples:
        hit = _hit(system, field, px, py, wl, defocus)
        if hit is None:
            continue
        xs.append(hit[0] - rx)
        ys.append(hit[1] - ry)
    n = len(xs)
    if n == 0:
        return {"field": field, "freqs": list(freqs), "tangential": [0.0] * len(freqs),
                "sagittal": [0.0] * len(freqs)}
    tang, sag = [], []
    two_pi = 2.0 * math.pi
    for f in freqs:
        st = ct = ss = cs = 0.0
        for i in range(n):
            a = -two_pi * f * ys[i]
            b = -two_pi * f * xs[i]
            ct += math.cos(a)
            st += math.sin(a)
            cs += math.cos(b)
            ss += math.sin(b)
        tang.append(math.hypot(ct, st) / n)
        sag.append(math.hypot(cs, ss) / n)
    return {
        "field": field,
        "freqs": list(freqs),
        "tangential": tang,
        "sagittal": sag,
        "rays": n,
        "wavelength": wl,
    }


def diffraction_mtf(
    f_number: float, wavelength_nm: float, freqs: Sequence[float]
) -> List[float]:
    """圆孔衍射极限 MTF（截止频率 1/(λF)）。"""
    if not math.isfinite(f_number) or f_number <= 0:
        return [0.0] * len(freqs)
    lam_mm = wavelength_nm * 1e-6
    cutoff = 1.0 / (lam_mm * f_number)
    out = []
    for f in freqs:
        x = abs(f) / cutoff
        if x >= 1.0:
            out.append(0.0)
        else:
            out.append((2.0 / math.pi) * (math.acos(x) - x * math.sqrt(1.0 - x * x)))
    return out


def mtf_analysis(
    system: OpticalSystem,
    fields: Optional[Sequence[float]] = None,
    max_freq: Optional[float] = None,
    points: int = 60,
    n_rays: int = 2500,
    defocus: float = 0.0,
) -> dict:
    """多视场 MTF 汇总 + 衍射极限参考曲线。"""
    fl = list(fields or system.fields)
    if max_freq is None:
        lam_mm = system.primary_wavelength * 1e-6
        fnum = system.f_number if math.isfinite(system.f_number) else 8.0
        max_freq = 1.0 / (lam_mm * fnum)
    freqs = [max_freq * i / (points - 1) for i in range(points)]
    curves = [geometric_mtf(system, f, freqs, n_rays, defocus, None) for f in fl]
    return {
        "freqs": freqs,
        "curves": curves,
        "diffraction_limit": diffraction_mtf(system.f_number, system.primary_wavelength, freqs),
        "cutoff": max_freq,
        "defocus": defocus,
    }


# --------------------------------------------------------------------------
# 能量集中度
# --------------------------------------------------------------------------
def encircled_energy(system: OpticalSystem, field: float, n_rays: int = 2000) -> dict:
    wl = system.primary_wavelength
    samples = pupil_disk(n_rays)
    ref = _hit(system, field, 0.0, 0.0, wl, 0.0)
    if ref is None:
        return {"radius_um": [], "fraction": []}
    r = []
    for px, py in samples:
        hit = _hit(system, field, px, py, wl, 0.0)
        if hit is None:
            continue
        r.append(math.hypot(hit[0] - ref[0], hit[1] - ref[1]) * 1000.0)
    if not r:
        return {"radius_um": [], "fraction": []}
    r.sort()
    n = len(r)
    idx = [int(round(k * (n - 1) / 49)) for k in range(50)]
    return {
        "radius_um": [r[i] for i in idx],
        "fraction": [(i + 1) / n for i in idx],
        "rms_um": math.sqrt(sum(v * v for v in r) / n),
    }


# --------------------------------------------------------------------------
# 布局
# --------------------------------------------------------------------------
def surface_profiles(system: OpticalSystem, samples: int = 60) -> List[dict]:
    """各面的二维剖面折线（全局坐标 z-y）。"""
    zpos = system.surface_positions()
    out = []
    for i, s in enumerate(system.surfaces):
        xs, ys = [], []
        r = s.semi_diameter
        for k in range(samples + 1):
            y = -r + 2.0 * r * k / samples
            try:
                sag = s.sag(y * y)
            except ValueError:
                continue
            xs.append(zpos[i] + sag)
            ys.append(y)
        out.append(
            {
                "index": i + 1,
                "z": zpos[i],
                "points": list(zip(xs, ys)),
                "material": s.material,
                "is_stop": s.is_stop,
                "semi_diameter": r,
                "radius": s.radius,
            }
        )
    return out


def layout(
    system: OpticalSystem,
    rays_per_field: Sequence[float] = (-1.0, -0.6, -0.2, 0.2, 0.6, 1.0),
    wavelengths: Optional[Sequence[float]] = None,
    fields: Optional[Sequence[float]] = None,
) -> dict:
    """2D 剖面布局：面型折线 + 各视场光线折线。"""
    profiles = surface_profiles(system)
    wls = list(wavelengths or [system.primary_wavelength])
    fl = list(fields or system.fields)
    ray_lines = []
    for f in fl:
        for py in rays_per_field:
            for wl in wls:
                res = system.trace_to_image(f, 0.0, py, wl)
                if not res.ok:
                    continue
                ray_lines.append(
                    {
                        "field": f,
                        "pupil": py,
                        "wavelength": wl,
                        "points": [(p[2], p[1]) for p in res.points],
                    }
                )
    z = system.surface_positions()
    # 纵向范围要同时容纳镜片外形和全部光线（像高可能远大于口径）
    ys: List[float] = [0.0]
    for prof in profiles:
        ys.extend(abs(p[1]) for p in prof.get("points", []))
    for line in ray_lines:
        ys.extend(abs(p[1]) for p in line["points"])
    y_max = max(ys)
    pad = 0.06 * max(y_max, 1e-6)
    return {
        "profiles": profiles,
        "rays": ray_lines,
        "image_plane": system.image_plane_z,
        "z_range": [z[0] - 5.0, system.image_plane_z + 5.0],
        "y_range": [-(y_max + pad), y_max + pad],
        "semi_diameter": max(s.semi_diameter for s in system.surfaces),
    }
