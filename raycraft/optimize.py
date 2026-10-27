"""镜头自动优化：阻尼最小二乘法（Damped Least Squares, Levenberg-Marquardt）。

目标函数
--------
取各视场、各波长点列图上光线相对质心的横向偏离作为残差向量 ``r``：

    merit = Σ r_i²  +  Σ w_c · max(0, 限值 - 约束量)²

像差残差衡量成像质量，约束项负责保证结构可实现（中心厚度、边缘厚度
不为负、曲率半径不越界）。

迭代
----
雅可比矩阵用前向差分构造，法方程加阻尼后求解::

    (JᵀJ + λ·diag(JᵀJ)) δ = -Jᵀ r

 merit 下降则接受步长并令 ``λ ← λ/3``，否则丢弃并令 ``λ ← λ×10``。
该策略在远离极小点时接近高斯-牛顿（快），在病态区域退化为小步长
梯度下降（稳），是镜头设计软件的标准做法。
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from . import analysis, materials
from .system import OpticalSystem

__all__ = ["Variable", "Optimizer", "solve_linear"]


@dataclass
class Variable:
    """一个可优化的结构量。"""

    surface: int          # 面序号（0 起）；param="focus" 时忽略
    param: str = "radius"  # radius | thickness | conic | coefficient | focus
    lower: float = -1e6
    upper: float = 1e6
    step: float = 1e-3     # 相对差分步长
    index: int = 0         # coefficient 时的非球面项序号
    min_step: float = 1e-3  # 绝对差分步长下限
    mirror: Optional[Tuple[int, float]] = None  # 对称联动：(面序号, 符号)

    def get(self, system: OpticalSystem) -> float:
        if self.param == "focus":
            return system.image_shift
        s = system.surfaces[self.surface]
        if self.param == "coefficient":
            return s.coefficients[self.index] if self.index < len(s.coefficients) else 0.0
        return getattr(s, self.param)

    def set(self, system: OpticalSystem, value: float) -> None:
        value = min(self.upper, max(self.lower, value))
        if self.param == "focus":
            system.image_shift = value
            return
        s = system.surfaces[self.surface]
        if self.param == "coefficient":
            co = list(s.coefficients)
            while len(co) <= self.index:
                co.append(0.0)
            co[self.index] = value
            s.coefficients = tuple(co)
        else:
            setattr(s, self.param, value)
        if self.mirror is not None and self.param == "radius":
            linked, sign = self.mirror
            setattr(system.surfaces[linked], "radius", sign * value)


def solve_linear(a: List[List[float]], b: List[float]) -> List[float]:
    """列主元高斯消元，A x = b。"""
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-14:
            return [0.0] * n
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        for r in range(col + 1, n):
            f = m[r][col] / pv
            if f == 0.0:
                continue
            for c in range(col, n + 1):
                m[r][c] -= f * m[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        acc = m[r][n] - sum(m[r][c] * x[c] for c in range(r + 1, n))
        x[r] = acc / m[r][r]
    return x


class Optimizer:
    """对 :class:`OpticalSystem` 做 RMS 点列图最小化。"""

    def __init__(
        self,
        system: OpticalSystem,
        variables: Sequence[Variable],
        fields: Optional[Sequence[float]] = None,
        rings: int = 2,
        wavelengths: Optional[Sequence[float]] = None,
        scale: float = 0.02,
        min_center_thickness: float = 1.0,
        min_edge_thickness: float = 0.5,
        target_efl: Optional[float] = None,
        efl_weight: float = 300.0,
        target_f_number: Optional[float] = None,
        fn_weight: float = 150.0,
        min_image_distance: float = 2.0,
    ) -> None:
        self.system = copy.deepcopy(system)
        self.variables = list(variables)
        self.fields = list(fields or system.fields)
        self.rings = rings
        self.wavelengths = list(wavelengths or system.wavelengths)
        self.scale = scale
        self.min_center_thickness = min_center_thickness
        self.min_edge_thickness = min_edge_thickness
        self.target_efl = target_efl
        self.efl_weight = efl_weight
        self.target_f_number = target_f_number
        self.fn_weight = fn_weight
        self.min_image_distance = min_image_distance
        self.history: List[float] = []
        self._samples = analysis.pupil_hexapolar(rings)

    # ------------------------------------------------------------------
    def residuals(self) -> List[float]:
        """残差向量。

        每个视场以**主波长**像斑质心为唯一参考点，再计入其余波长的
        落点，这样色差与单色像差会同时被压低（若各波长各自取质心，
        色差将被"抵消"掉，优化器就会放弃消色差）。
        """
        sys_ = self.system
        out: List[float] = []
        inv = 1.0 / self.scale
        for f in self.fields:
            data: List[Tuple[float, List[Tuple[float, float]]]] = []
            for wl in self.wavelengths:
                pts: List[Tuple[float, float]] = []
                for px, py in self._samples:
                    res = sys_.trace_to_image(f, px, py, wl)
                    if not res.ok:
                        pts.append((float("nan"), float("nan")))
                        continue
                    pts.append((res.image_point[0], res.image_point[1]))
                data.append((wl, pts))
            primary = min(data, key=lambda kv: abs(kv[0] - sys_.primary_wavelength))[1]
            good = [(x, y) for x, y in primary if x == x]
            if not good:
                continue
            cx = sum(p[0] for p in good) / len(good)
            cy = sum(p[1] for p in good) / len(good)
            for wl, pts in data:
                for x, y in pts:
                    if x != x:  # NaN：光线溢出
                        out.append(8.0)
                        out.append(8.0)
                    else:
                        out.append((x - cx) * inv)
                        out.append((y - cy) * inv)
        # 焦距约束：否则优化器会靠缩短焦距来缩小像斑
        if self.target_efl:
            efl = sys_.efl
            if not math.isfinite(efl):
                out.append(50.0)
            else:
                out.append(self.efl_weight * (efl - self.target_efl) / self.target_efl)
        if self.target_f_number:
            fn = sys_.f_number
            if not math.isfinite(fn):
                out.append(40.0)
            else:
                out.append(self.fn_weight * (fn - self.target_f_number) / self.target_f_number)
        # 像面必须落在系统之外（含调焦偏移）
        img_gap = sys_.image_distance + sys_.image_shift
        if img_gap < self.min_image_distance:
            out.append(60.0 * (self.min_image_distance - img_gap))
        # 结构约束
        ns = len(sys_.surfaces)
        for i, s in enumerate(sys_.surfaces[:-1]):
            gap = s.thickness
            if gap < self.min_center_thickness and s.material.lower() != "air":
                out.append(200.0 * (self.min_center_thickness - gap))
            elif gap < 0.2:
                out.append(150.0 * (0.2 - gap))
            nxt = sys_.surfaces[i + 1]
            if s.material.lower() != "air" and nxt.material.lower() != "air":
                r = min(s.semi_diameter, nxt.semi_diameter)
                try:
                    edge = gap - s.sag(r * r) + nxt.sag(r * r)
                except ValueError:
                    edge = -1.0
                if edge < self.min_edge_thickness:
                    out.append(200.0 * (self.min_edge_thickness - edge))
        return out

    def merit(self) -> float:
        return sum(v * v for v in self.residuals())

    # ------------------------------------------------------------------
    def _jacobian(self, r0: List[float]) -> List[List[float]]:
        cols: List[List[float]] = []
        for var in self.variables:
            orig = var.get(self.system)
            if not math.isfinite(orig):
                delta = 0.5
            else:
                delta = max(var.step * abs(orig), var.min_step)
            var.set(self.system, orig + delta)
            self.system.refocus()
            r1 = self.residuals()
            var.set(self.system, orig)
            self.system.refocus()
            col = []
            for i in range(len(r0)):
                a = r0[i]
                b = r1[i] if i < len(r1) else a
                col.append((b - a) / delta)
            cols.append(col)
        return cols

    def step(self, lam: float = 1e-3) -> Tuple[float, float]:
        """执行一次 LM 迭代，返回 (旧 merit, 新 merit)。"""
        r0 = self.residuals()
        m0 = sum(v * v for v in r0)
        cols = self._jacobian(r0)
        nv = len(self.variables)
        nr = len(r0)
        jtj = [[0.0] * nv for _ in range(nv)]
        jtr = [0.0] * nv
        for a in range(nv):
            ca = cols[a]
            for b in range(a, nv):
                cb = cols[b]
                acc = 0.0
                for i in range(nr):
                    acc += ca[i] * cb[i]
                jtj[a][b] = acc
                jtj[b][a] = acc
            acc = 0.0
            for i in range(nr):
                acc += ca[i] * r0[i]
            jtr[a] = acc
        diag = [jtj[i][i] if jtj[i][i] > 1e-12 else 1.0 for i in range(nv)]
        for i in range(nv):
            jtj[i][i] += lam * diag[i]
        delta = solve_linear(jtj, [-v for v in jtr])

        backup = [var.get(self.system) for var in self.variables]
        for var, d in zip(self.variables, delta):
            var.set(self.system, var.get(self.system) + d)
        self.system.refocus()
        m1 = self.merit()
        if not math.isfinite(m1) or m1 > m0:
            for var, v in zip(self.variables, backup):
                var.set(self.system, v)
            self.system.refocus()
            return m0, m0
        return m0, m1

    def run(self, iterations: int = 20, verbose: bool = False) -> dict:
        lam = 1e-3
        m = self.merit()
        self.history = [m]
        for it in range(iterations):
            m_old, m_new = self.step(lam)
            if m_new >= m_old:
                lam *= 10.0
                if lam > 1e8:
                    break
                continue
            lam = max(lam / 3.0, 1e-9)
            m = m_new
            self.history.append(m)
            if verbose:
                print(f"  iter {it + 1:3d}  merit={m:.6e}  λ={lam:.2e}")
            if len(self.history) > 3 and self.history[-2] - self.history[-1] < 1e-10 * self.history[-2]:
                break
        return {
            "merit": m,
            "history": self.history,
            "iterations": len(self.history) - 1,
            "system": self.system,
        }


# --------------------------------------------------------------------------
# 预设变量集
# --------------------------------------------------------------------------
def radius_variables(system: OpticalSystem, indices: Sequence[int], **kw) -> List[Variable]:
    out = []
    for i in indices:
        r = system.surfaces[i].radius
        lo = kw.get("lower", 8.0 if r > 0 else -3000.0)
        hi = kw.get("upper", 3000.0 if r > 0 else -8.0)
        out.append(Variable(i, "radius", lo, hi))
    return out


def thickness_variables(
    system: OpticalSystem, indices: Sequence[int], lower: float = 0.5, upper: float = 60.0
) -> List[Variable]:
    return [Variable(i, "thickness", lower, upper) for i in indices]


def default_variables(
    system: OpticalSystem,
    vary_glass_thickness: bool = False,
    vary_focus: bool = True,
) -> List[Variable]:
    """默认变量集：全部球面半径 + 空气间隔（+ 整组调焦量）。

    玻璃中心厚度一般由工艺与边缘厚度决定，默认不参与优化（这也是
    实际设计流程的做法：先定厚度，再优化面形与间隔）。需要时可打开
    ``vary_glass_thickness``，上限取通光口径的 2 倍以避免出现"厚砖头"。

    ``vary_focus`` 让像面可以整体前后移动——这是控制场曲的关键自由度：
    平场探测器的位置本身就是设计中要选的量。
    """
    vars_: List[Variable] = []
    if vary_focus:
        vars_.append(Variable(-1, "focus", -6.0, 6.0, step=0.0, min_step=0.05))
    for i, s in enumerate(system.surfaces):
        if math.isfinite(s.radius) and not s.is_stop:
            r = s.radius
            vars_.append(Variable(i, "radius", 8.0 if r > 0 else -3000.0,
                                  3000.0 if r > 0 else -8.0))
    for i, s in enumerate(system.surfaces[:-1]):
        if s.material.lower() != "air":
            if vary_glass_thickness:
                vars_.append(
                    Variable(i, "thickness", 1.2, max(3.0 * s.semi_diameter, 12.0))
                )
        elif not s.is_stop:
            vars_.append(Variable(i, "thickness", 0.5, 60.0))
    return vars_


def multi_start(
    system: OpticalSystem,
    starts: int = 5,
    seed: int = 2026,
    radius_jitter: float = 0.22,
    iterations: int = 60,
    variables_factory=default_variables,
    **kwargs,
) -> dict:
    """多起点优化：随机扰动面形后分别优化，取 merit 最小者。

    DLS 只保证局部收敛，镜头设计的 merit 面又常有多解（不同的弯曲
    方式都能满足一阶条件）。工业上常用「随机起点 + 局部优化」来跳出
    局部极小，这里给一个轻量实现。
    """
    import random

    rng = random.Random(seed)
    best: Optional[dict] = None
    trials: List[dict] = []
    for k in range(max(1, starts)):
        trial = copy.deepcopy(system)
        if k > 0:
            for surf in trial.surfaces:
                if math.isfinite(surf.radius) and not surf.is_stop:
                    surf.radius *= 1.0 + rng.uniform(-radius_jitter, radius_jitter)
            trial.refocus()
        opt = Optimizer(trial, variables_factory(trial), **kwargs)
        res = opt.run(iterations=iterations)
        res["start"] = k
        trials.append(
            {
                "start": k,
                "merit": res["merit"],
                "iterations": res["iterations"],
                "start_merit": res["history"][0],
            }
        )
        if best is None or res["merit"] < best["merit"]:
            best = res
    best["trials"] = trials  # type: ignore[index]
    return best  # type: ignore[return-value]
