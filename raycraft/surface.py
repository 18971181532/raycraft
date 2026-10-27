"""光学表面几何：球面、二次曲面（conic）与偶次非球面。

坐标约定（右手系，光沿 +z 传播）
--------------------------------
* 每个面的顶点位于自身局部坐标原点，局部 z 轴为光轴。
* ``radius`` 为顶点曲率半径，正号表示曲率中心位于顶点右侧（+z）。
* 子午面（tangential）取 y-z 面，弧矢面（sagittal）取 x-z 面。

面型矢高::

    z(r) = c r² / (1 + sqrt(1 - (1+k) c² r²)) + Σ A_i r^(2i+4)

其中 ``c = 1/R``，``k`` 为二次曲面常数（k=0 球面，k=-1 抛物面，
k<-1 双曲面，k>0 扁椭球）。

求交采用「解析二次方程（球面）+ 牛顿迭代（非球面）」两条路径：
球面走闭式解既快又不会失根；非球面用隐式方程 ``F = z - sag(r)``
做阻尼牛顿，初值取切平面交点。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

Vec3 = Tuple[float, float, float]

__all__ = ["Surface", "Vec3", "sag_from_curvature"]


def norm(v: Vec3) -> float:
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def normalize(v: Vec3) -> Vec3:
    n = norm(v)
    if n == 0.0:
        raise ValueError("零向量无法归一化")
    return (v[0] / n, v[1] / n, v[2] / n)


def sag_from_curvature(c: float, r2: float, conic: float = 0.0) -> float:
    """给定曲率与 r² 计算球面/二次曲面矢高。"""
    if c == 0.0:
        return 0.0
    arg = 1.0 - (1.0 + conic) * c * c * r2
    if arg < 0.0:
        raise ValueError("r 超出二次曲面有效半口径")
    return c * r2 / (1.0 + math.sqrt(arg))


@dataclass
class Surface:
    """单个折射/反射面。

    Attributes:
        radius: 顶点曲率半径 (mm)，``math.inf`` 表示平面。
        thickness: 到下一面的轴向距离 (mm)。
        material: 该面之后的介质名称（玻璃库键）。
        semi_diameter: 通光半口径 (mm)，超出即判定为拦光。
        conic: 二次曲面常数 k。
        coefficients: 偶次非球面系数 ``(A4, A6, A8, ...)``，单位 mm^(-(2i+3))。
        is_stop: 是否为孔径光阑。
        is_mirror: 反射面（仅布局与追迹支持，优化器会跳过）。
        label: 便于阅读的标注。
    """

    radius: float = math.inf
    thickness: float = 0.0
    material: str = "air"
    semi_diameter: float = 10.0
    conic: float = 0.0
    coefficients: Tuple[float, ...] = ()
    is_stop: bool = False
    is_mirror: bool = False
    label: str = ""

    # ---- 派生量 ---------------------------------------------------------
    @property
    def curvature(self) -> float:
        return 0.0 if not math.isfinite(self.radius) else 1.0 / self.radius

    @property
    def is_plane(self) -> bool:
        return self.curvature == 0.0 and not self.coefficients

    @property
    def is_asphere(self) -> bool:
        return bool(self.coefficients)

    # ---- 矢高与斜率 -----------------------------------------------------
    def sag(self, r2: float) -> float:
        """由 r² 求矢高；r 超出二次曲面有效范围时抛 ValueError。"""
        base = sag_from_curvature(self.curvature, r2, self.conic)
        if not self.coefficients:
            return base
        # A_i 乘以 r^(2i+4) = (r²)^(i+2)
        acc = 0.0
        p = r2 * r2  # r^4
        for a in self.coefficients:
            acc += a * p
            p *= r2
        return base + acc

    def sag_and_slope(self, r2: float) -> Tuple[float, float]:
        """返回 (矢高, dz/dr)。r²=0 处斜率为 0。"""
        c = self.curvature
        sag = 0.0
        slope = 0.0
        if c != 0.0:
            arg = 1.0 - (1.0 + self.conic) * c * c * r2
            if arg < 0.0:
                raise ValueError("r 超出二次曲面有效半口径")
            root = math.sqrt(arg)
            sag = c * r2 / (1.0 + root)
            slope = c * math.sqrt(r2) / root if r2 > 0.0 else 0.0
        if self.coefficients:
            r2k = r2 * r2  # r^4
            r = math.sqrt(r2)
            for i, a in enumerate(self.coefficients):
                order = 2 * i + 4
                sag += a * r2k
                if r2 > 0.0:
                    slope += order * a * (r ** (order - 1))
                r2k *= r2
        return sag, slope

    def point_and_normal(self, x: float, y: float) -> Tuple[Vec3, Vec3]:
        """给定横向坐标，返回面上点与单位法线（朝 +z 一侧）。"""
        r2 = x * x + y * y
        sag, slope = self.sag_and_slope(r2)
        if r2 > 0.0:
            r = math.sqrt(r2)
            gx, gy = -slope * x / r, -slope * y / r
        else:
            gx, gy = 0.0, 0.0
        # 隐式面 F = z - sag(x,y) = 0，梯度 (-∂sag/∂x, -∂sag/∂y, 1)
        inv = 1.0 / math.sqrt(gx * gx + gy * gy + 1.0)
        return (x, y, sag), (gx * inv, gy * inv, inv)

    # ---- 求交 -----------------------------------------------------------
    def intersect(
        self, origin: Vec3, direction: Vec3, t_min: float = -1e-9
    ) -> Optional[Tuple[float, Vec3, Vec3]]:
        """求射线 ``origin + t·direction`` 与面的第一个交点。

        Returns:
            ``(t, 交点, 单位法线)``；无交点返回 ``None``。
        """
        if not self.is_asphere:
            return self._intersect_quadric(origin, direction, t_min)
        return self._intersect_asphere(origin, direction, t_min)

    def _intersect_quadric(
        self, origin: Vec3, direction: Vec3, t_min: float
    ) -> Optional[Tuple[float, Vec3, Vec3]]:
        ox, oy, oz = origin
        dx, dy, dz = direction
        c = self.curvature
        if c == 0.0:
            if abs(dz) < 1e-14:
                return None
            t = -oz / dz
            if t < t_min:
                return None
            return t, (ox + t * dx, oy + t * dy, 0.0), (0.0, 0.0, 1.0)

        # 二次曲面隐式方程：x² + y² + (1+k)z² - 2z/c = 0
        # k=0 球面，k=-1 抛物面，(1+k)<0 双曲面，>0 扁椭球
        kk = 1.0 + self.conic
        t_plane = -oz / dz if abs(dz) > 1e-14 else 0.0
        a2 = dx * dx + dy * dy + kk * dz * dz
        b = 2.0 * (ox * dx + oy * dy + kk * oz * dz) - 2.0 * dz / c
        cc = ox * ox + oy * oy + kk * oz * oz - 2.0 * oz / c
        if abs(a2) < 1e-12:
            if abs(b) < 1e-14:
                return None
            cand = [-cc / b]
        else:
            disc = b * b - 4.0 * a2 * cc
            if disc < 0.0:
                return None
            root = math.sqrt(disc)
            cand = [(-b - root) / (2.0 * a2), (-b + root) / (2.0 * a2)]
        valid = [t for t in cand if t >= t_min]
        if not valid:
            return None
        # 取靠近切平面猜想的那个根，避免落到面背侧
        t = min(valid, key=lambda v: abs(v - t_plane))
        qx, qy, qz = ox + t * dx, oy + t * dy, oz + t * dz
        nx, ny, nz = 2.0 * qx, 2.0 * qy, 2.0 * kk * qz - 2.0 / c
        inv = 1.0 / math.sqrt(nx * nx + ny * ny + nz * nz)
        nx, ny, nz = nx * inv, ny * inv, nz * inv
        if nz < 0.0:
            nx, ny, nz = -nx, -ny, -nz
        return t, (qx, qy, qz), (nx, ny, nz)

    def _intersect_asphere(
        self, origin: Vec3, direction: Vec3, t_min: float
    ) -> Optional[Tuple[float, Vec3, Vec3]]:
        ox, oy, oz = origin
        dx, dy, dz = direction
        t = -oz / dz if abs(dz) > 1e-14 else 0.0
        if t < t_min:
            t = t_min
        for _ in range(60):
            qx, qy, qz = ox + t * dx, oy + t * dy, oz + t * dz
            r2 = qx * qx + qy * qy
            try:
                sag, slope = self.sag_and_slope(r2)
            except ValueError:
                return None
            f = qz - sag
            if abs(f) < 1e-11:
                break
            if r2 > 0.0:
                r = math.sqrt(r2)
                gx, gy = -slope * qx / r, -slope * qy / r
            else:
                gx, gy = 0.0, 0.0
            denom = gx * dx + gy * dy + dz
            if abs(denom) < 1e-14:
                return None
            step = f / denom
            # 阻尼：限制单步位移，避免强非球面处越界振荡
            limit = 0.5 * max(self.semi_diameter, 1.0)
            if step > limit:
                step = limit
            elif step < -limit:
                step = -limit
            t -= step
            if t < t_min and step < 0.0:
                t = t_min
        else:
            return None
        qx, qy, qz = ox + t * dx, oy + t * dy, oz + t * dz
        r2 = qx * qx + qy * qy
        try:
            sag, slope = self.sag_and_slope(r2)
        except ValueError:
            return None
        if abs(qz - sag) > 1e-6:
            return None
        r = math.sqrt(r2)
        if r > 0.0:
            gx, gy = -slope * qx / r, -slope * qy / r
        else:
            gx, gy = 0.0, 0.0
        inv = 1.0 / math.sqrt(gx * gx + gy * gy + 1.0)
        return t, (qx, qy, qz), (gx * inv, gy * inv, inv)

    # ---- 辅助 -----------------------------------------------------------
    def edge_sag(self) -> float:
        """口径边缘的矢高，用于绘制与边缘厚度约束。"""
        return self.sag(self.semi_diameter ** 2)

    def with_updates(self, **kwargs) -> "Surface":
        from dataclasses import replace

        return replace(self, **kwargs)
