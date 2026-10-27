"""光学系统：近轴矩阵求解、入瞳定位、三维实光线追迹、自动定口径。

近轴理论约定
------------
光线状态用 ``(y, u)`` 表示，``y`` 为高度、``u`` 为与光轴夹角。
球面折射写成矩阵形式（``Φ = (n'-n)c``）::

    [y']   [ 1     0 ] [y]
    [u'] = [-Φ/n'  n/n'] [u]

平移矩阵为 ``[[1, t], [0, 1]]``。系统矩阵 ``M`` 的 ``C`` 元素即 ``-1/f``，
因此 ``EFL = -1/M[1][0]``、``BFL = -M[0][0]/M[1][0]``，对有限共轭同样成立。

入瞳
----
孔径光阑之前的全部面构成矩阵 ``M_s``。轴上无穷远光线的入射高度与
光阑高度满足 ``y_stop = M_s[0][0] · y_0``，故入瞳半径
``R_ep = R_stop / M_s[0][0]``；主光线在物空间的延长线与光轴交点给出
入瞳位置 ``z_ep = M_s[0][1] / M_s[0][0]``（相对第一面顶点）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from . import materials
from .surface import Surface, Vec3, normalize

Matrix = Tuple[Tuple[float, float], Tuple[float, float]]

IDENTITY: Matrix = ((1.0, 0.0), (0.0, 1.0))


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def refract(d: Vec3, normal: Vec3, n1: float, n2: float) -> Optional[Vec3]:
    """矢量形式的斯涅尔定律。

    Args:
        d: 入射单位方向。
        normal: 面法线单位向量（朝向任意，内部会自动取向）。
        n1, n2: 入射侧/出射侧折射率。

    Returns:
        出射单位方向；全反射或无解时返回 ``None``。
    """
    cos_i = -(d[0] * normal[0] + d[1] * normal[1] + d[2] * normal[2])
    if cos_i < 0.0:
        normal = (-normal[0], -normal[1], -normal[2])
        cos_i = -cos_i
    mu = n1 / n2
    sin2_t = mu * mu * (1.0 - cos_i * cos_i)
    if sin2_t > 1.0:
        return None
    cos_t = math.sqrt(1.0 - sin2_t)
    k = mu * cos_i - cos_t
    return normalize(
        (
            mu * d[0] + k * normal[0],
            mu * d[1] + k * normal[1],
            mu * d[2] + k * normal[2],
        )
    )


@dataclass
class TraceResult:
    """单条光线的追迹结果。"""

    ok: bool
    points: List[Vec3] = field(default_factory=list)
    directions: List[Vec3] = field(default_factory=list)
    opl: float = 0.0
    image_point: Optional[Vec3] = None
    reason: str = ""

    @property
    def height(self) -> float:
        """像面上的横向高度 y (mm)。"""
        return self.image_point[1] if self.image_point else float("nan")

    @property
    def x(self) -> float:
        return self.image_point[0] if self.image_point else float("nan")

    @property
    def y(self) -> float:
        return self.image_point[1] if self.image_point else float("nan")


class OpticalSystem:
    """序列式（同轴）旋转对称光学系统。"""

    def __init__(
        self,
        name: str,
        surfaces: Sequence[Surface],
        fields: Sequence[float] = (0.0,),
        wavelengths: Sequence[float] = (materials.D_LINE,),
        object_distance: float = math.inf,
        object_height: float = 0.0,
        primary_wavelength: float = materials.D_LINE,
        image_distance: Optional[float] = None,
        description: str = "",
    ) -> None:
        self.name = name
        self.surfaces: List[Surface] = list(surfaces)
        self.fields: List[float] = list(fields)
        self.wavelengths: List[float] = list(wavelengths)
        self.object_distance = object_distance
        self.object_height = object_height
        self.primary_wavelength = primary_wavelength
        self.description = description
        self.object_index = 1.0
        self._image_distance_override = image_distance
        # 未标注光阑时约定第一面为孔径光阑，保证 F/# 与入瞳始终有定义
        if self.surfaces and not any(s.is_stop for s in self.surfaces):
            self.surfaces[0].is_stop = True
        # 近轴解算结果
        self.efl = 0.0
        self.bfl = 0.0
        self.image_distance = 0.0
        # 像面相对近轴焦点的整体偏移（"调焦"变量，用于平衡场曲）
        self.image_shift = 0.0
        self.entrance_pupil_z = 0.0
        self.entrance_pupil_radius = 0.0
        self.system_matrix: Matrix = IDENTITY
        self._solve()

    # ------------------------------------------------------------------
    # 基本属性
    # ------------------------------------------------------------------
    @property
    def is_infinite_object(self) -> bool:
        return not math.isfinite(self.object_distance)

    @property
    def stop_index(self) -> int:
        for i, s in enumerate(self.surfaces):
            if s.is_stop:
                return i
        # 未指定光阑时取第一面（常见于摄影镜头前光阑）
        return 0

    @property
    def f_number(self) -> float:
        if self.entrance_pupil_radius <= 0.0 or abs(self.efl) < 1e-12:
            return float("inf")
        return abs(self.efl) / (2.0 * self.entrance_pupil_radius)

    @property
    def numerical_aperture(self) -> float:
        """近轴数值孔径（像方）。"""
        if not math.isfinite(self.f_number):
            return 0.0
        return 1.0 / (2.0 * self.f_number)

    def surface_positions(self) -> List[float]:
        """各面顶点的轴向坐标 z (mm)，第一位为 0。"""
        z = [0.0]
        for s in self.surfaces[:-1]:
            z.append(z[-1] + s.thickness)
        return z

    @property
    def image_plane_z(self) -> float:
        return self.surface_positions()[-1] + self.image_distance + self.image_shift

    # ------------------------------------------------------------------
    # 近轴解算
    # ------------------------------------------------------------------
    def index_before(self, i: int, wavelength_nm: float) -> float:
        if i == 0:
            return self.object_index
        return materials.index_of(self.surfaces[i - 1].material, wavelength_nm)

    def index_after(self, i: int, wavelength_nm: float) -> float:
        return materials.index_of(self.surfaces[i].material, wavelength_nm)

    def refraction_matrix(self, i: int, wavelength_nm: float) -> Matrix:
        s = self.surfaces[i]
        n1 = self.index_before(i, wavelength_nm)
        n2 = self.index_after(i, wavelength_nm)
        phi = (n2 - n1) * s.curvature
        return ((1.0, 0.0), (-phi / n2, n1 / n2))

    def system_matrix_at(self, wavelength_nm: float, upto: Optional[int] = None) -> Matrix:
        """从第一面之前到第 ``upto`` 面之前的系统矩阵（不含该面折射）。

        ``upto=None`` 表示遍历全部面：包含每个面的折射，以及面间平移
        （不包含最后一面之后到像面的距离，像距由 ``image_distance`` 单独给出）。
        """
        m = IDENTITY
        limit = len(self.surfaces) if upto is None else upto
        for i in range(limit):
            m = mat_mul(self.refraction_matrix(i, wavelength_nm), m)
            if i < limit - 1:
                t = self.surfaces[i].thickness
                m = mat_mul(((1.0, t), (0.0, 1.0)), m)
        return m

    def _solve(self) -> None:
        wl = self.primary_wavelength
        m = self.system_matrix_at(wl)
        self.system_matrix = m

        # 有限共轭时把物面平移矩阵并入
        if not self.is_infinite_object:
            d = self.object_distance
            m_tot = mat_mul(m, ((1.0, d), (0.0, 1.0)))
        else:
            m_tot = m

        a, b = m_tot[0][0], m_tot[0][1]
        c, d_ = m_tot[1][0], m_tot[1][1]
        self.efl = -1.0 / c if abs(c) > 1e-15 else float("inf")
        self.bfl = -a / c if abs(c) > 1e-15 else float("inf")
        if self._image_distance_override is not None:
            self.image_distance = self._image_distance_override
        elif math.isfinite(self.bfl):
            self.image_distance = self.bfl
        else:
            # 无焦系统：保持名义后截距
            self.image_distance = self.surfaces[-1].thickness

        # 入瞳
        k = self.stop_index
        ms = self.system_matrix_at(wl, upto=k)
        m00, m01 = ms[0][0], ms[0][1]
        r_stop = self.surfaces[k].semi_diameter
        self.entrance_pupil_radius = abs(r_stop / m00) if abs(m00) > 1e-15 else r_stop
        self.entrance_pupil_z = m01 / m00 if abs(m00) > 1e-15 else 0.0
        self._zpos = self.surface_positions()

    def refocus(self) -> "OpticalSystem":
        """按当前参数重新解算（优化器每步调用）。"""
        self._solve()
        return self

    def paraxial_height(self, y0: float, u0: float, wavelength_nm: float) -> float:
        """近轴光线在像面上的高度，用于畸变等理想参考量。"""
        y, u = y0, u0
        n = self.object_index
        ns = len(self.surfaces)
        for i, s in enumerate(self.surfaces):
            n2 = self.index_after(i, wavelength_nm)
            phi = (n2 - n) * s.curvature
            u = (n * u - y * phi) / n2
            n = n2
            t = s.thickness if i < ns - 1 else self.image_distance + self.image_shift
            y = y + u * t
        return y

    def paraxial_focus(self, wavelength_nm: float) -> float:
        """近轴焦点相对最后一面的距离（负值表示焦点在面前方）。"""
        m = self.system_matrix_at(wavelength_nm)
        c = m[1][0]
        if abs(c) < 1e-15:
            return float("inf")
        return -m[0][0] / c

    # ------------------------------------------------------------------
    # 比例缩放：同种玻璃的相似变换，像差比例完全保持不变
    # ------------------------------------------------------------------
    def scale(self, factor: float) -> "OpticalSystem":
        """把整个系统按比例缩放（焦距随之缩放，像差比例不变）。"""
        for s in self.surfaces:
            if math.isfinite(s.radius):
                s.radius *= factor
            s.thickness *= factor
            s.semi_diameter *= factor
            if s.coefficients:
                s.coefficients = tuple(
                    a / factor ** (2 * i + 3) for i, a in enumerate(s.coefficients)
                )
        self._image_distance_override = None
        self.image_shift *= factor
        self._solve()
        return self

    def scale_to_efl(self, target_efl: float) -> "OpticalSystem":
        """缩放使近轴焦距等于目标值。"""
        if not math.isfinite(self.efl) or abs(self.efl) < 1e-9:
            return self
        return self.scale(target_efl / self.efl)

    # ------------------------------------------------------------------
    # 光线发射
    # ------------------------------------------------------------------
    def field_angle_rad(self, field: float) -> float:
        if self.is_infinite_object:
            return math.radians(field)
        return 0.0

    def launch(
        self, field: float, px: float, py: float, wavelength_nm: Optional[float] = None
    ) -> Tuple[Vec3, Vec3]:
        """按归一化入瞳坐标 (px, py) 发射光线（近轴瞄准）。"""
        r_ep = self.entrance_pupil_radius
        x_p = px * r_ep
        y_p = py * r_ep
        z_ep = self.entrance_pupil_z

        if self.is_infinite_object:
            theta = math.radians(field)
            d = (0.0, math.sin(theta), math.cos(theta))
            if z_ep < -1e-9:
                return (x_p, y_p, z_ep), d
            # 入瞳落在系统内部（罕见）：退化为在首面平面发射
            y0 = y_p + d[1] * (0.0 - z_ep)
            return (x_p, y0, 0.0), d

        d_obj = self.object_distance
        h = self.object_height
        # 物点 (0, h, -d_obj) 指向入瞳点 (x_p, y_p, z_ep)
        vx, vy, vz = x_p - 0.0, y_p - h, z_ep + d_obj
        ln = math.sqrt(vx * vx + vy * vy + vz * vz)
        return (0.0, h, -d_obj), (vx / ln, vy / ln, vz / ln)

    def chief_ray(self, field: float, wavelength_nm: Optional[float] = None) -> TraceResult:
        return self.trace_to_image(field, 0.0, 0.0, wavelength_nm)

    def marginal_ray(self, field: float, wavelength_nm: Optional[float] = None) -> TraceResult:
        return self.trace_to_image(field, 0.0, 1.0, wavelength_nm)

    # ------------------------------------------------------------------
    # 实光线追迹
    # ------------------------------------------------------------------
    def trace(
        self,
        origin: Vec3,
        direction: Vec3,
        wavelength_nm: float,
        check_aperture: bool = True,
        defocus: float = 0.0,
    ) -> TraceResult:
        """追迹一条光线直到像面（可附加离焦量）。"""
        p = origin
        d = normalize(direction)
        n = self.object_index
        points: List[Vec3] = [p]
        directions: List[Vec3] = [d]
        opl = 0.0
        ns = len(self.surfaces)
        zpos = self._zpos

        for i in range(ns):
            s = self.surfaces[i]
            # 各面在自身顶点局部坐标下求交，这里做 z 平移
            pl = (p[0], p[1], p[2] - zpos[i])
            hit = s.intersect(pl, d)
            if hit is None:
                return TraceResult(False, points, directions, opl, reason="未与面相交")
            t, ql, nrm = hit
            if t < -1e-9:
                return TraceResult(False, points, directions, opl, reason="交点位于光线后方")
            q = (ql[0], ql[1], ql[2] + zpos[i])
            opl += n * t
            p = q
            points.append(q)
            # 光阑口径只用于定义近轴入瞳，实光线不再二次拦光，
            # 这样出瞳采样与近轴瞄准保持自洽（详见 README）。
            if check_aperture and not s.is_stop:
                r2 = q[0] * q[0] + q[1] * q[1]
                lim = s.semi_diameter
                if r2 > lim * lim * (1.0 + 1e-9) + 1e-12:
                    return TraceResult(False, points, directions, opl, reason="超出通光口径")
            n2 = self.index_after(i, wavelength_nm)
            if s.is_mirror:
                cos_i = -(d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2])
                nn = nrm if cos_i > 0 else (-nrm[0], -nrm[1], -nrm[2])
                dd = (d[0] + 2 * nn[0], d[1] + 2 * nn[1], d[2] + 2 * nn[2])
                d = normalize(dd)
                n = n2
            else:
                nd = refract(d, nrm, n, n2)
                if nd is None:
                    return TraceResult(False, points, directions, opl, reason="全反射")
                d = nd
                n = n2
            directions.append(d)

        # 像面是垂直于光轴的平面 z = z_img，必须与平面求交，
        # 而不是沿光线方向平移“像距”长度（大视场下会低估像高）。
        z_img = self.image_plane_z + defocus
        if abs(d[2]) > 1e-12:
            travel = (z_img - p[2]) / d[2]
        else:
            travel = self.image_distance + defocus
        img = (p[0] + d[0] * travel, p[1] + d[1] * travel, z_img)
        opl += n * travel
        points.append(img)
        return TraceResult(True, points, directions, opl, img)

    def trace_to_image(
        self,
        field: float,
        px: float,
        py: float,
        wavelength_nm: Optional[float] = None,
        defocus: float = 0.0,
        check_aperture: bool = False,
    ) -> TraceResult:
        wl = self.primary_wavelength if wavelength_nm is None else wavelength_nm
        origin, direction = self.launch(field, px, py, wl)
        return self.trace(origin, direction, wl, check_aperture=check_aperture, defocus=defocus)

    # ------------------------------------------------------------------
    # 自动定口径
    # ------------------------------------------------------------------
    def autosize(self, margin: float = 1.06, samples: int = 16) -> "OpticalSystem":
        """按实际光线足迹设置各面通光半口径（避免非预期渐晕）。"""
        saved = [s.semi_diameter for s in self.surfaces]
        for s in self.surfaces:
            s.semi_diameter = 1e6
        heights = [0.0] * len(self.surfaces)
        for f in self.fields:
            for k in range(samples + 1):
                ang = math.pi * k / max(samples, 1) - math.pi / 2.0
                for rho in (1.0, 0.86, 0.5):
                    px, py = rho * math.cos(ang), rho * math.sin(ang)
                    res = self.trace_to_image(f, px, py)
                    if not res.ok:
                        continue
                    for i, q in enumerate(res.points[1:-1]):
                        h = math.hypot(q[0], q[1])
                        if h > heights[i]:
                            heights[i] = h
        for i, s in enumerate(self.surfaces):
            h = heights[i]
            if h <= 0.0:
                s.semi_diameter = saved[i]
            else:
                s.semi_diameter = round(h * margin, 3)
            if saved[i] > 0 and heights[i] <= 0:
                s.semi_diameter = saved[i]
        # 光阑口径保持不变（它定义了入瞳）
        k = self.stop_index
        if self.surfaces[k].is_stop and saved[k] > 0:
            self.surfaces[k].semi_diameter = saved[k]
        self._solve()
        return self

    # ------------------------------------------------------------------
    # 概要
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "efl": self.efl,
            "bfl": self.bfl,
            "image_distance": self.image_distance,
            "image_shift": self.image_shift,
            "image_plane_z": self.image_plane_z,
            "f_number": self.f_number,
            "numerical_aperture": self.numerical_aperture,
            "entrance_pupil_z": self.entrance_pupil_z,
            "entrance_pupil_radius": self.entrance_pupil_radius,
            "stop_index": self.stop_index,
            "fields": list(self.fields),
            "wavelengths": list(self.wavelengths),
            "surfaces": len(self.surfaces),
            "object_distance": self.object_distance,
        }

    def prescription_rows(self) -> List[dict]:
        z = self.surface_positions()
        rows = []
        for i, s in enumerate(self.surfaces):
            before = self.index_before(i, self.primary_wavelength)
            after = self.index_after(i, self.primary_wavelength)
            rows.append(
                {
                    "index": i + 1,
                    "radius": s.radius,
                    "thickness": s.thickness,
                    "material": s.material,
                    "semi_diameter": s.semi_diameter,
                    "conic": s.conic,
                    "z": z[i],
                    "n_before": before,
                    "n_after": after,
                    "is_stop": s.is_stop,
                    "label": s.label,
                }
            )
        return rows
