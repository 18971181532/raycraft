"""RayCraft — 光线追迹与镜头设计工作台。

零第三方依赖的序列式光学系统仿真库，覆盖：

* 几何：球面 / 二次曲面 / 偶次非球面，三维矢量斯涅尔折射
* 近轴：系统矩阵、焦距、后截距、入瞳、F/#
* 分析：点列图、光线扇形图、场曲、畸变、几何 MTF、能量集中度
* 优化：阻尼最小二乘（DLS）+ 多起点 + 对称约束

典型用法::

    from raycraft import library, analysis
    system = library.get("cooke_triplet")
    print(system.efl, system.f_number)
    spot = analysis.spot_diagram(system, 0.0)
    print(spot["rms"])
"""

from . import analysis, library, materials, optimize, surface, system  # noqa: F401
from .surface import Surface  # noqa: F401
from .system import OpticalSystem, refract, TraceResult  # noqa: F401

__version__ = "1.0.0"
__all__ = [
    "analysis",
    "library",
    "materials",
    "optimize",
    "surface",
    "system",
    "Surface",
    "OpticalSystem",
    "refract",
    "TraceResult",
    "__version__",
]
