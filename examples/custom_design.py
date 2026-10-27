#!/usr/bin/env python3
"""示例：从零写一个镜头处方，分析它，再用 DLS 优化器把它改好。

这个脚本用到的全部 API 都在 ``raycraft`` 包里，不需要任何第三方库::

    python examples/custom_design.py

输出：
  1. 自定义三片式的初始像质（RMS 点列图）
  2. 运行阻尼最小二乘优化
  3. 优化后像质
  4. 把处方导出成 examples/out/custom_triplet.json
"""

from __future__ import annotations

import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from raycraft import analysis, materials  # noqa: E402
from raycraft.console import force_utf8  # noqa: E402
from raycraft.optimize import Optimizer, Variable, default_variables  # noqa: E402
from raycraft.surface import Surface  # noqa: E402
from raycraft.system import OpticalSystem  # noqa: E402

force_utf8()

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def build() -> OpticalSystem:
    """一个刻意"没调好"的柯克三片式：半径随手填，像质很差。"""
    surfaces = [
        Surface(20.0, 3.5, "N-SK16", 8.0, label="正片 前"),
        Surface(-160.0, 4.0, "air", 7.5, label="正片 后"),
        Surface(-40.0, 1.2, "F2", 4.0, label="负片 前"),
        Surface(18.0, 3.0, "air", 4.0, label="负片 后"),
        Surface(math.inf, 4.5, "air", 4.5, is_stop=True, label="光阑"),
        Surface(30.0, 3.0, "N-SK16", 8.0, label="正片 前"),
        Surface(-80.0, 0.0, "air", 8.0, label="正片 后"),
    ]
    return OpticalSystem(
        "custom_triplet",
        surfaces,
        fields=(0.0, 8.0, 16.0),
        wavelengths=(materials.F_LINE, materials.D_LINE, materials.C_LINE),
        description="示例：手写初值的柯克三片式",
    )


def report(label: str, s: OpticalSystem) -> None:
    print(f"\n{label}")
    print(f"  EFL = {s.efl:8.3f} mm   F/# = {s.f_number:.2f}   像距 = {s.image_distance + s.image_shift:.3f} mm")
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=5)
        print(
            f"  视场 {f:>5.1f}°   RMS = {sp['rms'] * 1000:8.2f} µm   "
            f"几何半径 = {sp['geo_radius'] * 1000:8.2f} µm   拦光 {sp['vignetted']}"
        )


def export(s: OpticalSystem, path: str) -> None:
    data = {
        "name": s.name,
        "description": s.description,
        "fields": s.fields,
        "wavelengths": s.wavelengths,
        "primary_wavelength": s.primary_wavelength,
        "image_shift": s.image_shift,
        "surfaces": [
            {
                "radius": None if not math.isfinite(x.radius) else round(x.radius, 6),
                "thickness": round(x.thickness, 6),
                "material": x.material,
                "semi_diameter": round(x.semi_diameter, 6),
                "conic": x.conic,
                "coefficients": list(x.coefficients),
                "is_stop": x.is_stop,
                "label": x.label,
            }
            for x in s.surfaces
        ],
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    print(f"\n处方已写入 {os.path.relpath(path, ROOT)}")


def main() -> int:
    s = build()
    # 光阑口径决定 F/#：这里缩放到 F/5.6
    cur = s.f_number
    s.surfaces[s.stop_index].semi_diameter *= cur / 5.6
    s._solve()
    s.autosize()
    s.surfaces[s.stop_index].semi_diameter *= s.f_number / 5.6
    s._solve()

    report("初始像质", s)

    variables = default_variables(s)
    print(f"\n优化变量 {len(variables)} 个：" +
          ", ".join(
              ("整组调焦" if v.param == "focus" else f"面{v.surface + 1}.{v.param}")
              for v in variables
          ))

    opt = Optimizer(s, variables, rings=3, target_efl=s.efl, target_f_number=5.6)
    res = opt.run(iterations=80, verbose=True)
    best = res["system"]
    best.autosize()

    report("优化后像质", best)

    mtf = analysis.mtf_analysis(best, points=24, n_rays=1200)
    i40 = min(range(len(mtf["freqs"])), key=lambda i: abs(mtf["freqs"][i] - 40.0))
    print(f"\n  截止频率 {mtf['cutoff']:.0f} lp/mm")
    for c in mtf["curves"]:
        print(
            f"  视场 {c['field']:>5.1f}°   MTF(40 lp/mm) "
            f"T={c['tangential'][i40]:.3f}  S={c['sagittal'][i40]:.3f}"
        )

    export(best, os.path.join(OUT_DIR, "custom_triplet.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
