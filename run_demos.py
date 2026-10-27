#!/usr/bin/env python3
"""RayCraft 演示脚本：跑一遍引擎的关键能力，把结果打印成报告。

用法::

    python run_demos.py            # 全部演示
    python run_demos.py 1 3 5      # 只跑指定编号
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import raycraft as rc  # noqa: E402
from raycraft import analysis, library, materials  # noqa: E402
from raycraft.console import force_utf8  # noqa: E402
from raycraft.optimize import Optimizer, Variable, default_variables, multi_start  # noqa: E402
from raycraft.surface import Surface  # noqa: E402
from raycraft.system import OpticalSystem  # noqa: E402

force_utf8()

BAR = "═" * 74


def title(n: int, text: str) -> None:
    print()
    print(BAR)
    print(f"  {n}. {text}")
    print(BAR)


# --------------------------------------------------------------------------
def demo_materials() -> None:
    title(1, "玻璃库：Sellmeier 色散与阿贝数")
    print(f"{'玻璃':<10}{'nd':>9}{'Vd(库)':>9}{'Vd(算)':>9}{'nF-nC':>10}   备注")
    for g in materials.list_glasses():
        calc = g.abbe()
        print(
            f"{g.name:<10}{g.nd:>9.4f}{g.vd:>9.2f}{calc:>9.2f}"
            f"{g.dispersion():>10.5f}   {g.note}"
        )
    print()
    bk7 = materials.glass("N-BK7")
    print("N-BK7 在 400–700 nm 的折射率曲线（每 50 nm）：")
    print("  " + "  ".join(f"{w:.0f}:{bk7.index(w):.4f}" for w in range(400, 701, 50)))


# --------------------------------------------------------------------------
def demo_paraxial() -> None:
    title(2, "近轴矩阵：从处方直接算出焦距与主面")
    for name in library.names():
        s = library.get(name)
        m = s.system_matrix
        print(
            f"{name:<14} EFL={s.efl:>9.3f} mm   BFL={s.bfl:>9.3f} mm   "
            f"F/#={s.f_number:>5.2f}   入瞳 z={s.entrance_pupil_z:>+8.3f} mm"
        )
        print(
            f"{'':<14} 系统矩阵 [[{m[0][0]:>9.5f}, {m[0][1]:>9.5f}],"
            f" [{m[1][0]:>11.7f}, {m[1][1]:>9.5f}]]  det={m[0][0] * m[1][1] - m[0][1] * m[1][0]:.9f}"
        )


# --------------------------------------------------------------------------
def demo_spherical_aberration() -> None:
    title(3, "球差：同一镜头的近轴焦点 vs 边缘光线焦点")
    s = library.get("singlet")
    print("平凸单透镜 f=100 mm F/4，各光瞳带的光线在轴上的交点位置：")
    for py in (0.2, 0.4, 0.6, 0.8, 1.0):
        res = s.trace_to_image(0.0, 0.0, py)
        p, d = res.points[-2], res.directions[-1]
        z = p[2] - p[1] * d[2] / d[1]
        print(
            f"  光瞳 {py:>4.1f}  轴向交点 z={z:>9.4f} mm   "
            f"相对近轴焦点 {z - s.image_plane_z:>+8.4f} mm"
        )
    edge = s.trace_to_image(0.0, 0.0, 1.0)
    pe, de = edge.points[-2], edge.directions[-1]
    z_edge = pe[2] - pe[1] * de[2] / de[1]
    lsa = s.image_plane_z - z_edge
    print(f"\n  轴向球差 LSA（近轴−边缘）= {lsa:.4f} mm，占焦距的 {lsa / s.efl * 100:.2f}%")
    print("  这就是 F/4 单透镜不校正球差时的典型量级（约 f/60），也是必须用双胶合或多片结构的原因。")


# --------------------------------------------------------------------------
def demo_spot_and_mtf() -> None:
    title(4, "点列图与 MTF")
    s = library.get("cooke_triplet")
    print("柯克三片式 f=50 mm F/5.6")
    print(f"{'视场/°':>8}{'RMS/µm':>10}{'几何半径/µm':>14}{'像高/mm':>10}")
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=5)
        h = s.chief_ray(f).image_point[1]
        print(
            f"{f:>8.1f}{sp['rms'] * 1000:>10.2f}{sp['geo_radius'] * 1000:>14.2f}{h:>10.3f}"
        )
    print()
    m = analysis.mtf_analysis(s, points=41, n_rays=2000)
    probes = [0.1, 0.2, 0.3, 0.5, 0.7]
    idx = [min(range(len(m["freqs"])), key=lambda i: abs(m["freqs"][i] - q * m["cutoff"]))
           for q in probes]
    print(f"  截止频率 {m['cutoff']:.1f} lp/mm（λ={s.primary_wavelength:.1f} nm, F/{s.f_number:.1f}）")
    print(f"{'视场/°':>8}{'方向':>6}" + "".join(f"{m['freqs'][i]:>9.1f}" for i in idx))
    for c in m["curves"]:
        print(f"{c['field']:>8.1f}{'T':>6}" + "".join(f"{c['tangential'][i]:>9.3f}" for i in idx))
        print(f"{'':>8}{'S':>6}" + "".join(f"{c['sagittal'][i]:>9.3f}" for i in idx))
    print(f"{'衍射极限':>7}{'D':>6}" + "".join(f"{m['diffraction_limit'][i]:>9.3f}" for i in idx))


# --------------------------------------------------------------------------
def demo_field_curvature() -> None:
    title(5, "场曲与畸变")
    for name in ("cooke_triplet", "singlet"):
        s = library.get(name)
        fa = analysis.field_analysis(s)
        dist = {round(d["field"], 6): d for d in fa["distortion"]}
        print(f"\n{name}")
        print(f"{'视场/°':>8}{'T 焦移/mm':>13}{'S 焦移/mm':>13}{'T-S 像散/mm':>15}{'畸变/%':>10}")
        for c in fa["curves"]:
            d = dist.get(round(c["field"], 6), {})
            print(
                f"{c['field']:>8.1f}{c['tangential_focus']:>13.4f}"
                f"{c['sagittal_focus']:>13.4f}"
                f"{c['tangential_focus'] - c['sagittal_focus']:>15.4f}"
                f"{d.get('distortion', 0.0):>10.3f}"
            )


# --------------------------------------------------------------------------
def demo_asphere() -> None:
    title(6, "非球面：一片顶一组球面")
    a = library.get("asphere")
    print(f"非球面单透镜 f={a.efl:.2f} mm F/{a.f_number:.1f}，633 nm")
    for i, s in enumerate(a.surfaces, 1):
        co = ", ".join(f"A{2 * k + 4}={v:.3e}" for k, v in enumerate(s.coefficients))
        print(f"  面 {i}: R={s.radius:>9.3f}  k={s.conic:>+8.4f}  {co}")
    rms = [analysis.spot_diagram(a, f, rings=6)["rms"] * 1000 for f in a.fields]
    airy = 1.22 * 0.633 * a.f_number
    print(f"\n  艾里半径 {airy:.2f} µm；各视场 RMS：" +
          "  ".join(f"{f:.1f}°→{r:.2f} µm" for f, r in zip(a.fields, rms)))

    # 与同参数的球面版本对比
    sph = OpticalSystem(
        "spherical-ref",
        [
            Surface(32.7463, 6.0, "N-SK16", 5.75, is_stop=True, label="前"),
            Surface(-47.2228, 2.5, "air", 5.968, label="后"),
        ],
        fields=[0.0],
        wavelengths=[633.0],
    )
    library.set_f_number(sph, 2.8)
    rs = analysis.spot_diagram(sph, 0.0, rings=6)["rms"] * 1000
    print(f"  同样两个面去掉非球面项后轴上 RMS = {rs:.1f} µm，相差 {rs / rms[0]:.0f} 倍")


# --------------------------------------------------------------------------
def _bad_triplet() -> OpticalSystem:
    """一个刻意"没算过"的柯克三片式：半径随手填，像质很差。"""
    s = OpticalSystem(
        "bad_triplet",
        [
            Surface(20.0, 3.5, "N-SK16", 8.0, label="正片 前"),
            Surface(-160.0, 4.0, "air", 7.5, label="正片 后"),
            Surface(-40.0, 1.2, "F2", 4.0, label="负片 前"),
            Surface(18.0, 3.0, "air", 4.0, label="负片 后"),
            Surface(math.inf, 4.5, "air", 4.5, is_stop=True, label="光阑"),
            Surface(30.0, 3.0, "N-SK16", 8.0, label="正片 前"),
            Surface(-80.0, 0.0, "air", 8.0, label="正片 后"),
        ],
        fields=(0.0, 8.0, 16.0),
        wavelengths=(materials.F_LINE, materials.D_LINE, materials.C_LINE),
        description="演示用：半径随手填的柯克三片式",
    )
    library.set_f_number(s, 5.6)
    s.autosize()
    library.set_f_number(s, 5.6)
    return s


def _line(tag: str, s: OpticalSystem) -> str:
    rms = [analysis.spot_diagram(s, f, rings=4)["rms"] * 1000 for f in s.fields]
    m = analysis.mtf_analysis(s, points=20, n_rays=1200)
    i40 = min(range(len(m["freqs"])), key=lambda i: abs(m["freqs"][i] - 40.0))
    t40 = m["curves"][-1]["tangential"][i40]
    return (
        f"  {tag:<6} EFL={s.efl:>6.2f}  F/{s.f_number:.2f}  "
        f"RMS(µm)=" + "/".join(f"{v:7.1f}" for v in rms) +
        f"   边缘 MTF@40={t40:5.3f}"
    )


def demo_optimizer() -> None:
    title(7, "自动优化：把随手填半径的三片式压成可用镜头")
    s = _bad_triplet()
    print("  视场 " + " / ".join(f"{f:g}°" for f in s.fields) + "，目标 EFL=50 mm, F/5.6")
    print(_line("起始", s))

    variables = default_variables(s)
    print(f"\n  变量 {len(variables)} 个：" + ", ".join(
        "整组调焦" if v.param == "focus" else f"面{v.surface + 1}.{v.param}" for v in variables))

    res = Optimizer(s, variables, rings=3, target_efl=50.0, target_f_number=5.6).run(
        iterations=80, verbose=True
    )
    best = res["system"]
    library.set_f_number(best, 5.6)
    best.autosize()
    library.set_f_number(best, 5.6)
    print()
    print(_line("结果", best))
    print(f"\n  迭代 {res['iterations']} 次，merit {res['history'][0]:.3e} → {res['merit']:.3e}"
          f"（降低 {(1 - res['merit'] / res['history'][0]) * 100:.1f}%）")
    for i, sf in enumerate(best.surfaces, 1):
        r = "∞" if not math.isfinite(sf.radius) else f"{sf.radius:9.3f}"
        print(f"      面 {i}  R={r:>10}  t={sf.thickness:>6.3f}  {sf.material}")


# --------------------------------------------------------------------------
def demo_multi_start() -> None:
    title(8, "多起点优化：同一个坏初值，不同的收敛结果")
    s = _bad_triplet()
    plain = Optimizer(s, default_variables(s), rings=3, target_efl=50.0,
                      target_f_number=5.6).run(iterations=25)
    print(f"  单起点（未扰动）      merit = {plain['merit']:.4e}")
    print()
    multi = multi_start(s, starts=5, iterations=25, target_efl=50.0, target_f_number=5.6)
    print(f"\n  {'起点':<6}{'扰动后 merit':>16}{'优化后 merit':>16}{'迭代':>6}")
    for t in multi["trials"]:
        mark = "  ← 最优" if t["start"] == multi["start"] else ""
        print(f"  {t['start']:<6}{t['start_merit']:>16.3e}{t['merit']:>16.3e}"
              f"{t['iterations']:>6}{mark}")
    print(f"\n  多起点相比单起点提升 {(1 - multi['merit'] / plain['merit']) * 100:.1f}%")


# --------------------------------------------------------------------------
def demo_prescription_io() -> None:
    title(9, "配方导入导出：把设计存成 JSON")
    import json

    s = library.get("achromat")
    data = {
        "name": s.name,
        "fields": s.fields,
        "wavelengths": s.wavelengths,
        "image_shift": s.image_shift,
        "surfaces": [
            {
                "radius": None if not math.isfinite(x.radius) else x.radius,
                "thickness": x.thickness,
                "material": x.material,
                "semi_diameter": x.semi_diameter,
                "conic": x.conic,
                "coefficients": list(x.coefficients),
                "is_stop": x.is_stop,
                "label": x.label,
            }
            for x in s.surfaces
        ],
    }
    text = json.dumps(data, indent=2, ensure_ascii=False)
    print(text[:420] + "\n  ...")
    rebuilt = OpticalSystem(
        data["name"],
        [
            Surface(
                radius=(math.inf if r["radius"] is None else r["radius"]),
                thickness=r["thickness"],
                material=r["material"],
                semi_diameter=r["semi_diameter"],
                conic=r["conic"],
                coefficients=tuple(r["coefficients"]),
                is_stop=r["is_stop"],
                label=r["label"],
            )
            for r in data["surfaces"]
        ],
        fields=data["fields"],
        wavelengths=data["wavelengths"],
    )
    rebuilt.image_shift = data["image_shift"]
    print(f"\n  重新构建后 EFL = {rebuilt.efl:.6f} mm（原 {s.efl:.6f} mm）")


# --------------------------------------------------------------------------
DEMOS = {
    1: demo_materials,
    2: demo_paraxial,
    3: demo_spherical_aberration,
    4: demo_spot_and_mtf,
    5: demo_field_curvature,
    6: demo_asphere,
    7: demo_optimizer,
    8: demo_multi_start,
    9: demo_prescription_io,
}


def main(argv: list[str]) -> int:
    picks = [int(a) for a in argv if a.isdigit()] or sorted(DEMOS)
    print(f"RayCraft v{rc.__version__} — 演示")
    for n in picks:
        fn = DEMOS.get(n)
        if fn is None:
            print(f"没有编号 {n} 的演示")
            continue
        fn()
    print()
    print(BAR)
    print("  演示结束。启动交互式工作台： python run.py serve")
    print(BAR)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
