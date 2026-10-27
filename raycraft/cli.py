"""命令行界面：查看处方、像差分析、优化与启动 Web 服务。"""

from __future__ import annotations

import argparse
import math
import sys
from typing import List, Optional

from . import analysis, library, materials
from .optimize import Optimizer, Variable, default_variables, multi_start


def _fmt(v: float, width: int = 11, digits: int = 4) -> str:
    if isinstance(v, float) and not math.isfinite(v):
        return f"{'∞':>{width}}"
    return f"{v:>{width}.{digits}f}"


def _hr(title: str = "") -> str:
    line = "─" * 78
    return f"{title}\n{line}" if title else line


def cmd_systems(args: argparse.Namespace) -> int:
    print(_hr("RayCraft 镜头库"))
    print(f"{'名称':<14}{'焦距/mm':>10}{'F/#':>7}{'面数':>6}  说明")
    for item in library.catalog():
        print(
            f"{item['name']:<14}{item['efl']:>10.2f}{item['f_number']:>7.2f}"
            f"{item['surfaces']:>6}  {item['title']}"
        )
    print()
    print("玻璃目录：" + "、".join(g.name for g in materials.list_glasses()))
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    s = library.get(args.name)
    info = s.summary()
    print(_hr(f"{s.name} — {info['description']}"))
    print(
        f"焦距 EFL={info['efl']:.4f} mm   后焦距 BFL={info['bfl']:.4f} mm   "
        f"F/#={info['f_number']:.3f}   NA≈{info['numerical_aperture']:.4f}"
    )
    print(
        f"入瞳位置 z={info['entrance_pupil_z']:+.4f} mm   入瞳半径 "
        f"{info['entrance_pupil_radius']:.4f} mm   光阑=第 {info['stop_index'] + 1} 面"
    )
    print(f"视场 {s.fields} °    波长 {[round(w, 2) for w in s.wavelengths]} nm")
    print()
    print(_hr("处方"))
    print(
        f"{'#':>3}{'半径/mm':>12}{'厚度/mm':>10}{'材料':>10}"
        f"{'半口径':>9}{'k':>8}  位置 z"
    )
    for row in s.prescription_rows():
        r = "∞" if not math.isfinite(row["radius"]) else f"{row['radius']:.4f}"
        stop = " [光阑]" if row["is_stop"] else ""
        print(
            f"{row['index']:>3}{r:>12}{row['thickness']:>10.3f}"
            f"{row['material']:>10}{row['semi_diameter']:>9.3f}"
            f"{row['conic']:>8.3f}  {row['z']:>8.3f}{stop}"
        )
    print()
    if abs(s.image_shift) > 1e-9:
        print(
            f"像面位于 z = {s.image_plane_z:.4f} mm"
            f"（近轴焦点 {s.surface_positions()[-1] + s.image_distance:.4f} mm"
            f" + 调焦 {s.image_shift:+.4f} mm）"
        )
    else:
        print(f"像面位于 z = {s.image_plane_z:.4f} mm（近轴焦点）")
    return 0


def cmd_spot(args: argparse.Namespace) -> int:
    s = library.get(args.name)
    print(_hr(f"{s.name} 点列图"))
    print(f"{'视场/°':>8}{'RMS/µm':>10}{'几何半径/µm':>14}{'光线数':>8}{'拦光':>6}")
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=args.rings)
        print(
            f"{f:>8.2f}{sp['rms'] * 1000:>10.2f}{sp['geo_radius'] * 1000:>14.2f}"
            f"{sp['rays']:>8}{sp['vignetted']:>6}"
        )
    print()
    print("波长分解（RMS/µm）：")
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=args.rings)
        parts = []
        for g in sp["groups"]:
            n = len(g["x"])
            rms = math.sqrt(sum(x * x + y * y for x, y in zip(g["x"], g["y"])) / max(n, 1))
            parts.append(f"{g['wavelength']:.0f}nm:{rms * 1000:.2f}")
        print(f"  视场 {f:>5.1f}°  " + "  ".join(parts))
    return 0


def cmd_mtf(args: argparse.Namespace) -> int:
    s = library.get(args.name)
    m = analysis.mtf_analysis(s, points=args.points, n_rays=args.rays)
    print(_hr(f"{s.name} 调制传递函数（几何法，截止 {m['cutoff']:.1f} lp/mm）"))
    freqs = m["freqs"]
    probe = [0.1, 0.2, 0.3, 0.5, 0.7]
    idx = [min(range(len(freqs)), key=lambda i: abs(freqs[i] - q * m["cutoff"])) for q in probe]
    head = "".join(f"{freqs[i]:>9.1f}" for i in idx)
    print(f"{'视场/°':>8}{'方向':>6}{head}")
    for c in m["curves"]:
        print(f"{c['field']:>8.1f}{'T':>6}" + "".join(f"{c['tangential'][i]:>9.3f}" for i in idx))
        print(f"{'':>8}{'S':>6}" + "".join(f"{c['sagittal'][i]:>9.3f}" for i in idx))
    print(f"{'极限':>8}{'D':>6}" + "".join(f"{m['diffraction_limit'][i]:>9.3f}" for i in idx))
    print("（T 子午 / S 弧矢 / D 圆孔衍射极限）")
    return 0


def cmd_field(args: argparse.Namespace) -> int:
    s = library.get(args.name)
    fa = analysis.field_analysis(s)
    print(_hr(f"{s.name} 场曲与畸变"))
    print(f"{'视场/°':>8}{'子午焦移/mm':>14}{'弧矢焦移/mm':>14}{'实际像高':>12}{'畸变/%':>10}")
    dist = {round(d["field"], 6): d for d in fa["distortion"]}
    for c in fa["curves"]:
        d = dist.get(round(c["field"], 6), {})
        print(
            f"{c['field']:>8.2f}{c['tangential_focus']:>14.4f}"
            f"{c['sagittal_focus']:>14.4f}{d.get('real_height', 0.0):>12.4f}"
            f"{d.get('distortion', 0.0):>10.3f}"
        )
    return 0


def cmd_optimize(args: argparse.Namespace) -> int:
    s = library.get(args.name)
    f0, n0 = s.efl, s.f_number
    print(_hr(f"{s.name} 优化（阻尼最小二乘）"))
    print(f"起始：EFL={f0:.3f} mm  F/#={n0:.3f}  merit={Optimizer(s, default_variables(s)).merit():.4f}")
    res = multi_start(
        s,
        starts=args.starts,
        iterations=args.iterations,
        target_efl=args.target_efl or f0,
        target_f_number=n0,
    )
    best = res["system"]
    best.autosize()
    print(f"完成：{res['iterations']} 次迭代  merit={res['merit']:.4f}")
    print(f"结果：EFL={best.efl:.3f} mm  F/#={best.f_number:.3f}")
    for f in best.fields:
        sp = analysis.spot_diagram(best, f, rings=4)
        print(f"  视场 {f:>5.1f}°  RMS = {sp['rms'] * 1000:.2f} µm")
    print()
    print(_hr("优化后处方"))
    for i, row in enumerate(best.prescription_rows()):
        r = "∞" if not math.isfinite(row["radius"]) else f"{row['radius']:.4f}"
        co = ""
        if best.surfaces[i].coefficients:
            co = "  A4=%.3e" % best.surfaces[i].coefficients[0]
        print(
            f"{row['index']:>3}{r:>12}{row['thickness']:>10.4f}"
            f"{row['material']:>10}{row['semi_diameter']:>9.3f}"
            f"{row['conic']:>8.4f}{co}"
        )
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    from .server import serve

    serve(host=args.host, port=args.port)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="raycraft",
        description="RayCraft — 光线追迹与镜头设计工作台",
    )
    sub = p.add_subparsers(dest="command")

    sub.add_parser("systems", help="列出镜头库").set_defaults(func=cmd_systems)

    info = sub.add_parser("info", help="查看镜头处方与近轴参数")
    info.add_argument("--name", default="cooke_triplet")
    info.set_defaults(func=cmd_info)

    spot = sub.add_parser("spot", help="点列图统计")
    spot.add_argument("--name", default="cooke_triplet")
    spot.add_argument("--rings", type=int, default=6)
    spot.set_defaults(func=cmd_spot)

    mtf = sub.add_parser("mtf", help="MTF 曲线")
    mtf.add_argument("--name", default="cooke_triplet")
    mtf.add_argument("--points", type=int, default=41)
    mtf.add_argument("--rays", type=int, default=2000)
    mtf.set_defaults(func=cmd_mtf)

    field = sub.add_parser("field", help="场曲与畸变")
    field.add_argument("--name", default="cooke_triplet")
    field.set_defaults(func=cmd_field)

    opt = sub.add_parser("optimize", help="运行 DLS 优化")
    opt.add_argument("--name", default="cooke_triplet")
    opt.add_argument("--iterations", type=int, default=40)
    opt.add_argument("--starts", type=int, default=3)
    opt.add_argument("--target-efl", type=float, default=None)
    opt.set_defaults(func=cmd_optimize)

    srv = sub.add_parser("serve", help="启动 Web 工作台")
    srv.add_argument("--host", default="127.0.0.1")
    srv.add_argument("--port", type=int, default=8781)
    srv.set_defaults(func=cmd_serve)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
