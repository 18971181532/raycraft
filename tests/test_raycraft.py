"""RayCraft 测试套件（仅依赖标准库 unittest）。

运行::

    python -m unittest discover -s tests -v
    python run_tests.py
"""

from __future__ import annotations

import json
import math
import os
import sys
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import raycraft as rc  # noqa: E402
from raycraft import analysis, library, materials  # noqa: E402
from raycraft.optimize import (  # noqa: E402
    Optimizer,
    Variable,
    default_variables,
    multi_start,
    solve_linear,
)
from raycraft.surface import Surface  # noqa: E402
from raycraft.system import OpticalSystem, refract  # noqa: E402

D = materials.D_LINE


def _biconvex(radius=100.0, thickness=5.0, glass="N-BK7", semi=20.0) -> OpticalSystem:
    return OpticalSystem(
        "test-biconvex",
        [
            Surface(radius, thickness, glass, semi, is_stop=True, label="前"),
            Surface(-radius, 0.0, "air", semi, label="后"),
        ],
        fields=[0.0, 1.0],
        wavelengths=[D],
    )


# ==========================================================================
# 材料库
# ==========================================================================
class TestMaterials(unittest.TestCase):
    def test_catalog_not_empty(self):
        glasses = materials.list_glasses()
        self.assertGreaterEqual(len(glasses), 10)
        names = {g.name for g in glasses}
        self.assertIn("N-BK7", names)
        self.assertIn("F2", names)

    def test_bk7_reference_values(self):
        g = materials.glass("N-BK7")
        self.assertAlmostEqual(g.nd, 1.51680, places=4)
        self.assertAlmostEqual(g.abbe(), 64.17, places=1)

    def test_catalog_consistency(self):
        """库里每颗玻璃的 nd / Vd 必须与 Sellmeier 系数自洽。"""
        for g in materials.list_glasses():
            self.assertAlmostEqual(g.index(D), g.nd, places=4, msg=g.name)
            self.assertAlmostEqual(g.abbe(), g.vd, delta=0.35, msg=g.name)

    def test_index_dispersion_sign(self):
        """正常色散：短波折射率更高。"""
        g = materials.glass("N-SF6")
        self.assertGreater(g.index(materials.F_LINE), g.index(materials.C_LINE))
        self.assertGreater(g.dispersion(), 0.0)

    def test_air(self):
        self.assertAlmostEqual(materials.glass("air").index(500.0), 1.0, places=12)

    def test_unknown_glass(self):
        with self.assertRaises(KeyError):
            materials.glass("UNOBTAINIUM")


# ==========================================================================
# 面型几何
# ==========================================================================
class TestSurface(unittest.TestCase):
    def test_sphere_sag(self):
        s = Surface(radius=50.0)
        y = 10.0
        expect = 50.0 - math.sqrt(50.0 ** 2 - y * y)
        self.assertAlmostEqual(s.sag(y * y), expect, places=10)

    def test_plane_sag(self):
        self.assertEqual(Surface(radius=math.inf).sag(25.0), 0.0)

    def test_conic_parabola(self):
        """k=-1 应为严格抛物面 z = r²/(2R)。"""
        s = Surface(radius=80.0, conic=-1.0)
        y = 12.0
        self.assertAlmostEqual(s.sag(y * y), y * y / (2 * 80.0), places=10)

    def test_sag_matches_point(self):
        for surf in (
            Surface(radius=40.0),
            Surface(radius=-35.0, conic=-1.5),
            Surface(radius=60.0, conic=0.6),
            Surface(radius=50.0, conic=-2.0, coefficients=(2e-6, 3e-9)),
        ):
            for y in (2.0, 6.0, 11.0):
                p, nrm = surf.point_and_normal(y, 0.0)
                self.assertAlmostEqual(p[2], surf.sag(y * y), places=9)
                self.assertAlmostEqual(math.sqrt(sum(v * v for v in nrm)), 1.0, places=10)

    def test_intersect_on_sphere(self):
        s = Surface(radius=50.0)
        hit = s.intersect((0.0, 8.0, -20.0), (0.0, 0.0, 1.0))
        self.assertIsNotNone(hit)
        t, p, _ = hit
        self.assertAlmostEqual(p[1], 8.0, places=9)
        self.assertAlmostEqual(p[2], 50.0 - math.sqrt(2500.0 - 64.0), places=8)
        self.assertAlmostEqual(t, p[2] + 20.0, places=9)

    def test_intersect_asphere_consistency(self):
        s = Surface(radius=50.0, conic=-1.8, coefficients=(3e-6, -2e-9))
        y = 9.0
        hit = s.intersect((0.0, y, -30.0), (0.0, 0.0, 1.0))
        self.assertIsNotNone(hit)
        _, p, _ = hit
        self.assertAlmostEqual(p[0], 0.0, places=12)
        self.assertAlmostEqual(p[1], y, places=7)
        self.assertAlmostEqual(p[2], s.sag(y * y), places=7)

    def test_miss_on_steep_ray(self):
        """横向掠过球面之外应当判为不相交。"""
        s = Surface(radius=10.0, semi_diameter=3.0)
        self.assertIsNone(s.intersect((0.0, 40.0, -10.0), (0.0, 0.0, 1.0)))


# ==========================================================================
# 近轴矩阵与折射
# ==========================================================================
class TestRefraction(unittest.TestCase):
    def test_snell_flat_surface(self):
        """平面折射：n1 sinθ1 = n2 sinθ2。"""
        th = math.radians(30.0)
        d = (0.0, math.sin(th), math.cos(th))
        out = refract(d, (0.0, 0.0, -1.0), 1.0, 1.5)
        self.assertIsNotNone(out)
        self.assertAlmostEqual(math.sin(th), 1.5 * out[1], places=9)
        self.assertAlmostEqual(math.sqrt(sum(v * v for v in out)), 1.0, places=12)

    def test_total_internal_reflection(self):
        th = math.radians(60.0)  # > asin(1/1.5) ≈ 41.8°
        d = (0.0, math.sin(th), math.cos(th))
        self.assertIsNone(refract(d, (0.0, 0.0, -1.0), 1.5, 1.0))

    def test_no_deviation_at_normal_incidence(self):
        d = (0.0, 0.0, 1.0)
        out = refract(d, (0.0, 0.0, -1.0), 1.0, 1.7)
        self.assertIsNotNone(out)
        self.assertAlmostEqual(out[1], 0.0, places=12)
        self.assertAlmostEqual(out[2], 1.0, places=12)


class TestParaxial(unittest.TestCase):
    def test_thick_lens_efl(self):
        """EFL 必须与厚透镜解析式一致。"""
        n, r1, r2, t = 1.51680, 100.0, -100.0, 5.0
        phi = (n - 1.0) * (1.0 / r1 - 1.0 / r2 + (n - 1.0) * t / (n * r1 * r2))
        s = _biconvex(r1, t)
        self.assertAlmostEqual(s.efl, 1.0 / phi, places=4)

    def test_matrix_determinant(self):
        """系统矩阵行列式 = n_out / n_in。"""
        s = _biconvex()
        m = s.system_matrix
        self.assertAlmostEqual(m[0][0] * m[1][1] - m[0][1] * m[1][0], 1.0, places=10)

    def test_bfl_matches_paraxial_focus(self):
        s = _biconvex()
        self.assertAlmostEqual(s.bfl, s.paraxial_focus(D), places=8)

    def test_entrance_pupil_equals_stop_for_front_stop(self):
        """光阑在首面时，入瞳就是光阑本身。"""
        s = _biconvex(semi=20.0)
        self.assertAlmostEqual(s.entrance_pupil_z, 0.0, places=9)
        self.assertAlmostEqual(s.entrance_pupil_radius, 20.0, places=9)
        self.assertAlmostEqual(s.f_number, s.efl / 40.0, places=6)

    def test_scale_keeps_f_number(self):
        s = library.get("cooke_triplet")
        f0, n0 = s.efl, s.f_number
        s.scale(2.0)
        self.assertAlmostEqual(s.efl, 2.0 * f0, places=6)
        self.assertAlmostEqual(s.f_number, n0, places=6)

    def test_scale_to_efl(self):
        s = library.get("achromat")
        s.scale_to_efl(150.0)
        self.assertAlmostEqual(s.efl, 150.0, places=5)

    def test_paraxial_height_reference(self):
        """近轴像高应等于 f·tanθ（物在无穷远）。"""
        s = _biconvex()
        y = s.paraxial_height(0.0, math.tan(math.radians(1.0)), D)
        self.assertAlmostEqual(y, s.efl * math.tan(math.radians(1.0)), places=6)


# ==========================================================================
# 实光线追迹
# ==========================================================================
class TestTracing(unittest.TestCase):
    def test_axial_ray_lands_on_image_plane(self):
        s = _biconvex()
        res = s.trace_to_image(0.0, 0.0, 0.0)
        self.assertTrue(res.ok, res.reason)
        self.assertAlmostEqual(res.image_point[2], s.image_plane_z, places=9)
        self.assertAlmostEqual(res.image_point[1], 0.0, places=9)

    def test_all_rays_terminate_on_flat_plane(self):
        """像面是平面：任意视场、任意光瞳位置的光线都落在 z = image_plane_z。"""
        s = library.get("cooke_triplet")
        for f in s.fields:
            for py in (-1.0, -0.5, 0.0, 0.5, 1.0):
                res = s.trace_to_image(f, 0.0, py)
                if not res.ok:
                    continue
                self.assertAlmostEqual(res.image_point[2], s.image_plane_z, places=9)

    def test_image_height_matches_focal_length(self):
        """无限远物：像高 ≈ f·tanθ。"""
        s = library.get("cooke_triplet")
        f = 8.0
        chief = s.chief_ray(f)
        self.assertTrue(chief.ok)
        expect = s.efl * math.tan(math.radians(f))
        self.assertAlmostEqual(chief.image_point[1], expect, delta=0.02 * expect)

    def test_vignetting_detected(self):
        """故意把口径收到极小，远轴光线应被拦。"""
        s = library.get("cooke_triplet")
        for surf in s.surfaces:
            surf.semi_diameter = 0.4
        res = s.trace_to_image(16.0, 0.0, 1.0, check_aperture=True)
        self.assertFalse(res.ok)
        self.assertTrue(res.reason)

    def test_defocus_moves_plane(self):
        s = library.get("achromat")
        base = s.trace_to_image(0.0, 0.0, 0.3)
        shifted = s.trace_to_image(0.0, 0.0, 0.3, defocus=0.5)
        self.assertAlmostEqual(
            shifted.image_point[2] - base.image_point[2], 0.5, places=9
        )

    def test_asphere_intersection_is_smooth(self):
        """非球面系统的像高随光瞳坐标应连续（间接验证牛顿求交的稳定性）。"""
        s = library.get("asphere")
        ys = []
        for k in range(41):
            py = -1.0 + k / 20.0
            res = s.trace_to_image(0.0, 0.0, py)
            self.assertTrue(res.ok, res.reason)
            ys.append(res.image_point[1])
        for a, b in zip(ys, ys[1:]):
            self.assertLess(abs(b - a), 0.05)

    def test_autosize_covers_all_rays(self):
        for name in library.names():
            s = library.get(name)
            for f in s.fields:
                for py in (-1.0, -0.7, 0.0, 0.7, 1.0):
                    res = s.trace_to_image(f, 0.0, py, check_aperture=True)
                    self.assertTrue(res.ok, f"{name} f={f} py={py}: {res.reason}")


# ==========================================================================
# 像差分析
# ==========================================================================
class TestAnalysis(unittest.TestCase):
    def test_spot_on_axis_centroid(self):
        s = library.get("cooke_triplet")
        sp = analysis.spot_diagram(s, 0.0, rings=5)
        self.assertAlmostEqual(sum(sp["groups"][1]["x"]) / len(sp["groups"][1]["x"]), 0.0, places=6)
        self.assertLess(sp["rms"], 0.02)  # < 20 µm
        self.assertEqual(sp["vignetted"], 0)

    def test_spot_rms_monotonic_with_field(self):
        s = library.get("cooke_triplet")
        rms = [analysis.spot_diagram(s, f, rings=5)["rms"] for f in s.fields]
        self.assertLess(rms[0], rms[-1])

    def test_library_regression(self):
        """镜头库像质回归：与 library.REFERENCE 偏差不得超过 8%。"""
        for name in library.names():
            s = library.get(name)
            for f, ref in zip(s.fields, library.REFERENCE[name]):
                got = analysis.spot_diagram(s, f, rings=5)["rms"] * 1000.0
                self.assertAlmostEqual(
                    got, ref, delta=max(0.08 * ref, 0.2),
                    msg=f"{name} 视场 {f}°: {got:.2f} µm vs 基准 {ref} µm",
                )

    def test_achromat_reduces_chromatic_spread(self):
        """消色差物镜在 F 与 C 光的焦移应远小于同口径单透镜。"""
        ach = library.get("achromat")
        sin_ = library.get("singlet")
        dach = abs(ach.paraxial_focus(materials.F_LINE) - ach.paraxial_focus(materials.C_LINE))
        dsin = abs(sin_.paraxial_focus(materials.F_LINE) - sin_.paraxial_focus(materials.C_LINE))
        self.assertLess(dach, 0.2 * dsin)

    def test_distortion_near_zero_for_symmetric(self):
        s = library.get("achromat")
        for row in analysis.distortion(s):
            self.assertLess(abs(row["distortion"]), 0.5)  # %

    def test_ray_fan_odd_symmetry_on_axis(self):
        """轴上无彗差时，子午扇形图对光瞳坐标应为奇函数。"""
        s = library.get("cooke_triplet")
        fan = analysis.ray_fan(s, 0.0, n=21)
        vals = fan["tangential"][1]["value"]  # 主波长
        for a, b in zip(vals, reversed(vals)):
            self.assertAlmostEqual(a, -b, places=3)

    def test_mtf_normalisation(self):
        s = library.get("cooke_triplet")
        m = analysis.mtf_analysis(s, [0.0], points=16, n_rays=1200)
        self.assertAlmostEqual(m["curves"][0]["tangential"][0], 1.0, places=6)
        self.assertAlmostEqual(m["curves"][0]["sagittal"][0], 1.0, places=6)
        for v in m["curves"][0]["tangential"]:
            self.assertLessEqual(v, 1.0 + 1e-9)

    def test_mtf_axis_t_equals_s(self):
        s = library.get("cooke_triplet")
        m = analysis.mtf_analysis(s, [0.0], points=12, n_rays=1500)
        # 六角极坐标采样并非 90° 旋转对称，允许 2e-3 的采样残差
        for t, sg in zip(m["curves"][0]["tangential"], m["curves"][0]["sagittal"]):
            self.assertLess(abs(t - sg), 2e-3)

    def test_mtf_decreases_with_frequency(self):
        s = library.get("asphere")
        m = analysis.mtf_analysis(s, [0.0], points=24, n_rays=2500)
        t = m["curves"][0]["tangential"]
        self.assertGreater(t[1], t[-1])

    def test_diffraction_mtf_cutoff(self):
        lam, fno = 0.55e-3, 4.0
        f = [0.0, 1.0 / (lam * fno)]
        mtf = analysis.diffraction_mtf(fno, 550.0, f)
        self.assertAlmostEqual(mtf[0], 1.0, places=9)
        self.assertAlmostEqual(mtf[1], 0.0, places=2)

    def test_encircled_energy_monotonic(self):
        s = library.get("cooke_triplet")
        ee = analysis.encircled_energy(s, 0.0, n_rays=1500)
        vals = ee["fraction"]
        for a, b in zip(vals, vals[1:]):
            self.assertGreaterEqual(b, a - 1e-12)
        self.assertLessEqual(max(vals), 1.0 + 1e-9)

    def test_layout_touches_image_plane(self):
        s = library.get("asphere")
        lay = analysis.layout(s)
        for ray in lay["rays"]:
            self.assertAlmostEqual(ray["points"][-1][0], lay["image_plane"], places=9)

    def test_field_analysis_shapes(self):
        s = library.get("cooke_triplet")
        fa = analysis.field_analysis(s)
        self.assertEqual(len(fa["curves"]), len(s.fields))
        for c in fa["curves"]:
            self.assertIn("tangential_focus", c)
            self.assertIn("sagittal_focus", c)


# ==========================================================================
# 优化器
# ==========================================================================
class TestOptimizer(unittest.TestCase):
    def test_solve_linear(self):
        a = [[4.0, 1.0], [1.0, 3.0]]
        b = [1.0, 2.0]
        x = solve_linear(a, b)
        self.assertAlmostEqual(4 * x[0] + x[1], 1.0, places=10)
        self.assertAlmostEqual(x[0] + 3 * x[1], 2.0, places=10)

    def test_singular_matrix(self):
        self.assertEqual(solve_linear([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0]), [0.0, 0.0])

    def test_merit_decreases(self):
        s = library.get("cooke_triplet")
        opt = Optimizer(s, default_variables(s), rings=2)
        m0, m1 = opt.step()
        self.assertLessEqual(m1, m0 + 1e-12)

    def test_optimiser_improves_spot(self):
        s = _biconvex(radius=80.0)
        before = analysis.spot_diagram(s, 0.0, rings=3)["rms"]
        res = Optimizer(s, default_variables(s), rings=3, target_efl=s.efl).run(iterations=25)
        after = analysis.spot_diagram(res["system"], 0.0, rings=3)["rms"]
        self.assertLessEqual(after, before * 1.001)

    def test_focus_variable_moves_image_plane(self):
        s = library.get("cooke_triplet")
        base = s.surface_positions()[-1] + s.image_distance
        var = Variable(-1, "focus", -2.0, 2.0)
        var.set(s, -0.5)
        self.assertAlmostEqual(s.image_shift, -0.5, places=12)
        self.assertAlmostEqual(s.image_plane_z, base - 0.5, places=9)
        self.assertAlmostEqual(var.get(s), -0.5, places=12)

    def test_variable_clamping(self):
        s = library.get("cooke_triplet")
        var = Variable(0, "radius", 10.0, 50.0)
        var.set(s, 999.0)
        self.assertAlmostEqual(s.surfaces[0].radius, 50.0, places=12)

    def test_mirror_variable(self):
        s = library.get("cooke_triplet")
        var = Variable(0, "radius", 1.0, 100.0, mirror=(6, -1.0))
        var.set(s, 20.0)
        self.assertAlmostEqual(s.surfaces[0].radius, 20.0, places=12)
        self.assertAlmostEqual(s.surfaces[6].radius, -20.0, places=12)

    def test_multi_start_returns_best(self):
        s = library.get("singlet")
        res = multi_start(
            s, starts=2, iterations=6, target_efl=s.efl, target_f_number=s.f_number
        )
        self.assertIn("system", res)
        self.assertLessEqual(res["merit"], Optimizer(s, default_variables(s), rings=2).merit() * 1.5)

    def test_constraints_protect_structure(self):
        """优化后中心厚度/空气间隔不得越界。"""
        s = library.get("cooke_triplet")
        res = Optimizer(s, default_variables(s), rings=2, target_efl=s.efl).run(iterations=15)
        out = res["system"]
        for i, surf in enumerate(out.surfaces[:-1]):
            self.assertGreater(surf.thickness, 0.2, f"第 {i + 1} 面间距非正")
        self.assertGreater(out.image_distance + out.image_shift, 1.0)


# ==========================================================================
# 镜头库
# ==========================================================================
class TestLibrary(unittest.TestCase):
    def test_names_and_catalog(self):
        names = library.names()
        self.assertGreaterEqual(len(names), 4)
        cat = library.catalog()
        self.assertEqual(len(cat), len(names))
        for item in cat:
            self.assertGreater(item["efl"], 0)
            self.assertGreater(item["f_number"], 0.5)

    def test_build_all(self):
        for name in library.names():
            s = library.get(name)
            self.assertEqual(len(s.surfaces), len(s.prescription_rows()))
            self.assertTrue(math.isfinite(s.efl))
            self.assertTrue(math.isfinite(s.image_distance))

    def test_get_returns_deep_copy(self):
        a = library.get("singlet")
        a.surfaces[0].radius = 1.0
        b = library.get("singlet")
        self.assertNotAlmostEqual(b.surfaces[0].radius, 1.0, places=3)

    def test_unknown_name(self):
        with self.assertRaises(KeyError):
            library.get("no_such_lens")

    def test_asphere_has_terms(self):
        s = library.get("asphere")
        self.assertTrue(s.surfaces[0].coefficients)
        self.assertNotEqual(s.surfaces[0].conic, 0.0)

    def test_asphere_beats_spherical_equivalent(self):
        """同样 F/2.8 下非球面的轴上光斑应显著小于球面单透镜。"""
        asp = library.get("asphere")
        sph = _biconvex(radius=2.0 * 32.201 * 0.5168, thickness=6.0, semi=5.75)
        sphere_rms = analysis.spot_diagram(sph, 0.0, rings=4)["rms"]
        asphere_rms = analysis.spot_diagram(asp, 0.0, rings=4)["rms"]
        self.assertLess(asphere_rms * 6.0, sphere_rms)

    def test_prescription_rows_complete(self):
        s = library.get("cooke_triplet")
        rows = s.prescription_rows()
        self.assertEqual(rows[0]["z"], 0.0)
        for a, b in zip(rows, rows[1:]):
            self.assertGreater(b["z"], a["z"] - 1e-9)


# ==========================================================================
# HTTP 服务
# ==========================================================================
class TestServer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from raycraft import server

        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.httpd.daemon_threads = True
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()
        cls.thread.join(timeout=3)

    def _get(self, path):
        url = f"http://127.0.0.1:{self.port}{path}"
        with urllib.request.urlopen(url, timeout=20) as resp:
            return resp.status, resp.headers.get("Content-Type", ""), resp.read()

    def _json(self, path):
        status, ctype, body = self._get(path)
        self.assertEqual(status, 200)
        self.assertIn("application/json", ctype)
        return json.loads(body.decode("utf-8"))

    def test_index_served(self):
        status, ctype, body = self._get("/")
        self.assertEqual(status, 200)
        self.assertIn("text/html", ctype)
        self.assertIn(b"RayCraft", body)

    def test_static_assets(self):
        for path in ("/app.js", "/style.css"):
            status, _, body = self._get(path)
            self.assertEqual(status, 200)
            self.assertGreater(len(body), 100)

    def test_api_systems(self):
        data = self._json("/api/systems")
        self.assertEqual(len(data["systems"]), len(library.names()))

    def test_api_design(self):
        data = self._json("/api/design?name=cooke_triplet")
        for key in ("summary", "prescription", "layout", "spots", "distortion", "image_plane_z"):
            self.assertIn(key, data)
        self.assertGreater(len(data["prescription"]), 3)

    def test_api_spots(self):
        data = self._json("/api/spots?name=achromat")
        self.assertEqual(len(data["fields"]), len(library.get("achromat").fields))

    def test_api_fans_has_vmax(self):
        data = self._json("/api/fans?name=cooke_triplet")
        self.assertIn("vmax", data)
        self.assertGreater(data["vmax"], 0.0)

    def test_api_mtf(self):
        data = self._json("/api/mtf?name=asphere&rays=800")
        self.assertIn("cutoff", data)
        self.assertGreater(data["cutoff"], 100.0)
        self.assertEqual(len(data["curves"]), len(library.get("asphere").fields))

    def test_api_field(self):
        data = self._json("/api/field?name=cooke_triplet")
        self.assertIn("curves", data)
        self.assertIn("distortion", data)

    def test_api_glasses(self):
        data = self._json("/api/glasses")
        self.assertGreaterEqual(len(data["glasses"]), 10)

    def test_api_unknown_lens_404(self):
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            self._get("/api/design?name=nope")
        self.assertEqual(ctx.exception.code, 404)

    def test_api_unknown_route_404(self):
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            self._get("/api/nope")
        self.assertEqual(ctx.exception.code, 404)

    def test_path_traversal_blocked(self):
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            self._get("/../raycraft/__init__.py")
        self.assertEqual(ctx.exception.code, 404)


# ==========================================================================
# 命令行
# ==========================================================================
class TestCLI(unittest.TestCase):
    def _run(self, argv):
        import contextlib
        import io

        from raycraft import cli

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = cli.main(argv)
        return code, buf.getvalue()

    def test_systems(self):
        code, out = self._run(["systems"])
        self.assertEqual(code, 0)
        self.assertIn("RayCraft", out)
        for name in library.names():
            self.assertIn(name, out)

    def test_info(self):
        code, out = self._run(["info", "--name", "cooke_triplet"])
        self.assertEqual(code, 0)
        self.assertIn("EFL", out)
        self.assertIn("光阑", out)

    def test_spot(self):
        code, out = self._run(["spot", "--name", "achromat", "--rings", "3"])
        self.assertEqual(code, 0)
        self.assertIn("RMS", out)

    def test_mtf(self):
        code, out = self._run(["mtf", "--name", "asphere", "--points", "12", "--rays", "600"])
        self.assertEqual(code, 0)
        self.assertIn("lp/mm", out)

    def test_field(self):
        code, out = self._run(["field", "--name", "singlet"])
        self.assertEqual(code, 0)
        self.assertIn("畸变", out)

    def test_no_args_prints_help(self):
        code, out = self._run([])
        self.assertEqual(code, 0)
        self.assertIn("usage", out.lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
