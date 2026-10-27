"""零依赖 HTTP 服务：静态页面 + JSON 分析接口。

接口一览
--------
``GET /api/systems``            镜头库清单
``GET /api/design?name=&defocus=`` 近轴参数、处方表、2D 布局、点列摘要
``GET /api/spots?name=&defocus=``  多视场点列图（分波长）
``GET /api/fans?name=&defocus=``   子午/弧矢光线扇形图
``GET /api/mtf?name=&defocus=``    几何 MTF + 衍射极限
``GET /api/field?name=``           场曲与畸变
``GET /api/glasses``               玻璃目录
"""

from __future__ import annotations

import json
import math
import mimetypes
import os
import posixpath
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, Optional, Tuple

from .console import force_utf8

from . import analysis, library, materials

WEB_ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "web")

_CACHE: Dict[Tuple[str, str, str], Any] = {}
_LOCK = threading.Lock()


def _clean(obj: Any) -> Any:
    """把 inf/nan 换成 None，保证是合法 JSON。"""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


def design_payload(name: str, defocus: float) -> dict:
    s = library.get(name)
    spots = []
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=4, defocus=defocus)
        spots.append({"field": f, "rms": sp["rms"], "geo_radius": sp["geo_radius"]})
    rows = []
    for row in s.prescription_rows():
        sf = s.surfaces[row["index"] - 1]
        rows.append(
            {
                **row,
                "coefficients": list(sf.coefficients),
            }
        )
    return {
        "summary": s.summary(),
        "prescription": rows,
        "layout": analysis.layout(s, fields=s.fields),
        "spots": spots,
        "distortion": analysis.distortion(s),
        "image_plane_z": s.image_plane_z,
        "defocus": defocus,
    }


def spots_payload(name: str, defocus: float) -> dict:
    s = library.get(name)
    fields = []
    for f in s.fields:
        sp = analysis.spot_diagram(s, f, rings=6, defocus=defocus)
        fields.append(sp)
    radius = max((sp["geo_radius"] for sp in fields), default=1e-3)
    return {
        "fields": fields,
        "max_radius": max(radius, 1e-4),
        "defocus": defocus,
        "wavelengths": list(s.wavelengths),
    }


def fans_payload(name: str, defocus: float) -> dict:
    s = library.get(name)
    fields = [analysis.ray_fan(s, f, n=41, defocus=defocus) for f in s.fields]
    vmax = 1e-6
    for fd in fields:
        for key in ("tangential", "sagittal"):
            for g in fd[key]:
                for v in g["value"]:
                    vmax = max(vmax, abs(v))
    return {"fields": fields, "vmax": vmax, "defocus": defocus}


def mtf_payload(name: str, defocus: float, rays: int = 2200) -> dict:
    s = library.get(name)
    return analysis.mtf_analysis(s, points=61, n_rays=rays, defocus=defocus)


def field_payload(name: str) -> dict:
    s = library.get(name)
    return analysis.field_analysis(s)


def glasses_payload() -> dict:
    out = []
    for g in materials.list_glasses():
        out.append(
            {
                "name": g.name,
                "nd": g.index(materials.D_LINE),
                "vd": g.abbe(),
                "note": g.note,
            }
        )
    return {"glasses": out, "lines": {"F": materials.F_LINE, "d": materials.D_LINE, "C": materials.C_LINE}}


ROUTES = {
    "/api/systems": lambda q: {"systems": library.catalog()},
    "/api/design": lambda q: design_payload(q.get("name", ["cooke_triplet"])[0],
                                            float(q.get("defocus", ["0"])[0])),
    "/api/spots": lambda q: spots_payload(q.get("name", ["cooke_triplet"])[0],
                                          float(q.get("defocus", ["0"])[0])),
    "/api/fans": lambda q: fans_payload(q.get("name", ["cooke_triplet"])[0],
                                        float(q.get("defocus", ["0"])[0])),
    "/api/mtf": lambda q: mtf_payload(q.get("name", ["cooke_triplet"])[0],
                                      float(q.get("defocus", ["0"])[0]),
                                      int(q.get("rays", ["2200"])[0])),
    "/api/field": lambda q: field_payload(q.get("name", ["cooke_triplet"])[0]),
    "/api/glasses": lambda q: glasses_payload(),
}


class Handler(BaseHTTPRequestHandler):
    server_version = "RayCraft/1.0"

    def log_message(self, fmt: str, *args: Any) -> None:  # 静默
        pass

    # ---------------------------------------------------------------
    def _send(self, code: int, body: bytes, content_type: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def _json(self, obj: Any, code: int = 200) -> None:
        payload = json.dumps(_clean(obj), ensure_ascii=False).encode("utf-8")
        self._send(code, payload, "application/json; charset=utf-8")

    def do_GET(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        if path.startswith("/api/"):
            route = ROUTES.get(path)
            if route is None:
                self._json({"error": "未知接口", "path": path}, 404)
                return
            try:
                key = (path, query.get("name", [""])[0], json.dumps(query, sort_keys=True))
                with _LOCK:
                    cached = _CACHE.get(key)
                if cached is not None:
                    self._json(cached)
                    return
                data = route(query)
                with _LOCK:
                    if len(_CACHE) > 256:
                        _CACHE.clear()
                    _CACHE[key] = data
                self._json(data)
            except KeyError as exc:
                self._json({"error": str(exc)}, 404)
            except Exception as exc:  # noqa: BLE001
                self._json({"error": f"{type(exc).__name__}: {exc}"}, 500)
            return

        # 静态文件
        if path in ("/", "/index.html"):
            rel = "index.html"
        else:
            rel = posixpath.normpath(path.lstrip("/"))
        full = os.path.join(WEB_ROOT, rel)
        if not os.path.abspath(full).startswith(os.path.abspath(WEB_ROOT)) or not os.path.isfile(full):
            self._send(404, b"not found", "text/plain; charset=utf-8")
            return
        ctype = mimetypes.guess_type(full)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/javascript", "application/json"):
            ctype += "; charset=utf-8"
        with open(full, "rb") as fh:
            self._send(200, fh.read(), ctype)


def serve(host: str = "127.0.0.1", port: int = 8781) -> None:
    force_utf8()
    httpd = ThreadingHTTPServer((host, port), Handler)
    httpd.daemon_threads = True
    print(f"RayCraft 工作台已启动： http://{host}:{port}")
    print("按 Ctrl+C 停止")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止")
    finally:
        httpd.server_close()
