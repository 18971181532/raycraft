/* RayCraft 前端：手写 Canvas 绘图，无任何第三方库 */

const FIELD_COLORS = ["#1f5fa8", "#c05621", "#2f7d5a"];
const WL_COLORS = ["#2b6cb0", "#b7791f", "#c53030"];
const INK = "#17202b";
const MUTED = "#5d6b7a";
const LINE = "#dfe4ea";
const GRID = "#eef2f6";

const state = {
  name: "cooke_triplet",
  defocus: 0,
  tab: "layout",
  cache: new Map(),
  systems: [],
};

const $ = (id) => document.getElementById(id);

/* ------------------------------------------------------------------ */
/* 基础工具                                                             */
/* ------------------------------------------------------------------ */
function canvas2d(canvas, height) {
  const dpr = window.devicePixelRatio || 1;
  const w = Math.max(320, canvas.clientWidth || canvas.parentElement.clientWidth);
  const h = height || canvas.getAttribute("data-h") || 380;
  canvas.style.height = h + "px";
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  ctx.lineJoin = "round";
  ctx.lineCap = "round";
  return { ctx, w, h };
}

function line(ctx, x1, y1, x2, y2, color, width = 1, dash = null) {
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  if (dash) ctx.setLineDash(dash);
  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(x2, y2);
  ctx.stroke();
  ctx.restore();
}

function text(ctx, str, x, y, { color = MUTED, size = 11, align = "left", weight = 400 } = {}) {
  ctx.save();
  ctx.fillStyle = color;
  ctx.font = `${weight} ${size}px "Inter","Segoe UI","PingFang SC",system-ui,sans-serif`;
  ctx.textAlign = align;
  ctx.textBaseline = "middle";
  ctx.fillText(str, x, y);
  ctx.restore();
}

function fmt(v, digits = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  return Number(v).toFixed(digits);
}

async function api(path, params) {
  const qs = new URLSearchParams({ name: state.name, ...params });
  const res = await fetch(`/api/${path}?${qs}`);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return res.json();
}

function setStatus(msg, busy) {
  const el = $("status");
  el.textContent = msg;
  el.className = busy ? "busy" : "";
}

/* ------------------------------------------------------------------ */
/* 光学布局                                                             */
/* ------------------------------------------------------------------ */
function drawLayout(data) {
  const { ctx, w, h } = canvas2d($("layout-canvas"), 520);
  const L = data.layout;
  const pad = { l: 54, r: 26, t: 30, b: 42 };

  let ymax = 0;
  L.profiles.forEach((p) => p.points.forEach((pt) => (ymax = Math.max(ymax, Math.abs(pt[1])))));
  // 光线（尤其大视场主光线）可以远超镜片口径，纵向范围必须一起考虑
  L.rays.forEach((r) => r.points.forEach((pt) => (ymax = Math.max(ymax, Math.abs(pt[1])))));
  if (L.y_range) ymax = Math.max(ymax, Math.abs(L.y_range[0]), Math.abs(L.y_range[1]));
  ymax = Math.max(ymax * 1.06, 1e-3);
  const z0 = L.z_range[0], z1 = L.z_range[1];
  const availW = w - pad.l - pad.r;
  const availH = h - pad.t - pad.b;
  const scale = Math.min(availW / (z1 - z0), availH / (2 * ymax));

  const cx = pad.l + availW / 2;
  const cy = pad.t + availH / 2;
  const zmid = (z0 + z1) / 2;
  const P = (z, y) => [cx + (z - zmid) * scale, cy - y * scale];

  // 网格
  ctx.save();
  ctx.strokeStyle = GRID;
  ctx.lineWidth = 1;
  for (let y = -Math.floor(ymax / 5) * 5; y <= ymax; y += 5) {
    const [, py] = P(z0, y);
    if (py < pad.t || py > h - pad.b) continue;
    ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(w - pad.r, py); ctx.stroke();
  }
  ctx.restore();

  // 光轴
  const [, ay] = P(z0, 0);
  line(ctx, pad.l, ay, w - pad.r, ay, "#b9c3cd", 1, [5, 5]);
  text(ctx, "光轴", w - pad.r - 2, ay - 10, { align: "right", size: 10 });

  // 玻璃元件
  const prof = L.profiles;
  for (let i = 0; i < prof.length - 1; i++) {
    const mat = prof[i].material;
    if (!mat || mat.toLowerCase() === "air") continue;
    const a = prof[i].points, b = prof[i + 1].points;
    ctx.beginPath();
    a.forEach(([z, y], k) => { const [px, py] = P(z, y); k ? ctx.lineTo(px, py) : ctx.moveTo(px, py); });
    for (let k = b.length - 1; k >= 0; k--) { const [px, py] = P(b[k][0], b[k][1]); ctx.lineTo(px, py); }
    ctx.closePath();
    const grad = ctx.createLinearGradient(0, cy - ymax * scale, 0, cy + ymax * scale);
    grad.addColorStop(0, "rgba(168, 205, 240, .55)");
    grad.addColorStop(0.5, "rgba(206, 228, 248, .40)");
    grad.addColorStop(1, "rgba(168, 205, 240, .55)");
    ctx.fillStyle = grad;
    ctx.fill();
    ctx.strokeStyle = "rgba(31, 95, 168, .85)";
    ctx.lineWidth = 1.6;
    ctx.stroke();
  }

  // 光线
  const fields = [...new Set(L.rays.map((r) => r.field))].sort((a, b) => a - b);
  L.rays.forEach((ray) => {
    const idx = fields.indexOf(ray.field);
    ctx.beginPath();
    ray.points.forEach(([z, y], k) => {
      const [px, py] = P(z, y);
      k ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
    });
    ctx.strokeStyle = FIELD_COLORS[idx % FIELD_COLORS.length];
    ctx.globalAlpha = 0.62;
    ctx.lineWidth = 1.05;
    ctx.stroke();
    ctx.globalAlpha = 1;
  });

  // 像面
  const [ix] = P(L.image_plane, 0);
  line(ctx, ix, pad.t - 6, ix, h - pad.b + 6, "#c53030", 1.4, [6, 4]);
  text(ctx, "像面", ix + 6, pad.t - 12, { color: "#c53030", size: 11 });

  // 光阑
  const stop = prof.find((p) => p.is_stop);
  if (stop) {
    const [sx] = P(stop.z, 0);
    const [, ytop] = P(stop.z, stop.semi_diameter);
    const [, ybot] = P(stop.z, -stop.semi_diameter);
    ctx.save();
    ctx.strokeStyle = "#8a6d3b";
    ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(sx, ytop - 8); ctx.lineTo(sx, ytop + 16); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(sx, ybot + 8); ctx.lineTo(sx, ybot - 16); ctx.stroke();
    ctx.restore();
    text(ctx, "光阑", sx + 6, ybot - 12, { color: "#8a6d3b", size: 10 });
  }

  // 标尺
  const bar = 10 * scale;
  const by = h - 18;
  line(ctx, pad.l, by, pad.l + bar, by, INK, 1.6);
  line(ctx, pad.l, by - 3, pad.l, by + 3, INK, 1.6);
  line(ctx, pad.l + bar, by - 3, pad.l + bar, by + 3, INK, 1.6);
  text(ctx, "10 mm", pad.l + bar + 8, by, { size: 10.5 });

  // 视场图例
  $("layout-legend").innerHTML = fields
    .map((f, i) => `<span style="color:${FIELD_COLORS[i % 3]}"><i></i>视场 ${fmt(f, 1)}°</span>`)
    .join("") + `<span style="color:#8a6d3b"><i style="height:9px;width:3px"></i>孔径光阑</span>`;

  $("layout-note").textContent =
    `共 ${data.prescription.length} 面 · 像面位于 z=${fmt(data.image_plane_z, 2)} mm · ` +
    `全口径 ${fmt(2 * Math.max(...L.profiles.map(p => p.semi_diameter)), 1)} mm`;
}

/* ------------------------------------------------------------------ */
/* 点列图                                                               */
/* ------------------------------------------------------------------ */
function drawSpot(data) {
  const cv = $("spot-canvas");
  const { ctx, w, h } = canvas2d(cv, 330);
  const n = data.fields.length;
  const gap = 14;
  const pw = (w - gap * (n - 1)) / n;
  const rmax = data.max_radius * 1000; // µm
  const wls = data.wavelengths || [];

  data.fields.forEach((f, i) => {
    const x0 = i * (pw + gap);
    const size = Math.min(pw, h - 46);
    const cx = x0 + pw / 2;
    const cy = 26 + size / 2;
    const R = size / 2 - 8;
    const sc = R / (rmax * 1.08);

    ctx.save();
    ctx.strokeStyle = LINE;
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, Math.PI * 2); ctx.stroke();
    ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.arc(cx, cy, R * 0.5, 0, Math.PI * 2); ctx.stroke();
    ctx.restore();
    line(ctx, cx - R, cy, cx + R, cy, "#dfe4ea", 1);
    line(ctx, cx, cy - R, cx, cy + R, "#dfe4ea", 1);

    f.groups.forEach((g, gi) => {
      const idx = wls.length > 1 ? Math.min(gi, 2) : 2;
      ctx.fillStyle = WL_COLORS[idx % 3];
      ctx.globalAlpha = 0.72;
      for (let k = 0; k < g.x.length; k++) {
        ctx.beginPath();
        ctx.arc(cx + g.x[k] * 1000 * sc, cy - g.y[k] * 1000 * sc, 1.35, 0, Math.PI * 2);
        ctx.fill();
      }
      ctx.globalAlpha = 1;
    });

    text(ctx, `${fmt(f.field, 1)}°`, cx, 13, { align: "center", size: 12, color: INK, weight: 600 });
    text(ctx, `RMS ${fmt(f.rms * 1000, 2)} µm`, cx, h - 12, { align: "center", size: 11 });
  });

  $("spot-legend").innerHTML =
    wls.map((wl, i) => `<span style="color:${WL_COLORS[(wls.length > 1 ? i : 2) % 3]}">
      <i class="dot"></i>${fmt(wl, 1)} nm</span>`).join("") +
    `<span>外圈 = 最大弥散半径 ${fmt(rmax, 1)} µm</span>`;
}

/* ------------------------------------------------------------------ */
/* 光线扇形图                                                           */
/* ------------------------------------------------------------------ */
function fanPanel(ctx, data, key, x0, y0, pw, ph, vmaxIn, title) {
  const pad = { l: 46, r: 14, t: 26, b: 30 };
  const iw = pw - pad.l - pad.r;
  const ih = ph - pad.t - pad.b;
  const vmax = Math.max(vmaxIn || 0, 1e-3);
  const X = (p) => pad.l + x0 + (p + 1) / 2 * iw;
  const Y = (v) => pad.t + y0 + ih / 2 - (v / vmax) * (ih / 2) * 0.94;

  ctx.save();
  ctx.fillStyle = "#fcfdfe";
  ctx.fillRect(x0 + pad.l, y0 + pad.t, iw, ih);
  ctx.strokeStyle = LINE;
  ctx.lineWidth = 1;
  ctx.strokeRect(x0 + pad.l, y0 + pad.t, iw, ih);
  ctx.restore();

  line(ctx, X(-1), Y(0), X(1), Y(0), "#c9d2dc", 1);
  line(ctx, X(0), Y(-vmax), X(0), Y(vmax), "#eef2f6", 1);

  text(ctx, title, x0 + pad.l, y0 + 13, { size: 11.5, color: INK, weight: 600 });
  text(ctx, `+${fmt(vmax, 1)}`, x0 + 4, Y(vmax), { size: 9.5 });
  text(ctx, `-${fmt(vmax, 1)}`, x0 + 4, Y(-vmax), { size: 9.5 });
  text(ctx, "-1", X(-1), y0 + pad.t + ih + 12, { align: "center", size: 9.5 });
  text(ctx, "0", X(0), y0 + pad.t + ih + 12, { align: "center", size: 9.5 });
  text(ctx, "+1", X(1), y0 + pad.t + ih + 12, { align: "center", size: 9.5 });
  text(ctx, "光瞳坐标", x0 + pad.l + iw, y0 + pad.t + ih + 12, { align: "right", size: 9.5 });

  data[key].forEach((g, gi) => {
    ctx.beginPath();
    g.pupil.forEach((p, k) => {
      const px = X(p), py = Y(g.value[k]);
      k ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
    });
    ctx.strokeStyle = WL_COLORS[(data[key].length > 1 ? gi : 2) % 3];
    ctx.lineWidth = 1.5;
    ctx.stroke();
  });
}

function drawFan(data) {
  const { ctx, w, h } = canvas2d($("fan-canvas"), 420);
  const rows = data.fields.length;
  const ph = h / rows;
  const vmax = data.vmax;
  data.fields.forEach((fd, i) => {
    const y0 = i * ph;
    // 每个视场单独定标，否则小视场的曲线会被大视场的边缘光线压成一条直线
    let vm = 0;
    ["tangential", "sagittal"].forEach((k) =>
      fd[k].forEach((g) => g.value.forEach((v) => (vm = Math.max(vm, Math.abs(v))))));
    vm = Math.max(vm * 1.08, 1, (vmax || 0) * 0.12);
    fanPanel(ctx, fd, "tangential", 0, y0, w / 2, ph, vm, `子午 T · 视场 ${fmt(fd.field, 1)}°`);
    fanPanel(ctx, fd, "sagittal", w / 2, y0, w / 2, ph, vm, `弧矢 S · 视场 ${fmt(fd.field, 1)}°`);
    if (i > 0) line(ctx, 0, y0, w, y0, "#eef2f6", 1);
  });
  const wls = data.fields[0].tangential;
  $("fan-legend").innerHTML =
    wls.map((g, i) => `<span style="color:${WL_COLORS[(wls.length > 1 ? i : 2) % 3]}">
      <i></i>${fmt(g.wavelength, 1)} nm</span>`).join("") +
    `<span>纵轴单位 µm，离焦 ${fmt(data.defocus, 2)} mm</span>`;
}

/* ------------------------------------------------------------------ */
/* MTF                                                                 */
/* ------------------------------------------------------------------ */
function drawMtf(data) {
  const canvas = $("mtf-canvas");
  const { ctx, w, h } = canvas2d(canvas, 440);
  const pad = { l: 58, r: 132, t: 24, b: 48 };
  const iw = w - pad.l - pad.r;
  const ih = h - pad.t - pad.b;
  const X = (f) => pad.l + (f / data.cutoff) * iw;
  const Y = (v) => pad.t + ih - v * ih;

  // 网格与刻度
  ctx.save();
  ctx.strokeStyle = GRID;
  ctx.lineWidth = 1;
  for (let k = 0; k <= 5; k++) {
    const y = Y(k / 5);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + iw, y); ctx.stroke();
    text(ctx, fmt(k / 5, 1), pad.l - 10, y, { align: "right", size: 10 });
  }
  for (let k = 0; k <= 5; k++) {
    const x = pad.l + (k / 5) * iw;
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ih); ctx.stroke();
    text(ctx, fmt((k / 5) * data.cutoff, 0), x, pad.t + ih + 16, { align: "center", size: 10 });
  }
  ctx.restore();
  line(ctx, pad.l, Y(0), pad.l + iw, Y(0), "#c9d2dc", 1);
  line(ctx, pad.l, pad.t, pad.l, Y(0), "#c9d2dc", 1);
  text(ctx, "空间频率 lp/mm", pad.l + iw, pad.t + ih + 34, { align: "right", size: 11 });
  text(ctx, "MTF", pad.l - 10, pad.t - 12, { align: "right", size: 11 });

  // 衍射极限
  ctx.save();
  ctx.beginPath();
  data.diffraction_limit.forEach((v, i) => {
    const x = X(data.freqs[i]), y = Y(v);
    i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
  });
  ctx.strokeStyle = "#8b98a6";
  ctx.lineWidth = 1.4;
  ctx.setLineDash([5, 4]);
  ctx.stroke();
  ctx.restore();

  // 各视场
  data.curves.forEach((c, i) => {
    const color = FIELD_COLORS[i % 3];
    [["tangential", []], ["sagittal", [5, 4]]].forEach(([key, dash]) => {
      ctx.save();
      ctx.beginPath();
      c[key].forEach((v, k) => {
        const x = X(data.freqs[k]), y = Y(v);
        k ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      });
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.8;
      if (dash.length) ctx.setLineDash(dash);
      ctx.stroke();
      ctx.restore();
    });
  });

  // 右侧标注
  const bx = pad.l + iw + 18;
  data.curves.forEach((c, i) => {
    const y = pad.t + 12 + i * 54;
    text(ctx, `视场 ${fmt(c.field, 1)}°`, bx, y, { color: FIELD_COLORS[i % 3], size: 12, weight: 600 });
    text(ctx, `T  ${fmt(c.tangential[Math.floor(c.tangential.length * 0.3)], 3)}`, bx, y + 16, { size: 10.5 });
    text(ctx, `S  ${fmt(c.sagittal[Math.floor(c.sagittal.length * 0.3)], 3)}`, bx, y + 31, { size: 10.5 });
  });
  text(ctx, `截止 ${fmt(data.cutoff, 0)} lp/mm`, bx, pad.t + 12 + data.curves.length * 54 + 4, { size: 10.5 });
  text(ctx, `离焦 ${fmt(data.defocus, 2)} mm`, bx, pad.t + 12 + data.curves.length * 54 + 20, { size: 10.5 });

  $("mtf-legend").innerHTML =
    data.curves.map((c, i) => `<span style="color:${FIELD_COLORS[i % 3]}"><i></i>视场 ${fmt(c.field, 1)}°</span>`).join("") +
    `<span style="color:#8b98a6"><i style="border-top:2px dashed #8b98a6;background:none;height:0"></i>衍射极限</span>` +
    `<span>实线 T（子午）/ 虚线 S（弧矢）</span>`;
}

/* ------------------------------------------------------------------ */
/* 场曲与畸变                                                           */
/* ------------------------------------------------------------------ */
function drawField(data) {
  const curves = data.curves;
  const dist = data.distortion;
  const fmax = Math.max(...curves.map((c) => Math.abs(c.field))) || 1;

  // ---- 场曲 ----
  {
    const { ctx, w, h } = canvas2d($("fc-canvas"), 380);
    const pad = { l: 62, r: 26, t: 26, b: 46 };
    const iw = w - pad.l - pad.r, ih = h - pad.t - pad.b;
    const vals = curves.flatMap((c) => [c.tangential_focus, c.sagittal_focus]);
    const fmin = Math.min(...vals), fmaxv = Math.max(...vals);
    const span = Math.max(Math.abs(fmin), Math.abs(fmaxv), 0.02) * 1.15;
    const X = (v) => pad.l + iw / 2 + (v / span) * (iw / 2) * 0.9;
    const Y = (f) => pad.t + ih - (f / fmax) * ih;

    ctx.save(); ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (let k = 0; k <= 4; k++) {
      const y = pad.t + (k / 4) * ih;
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + iw, y); ctx.stroke();
    }
    ctx.restore();
    line(ctx, X(0), pad.t, X(0), pad.t + ih, "#c9d2dc", 1);
    line(ctx, pad.l, Y(fmax), pad.l + iw, Y(fmax), "#c9d2dc", 1);

    [["tangential_focus", "T 子午"], ["sagittal_focus", "S 弧矢"]].forEach(([key, label], i) => {
      ctx.beginPath();
      curves.forEach((c, k) => {
        const x = X(c[key]), y = Y(c.field);
        k ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      });
      ctx.strokeStyle = i ? "#2f7d5a" : "#c05621";
      ctx.lineWidth = 2;
      ctx.stroke();
      curves.forEach((c) => {
        ctx.beginPath(); ctx.arc(X(c[key]), Y(c.field), 3, 0, Math.PI * 2);
        ctx.fillStyle = i ? "#2f7d5a" : "#c05621"; ctx.fill();
      });
      text(ctx, label, pad.l + 6 + i * 78, pad.t - 12, { color: i ? "#2f7d5a" : "#c05621", size: 11, weight: 600 });
    });

    for (let k = 0; k <= 4; k++) {
      const f = (k / 4) * fmax;
      text(ctx, fmt(f, 1) + "°", pad.l - 10, Y(f), { align: "right", size: 10 });
    }
    ["-" + fmt(span, 2), "0", "+" + fmt(span, 2)].forEach((s, i) => {
      text(ctx, s, pad.l + (i / 2) * iw, pad.t + ih + 18, { align: "center", size: 10 });
    });
    text(ctx, "焦面偏移 / mm（相对近轴像面）", pad.l + iw / 2, h - 10, { align: "center", size: 11 });
  }

  // ---- 畸变 ----
  {
    const { ctx, w, h } = canvas2d($("dist-canvas"), 380);
    const pad = { l: 62, r: 26, t: 26, b: 46 };
    const iw = w - pad.l - pad.r, ih = h - pad.t - pad.b;
    const dmax = Math.max(0.5, ...dist.map((d) => Math.abs(d.distortion))) * 1.15;
    const X = (v) => pad.l + iw / 2 + (v / dmax) * (iw / 2) * 0.9;
    const Y = (f) => pad.t + ih - (Math.abs(f) / fmax) * ih;

    ctx.save(); ctx.strokeStyle = GRID; ctx.lineWidth = 1;
    for (let k = 0; k <= 4; k++) {
      const y = pad.t + (k / 4) * ih;
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + iw, y); ctx.stroke();
    }
    ctx.restore();
    line(ctx, X(0), pad.t, X(0), pad.t + ih, "#c9d2dc", 1);

    ctx.beginPath();
    dist.forEach((d, k) => {
      const x = X(d.distortion), y = Y(d.field);
      k ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
    });
    ctx.strokeStyle = "#1f5fa8";
    ctx.lineWidth = 2;
    ctx.stroke();
    dist.forEach((d) => {
      ctx.beginPath(); ctx.arc(X(d.distortion), Y(d.field), 3, 0, Math.PI * 2);
      ctx.fillStyle = "#1f5fa8"; ctx.fill();
    });

    dist.forEach((d) => {
      text(ctx, fmt(Math.abs(d.field), 1) + "°", pad.l - 10, Y(d.field), { align: "right", size: 10 });
      text(ctx, fmt(d.distortion, 2) + "%", X(d.distortion), Y(d.field) - 12,
        { align: "center", size: 10, color: "#1f5fa8" });
    });
    ["-" + fmt(dmax, 2) + "%", "0", "+" + fmt(dmax, 2) + "%"].forEach((s, i) => {
      text(ctx, s, pad.l + (i / 2) * iw, pad.t + ih + 18, { align: "center", size: 10 });
    });
    text(ctx, "相对畸变", pad.l + iw / 2, h - 10, { align: "center", size: 11 });
  }
}

/* ------------------------------------------------------------------ */
/* 处方表                                                               */
/* ------------------------------------------------------------------ */
function drawPrescription(data) {
  const rows = data.prescription;
  const head = ["#", "半径", "厚度", "材料", "半口径", "k", "z 位置", "说明"];
  let html = "<thead><tr>" + head.map((t) => `<th>${t}</th>`).join("") + "</tr></thead><tbody>";
  rows.forEach((r) => {
    const radius = r.radius === null ? "∞" : fmt(r.radius, 3);
    const tag = r.is_stop ? ' <span class="tag">光阑</span>' : "";
    const k = r.conic ? fmt(r.conic, 3) : "—";
    html += `<tr>
      <td>${r.index}</td>
      <td>${radius}</td>
      <td>${fmt(r.thickness, 3)}</td>
      <td class="${r.is_stop ? "stop" : ""}">${r.material}${tag}</td>
      <td>${fmt(r.semi_diameter, 3)}</td>
      <td>${k}</td>
      <td>${fmt(r.z, 3)}</td>
      <td style="text-align:left;color:#5d6b7a">${r.label || ""}</td>
    </tr>`;
  });
  const s = data.summary;
  const gap = data.image_plane_z - (rows[rows.length - 1].z || 0);
  const maxH = Math.max(0, ...(data.distortion || []).map((d) => Math.abs(d.real_height)));
  html += `<tr class="img-row">
      <td>像面</td>
      <td>—</td>
      <td>${fmt(gap, 3)}</td>
      <td>air</td>
      <td>${fmt(maxH * 1.05, 3)}</td>
      <td>—</td>
      <td>${fmt(data.image_plane_z, 3)}</td>
      <td style="text-align:left;color:#5d6b7a">探测器</td>
    </tr>`;
  html += "</tbody>";
  $("rx-table").innerHTML = html;
  $("rx-note").textContent =
    `半径与厚度单位 mm；EFL=${fmt(s.efl, 2)} mm，BFL=${fmt(s.bfl, 2)} mm，` +
    `入瞳半径 ${fmt(s.entrance_pupil_radius, 2)} mm，` +
    `像距 ${fmt(s.image_distance + s.image_shift, 2)} mm（含 ${fmt(s.image_shift, 3)} mm 调焦）。`;
}

/* ------------------------------------------------------------------ */
/* 指标条与控制                                                         */
/* ------------------------------------------------------------------ */
function renderMetrics(design) {
  const s = design.summary;
  const maxRms = Math.max(...design.spots.map((x) => x.rms)) * 1000;
  const items = [
    ["焦距 EFL", fmt(s.efl, 2), "mm"],
    ["相对孔径", "F/" + fmt(s.f_number, 1), ""],
    ["后截距 BFL", fmt(s.bfl, 2), "mm"],
    ["最大弥散", fmt(maxRms, 1), "µm"],
    ["视场", `${fmt(Math.min(...s.fields), 0)}–${fmt(Math.max(...s.fields), 0)}`, "°"],
    ["面数", String(s.surfaces), ""],
  ];
  $("metrics").innerHTML = items
    .map(([k, v, u]) => `<div class="metric"><div class="k">${k}</div><div class="v">${v}<small>${u}</small></div></div>`)
    .join("");
}

async function renderTab(tab) {
  const key = `${tab}|${state.name}|${state.defocus}`;
  if (state.cache.has(key)) {
    applyTab(tab, state.cache.get(key));
    return;
  }
  setStatus("计算中…", true);
  try {
    let data;
    if (tab === "layout") data = await api("design", { defocus: state.defocus });
    else if (tab === "spot") data = await api("spots", { defocus: state.defocus });
    else if (tab === "fan") data = await api("fans", { defocus: state.defocus });
    else if (tab === "mtf") data = await api("mtf", { defocus: state.defocus });
    else if (tab === "field") data = await api("field", {});
    else if (tab === "rx") data = await api("design", { defocus: state.defocus });
    state.cache.set(key, data);
    if (state.cache.size > 24) state.cache.clear();
    applyTab(tab, data);
    setStatus("就绪", false);
  } catch (err) {
    setStatus("计算失败：" + err.message, false);
  }
}

function applyTab(tab, data) {
  if (tab === "layout") { drawLayout(data); renderMetrics(data); }
  else if (tab === "spot") drawSpot(data);
  else if (tab === "fan") drawFan(data);
  else if (tab === "mtf") drawMtf(data);
  else if (tab === "field") drawField(data);
  else if (tab === "rx") drawPrescription(data);
}

/* ------------------------------------------------------------------ */
/* 初始化                                                               */
/* ------------------------------------------------------------------ */
async function boot() {
  const res = await fetch("/api/systems");
  const data = await res.json();
  state.systems = data.systems;
  $("system-select").innerHTML = data.systems
    .map((s) => `<option value="${s.name}">${s.title}</option>`)
    .join("");
  state.name = data.systems.find((s) => s.name === "cooke_triplet")?.name || data.systems[0].name;
  $("system-select").value = state.name;
  renderTab(state.tab);
}

document.querySelectorAll(".tab").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".view").forEach((v) => v.classList.remove("active"));
    btn.classList.add("active");
    state.tab = btn.dataset.tab;
    $("view-" + state.tab).classList.add("active");
    renderTab(state.tab);
  });
});

$("system-select").addEventListener("change", (e) => {
  state.name = e.target.value;
  state.cache.clear();
  renderTab(state.tab);
});

$("defocus").addEventListener("input", (e) => {
  state.defocus = parseFloat(e.target.value);
  $("defocus-value").textContent = state.defocus.toFixed(2);
  if (["spot", "fan", "mtf", "layout", "rx"].includes(state.tab)) renderTab(state.tab);
});

$("reset").addEventListener("click", () => {
  $("defocus").value = 0;
  state.defocus = 0;
  $("defocus-value").textContent = "0.00";
  state.cache.clear();
  renderTab(state.tab);
});

window.addEventListener("resize", () => {
  const key = `${state.tab}|${state.name}|${state.defocus}`;
  if (state.cache.has(key)) applyTab(state.tab, state.cache.get(key));
});

boot();
