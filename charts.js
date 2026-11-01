// Small canvas charts for the evidence panel. Each function redraws the whole canvas.
const PAD = { l: 46, r: 12, t: 36, b: 30 };
const DARK = { grid: "rgba(255,255,255,.07)", axis: "#8a8a8a", text: "#b3b3b3", limit: "#f87171", highlight: "#fbbf24" };
const LIGHT = { grid: "rgba(0,0,0,.08)", axis: "#6b7280", text: "#4b5563", limit: "#dc2626", highlight: "#d97706" };
const light = matchMedia("(prefers-color-scheme: light)");
let C = light.matches ? LIGHT : DARK;
const FONT = "11px system-ui, sans-serif";

export function setup(canvas, height) {
  C = light.matches ? LIGHT : DARK;
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth || 600;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(height * dpr);
  canvas.style.height = `${height}px`;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, height);
  ctx.font = FONT;
  return { ctx, w, h: height };
}

function ticks(lo, hi, n = 5) {
  const span = hi - lo || 1;
  const raw = span / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

const fmt = (v) => (Math.abs(v) >= 100 || Number.isInteger(v) ? String(Math.round(v)) : v.toFixed(Math.abs(v) < 1 ? 2 : 1));

// Plot frame with a y grid; returns scale functions.
function frame(ctx, w, h, x0, x1, y0, y1, o = {}) {
  const X = (v) => PAD.l + ((v - x0) / (x1 - x0 || 1)) * (w - PAD.l - PAD.r);
  const Y = (v) => h - PAD.b - ((v - y0) / (y1 - y0 || 1)) * (h - PAD.t - PAD.b);
  ctx.fillStyle = C.text;
  ctx.strokeStyle = C.grid;
  ctx.lineWidth = 1;
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const t of ticks(y0, y1, 4)) {
    ctx.beginPath();
    ctx.moveTo(PAD.l, Y(t) + 0.5);
    ctx.lineTo(w - PAD.r, Y(t) + 0.5);
    ctx.stroke();
    ctx.fillText(fmt(t), PAD.l - 6, Y(t));
  }
  if (!o.noXTicks) {
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const t of ticks(x0, x1, Math.max(2, Math.floor((w - PAD.l) / 70)))) ctx.fillText(fmt(t), X(t), h - PAD.b + 5);
  }
  if (o.xLabel) {
    ctx.textAlign = "right";
    ctx.textBaseline = "bottom";
    ctx.fillText(o.xLabel, w - PAD.r, h - 1);
  }
  if (o.yLabel) {
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillText(o.yLabel, 2, 14);
  }
  return { X, Y };
}

function limitLine(ctx, w, Y, value, label) {
  if (value == null) return;
  ctx.save();
  ctx.strokeStyle = C.limit;
  ctx.setLineDash([5, 4]);
  ctx.beginPath();
  ctx.moveTo(PAD.l, Y(value) + 0.5);
  ctx.lineTo(w - PAD.r, Y(value) + 0.5);
  ctx.stroke();
  ctx.fillStyle = C.limit;
  ctx.textAlign = "right";
  ctx.textBaseline = "bottom";
  ctx.fillText(label, w - PAD.r, Y(value) - 2);
  ctx.restore();
}

function legend(ctx, w, series) {
  let x = w - PAD.r;
  ctx.textAlign = "right";
  ctx.textBaseline = "top";
  for (const s of [...series].reverse()) {
    if (!s.name) continue;
    ctx.fillStyle = C.text;
    ctx.fillText(s.name, x, 0);
    x -= ctx.measureText(s.name).width + 4;
    ctx.fillStyle = s.color;
    ctx.fillRect(x - 9, 2, 9, 9);
    x -= 22;
  }
}

/** Grouped bars over categorical labels. */
export function bars(canvas, { labels, series, threshold, thresholdLabel, yMin = 0, yMax, yLabel, xLabel, height = 190 }) {
  const { ctx, w, h } = setup(canvas, height);
  const values = series.flatMap((s) => s.values.filter((v) => v != null));
  const top = yMax ?? Math.max(1, threshold ?? 0, ...values) * 1.1;
  const lo = Math.min(yMin, ...values);
  const { Y } = frame(ctx, w, h, 0, 1, lo, top, { noXTicks: true, yLabel, xLabel });
  const slot = (w - PAD.l - PAD.r) / labels.length;
  const bw = Math.max(2, (slot * 0.75) / series.length);
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  labels.forEach((label, i) => {
    const cx = PAD.l + slot * (i + 0.5);
    series.forEach((s, j) => {
      const v = s.values[i];
      if (v == null) return;
      const x = cx - (bw * series.length) / 2 + j * bw;
      ctx.fillStyle = typeof s.color === "function" ? s.color(v, i) : s.color;
      const y0 = Y(Math.max(0, lo)), y1 = Y(v);
      ctx.fillRect(x + 0.5, Math.min(y0, y1), bw - 1, Math.max(1, Math.abs(y1 - y0)));
    });
    ctx.fillStyle = C.text;
    ctx.fillText(label, cx, h - PAD.b + 5);
  });
  limitLine(ctx, w, Y, threshold, thresholdLabel);
  legend(ctx, w, series);
}

/** Values on a regular x grid (x0 + i*dx), drawn as columns or a line. */
export function series(canvas, { x0 = 0, dx = 1, series: list, threshold, thresholdLabel, yMax, xMax, yLabel, xLabel,
                                 bands = [], height = 170, columns = true }) {
  const { ctx, w, h } = setup(canvas, height);
  const n = Math.max(...list.map((s) => s.values.length), 1);
  const x1 = xMax ?? x0 + n * dx;
  const top = yMax ?? Math.max(1, threshold ?? 0, ...list.flatMap((s) => s.values)) * 1.1;
  const { X, Y } = frame(ctx, w, h, x0, x1, 0, top, { yLabel, xLabel });
  for (const b of bands) {
    ctx.fillStyle = b.color;
    ctx.fillRect(X(b.from), PAD.t, Math.max(2, X(b.to) - X(b.from)), h - PAD.t - PAD.b);
  }
  for (const s of list) {
    ctx.fillStyle = s.color;
    ctx.strokeStyle = s.color;
    if (columns && !s.line) {
      const cw = Math.max(1, X(x0 + dx) - X(x0) - (list.length > 1 ? 0 : 1));
      s.values.forEach((v, i) => {
        if (v > 0) ctx.fillRect(X(x0 + i * dx), Y(Math.min(v, top)), cw, Y(0) - Y(Math.min(v, top)));
      });
    } else {
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      s.values.forEach((v, i) => {
        const x = X(x0 + (i + 0.5) * dx), y = Y(v == null ? 0 : Math.min(v, top));
        i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
      });
      ctx.stroke();
    }
  }
  limitLine(ctx, w, Y, threshold, thresholdLabel);
  legend(ctx, w, list);
}

/** Matched hash pairs, query time across, reference time up. */
export function scatter(canvas, { xs, ys, xMax, yMax, color, highlight, xLabel, yLabel, height = 260 }) {
  const { ctx, w, h } = setup(canvas, height);
  const { X, Y } = frame(ctx, w, h, 0, xMax, 0, yMax, { xLabel, yLabel });
  ctx.fillStyle = color;
  for (let i = 0; i < xs.length; i++) {
    ctx.globalAlpha = highlight && highlight(xs[i], ys[i]) ? 1 : 0.35;
    ctx.fillStyle = highlight && highlight(xs[i], ys[i]) ? C.highlight : color;
    ctx.fillRect(X(xs[i]) - 1, Y(ys[i]) - 1, 2.5, 2.5);
  }
  ctx.globalAlpha = 1;
}
