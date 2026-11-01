// Display-only analysis for the page: waveform, spectrogram and chroma images computed in
// JavaScript from the decoded samples. The comparison itself is done by the WebAssembly engine.
export const PITCHES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"];

const N = 4096;
const BITS = 12;
const rev = new Uint32Array(N);
for (let i = 0; i < N; i++) {
  let r = 0;
  for (let b = 0; b < BITS; b++) r |= ((i >> b) & 1) << (BITS - 1 - b);
  rev[i] = r;
}
const cosT = new Float64Array(N / 2);
const sinT = new Float64Array(N / 2);
for (let i = 0; i < N / 2; i++) {
  cosT[i] = Math.cos((2 * Math.PI * i) / N);
  sinT[i] = -Math.sin((2 * Math.PI * i) / N);
}
const hann = new Float64Array(N);
for (let i = 0; i < N; i++) hann[i] = 0.5 - 0.5 * Math.cos((2 * Math.PI * i) / N);

// power[0..N/2] = |FFT(x[start:start+N] * hann)|^2, zero-padded outside the signal.
function powerSpectrum(x, start, re, im, power) {
  for (let i = 0; i < N; i++) {
    const j = start + i;
    re[rev[i]] = j >= 0 && j < x.length ? x[j] * hann[i] : 0;
    im[rev[i]] = 0;
  }
  for (let size = 2; size <= N; size <<= 1) {
    const half = size >> 1;
    const step = N / size;
    for (let s = 0; s < N; s += size) {
      for (let k = 0; k < half; k++) {
        const tr = cosT[k * step];
        const ti = sinT[k * step];
        const a = s + k;
        const b = a + half;
        const xr = re[b] * tr - im[b] * ti;
        const xi = re[b] * ti + im[b] * tr;
        re[b] = re[a] - xr;
        im[b] = im[a] - xi;
        re[a] += xr;
        im[a] += xi;
      }
    }
  }
  for (let k = 0; k <= N / 2; k++) power[k] = re[k] * re[k] + im[k] * im[k];
}

/** Summarise a mono signal into `cols` time columns for drawing. */
export function analyse(x, rate, cols = 1000, rows = 160) {
  const df = rate / N;
  const fMin = 50;
  const fMax = rate / 2;
  const lo = new Int32Array(rows);
  const hi = new Int32Array(rows);
  for (let r = 0; r < rows; r++) {
    const f0 = fMin * Math.pow(fMax / fMin, r / rows);
    const f1 = fMin * Math.pow(fMax / fMin, (r + 1) / rows);
    lo[r] = Math.min(N / 2, Math.round(f0 / df));
    hi[r] = Math.min(N / 2, Math.max(lo[r], Math.round(f1 / df) - 1));
  }
  // Nearest pitch class of each bin from 65 to 2100 Hz, as the engine's chroma does.
  const pc = new Int8Array(N / 2 + 1).fill(-1);
  for (let k = 1; k <= N / 2; k++) {
    const f = k * df;
    if (f >= 65 && f <= 2100) pc[k] = (((Math.round(69 + 12 * Math.log2(f / 440)) % 12) + 12) % 12);
  }

  const spec = new Float32Array(rows * cols);
  const chroma = new Float32Array(12 * cols);
  const energy = new Float64Array(cols);
  const waveMin = new Float32Array(cols);
  const waveMax = new Float32Array(cols);
  const re = new Float64Array(N);
  const im = new Float64Array(N);
  const power = new Float64Array(N / 2 + 1);
  let peak = 0;

  for (let c = 0; c < cols; c++) {
    const a = Math.floor((c * x.length) / cols);
    const b = Math.max(a + 1, Math.floor(((c + 1) * x.length) / cols));
    let mn = 0;
    let mx = 0;
    for (let i = a; i < b && i < x.length; i++) {
      const v = x[i];
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
    waveMin[c] = mn;
    waveMax[c] = mx;
    peak = Math.max(peak, mx, -mn);

    powerSpectrum(x, Math.round(((c + 0.5) * x.length) / cols) - N / 2, re, im, power);
    for (let r = 0; r < rows; r++) {
      let m = 0;
      for (let k = lo[r]; k <= hi[r]; k++) if (power[k] > m) m = power[k];
      spec[r * cols + c] = 10 * Math.log10(m + 1e-12);
    }
    for (let k = 1; k <= N / 2; k++) {
      if (pc[k] >= 0) {
        chroma[pc[k] * cols + c] += power[k];
        energy[c] += power[k];
      }
    }
  }

  // Per column: pitch classes relative to the strongest one, faded where the music is quiet.
  let eMax = 0;
  for (const e of energy) if (e > eMax) eMax = e;
  const chromaN = new Float32Array(12 * cols);
  for (let c = 0; c < cols; c++) {
    if (!(energy[c] > 0) || !(eMax > 0)) continue;
    const loud = Math.min(1, Math.max(0, (10 * Math.log10(energy[c] / eMax) + 40) / 40));
    let colMax = 0;
    for (let k = 0; k < 12; k++) colMax = Math.max(colMax, chroma[k * cols + c]);
    for (let k = 0; k < 12; k++) {
      chromaN[k * cols + c] = (Math.log1p((100 * chroma[k * cols + c]) / colMax) / Math.log1p(100)) * loud;
    }
  }

  return { cols, rows, spec, chroma: chromaN, waveMin, waveMax, peak, duration: x.length / rate };
}

// An inferno-like colour ramp.
const STOPS = [
  [0, 0, 4], [31, 12, 72], [85, 15, 109], [136, 34, 106], [186, 54, 85],
  [227, 89, 51], [249, 140, 10], [249, 201, 50], [252, 255, 164],
];
const LUT = new Uint8ClampedArray(256 * 3);
for (let i = 0; i < 256; i++) {
  const t = (i / 255) * (STOPS.length - 1);
  const j = Math.floor(t);
  const k = Math.min(j + 1, STOPS.length - 1);
  for (let ch = 0; ch < 3; ch++) LUT[i * 3 + ch] = STOPS[j][ch] + (STOPS[k][ch] - STOPS[j][ch]) * (t - j);
}

function image(width, height, value) {
  const canvas = document.createElement("canvas");
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext("2d");
  const img = ctx.createImageData(width, height);
  const px = img.data;
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const t = Math.min(1, Math.max(0, value(x, height - 1 - y)));
      const i = Math.round(t * 255) * 3;
      const o = (y * width + x) * 4;
      px[o] = LUT[i];
      px[o + 1] = LUT[i + 1];
      px[o + 2] = LUT[i + 2];
      px[o + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
  return canvas;
}

/** Log-frequency spectrogram, low pitches at the bottom, 80 dB range. */
export function spectrogramImage(v) {
  let max = -Infinity;
  for (const d of v.spec) if (d > max) max = d;
  return image(v.cols, v.rows, (x, r) => (v.spec[r * v.cols + x] - (max - 80)) / 80);
}

/** Chromagram, C at the bottom. */
export function chromaImage(v) {
  return image(v.cols, 12, (x, k) => v.chroma[k * v.cols + x]);
}

/** Mean pitch-class profile, scaled so the strongest class is 1. */
export function pitchProfile(v) {
  const p = new Array(12).fill(0);
  for (let k = 0; k < 12; k++) for (let c = 0; c < v.cols; c++) p[k] += v.chroma[k * v.cols + c];
  const max = Math.max(...p);
  return p.map((x) => (max > 0 ? x / max : 0));
}
