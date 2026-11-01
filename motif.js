// JavaScript wrapper around the Motif WebAssembly engine. Works in browsers and Node.
import createMotifEngine from "./motif_engine.js";

export const SAMPLE_RATE = 22050;

let enginePromise = null;

function engine() {
  enginePromise ??= createMotifEngine();
  return enginePromise;
}

function copyIn(m, samples) {
  const ptr = m._motif_malloc(samples.length * 4);
  if (!ptr) throw new Error("WebAssembly engine is out of memory");
  m.HEAPF32.set(samples, ptr >>> 2);
  return ptr;
}

function cString(m, value) {
  return value == null ? 0 : m.stringToNewUTF8(value);
}

/**
 * Compare two mono Float32Array signals at 22050 Hz.
 * `query` and `reference` are {samples, name, sha256} where sha256 is the hex digest of the
 * original file bytes, or null when unknown. Returns the parsed `motif compare` report.
 */
export async function compare(query, reference, decoder) {
  const m = await engine();
  const owned = [];
  try {
    const q = copyIn(m, query.samples);
    owned.push(q);
    const r = copyIn(m, reference.samples);
    owned.push(r);
    const strings = [query.name, query.sha256, reference.name, reference.sha256, decoder].map((s) => cString(m, s));
    owned.push(...strings.filter(Boolean));
    const out = m._motif_compare(q, query.samples.length, strings[0], strings[1],
                                 r, reference.samples.length, strings[2], strings[3], strings[4]);
    if (!out) throw new Error("WebAssembly engine ran out of memory");
    owned.push(out);
    return JSON.parse(m.UTF8ToString(out));
  } finally {
    for (const p of owned) m._motif_free(p);
  }
}
