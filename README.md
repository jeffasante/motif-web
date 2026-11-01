# Motif (web)

Compare two audio files in your browser. Motif tells you whether they share audio
(a sample, a re-release, the same recording) or are musically similar (a cover:
same melody and chords, possibly in another key or tempo).

Everything runs in the page. Your files are not uploaded. The page follows your
system's light or dark setting.

**Use it:** https://jeffasante.github.io/motif-web/

Choose a query file and a reference file, then press **Compare**. You get:

- two track cards (duration, key, tempo) joined by a connection badge, and a plain
  answer: a shared sample, a musical resemblance only (not a sample), or nothing found;
- linked timelines, each with a waveform, a spectrogram and a chroma lane; matched
  passages are highlighted and connected across the two tracks;
- **Matched passages**: a zoomed waveform of where each passage sits in both tracks
  (click it to play from there); play a passage in both tracks together or one after the other,
  trim or split it, and save either side as a WAV clip (cut from the full-quality
  decode, all channels kept, 16-bit at the browser's audio rate);
- **Why this result**: graphs of the engine's own evidence against its thresholds,
  including how short, pitch-shifted, heavily processed samples and dense mixes
  fall short, plus plain-language hints;
- the full JSON report.

## How it is built

The analysis engine is written in C (`src/`) and compiled to WebAssembly with
[Emscripten](https://emscripten.org/). It uses only the C standard and math libraries.

```
audio file
  -> browser decodes it (Web Audio), mixes to mono, resamples to 22050 Hz
  -> motif_engine.wasm:
       spectral peaks -> landmark hashes -> offset histogram  => shared audio
       chroma (12 pitch classes) -> 8 s windows x 12 keys x 5 tempos => musical similarity
       key and tempo estimates
  -> JSON report, with a "diagnostics" section: matched hash pairs, the offset
     histogram, per-second matches and peak density, a pitch/speed sweep of
     the fingerprints, and chroma similarity per transposition and tempo
```

The graphs in the page are drawn from those diagnostics. The pitch sweep is
informational only: the detector itself compares at the original pitch.

| File | Role |
|---|---|
| `src/dsp.c` | Real FFT, Hann window, memory arena |
| `src/fingerprint.c` | Spectral peak picking and landmark hashes |
| `src/musical.c` | Chroma, key and tempo estimation |
| `src/matching.c` | Shared-audio segments and musical-similarity search |
| `src/api.c` | Input checks, assessment, JSON output; exports `motif_compare` |
| `motif.js` | Loads the engine, copies audio into WebAssembly memory, returns the report |
| `index.html` | The page: decoding, SHA-256 of each file, display |
| `viz.js` | Display-only spectrogram, chroma and pitch profile (not used for the result) |
| `charts.js` | Small dependency-free canvas charts for the evidence graphs |
| `motif_engine.js` / `.wasm` | Build output (committed so GitHub Pages can serve it) |

### Rebuild

```sh
sh build.sh            # needs emcc on PATH
python3 -m http.server 8765 --bind 127.0.0.1
# open http://127.0.0.1:8765/
```

The page must be served over HTTP; browsers do not load WebAssembly modules from
`file://`.

## Limitations

- Which formats open depends on your browser's decoders (WAV, MP3 and FLAC work
  widely; AAC/M4A does not decode in every browser).
- Settings are fixed at Motif's defaults.
- Files up to 30 minutes. The whole file is held in memory.
- Results are evidence, not proof. `no_evidence` does not mean the recordings are
  unrelated; short, pitch-shifted or heavily processed samples can fall below the
  thresholds.
