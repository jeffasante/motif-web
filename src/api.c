/* motif_compare: the exported entry point. Validates two mono PCM buffers, extracts
 * features, matches them and returns the `motif compare` JSON report. */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "motif.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

static const double SPEED_RATIOS[] = {1.0};
static const double TEMPO_RATIOS[] = {0.8, 0.9, 1.0, 1.1, 1.25};

static const FingerprintConfig FP = {2048, 256, 100.0, 5000.0, -60.0, 10, 10, 50, 10, 1, 64, 100};
static const MusicalConfig MC = {4096, 1024, 65.0, 2100.0, 0.25, 0.5, 60.0, 200.0, 0.1};
static const MatchingConfig MATCH = {40, 2, 5, SPEED_RATIOS, 1, 8.0, 2.0, 0.8, 3, TEMPO_RATIOS, 5};

#define MIN_DURATION_SECONDS 1.0
#define MAX_DURATION_SECONDS 1800.0
#define SILENCE_DBFS -90.0

static const char *const ASSESS_SHARED =
    "The same recorded audio occurs in both files. This is consistent with a sample, a re-release or an edit "
    "of one recording. It does not show which came first or whether the use was licensed.";
static const char *const ASSESS_MUSICAL =
    "No shared recorded audio was found, but the pitch-class (chroma) sequences agree for a sustained stretch. "
    "This is a lead for a listener, not a finding: it is consistent with a cover, an interpolation or a "
    "re-recording, and also with two pieces that share a common chord progression.";
static const char *const ASSESS_NONE =
    "Neither shared audio nor sustained musical similarity reached the configured thresholds. This does not "
    "show that the recordings are unrelated: short, pitch-shifted or heavily processed samples and very dense "
    "mixes can fall below them.";
static const char *const ASSESS_IDENTICAL = "The two inputs are the same file content (equal SHA-256).";
static const char *const ASSESS_NOT_COMPARED = "At least one input could not be analysed; see its status and error.";

/* ---- JSON output ------------------------------------------------------- */

typedef struct {
    char *data;
    size_t len, cap;
    int failed;
} Buf;

static void put(Buf *b, const char *s, size_t n) {
    if (b->failed) return;
    if (b->len + n + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        while (b->len + n + 1 > cap) cap *= 2;
        char *d = realloc(b->data, cap);
        if (!d) { b->failed = 1; return; }
        b->data = d;
        b->cap = cap;
    }
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = 0;
}

static void puts_(Buf *b, const char *s) { put(b, s, strlen(s)); }

static void printf_(Buf *b, const char *fmt, ...) {
    char tmp[512];
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    if (n < 0 || (size_t)n >= sizeof tmp) { b->failed = 1; return; }
    put(b, tmp, (size_t)n);
}

static void str(Buf *b, const char *s) {
    if (!s) { puts_(b, "null"); return; }
    puts_(b, "\"");
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (*p == '"' || *p == '\\') printf_(b, "\\%c", *p);
        else if (*p < 0x20) printf_(b, "\\u%04x", *p);
        else put(b, (const char *)p, 1);
    }
    puts_(b, "\"");
}

/* A float rounded to `digits` decimals, printed as Python prints the rounded value. */
static void num(Buf *b, double v, int digits) {
    char tmp[64];
    snprintf(tmp, sizeof tmp, "%.*f", digits, v);
    if (strcmp(tmp, "-0") == 0 || strncmp(tmp, "-0.", 3) == 0) {
        int zero = 1;
        for (const char *p = tmp + 1; *p; p++) if (*p != '0' && *p != '.') zero = 0;
        if (zero) memmove(tmp, tmp + 1, strlen(tmp));
    }
    char *dot = strchr(tmp, '.');
    if (dot) {
        char *end = tmp + strlen(tmp) - 1;
        while (end > dot + 1 && *end == '0') *end-- = 0;
    } else {
        strcat(tmp, ".0");
    }
    puts_(b, tmp);
}

static void key(Buf *b, const char *k) {
    str(b, k);
    puts_(b, ": ");
}

/* ---- per-input analysis ------------------------------------------------ */

typedef struct {
    const char *name, *sha256;
    const char *status;
    char error[256];
    double duration;
    int completed;
    Peaks peaks;
    Hashes hashes, sorted;
    Chroma chroma;
    KeyEstimate key;
    TempoEstimate tempo;
} Input;

/* The checks of motif.audio.validation that apply to host-decoded mono PCM. */
static void validate(Input *in, const float *x, size_t n) {
    in->duration = (double)n / SAMPLE_RATE;
    size_t non_finite = 0;
    float peak = 0;
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(x[i])) { non_finite++; continue; }
        float v = fabsf(x[i]);
        if (v > peak) peak = v;
    }
    char failed[3][128];
    int n_failed = 0;
    const char *status = NULL;
    if (n == 0) {
        snprintf(failed[n_failed++], 128, "duration: decoder produced no audio samples");
        status = "insufficient_audio";
    } else if (in->duration < MIN_DURATION_SECONDS) {
        snprintf(failed[n_failed++], 128, "duration: %.3f s is shorter than the minimum %g s", in->duration,
                 MIN_DURATION_SECONDS);
        status = "insufficient_audio";
    } else if (in->duration > MAX_DURATION_SECONDS) {
        snprintf(failed[n_failed++], 128, "duration: longer than the maximum %g s", MAX_DURATION_SECONDS);
        status = "failed";
    }
    if (non_finite) {
        snprintf(failed[n_failed++], 128, "finite_samples: %zu NaN or infinite sample values", non_finite);
        if (!status) status = "failed";
    }
    if (n == 0) {
        snprintf(failed[n_failed++], 128, "signal_present: no samples to measure");
    } else if (non_finite == n) {
        snprintf(failed[n_failed++], 128, "signal_present: no finite samples to measure");
        if (!status) status = "insufficient_audio";
    } else if (peak <= 0 || 20 * log10(peak) < SILENCE_DBFS) {
        snprintf(failed[n_failed++], 128, "signal_present: peak level below %g dBFS; the audio is silent",
                 SILENCE_DBFS);
        if (!status) status = "insufficient_audio";
    }
    in->error[0] = 0;
    for (int i = 0; i < n_failed; i++) {
        if (i) strncat(in->error, "; ", sizeof in->error - strlen(in->error) - 1);
        strncat(in->error, failed[i], sizeof in->error - strlen(in->error) - 1);
    }
    in->status = status ? status : "completed";
    in->completed = status == NULL;
}

static int analyse(Arena *a, Input *in, const float *x, size_t n) {
    validate(in, x, n);
    if (!in->completed) return 0;
    if (find_peaks(a, x, n, &FP, &in->peaks) || hash_peaks(a, &in->peaks, &FP, 1.0, &in->hashes) ||
        sort_hashes(a, &in->hashes, &in->sorted) || chroma(a, x, n, &MC, &in->chroma) ||
        estimate_tempo(a, x, n, &MC, &in->tempo))
        return -1;
    estimate_key(&in->chroma, &MC, &in->key);
    return 0;
}

static void key_label(char *out, size_t size, int k) {
    snprintf(out, size, "%s %s", PITCH_CLASSES[k % 12], k < 12 ? "major" : "minor");
}

static void write_key(Buf *b, const KeyEstimate *k) {
    if (!k->has_scores) {
        puts_(b, "{\"key\": null, \"correlation\": null, \"runner_up\": null, \"runner_up_correlation\": null, "
                 "\"reason\": \"chroma has no pitch-class contrast\"}");
        return;
    }
    char best[16], second[16], reason[160];
    key_label(best, sizeof best, k->best);
    key_label(second, sizeof second, k->second);
    puts_(b, "{");
    key(b, "key");
    str(b, k->known ? best : NULL);
    puts_(b, ", ");
    key(b, "correlation");
    num(b, k->best_score, 3);
    puts_(b, ", ");
    key(b, "runner_up");
    str(b, second);
    puts_(b, ", ");
    key(b, "runner_up_correlation");
    num(b, k->second_score, 3);
    puts_(b, ", ");
    key(b, "reason");
    snprintf(reason, sizeof reason, "best key profile (%s) correlates %.3f, below musical.key_min_correlation=%g",
             best, k->best_score, MC.key_min_correlation);
    str(b, k->known ? NULL : reason);
    puts_(b, "}");
}

static void write_tempo(Buf *b, const TempoEstimate *t) {
    if (t->state == TEMPO_NO_ONSETS) {
        puts_(b, "{\"bpm\": null, \"strength\": null, \"alternatives\": [], "
                 "\"reason\": \"audio too short or without onsets for the tempo range\"}");
        return;
    }
    puts_(b, "{");
    key(b, "bpm");
    if (t->state == TEMPO_KNOWN) num(b, t->bpm, 1);
    else puts_(b, "null");
    puts_(b, ", ");
    key(b, "strength");
    num(b, t->strength, 3);
    puts_(b, ", ");
    key(b, "alternatives");
    puts_(b, "[");
    for (int i = 0; i < t->n_alternatives; i++) {
        if (i) puts_(b, ", ");
        num(b, t->alternatives[i], 1);
    }
    puts_(b, "], ");
    key(b, "reason");
    if (t->state == TEMPO_WEAK) {
        char reason[128];
        snprintf(reason, sizeof reason, "autocorrelation strength %.3f is below musical.tempo_min_strength=%g",
                 t->strength, MC.tempo_min_strength);
        str(b, reason);
    } else {
        puts_(b, "null");
    }
    puts_(b, "}");
}

static void write_summary(Buf *b, const Input *in) {
    puts_(b, "{");
    key(b, "input");
    str(b, in->name);
    puts_(b, ", ");
    key(b, "status");
    str(b, in->status);
    puts_(b, ", ");
    key(b, "sha256");
    str(b, in->sha256);
    puts_(b, ", ");
    key(b, "duration_seconds");
    num(b, in->duration, 3);
    puts_(b, ", ");
    key(b, "key");
    if (in->completed) write_key(b, &in->key);
    else puts_(b, "null");
    puts_(b, ", ");
    key(b, "tempo");
    if (in->completed) write_tempo(b, &in->tempo);
    else puts_(b, "null");
    puts_(b, ", ");
    key(b, "error");
    if (in->completed) {
        puts_(b, "null");
    } else {
        puts_(b, "{\"type\": \"ValidationFailed\", \"message\": ");
        str(b, in->error);
        puts_(b, "}");
    }
    puts_(b, "}");
}

static void write_ratios(Buf *b, const double *r, int n) {
    puts_(b, "[");
    for (int i = 0; i < n; i++) {
        if (i) puts_(b, ", ");
        num(b, r[i], 6);
    }
    puts_(b, "]");
}

static void write_envelope(Buf *b, const char *decoder) {
    puts_(b, "{\"schema_version\": \"1.0\", \"motif_version\": \"" MOTIF_VERSION "\", \"command\": \"compare\", ");
    puts_(b, "\"analysis\": {\"engine\": \"motif C/WebAssembly port\", ");
    key(b, "decoder");
    str(b, decoder);
    printf_(b, ", \"features\": {\"sample_rate\": %d, \"downmix_policy\": \"mean_of_channels\", "
               "\"fingerprint_version\": \"fingerprint-1\", \"musical_version\": \"musical-1\", ", SAMPLE_RATE);
    printf_(b, "\"fingerprint\": {\"n_fft\": %d, \"hop_length\": %d, \"f_min\": ", FP.n_fft, FP.hop_length);
    num(b, FP.f_min, 6);
    puts_(b, ", \"f_max\": ");
    num(b, FP.f_max, 6);
    puts_(b, ", \"min_peak_dbfs\": ");
    num(b, FP.min_peak_dbfs, 6);
    printf_(b, ", \"neighborhood_frames\": %d, \"neighborhood_bins\": %d, \"peaks_per_second\": %d, "
               "\"fan_out\": %d, \"target_min_frames\": %d, \"target_max_frames\": %d, \"target_max_bins\": %d}, ",
            FP.neighborhood_frames, FP.neighborhood_bins, FP.peaks_per_second, FP.fan_out, FP.target_min_frames,
            FP.target_max_frames, FP.target_max_bins);
    printf_(b, "\"musical\": {\"chroma_n_fft\": %d, \"chroma_hop_length\": %d, \"chroma_f_min\": ", MC.chroma_n_fft,
            MC.chroma_hop_length);
    num(b, MC.chroma_f_min, 6);
    puts_(b, ", \"chroma_f_max\": ");
    num(b, MC.chroma_f_max, 6);
    puts_(b, ", \"chroma_frame_seconds\": ");
    num(b, MC.chroma_frame_seconds, 6);
    puts_(b, ", \"key_min_correlation\": ");
    num(b, MC.key_min_correlation, 6);
    puts_(b, ", \"tempo_min_bpm\": ");
    num(b, MC.tempo_min_bpm, 6);
    puts_(b, ", \"tempo_max_bpm\": ");
    num(b, MC.tempo_max_bpm, 6);
    puts_(b, ", \"tempo_min_strength\": ");
    num(b, MC.tempo_min_strength, 6);
    printf_(b, "}}, \"thresholds\": {\"min_aligned_hashes\": %d, \"speed_ratios\": ", MATCH.min_aligned_hashes);
    write_ratios(b, MATCH.speed_ratios, MATCH.n_speed_ratios);
    puts_(b, ", \"melodic_min_similarity\": ");
    num(b, MATCH.melodic_min_similarity, 6);
    printf_(b, ", \"melodic_min_run_windows\": %d, \"melodic_window_seconds\": ", MATCH.melodic_min_run_windows);
    num(b, MATCH.melodic_window_seconds, 6);
    puts_(b, ", \"melodic_tempo_ratios\": ");
    write_ratios(b, MATCH.melodic_tempo_ratios, MATCH.n_melodic_tempo_ratios);
    puts_(b, "}}, ");
}

static void write_shared(Buf *b, const SharedAudio *s) {
    puts_(b, "{\"detected\": ");
    puts_(b, s->detected ? "true" : "false");
    puts_(b, ", \"speed_ratio\": ");
    num(b, s->speed_ratio, 6);
    printf_(b, ", \"matching_hashes\": %zu, \"query_hashes\": %zu, \"segments\": [", s->matching_hashes,
            s->query_hashes);
    for (int i = 0; i < s->n_segments; i++) {
        const Segment *g = &s->segments[i];
        if (i) puts_(b, ", ");
        puts_(b, "{\"query_start_seconds\": ");
        num(b, g->query_start, 2);
        puts_(b, ", \"query_end_seconds\": ");
        num(b, g->query_end, 2);
        puts_(b, ", \"reference_start_seconds\": ");
        num(b, g->reference_start, 2);
        puts_(b, ", \"reference_end_seconds\": ");
        num(b, g->reference_end, 2);
        printf_(b, ", \"aligned_hashes\": %d, \"offset_seconds\": ", g->aligned);
        num(b, g->offset, 3);
        puts_(b, "}");
    }
    puts_(b, "]}");
}

static void write_musical(Buf *b, const Musical *m) {
    if (!m->detected) {
        char reason[128];
        if (m->short_reference)
            snprintf(reason, sizeof reason, "reference is shorter than one %g s window", MATCH.melodic_window_seconds);
        else
            snprintf(reason, sizeof reason, "no run of %d windows reached similarity %g", MATCH.melodic_min_run_windows,
                     MATCH.melodic_min_similarity);
        puts_(b, "{\"detected\": false, \"best\": null, \"reason\": ");
        str(b, reason);
        puts_(b, "}");
        return;
    }
    printf_(b, "{\"detected\": true, \"best\": {\"windows\": %d, \"mean_similarity\": ", m->windows);
    num(b, m->mean_similarity, 3);
    printf_(b, ", \"transposition_semitones\": %d, \"tempo_ratio\": ", m->transposition);
    num(b, m->tempo_ratio, 6);
    puts_(b, ", \"query_start_seconds\": ");
    num(b, m->query_start, 2);
    puts_(b, ", \"query_end_seconds\": ");
    num(b, m->query_end, 2);
    puts_(b, ", \"reference_start_seconds\": ");
    num(b, m->reference_start, 2);
    puts_(b, ", \"reference_end_seconds\": ");
    num(b, m->reference_end, 2);
    puts_(b, "}, \"reason\": null}");
}

/* ---- diagnostics: the data behind the page's graphs -------------------- */

#define MAX_POINTS 4000
#define MAX_OFFSET_BINS 1500
#define SWEEP_SEMITONES 6

static void write_ints(Buf *b, const int *v, size_t n) {
    puts_(b, "[");
    for (size_t i = 0; i < n; i++) printf_(b, i ? ", %d" : "%d", v[i]);
    puts_(b, "]");
}

static void write_score(Buf *b, double v) {
    if (isfinite(v)) num(b, v, 3);
    else puts_(b, "null");
}

static void write_density(Buf *b, const Input *in) {
    printf_(b, "{\"hashes\": %zu, \"peaks_found_per_second\": ", in->hashes.n);
    write_ints(b, in->peaks.found_per_bucket, (size_t)in->peaks.n_buckets);
    puts_(b, "}");
}

static int write_diagnostics(Arena *a, Buf *b, const Input *q, const Input *r, const Musical *m) {
    double fs = (double)FP.hop_length / SAMPLE_RATE;
    int tol = MATCH.offset_tolerance_frames, fpb = q->peaks.frames_per_bucket;
    Alignment al;
    if (align_hashes(a, &q->hashes, &r->sorted, tol, &al)) return -1;

    puts_(b, "{\"frame_seconds\": ");
    num(b, fs, 6);
    puts_(b, ", \"bucket_seconds\": ");
    num(b, fpb * fs, 6);
    printf_(b, ", \"min_aligned_hashes\": %d, \"peaks_per_second_cap\": %d, \"melodic_min_similarity\": ",
            MATCH.min_aligned_hashes, FP.peaks_per_second);
    num(b, MATCH.melodic_min_similarity, 3);
    puts_(b, ", \"query\": ");
    write_density(b, q);
    puts_(b, ", \"reference\": ");
    write_density(b, r);

    printf_(b, ", \"alignment\": {\"query_hashes_matched\": %zu, \"pairs\": %zu, \"best_aligned\": %zu, "
               "\"best_offset_seconds\": ", al.matching, al.n_pairs, al.best_count);
    num(b, (double)al.best_offset * fs, 3);

    /* Matched hash pairs as (query time, reference time); a diagonal line is shared audio. */
    size_t stride = al.n_pairs > MAX_POINTS ? (al.n_pairs + MAX_POINTS - 1) / MAX_POINTS : 1;
    printf_(b, ", \"pairs_stride\": %zu, \"pair_query_seconds\": [", stride);
    for (size_t i = 0; i < al.n_pairs; i += stride) {
        if (i) puts_(b, ", ");
        num(b, al.q[i] * fs, 2);
    }
    puts_(b, "], \"pair_reference_seconds\": [");
    for (size_t i = 0; i < al.n_pairs; i += stride) {
        if (i) puts_(b, ", ");
        num(b, (double)(al.q[i] + al.offset[i]) * fs, 2);
    }
    puts_(b, "]");

    /* Offset histogram: pairs within +-tol frames of each offset, as the detector counts them,
     * then the maximum over groups of offsets so it fits in MAX_OFFSET_BINS. */
    int64_t lo = -(int64_t)ceil(q->duration / fs) - tol, hi = (int64_t)ceil(r->duration / fs) + tol;
    size_t span = (size_t)(hi - lo + 1);
    int *per = arena_alloc(a, span, sizeof(int)), *win = arena_alloc(a, span, sizeof(int));
    if (!per || !win) return -1;
    for (size_t i = 0; i < al.n_pairs; i++)
        if (al.offset[i] >= lo && al.offset[i] <= hi) per[al.offset[i] - lo]++;
    int sum = 0;
    for (size_t i = 0; i < span + (size_t)tol; i++) {
        if (i < span) sum += per[i];
        if (i >= (size_t)(2 * tol + 1)) sum -= per[i - (size_t)(2 * tol + 1)];
        if (i >= (size_t)tol) win[i - (size_t)tol] = sum;
    }
    size_t group = (span + MAX_OFFSET_BINS - 1) / MAX_OFFSET_BINS, n_bins = (span + group - 1) / group;
    for (size_t g = 0; g < n_bins; g++) {
        int best = 0;
        for (size_t i = g * group; i < (g + 1) * group && i < span; i++)
            if (win[i] > best) best = win[i];
        win[g] = best;
    }
    puts_(b, ", \"offset_histogram\": {\"start_seconds\": ");
    num(b, (double)lo * fs, 3);
    puts_(b, ", \"bin_seconds\": ");
    num(b, (double)group * fs, 6);
    puts_(b, ", \"counts\": ");
    write_ints(b, win, n_bins);
    puts_(b, "}");

    /* Per query second: hashes, hashes with any reference match, and pairs at the best offset. */
    size_t n_sec = (size_t)q->peaks.n_buckets;
    int *hashes = arena_alloc(a, n_sec + 1, sizeof(int)), *matched = arena_alloc(a, n_sec + 1, sizeof(int));
    int *aligned = arena_alloc(a, n_sec + 1, sizeof(int));
    if (!hashes || !matched || !aligned) return -1;
    for (size_t i = 0; i < q->hashes.n; i++) {
        size_t s = (size_t)(q->hashes.items[i].frame / fpb);
        if (s >= n_sec) continue;
        hashes[s]++;
        matched[s] += al.matched[i];
    }
    for (size_t i = 0; i < al.n_pairs; i++) {
        int64_t d = al.offset[i] - al.best_offset;
        size_t s = (size_t)(al.q[i] / fpb);
        if (s < n_sec && (d < 0 ? -d : d) <= tol) aligned[s]++;
    }
    puts_(b, ", \"query_hashes_per_second\": ");
    write_ints(b, hashes, n_sec);
    puts_(b, ", \"matched_per_second\": ");
    write_ints(b, matched, n_sec);
    puts_(b, ", \"aligned_per_second\": ");
    write_ints(b, aligned, n_sec);
    puts_(b, "}");

    /* Best aligned count if the query were pitch-shifted (time kept) or sped up (pitch and time). */
    int sweep[2][2 * SWEEP_SEMITONES + 1];
    for (int mode = 0; mode < 2; mode++) {
        for (int s = -SWEEP_SEMITONES; s <= SWEEP_SEMITONES; s++) {
            size_t mark = arena_mark(a);
            double ratio = pow(2.0, s / 12.0);
            Hashes h;
            Alignment x;
            if (hash_peaks_scaled(a, &q->peaks, &FP, ratio, mode ? ratio : 1.0, &h) ||
                align_hashes(a, &h, &r->sorted, tol, &x))
                return -1;
            sweep[mode][s + SWEEP_SEMITONES] = (int)x.best_count;
            arena_release(a, mark);
        }
    }
    puts_(b, ", \"pitch_sweep\": {\"semitones\": [");
    for (int s = -SWEEP_SEMITONES; s <= SWEEP_SEMITONES; s++) printf_(b, s > -SWEEP_SEMITONES ? ", %d" : "%d", s);
    puts_(b, "], \"pitch_only\": ");
    write_ints(b, sweep[0], 2 * SWEEP_SEMITONES + 1);
    puts_(b, ", \"speed\": ");
    write_ints(b, sweep[1], 2 * SWEEP_SEMITONES + 1);
    puts_(b, "}");

    puts_(b, ", \"chroma_transposition\": {\"semitones\": [");
    for (int s = -5; s <= 6; s++) printf_(b, s > -5 ? ", %d" : "%d", s);
    puts_(b, "], \"best_window_similarity\": [");
    for (int s = -5; s <= 6; s++) {
        if (s > -5) puts_(b, ", ");
        write_score(b, m->shift_best[(s + 12) % 12]);
    }
    puts_(b, "]}, \"chroma_tempo\": {\"tempo_ratios\": ");
    write_ratios(b, MATCH.melodic_tempo_ratios, MATCH.n_melodic_tempo_ratios);
    puts_(b, ", \"best_window_similarity\": [");
    for (int i = 0; i < MATCH.n_melodic_tempo_ratios && i < MAX_TEMPO_RATIOS; i++) {
        if (i) puts_(b, ", ");
        write_score(b, m->ratio_best[i]);
    }
    puts_(b, "]}}");
    return 0;
}

/* ---- exports ----------------------------------------------------------- */

EXPORT void *motif_malloc(size_t bytes) { return malloc(bytes ? bytes : 1); }

EXPORT void motif_free(void *p) { free(p); }

/* Compare two mono float32 signals at 22050 Hz. Names are shown as given; sha256 values are
 * the hex digests of the original file bytes computed by the host, or NULL when unknown.
 * Returns a NUL-terminated JSON report to be released with motif_free, or NULL when out of memory. */
EXPORT char *motif_compare(const float *query, size_t n_query, const char *query_name, const char *query_sha256,
                           const float *reference, size_t n_reference, const char *reference_name,
                           const char *reference_sha256, const char *decoder) {
    Arena arena = {0};
    Input q = {.name = query_name, .sha256 = query_sha256};
    Input r = {.name = reference_name, .sha256 = reference_sha256};
    SharedAudio shared;
    Musical musical;
    int identical = query_sha256 && reference_sha256 && query_sha256[0] && strcmp(query_sha256, reference_sha256) == 0;
    int compared = 0;
    int err = analyse(&arena, &q, query, n_query) || analyse(&arena, &r, reference, n_reference);
    if (!err && q.completed && r.completed && !identical) {
        err = shared_audio_any_speed(&arena, &q.peaks, &r.sorted, &FP, &MATCH, &shared) ||
              musical_similarity(&arena, &q.chroma, &r.chroma, &MC, &MATCH, &musical);
        compared = 1;
    }
    if (err || arena.oom) {
        arena_free_all(&arena);
        return NULL;
    }

    const char *label, *statement;
    if (identical) label = "identical_file", statement = ASSESS_IDENTICAL;
    else if (!compared) label = "not_compared", statement = ASSESS_NOT_COMPARED;
    else if (shared.detected) label = "shared_audio", statement = ASSESS_SHARED;
    else if (musical.detected) label = "musical_similarity", statement = ASSESS_MUSICAL;
    else label = "no_evidence", statement = ASSESS_NONE;

    Buf b = {0};
    write_envelope(&b, decoder);
    key(&b, "query");
    write_summary(&b, &q);
    puts_(&b, ", ");
    key(&b, "reference");
    write_summary(&b, &r);
    puts_(&b, ", \"relationship\": {\"assessment\": {\"label\": ");
    str(&b, label);
    puts_(&b, ", \"statement\": ");
    str(&b, statement);
    puts_(&b, "}, \"shared_audio\": ");
    if (compared) write_shared(&b, &shared);
    else puts_(&b, "null");
    puts_(&b, ", \"musical_similarity\": ");
    if (compared) write_musical(&b, &musical);
    else puts_(&b, "null");
    puts_(&b, "}, \"diagnostics\": ");
    if (compared) err = write_diagnostics(&arena, &b, &q, &r, &musical);
    else puts_(&b, "null");
    puts_(&b, "}\n");
    arena_free_all(&arena);
    if (err || b.failed) {
        free(b.data);
        return NULL;
    }
    return b.data;
}
