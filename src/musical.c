/* Port of motif/features/musical.py: chroma, key and tempo. */
#include <math.h>
#include <string.h>

#include "motif.h"

#define ONSET_N_FFT 1024
#define ONSET_HOP 256
#define ONSET_SMOOTH_SECONDS 0.5

const char *const PITCH_CLASSES[12] = {"C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"};

static const double MAJOR_PROFILE[12] = {6.35, 2.23, 3.48, 2.33, 4.38, 4.09, 2.52, 5.19, 2.39, 3.66, 2.29, 2.88};
static const double MINOR_PROFILE[12] = {6.33, 2.68, 3.52, 5.38, 2.60, 3.53, 2.54, 4.75, 3.98, 2.69, 3.34, 3.17};

int chroma(Arena *a, const float *x, size_t n, const MusicalConfig *c, Chroma *out) {
    int n_fft = c->chroma_n_fft, hop = c->chroma_hop_length, half = n_fft / 2;
    int n_stft = frame_count(n, n_fft, hop);
    if (n_stft == 0) return -1;
    Rfft fft;
    int *pc = arena_alloc(a, (size_t)half + 1, sizeof(int));
    double *power = arena_alloc(a, (size_t)half + 1, sizeof(double));
    if (rfft_init(a, &fft, n_fft) || !pc || !power) return -1;
    /* Nearest pitch class of each bin in f_min..f_max, -1 outside; freqs as np.linspace. */
    double step = (SAMPLE_RATE / 2.0) / half;
    for (int k = 0; k <= half; k++) {
        double f = k * step;
        pc[k] = -1;
        if (f >= c->chroma_f_min && f <= c->chroma_f_max) {
            double midi = 69 + 12 * log2(f / 440.0);
            int r = (int)rint(midi) % 12;
            pc[k] = r < 0 ? r + 12 : r;
        }
    }
    int last_group = (int)floor((double)(n_stft - 1) * hop / SAMPLE_RATE / c->chroma_frame_seconds + 1e-9);
    int n_out = last_group + 1;
    double *total = arena_alloc(a, (size_t)12 * n_out, sizeof(double));
    int *counts = arena_alloc(a, (size_t)n_out, sizeof(int));
    float *data = arena_alloc(a, (size_t)12 * n_out, sizeof(float));
    if (!total || !counts || !data) return -1;
    for (int t = 0; t < n_stft; t++) {
        int g = (int)floor((double)t * hop / SAMPLE_RATE / c->chroma_frame_seconds + 1e-9);
        frame_power(&fft, x, (size_t)t * hop, power);
        float pcs[12] = {0};
        for (int k = 0; k <= half; k++)
            if (pc[k] >= 0) pcs[pc[k]] += (float)power[k];
        for (int k = 0; k < 12; k++) total[k * n_out + g] += pcs[k];
        counts[g]++;
    }
    double peak = 0;
    for (int i = 0; i < 12 * n_out; i++) {
        int g = i % n_out;
        total[i] /= counts[g] > 0 ? counts[g] : 1;
        if (total[i] > peak) peak = total[i];
    }
    out->data = data;
    out->frames = n_out;
    if (peak <= 0) return 0;
    for (int g = 0; g < n_out; g++) {
        double v[12], norm = 0;
        for (int k = 0; k < 12; k++) {
            v[k] = log1p(1000.0 * total[k * n_out + g] / peak);
            norm += v[k] * v[k];
        }
        norm = sqrt(norm);
        if (norm > 0)
            for (int k = 0; k < 12; k++) data[k * n_out + g] = (float)(v[k] / norm);
    }
    return 0;
}

void estimate_key(const Chroma *ch, const MusicalConfig *c, KeyEstimate *out) {
    double mean[12], lo = INFINITY, hi = -INFINITY;
    int any = 0;
    memset(out, 0, sizeof *out);
    for (int k = 0; k < 12; k++) {
        double s = 0;
        for (int t = 0; t < ch->frames; t++) s += ch->data[k * ch->frames + t];
        mean[k] = ch->frames ? s / ch->frames : 0;
        if (mean[k] > 0) any = 1;
        if (mean[k] < lo) lo = mean[k];
        if (mean[k] > hi) hi = mean[k];
    }
    if (!any || hi == lo) return;
    double vm = 0, vn = 0, scores[24];
    for (int k = 0; k < 12; k++) vm += mean[k] / 12;
    for (int k = 0; k < 12; k++) vn += (mean[k] - vm) * (mean[k] - vm);
    vn = sqrt(vn);
    for (int key = 0; key < 24; key++) {
        const double *base = key < 12 ? MAJOR_PROFILE : MINOR_PROFILE;
        int tonic = key % 12;
        double p[12], pm = 0, pn = 0, dot = 0;
        for (int k = 0; k < 12; k++) p[k] = base[(k - tonic + 12) % 12]; /* np.roll(profile, tonic) */
        for (int k = 0; k < 12; k++) pm += p[k] / 12;
        for (int k = 0; k < 12; k++) {
            pn += (p[k] - pm) * (p[k] - pm);
            dot += (p[k] - pm) * (mean[k] - vm);
        }
        scores[key] = dot / (sqrt(pn) * vn);
    }
    int best = 0, second = -1;
    for (int key = 1; key < 24; key++)
        if (scores[key] > scores[best]) best = key;
    for (int key = 0; key < 24; key++)
        if (key != best && (second < 0 || scores[key] > scores[second])) second = key;
    out->has_scores = 1;
    out->best = best;
    out->second = second;
    out->best_score = scores[best];
    out->second_score = scores[second];
    out->known = scores[best] >= c->key_min_correlation;
}

/* Half-wave-rectified log-spectral flux minus its local mean; one value per ONSET_HOP. */
static double *onset_envelope(Arena *a, const float *x, size_t n, int *size) {
    int frames = frame_count(n, ONSET_N_FFT, ONSET_HOP), bins = ONSET_N_FFT / 2 + 1;
    Rfft fft;
    double *power = arena_alloc(a, (size_t)bins, sizeof(double));
    float *prev = arena_alloc(a, (size_t)bins, sizeof(float));
    float *cur = arena_alloc(a, (size_t)bins, sizeof(float));
    double *flux = arena_alloc(a, (size_t)frames + 1, sizeof(double));
    double *env = arena_alloc(a, (size_t)frames + 1, sizeof(double));
    if (rfft_init(a, &fft, ONSET_N_FFT) || !power || !prev || !cur || !flux || !env) return NULL;
    for (int t = 0; t < frames; t++) {
        frame_power(&fft, x, (size_t)t * ONSET_HOP, power);
        double sum = 0; /* NumPy sums the float32 diffs pairwise; a double sum rounds closer to that */
        for (int k = 0; k < bins; k++) {
            cur[k] = log1pf(1000.0f * (float)power[k]);
            float d = t ? cur[k] - prev[k] : 0.0f;
            if (d > 0) sum += d;
        }
        flux[t] = (float)sum;
        float *swap = prev;
        prev = cur;
        cur = swap;
    }
    /* np.convolve(flux, ones(w)/w, "same") */
    int w = (int)rint(ONSET_SMOOTH_SECONDS * SAMPLE_RATE / ONSET_HOP);
    if (w < 1) w = 1;
    int shift = (w - 1) / 2;
    for (int i = 0; i < frames; i++) {
        double s = 0;
        for (int j = 0; j < w; j++) {
            int m = i + shift - j;
            if (m >= 0 && m < frames) s += flux[m] / w;
        }
        double v = flux[i] - s;
        env[i] = v > 0 ? v : 0;
    }
    *size = frames;
    return env;
}

int estimate_tempo(Arena *a, const float *x, size_t n, const MusicalConfig *c, TempoEstimate *out) {
    memset(out, 0, sizeof *out);
    out->state = TEMPO_NO_ONSETS;
    int size = 0;
    double *env = onset_envelope(a, x, n, &size);
    if (!env) return -1;
    double fps = (double)SAMPLE_RATE / ONSET_HOP, mean = 0;
    int min_lag = (int)floor(60.0 / c->tempo_max_bpm * fps);
    int max_lag = (int)ceil(60.0 / c->tempo_min_bpm * fps);
    for (int i = 0; i < size; i++) mean += env[i];
    mean = size ? mean / size : 0;
    int nonzero = 0;
    for (int i = 0; i < size; i++) {
        env[i] -= mean;
        if (env[i] != 0) nonzero = 1;
    }
    if (size < 2 * max_lag + 2 || !nonzero) return 0;
    /* Direct autocorrelation over the lags used; equal to the FFT form up to rounding. */
    int first = min_lag > 1 ? min_lag : 1, hi = max_lag + 1 < size ? max_lag + 1 : size - 1;
    double zero = 0;
    for (int i = 0; i < size; i++) zero += env[i] * env[i];
    double acf[512];
    if (hi >= 512) return -1;
    for (int lag = first - 1; lag <= hi; lag++) {
        double s = 0;
        for (int i = 0; i + lag < size; i++) s += env[i] * env[i + lag];
        acf[lag] = s / zero;
    }
    int best = first;
    for (int lag = first + 1; lag <= max_lag; lag++)
        if (acf[lag] > acf[best]) best = lag;
    double strength = acf[best], refined = best;
    if (best > 0 && best < size - 1) {
        double p = acf[best - 1], q = acf[best], r = acf[best + 1], denom = p - 2 * q + r;
        if (denom != 0) refined = best + 0.5 * (p - r) / denom;
    }
    double bpm = 60.0 * fps / refined;
    static const double factors[2] = {0.5, 2.0};
    for (int i = 0; i < 2; i++) {
        double alt = bpm * factors[i];
        if (alt >= c->tempo_min_bpm && alt <= c->tempo_max_bpm) out->alternatives[out->n_alternatives++] = alt;
    }
    out->bpm = bpm;
    out->strength = strength;
    out->state = strength < c->tempo_min_strength ? TEMPO_WEAK : TEMPO_KNOWN;
    return 0;
}
