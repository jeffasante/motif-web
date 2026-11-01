/* Port of motif/features/fingerprint.py. */
#include <math.h>
#include <stdlib.h>

#include "motif.h"

static int by_strength(const void *pa, const void *pb) {
    const PeakRow *a = pa, *b = pb;
    if (a->level != b->level) return a->level > b->level ? -1 : 1;
    if (a->frame != b->frame) return a->frame < b->frame ? -1 : 1;
    return (a->bin > b->bin) - (a->bin < b->bin);
}

static int by_time(const void *pa, const void *pb) {
    const PeakRow *a = pa, *b = pb;
    if (a->frame != b->frame) return a->frame < b->frame ? -1 : 1;
    return (a->bin > b->bin) - (a->bin < b->bin);
}

/* Keep the strongest peaks_per_second peaks in each one-second bucket, in (frame, bin) order. */
static void cap_density(Peaks *p, const FingerprintConfig *c) {
    size_t kept = 0;
    for (size_t first = 0; first < p->n;) {
        size_t last = first;
        while (last < p->n && p->rows[last].bucket == p->rows[first].bucket) last++;
        qsort(p->rows + first, last - first, sizeof *p->rows, by_strength);
        size_t keep = last - first < (size_t)c->peaks_per_second ? last - first : (size_t)c->peaks_per_second;
        qsort(p->rows + first, keep, sizeof *p->rows, by_time);
        for (size_t i = 0; i < keep; i++) p->rows[kept++] = p->rows[first + i];
        first = last;
    }
    p->n = kept;
}

int find_peaks(Arena *a, const float *x, size_t n, const FingerprintConfig *c, Peaks *out) {
    int n_frames = frame_count(n, c->n_fft, c->hop_length);
    double hz_per_bin = (double)SAMPLE_RATE / c->n_fft;
    int lo = (int)ceil(c->f_min / hz_per_bin), hi = (int)floor(c->f_max / hz_per_bin);
    int nb = hi - lo + 1, radius = c->neighborhood_frames, rb = c->neighborhood_bins;
    int ring = 2 * radius + 1;
    int frames_per_bucket = (int)nearbyint((double)SAMPLE_RATE / c->hop_length);
    if (frames_per_bucket < 1) frames_per_bucket = 1;

    Rfft fft;
    double *power = arena_alloc(a, (size_t)c->n_fft / 2 + 1, sizeof(double));
    double *level = arena_alloc(a, (size_t)ring * nb, sizeof(double));
    if (rfft_init(a, &fft, c->n_fft) || !power || !level) return -1;
    /* |X| of a full-scale sine is sum(window)/2; this offset makes that 0 dB. */
    double db_offset = 10 * log10(4.0 / (fft.window_sum * fft.window_sum));

    PeakRow *rows = NULL;
    size_t count = 0, cap = 0;
    for (int f = 0; f < n_frames + radius; f++) {
        if (f < n_frames) {
            frame_power(&fft, x, (size_t)f * c->hop_length, power);
            double *row = level + (size_t)(f % ring) * nb;
            for (int b = 0; b < nb; b++) {
                /* float32 power and log, then float64, as the NumPy dtypes go */
                float pw = (float)power[lo + b];
                float db = 10.0f * (float)log10((double)(pw > 1e-20f ? pw : 1e-20f));
                row[b] = (double)db + db_offset;
            }
        }
        int t = f - radius;
        if (t < 0) continue;
        const double *row = level + (size_t)(t % ring) * nb;
        int t0 = t - radius < 0 ? 0 : t - radius, t1 = t + radius >= n_frames ? n_frames - 1 : t + radius;
        for (int b = 0; b < nb; b++) {
            double v = row[b];
            if (v < c->min_peak_dbfs) continue;
            int b0 = b - rb < 0 ? 0 : b - rb, b1 = b + rb >= nb ? nb - 1 : b + rb, peak = 1;
            for (int u = t0; u <= t1 && peak; u++) {
                const double *other = level + (size_t)(u % ring) * nb;
                for (int k = b0; k <= b1; k++)
                    if (other[k] > v) { peak = 0; break; }
            }
            if (!peak) continue;
            rows = vec_reserve(a, rows, &cap, count + 1, sizeof *rows);
            if (!rows) return -1;
            rows[count++] = (PeakRow){t, lo + b, t / frames_per_bucket, (float)v};
        }
    }
    out->rows = rows;
    out->n = count;
    out->frames_per_bucket = frames_per_bucket;
    out->n_buckets = n_frames > 0 ? (n_frames - 1) / frames_per_bucket + 1 : 0;
    out->found_per_bucket = arena_alloc(a, (size_t)out->n_buckets + 1, sizeof(int));
    if (!out->found_per_bucket) return -1;
    for (size_t i = 0; i < count; i++) out->found_per_bucket[rows[i].bucket]++;
    cap_density(out, c);
    return 0;
}

int hash_peaks(Arena *a, const Peaks *peaks, const FingerprintConfig *c, double speed_ratio, Hashes *out) {
    return hash_peaks_scaled(a, peaks, c, speed_ratio, speed_ratio, out);
}

int hash_peaks_scaled(Arena *a, const Peaks *peaks, const FingerprintConfig *c, double freq_ratio, double time_ratio,
                      Hashes *out) {
    const PeakRow *p = peaks->rows;
    Hash *items = NULL;
    size_t count = 0, cap = 0;
    for (size_t i = 0; i < peaks->n; i++) {
        int taken = 0;
        for (size_t j = i + 1; j < peaks->n && taken < c->fan_out; j++) {
            int gap = p[j].frame - p[i].frame;
            if (gap > c->target_max_frames) break;
            if (gap < c->target_min_frames || abs(p[j].bin - p[i].bin) > c->target_max_bins) continue;
            taken++;
            double f1 = rint(p[i].bin / freq_ratio), f2 = rint(p[j].bin / freq_ratio);
            double g = rint(gap * time_ratio);
            if (f1 < 0 || f1 >= 1024 || f2 < 0 || f2 >= 1024 || g < 1 || g > 255) continue;
            items = vec_reserve(a, items, &cap, count + 1, sizeof *items);
            if (!items) return -1;
            items[count++] = (Hash){(uint32_t)f1 << 18 | (uint32_t)f2 << 8 | (uint32_t)g,
                                    (int32_t)rint(p[i].frame * time_ratio)};
        }
    }
    out->items = items;
    out->n = count;
    return 0;
}
